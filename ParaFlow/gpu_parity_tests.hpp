// gpu_parity_tests.hpp — CPU-vs-GPU RK4 parity / smoke tests (DEBUG-ONLY)
//
// NOT a standalone header. It is textually #included exactly once by
// ParaFlow.cpp, from *inside* its `#ifdef OSUFLOW_ENABLE_CUDA` block, right
// after the gpu_smoke_env_* helpers. It depends on that surrounding context
// (ParaFlow.hpp, block.hpp, OSUFlow + GPU headers, and gpu_smoke_env_int/
// gpu_smoke_env_double) and must not be #included anywhere else.
//
// Every function early-returns unless the OSUFLOW_GPU_SMOKE env var is set —
// they exist only for CPU/GPU parity verification and never run in production.

// Host copy of FieldLine.C:geodesic_step — needed because the CPU version is
// a protected member of vtCFieldLine. Kept numerically identical (same order,
// same rotation matrix construction).
static bool host_geodesic_step(VECTOR3 pt_src, double r_src,
                               VECTOR3 h_vel,  double v_vel,
                               double fdt,
                               VECTOR3& pt_dst, double& r_dst)
{
    double vel_mag = h_vel.GetMag() * fdt;
    if (vel_mag == 0.0) { pt_dst = pt_src; r_dst = r_src; return true; }
    double r_new = r_src + v_vel * fdt;
    if (r_new <= 0.0) return false;
    double  omega    = vel_mag / r_src;
    VECTOR3 normal   = cross(pt_src, h_vel);
    MATRIX3 rotate_m = rotate_matrix_axis(normal, omega);
    pt_dst = rotate_m * pt_src;
    pt_dst.Normalize();
    pt_dst.scale(r_new);
    r_dst = r_new;
    return true;
}

// Host RK4 one-step, mirroring vtCFieldLine::MPASO_rk4 (FieldLine.C:247-316).
// Returns 0 on success, nonzero stage code on failure (1-4 = at_phys failed,
// 5-7 = geodesic failed, 8 = bad radius).
static int host_rk4_step(CVectorField* field,
                         VECTOR3& pt, int& /*fromCell*/, double t, double dt)
{
    // Force a fresh full-scan phys_to_cell on every at_phys call by passing
    // fromCell = -1. Trades speed for correctness — gives a ground-truth CPU
    // trajectory that the GPU (which uses a neighbor-ring walk) must match.
    VECTOR4 vel;
    PointInfo ci; ci.phyCoord = pt;
    double r0 = pt.GetMag();
    if (r0 <= 0.0) return 8;

    PointInfo ci_tmp = ci;
    if (field->at_phys(-1, pt, ci_tmp, t, vel, nullptr) != 1) return 1;
    VECTOR3 k1_h(vel[0], vel[1], vel[2]); double k1_v = vel[3];

    VECTOR3 pt1; double r1;
    if (!host_geodesic_step(pt, r0, k1_h, k1_v, 0.5 * dt, pt1, r1)) return 5;
    ci_tmp = ci; ci_tmp.phyCoord = pt1;
    if (field->at_phys(-1, pt1, ci_tmp, t, vel, nullptr) != 1) return 2;
    VECTOR3 k2_h(vel[0], vel[1], vel[2]); double k2_v = vel[3];

    VECTOR3 pt2; double r2;
    if (!host_geodesic_step(pt, r0, k2_h, k2_v, 0.5 * dt, pt2, r2)) return 6;
    ci_tmp = ci; ci_tmp.phyCoord = pt2;
    if (field->at_phys(-1, pt2, ci_tmp, t, vel, nullptr) != 1) return 3;
    VECTOR3 k3_h(vel[0], vel[1], vel[2]); double k3_v = vel[3];

    VECTOR3 pt3; double r3;
    if (!host_geodesic_step(pt, r0, k3_h, k3_v, dt, pt3, r3)) return 7;
    ci_tmp = ci; ci_tmp.phyCoord = pt3;
    if (field->at_phys(-1, pt3, ci_tmp, t, vel, nullptr) != 1) return 4;
    VECTOR3 k4_h(vel[0], vel[1], vel[2]); double k4_v = vel[3];

    VECTOR3 v_avg_h;
    for (int i = 0; i < 3; ++i)
        v_avg_h[i] = (k1_h[i] + 2.0*k2_h[i] + 2.0*k3_h[i] + k4_h[i]) / 6.0;
    double v_avg_v = (k1_v + 2.0*k2_v + 2.0*k3_v + k4_v) / 6.0;

    VECTOR3 pt_final; double r_final;
    if (!host_geodesic_step(pt, r0, v_avg_h, v_avg_v, dt, pt_final, r_final)) return 9;

    pt = pt_final;
    return 0;
}

static void run_gpu_smoke_test(Block* b, int rank)
{
    if (std::getenv("OSUFLOW_GPU_SMOKE") == nullptr) return;
    if (rank != 0) return;
    if (b == nullptr || b->currentSeeds.empty()) {
        std::fprintf(stderr, "[GPU parity] skipped: no seeds on rank 0 block\n");
        return;
    }

    CVectorField* field = b->osuflow->GetFlowField();
    MPASOGrid* grid = dynamic_cast<MPASOGrid*>(field->GetGrid());
    Solution*  pSol = field->GetHorizontalSolution();
    Solution*  vSol = field->GetVerticalSolution();
    if (!grid || !pSol || !vSol) {
        std::fprintf(stderr, "[GPU parity] skipped: grid/solution not available\n");
        return;
    }

    const int n_particles = std::min<int>((int)b->currentSeeds.size(),
                                          gpu_smoke_env_int("OSUFLOW_GPU_SMOKE_PARTICLES", 8));
    const int n_steps     = gpu_smoke_env_int("OSUFLOW_GPU_SMOKE_STEPS", 10);
    const double dt       = gpu_smoke_env_double("OSUFLOW_GPU_SMOKE_DT", 60.0);
    const int nVertLevels = grid->getNVertLevels();

    // Nudge seeds slightly below earth_radius — raw surface seeds (|r| = R)
    // hit Grid.C's `if (|r| > earth_radius) return -1;` check as soon as RK4
    // averaging produces any net upward vertical component.
    const double earth_r = grid->getEarthRadius();
    const double shrink  = (earth_r - 1.0) / earth_r;    // 1m below surface — enough to clear the |r|>R boundary check without diving below shallow bathymetry
    std::vector<VECTOR3> seeds(n_particles);
    std::vector<int>     seed_cells(n_particles);
    for (int i = 0; i < n_particles; ++i) {
        VECTOR3 s = b->currentSeeds[i];
        s[0] *= shrink; s[1] *= shrink; s[2] *= shrink;
        seeds[i] = s;
        seed_cells[i] = 0;  // overwritten below once we have a CPU-locate result
    }

    // ---- CPU reference trajectory ----
    std::vector<VECTOR3> cpu_trace(n_particles * (n_steps + 1));
    std::vector<int>     cpu_alive(n_particles, 1);
    std::vector<int>     cpu_steps_taken(n_particles, 0);
    for (int p = 0; p < n_particles; ++p) {
        VECTOR3 pt = seeds[p];
        int fromCell = (p < (int)b->fromCells.size()) ? b->fromCells[p] : -1;
        cpu_trace[p * (n_steps + 1)] = pt;

        // One-shot probe: sample velocity at the seed so we can report it,
        // and capture ci.inCell so we can give the GPU a valid cell hint.
        {
            VECTOR4 v4;
            PointInfo ci; ci.phyCoord = pt; ci.fromCell = fromCell; ci.inCell = fromCell;
            int ok = field->at_phys(fromCell, pt, ci, 0.0, v4, nullptr);
            // ci.inCell is encoded as localCell * nVertLevels + vLevel
            if (ok == 1 && ci.inCell >= 0) {
                seed_cells[p] = ci.inCell / nVertLevels;
                // also refresh fromCell so CPU trace uses the valid hint
                fromCell = ci.inCell;
            }
            std::fprintf(stderr,
                "[GPU parity] seed %d pos=(%.3e,%.3e,%.3e) fromCell=%d at_phys=%d v=(%.3e,%.3e,%.3e,%.3e) seed_cell=%d\n",
                p, pt[0],pt[1],pt[2], fromCell, ok, v4[0],v4[1],v4[2],v4[3], seed_cells[p]);
        }

        double t = 0.0;
        for (int s = 0; s < n_steps; ++s) {
            int fail = host_rk4_step(field, pt, fromCell, t, dt);
            if (fail != 0) {
                std::fprintf(stderr, "[GPU parity] seed %d step %d FAILED stage=%d\n",
                             p, s, fail);
                cpu_alive[p] = 0;
                for (int k = s + 1; k <= n_steps; ++k)
                    cpu_trace[p * (n_steps + 1) + k] = pt;
                break;
            }
            cpu_steps_taken[p] = s + 1;
            // Streamline semantics: t stays at the initial timestamp (mirrors
            // vtCFieldLine::MPASO_rk4, which only advances t under UNSTEADY).
            cpu_trace[p * (n_steps + 1) + s + 1] = pt;
        }
        std::fprintf(stderr, "[GPU parity] seed %d cpu_steps_taken=%d alive=%d\n",
                     p, cpu_steps_taken[p], cpu_alive[p]);
    }

    // ---- GPU trajectory ----
    std::vector<VECTOR3> gpu_trace(n_particles * (n_steps + 1));
    std::fprintf(stderr, "[GPU parity] launching: P=%d, steps=%d, dt=%.1f\n",
                 n_particles, n_steps, dt);
    mpaso_gpu_host::GPUBlockContext* gpu_context =
        mpaso_gpu_host::CreateGPUBlockContext(grid);
    if (gpu_context == nullptr) {
        std::fprintf(stderr, "[GPU parity] FAILED: could not create GPU context\n");
        return;
    }
    mpaso_gpu_host::UploadGPUVelocityWindow(
        gpu_context, grid, pSol, vSol, /*use_real_timestamps=*/false);
    mpaso_gpu_host::TraceParticlesOnGPUContext(
        gpu_context,
        seeds.data(), seed_cells.data(), n_particles,
        n_steps, dt, /*t_start=*/0.0, /*use_euler=*/false,
        gpu_trace.data());
    mpaso_gpu_host::DestroyGPUBlockContext(gpu_context);

    // ---- compare ----
    const double earth = earth_r;
    double max_abs = 0.0;
    double max_rel = 0.0;
    int    max_p = -1, max_s = -1;
    int    compared_count = 0;
    double first_seed_mag = seeds[0].GetMag();
    for (int p = 0; p < n_particles; ++p) {
        if (!cpu_alive[p]) continue;
        for (int s = 1; s <= n_steps; ++s) {
            const VECTOR3& a = cpu_trace[p * (n_steps + 1) + s];
            const VECTOR3& g = gpu_trace[p * (n_steps + 1) + s];
            double dx = a[0]-g[0], dy = a[1]-g[1], dz = a[2]-g[2];
            double d  = std::sqrt(dx*dx + dy*dy + dz*dz);
            double mag = std::sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
            double rel = (mag > 0.0) ? d / mag : d;
            ++compared_count;
            if (d > max_abs)   { max_abs = d; max_p = p; max_s = s; }
            if (rel > max_rel) { max_rel = rel; }
        }
    }

    const double rel_thresh = gpu_smoke_env_double("OSUFLOW_GPU_SMOKE_REL_THRESH", 1e-5);
    std::fprintf(stderr, "[GPU parity] compared=%d max_abs=%.3e m  max_rel=%.3e  (earth=%.3e, seed0|r|=%.3e)\n",
                 compared_count, max_abs, max_rel, earth, first_seed_mag);
    if (compared_count <= 0) {
        std::fprintf(stderr, "[GPU parity] FAIL: no comparable CPU/GPU samples\n");
    } else if (max_rel < rel_thresh) {
        std::fprintf(stderr, "[GPU parity] PASS (rel err < %.1e)\n", rel_thresh);
    } else {
        std::fprintf(stderr, "[GPU parity] FAIL at p=%d step=%d\n", max_p, max_s);
        if (max_p >= 0) {
            const VECTOR3& a = cpu_trace[max_p * (n_steps + 1) + max_s];
            const VECTOR3& g = gpu_trace[max_p * (n_steps + 1) + max_s];
            std::fprintf(stderr, "  cpu=(%.9e,%.9e,%.9e)\n", a[0],a[1],a[2]);
            std::fprintf(stderr, "  gpu=(%.9e,%.9e,%.9e)\n", g[0],g[1],g[2]);
        }
    }
}

// Host RK4 one-step for pathlines — like host_rk4_step but advances sample
// time per stage, mirroring vtCFieldLine::MPASO_rk4 UNSTEADY branch:
//   stage1 at t, stages 2/3 at t+dt/2, stage 4 at t+dt.
static int host_rk4_step_pathline(CVectorField* field,
                                  VECTOR3& pt, double t, double dt)
{
    VECTOR4 vel;
    PointInfo ci; ci.phyCoord = pt;
    double r0 = pt.GetMag();
    if (r0 <= 0.0) return 8;

    PointInfo ci_tmp = ci;
    if (field->at_phys(-1, pt, ci_tmp, t, vel, nullptr) != 1) return 1;
    VECTOR3 k1_h(vel[0], vel[1], vel[2]); double k1_v = vel[3];

    VECTOR3 pt1; double r1;
    if (!host_geodesic_step(pt, r0, k1_h, k1_v, 0.5 * dt, pt1, r1)) return 5;
    ci_tmp = ci; ci_tmp.phyCoord = pt1;
    if (field->at_phys(-1, pt1, ci_tmp, t + 0.5 * dt, vel, nullptr) != 1) return 2;
    VECTOR3 k2_h(vel[0], vel[1], vel[2]); double k2_v = vel[3];

    VECTOR3 pt2; double r2;
    if (!host_geodesic_step(pt, r0, k2_h, k2_v, 0.5 * dt, pt2, r2)) return 6;
    ci_tmp = ci; ci_tmp.phyCoord = pt2;
    if (field->at_phys(-1, pt2, ci_tmp, t + 0.5 * dt, vel, nullptr) != 1) return 3;
    VECTOR3 k3_h(vel[0], vel[1], vel[2]); double k3_v = vel[3];

    VECTOR3 pt3; double r3;
    if (!host_geodesic_step(pt, r0, k3_h, k3_v, dt, pt3, r3)) return 7;
    ci_tmp = ci; ci_tmp.phyCoord = pt3;
    if (field->at_phys(-1, pt3, ci_tmp, t + dt, vel, nullptr) != 1) return 4;
    VECTOR3 k4_h(vel[0], vel[1], vel[2]); double k4_v = vel[3];

    VECTOR3 v_avg_h;
    for (int i = 0; i < 3; ++i)
        v_avg_h[i] = (k1_h[i] + 2.0*k2_h[i] + 2.0*k3_h[i] + k4_h[i]) / 6.0;
    double v_avg_v = (k1_v + 2.0*k2_v + 2.0*k3_v + k4_v) / 6.0;

    VECTOR3 pt_final; double r_final;
    if (!host_geodesic_step(pt, r0, v_avg_h, v_avg_v, dt, pt_final, r_final)) return 9;

    pt = pt_final;
    return 0;
}

static void run_gpu_pathline_parity_test(Block* b, int rank)
{
    if (std::getenv("OSUFLOW_GPU_SMOKE") == nullptr) return;
    if (rank != 0) return;
    if (b == nullptr || b->currentSeeds.empty()) {
        std::fprintf(stderr, "[GPU pathline parity] skipped: no seeds on rank 0 block\n");
        return;
    }

    CVectorField* field = b->osuflow->GetFlowField();
    MPASOGrid* grid = dynamic_cast<MPASOGrid*>(field->GetGrid());
    Solution*  pSol = field->GetHorizontalSolution();
    Solution*  vSol = field->GetVerticalSolution();
    if (!grid || !pSol || !vSol) {
        std::fprintf(stderr, "[GPU pathline parity] skipped: grid/solution not available\n");
        return;
    }

    const int n_particles = std::min<int>((int)b->currentSeeds.size(),
                                          gpu_smoke_env_int("OSUFLOW_GPU_SMOKE_PARTICLES", 8));
    const int n_steps     = gpu_smoke_env_int("OSUFLOW_GPU_SMOKE_STEPS", 10);
    const double dt       = gpu_smoke_env_double("OSUFLOW_GPU_SMOKE_DT", 60.0);
    const int nVertLevels = grid->getNVertLevels();

    // Pick t_start from the loaded window — use timestamps[0] so at_phys has
    // both endpoints of the first interval available.
    const std::vector<double>& tsWin = grid->getZTopTimestamps();
    double t_start = tsWin.empty() ? 0.0 : tsWin[0];

    const double earth_r = grid->getEarthRadius();
    const double shrink  = (earth_r - 1.0) / earth_r;
    std::vector<VECTOR3> seeds(n_particles);
    std::vector<int>     seed_cells(n_particles);
    std::vector<double>  seed_t_start(n_particles, t_start);
    for (int i = 0; i < n_particles; ++i) {
        VECTOR3 s = b->currentSeeds[i];
        s[0] *= shrink; s[1] *= shrink; s[2] *= shrink;
        seeds[i] = s;
        seed_cells[i] = 0;
    }

    // ---- CPU reference trajectory (time-varying) ----
    std::vector<VECTOR3> cpu_trace(n_particles * (n_steps + 1));
    std::vector<int>     cpu_alive(n_particles, 1);
    std::vector<int>     cpu_steps_taken(n_particles, 0);
    for (int p = 0; p < n_particles; ++p) {
        VECTOR3 pt = seeds[p];
        int fromCell = (p < (int)b->fromCells.size()) ? b->fromCells[p] : -1;
        cpu_trace[p * (n_steps + 1)] = pt;

        {
            VECTOR4 v4;
            PointInfo ci; ci.phyCoord = pt; ci.fromCell = fromCell; ci.inCell = fromCell;
            int ok = field->at_phys(fromCell, pt, ci, t_start, v4, nullptr);
            if (ok == 1 && ci.inCell >= 0) {
                seed_cells[p] = ci.inCell / nVertLevels;
            }
            std::fprintf(stderr,
                "[GPU pathline parity] seed %d pos=(%.3e,%.3e,%.3e) at_phys=%d v=(%.3e,%.3e,%.3e,%.3e) seed_cell=%d t_start=%.3e\n",
                p, pt[0],pt[1],pt[2], ok, v4[0],v4[1],v4[2],v4[3], seed_cells[p], t_start);
        }

        double t = t_start;
        for (int s = 0; s < n_steps; ++s) {
            int fail = host_rk4_step_pathline(field, pt, t, dt);
            if (fail != 0) {
                std::fprintf(stderr, "[GPU pathline parity] seed %d step %d FAILED stage=%d\n",
                             p, s, fail);
                cpu_alive[p] = 0;
                for (int k = s + 1; k <= n_steps; ++k)
                    cpu_trace[p * (n_steps + 1) + k] = pt;
                break;
            }
            t += dt;
            cpu_steps_taken[p] = s + 1;
            cpu_trace[p * (n_steps + 1) + s + 1] = pt;
        }
        std::fprintf(stderr, "[GPU pathline parity] seed %d cpu_steps_taken=%d alive=%d final_t=%.3e\n",
                     p, cpu_steps_taken[p], cpu_alive[p], t_start + cpu_steps_taken[p] * dt);
    }

    // ---- GPU trajectory ----
    std::vector<VECTOR3> gpu_trace(n_particles * (n_steps + 1));
    std::vector<double>  gpu_final_time(n_particles, 0.0);
    std::vector<int>     gpu_steps_taken(n_particles, 0);
    std::vector<int>     gpu_final_cell(n_particles, -1);
    std::vector<int>     gpu_seed_max_steps(n_particles, n_steps);
    std::fprintf(stderr, "[GPU pathline parity] launching: P=%d, steps=%d, dt=%.1f, t_start=%.3e\n",
                 n_particles, n_steps, dt, t_start);
    // save_interval=1, max_saved_points=n_steps+1: store every step so the
    // step-by-step CPU/GPU comparison below works with the same row stride.
    mpaso_gpu_host::GPUBlockContext* gpu_context =
        mpaso_gpu_host::CreateGPUBlockContext(grid);
    if (gpu_context == nullptr) {
        std::fprintf(stderr, "[GPU pathline parity] FAILED: could not create GPU context\n");
        return;
    }
    mpaso_gpu_host::UploadGPUVelocityWindow(
        gpu_context, grid, pSol, vSol, /*use_real_timestamps=*/true);
    mpaso_gpu_host::TracePathlineBatchOnGPUContext(
        gpu_context,
        seeds.data(), seed_cells.data(), seed_t_start.data(),
        gpu_seed_max_steps.data(),
        n_particles, n_steps, dt,
        gpu_trace.data(),
        gpu_final_time.data(),
        gpu_steps_taken.data(),
        gpu_final_cell.data(),
        /*save_interval=*/1,
        /*max_saved_points=*/n_steps + 1);
    mpaso_gpu_host::DestroyGPUBlockContext(gpu_context);

    // ---- compare ----
    double max_abs = 0.0, max_rel = 0.0;
    int    max_p = -1, max_s = -1;
    int    compared_count = 0;
    bool   step_mismatch = false;
    for (int p = 0; p < n_particles; ++p) {
        if (!cpu_alive[p]) continue;
        int steps_cmp = std::min(cpu_steps_taken[p], gpu_steps_taken[p]);
        std::fprintf(stderr, "[GPU pathline parity] seed %d cpu_steps=%d gpu_steps=%d gpu_final_t=%.3e\n",
                     p, cpu_steps_taken[p], gpu_steps_taken[p], gpu_final_time[p]);
        if (gpu_steps_taken[p] != cpu_steps_taken[p]) step_mismatch = true;
        for (int s = 1; s <= steps_cmp; ++s) {
            const VECTOR3& a = cpu_trace[p * (n_steps + 1) + s];
            const VECTOR3& g = gpu_trace[p * (n_steps + 1) + s];
            double dx = a[0]-g[0], dy = a[1]-g[1], dz = a[2]-g[2];
            double d  = std::sqrt(dx*dx + dy*dy + dz*dz);
            double mag = std::sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
            double rel = (mag > 0.0) ? d / mag : d;
            ++compared_count;
            if (d > max_abs)   { max_abs = d; max_p = p; max_s = s; }
            if (rel > max_rel) { max_rel = rel; }
        }
    }

    const double rel_thresh = gpu_smoke_env_double("OSUFLOW_GPU_SMOKE_REL_THRESH", 1e-5);
    std::fprintf(stderr, "[GPU pathline parity] compared=%d max_abs=%.3e m  max_rel=%.3e\n",
                 compared_count, max_abs, max_rel);
    if (compared_count <= 0) {
        std::fprintf(stderr, "[GPU pathline parity] FAIL: no comparable CPU/GPU samples\n");
    } else if (step_mismatch) {
        std::fprintf(stderr, "[GPU pathline parity] FAIL: CPU/GPU step counts differ\n");
    } else if (max_rel < rel_thresh) {
        std::fprintf(stderr, "[GPU pathline parity] PASS (rel err < %.1e)\n", rel_thresh);
    } else {
        std::fprintf(stderr, "[GPU pathline parity] FAIL at p=%d step=%d\n", max_p, max_s);
        if (max_p >= 0) {
            const VECTOR3& a = cpu_trace[max_p * (n_steps + 1) + max_s];
            const VECTOR3& g = gpu_trace[max_p * (n_steps + 1) + max_s];
            std::fprintf(stderr, "  cpu=(%.9e,%.9e,%.9e)\n", a[0],a[1],a[2]);
            std::fprintf(stderr, "  gpu=(%.9e,%.9e,%.9e)\n", g[0],g[1],g[2]);
        }
    }
}
