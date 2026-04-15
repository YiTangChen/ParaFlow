#include "MPASOGPUTracer.h"

#include "Grid.h"
#include "Solution.h"
#include "MPASODeviceField.h"
#include "Kernel/MPASOTracerKernels.h"

#include <cuda_runtime.h>
#include <cstring>
#include <vector>

// VECTOR3 stores a single `double vec[3]` and no virtuals, so it is
// standard-layout with size == sizeof(mpaso_vec3) == 24 bytes. We exploit that
// to reinterpret_cast arrays of VECTOR3 to arrays of mpaso_vec3 without any
// per-element copying. Any future change to VECTOR3's layout would break this
// assumption — the static_assert below catches it at compile time.
static_assert(sizeof(VECTOR3) == sizeof(mpaso_vec3),
              "VECTOR3 and mpaso_vec3 must be bit-compatible for reinterpret_cast");

namespace {

// Pack Solution's jagged [ts][node] layout into a single flat [ts * node]
// buffer of mpaso_vec3 that uploadVelocityWindow expects.
std::vector<mpaso_vec3> flattenSolution(Solution* sol, int n_vert_nodes, int n_timesteps)
{
    std::vector<mpaso_vec3> out(static_cast<size_t>(n_timesteps) * n_vert_nodes);
    for (int ts = 0; ts < n_timesteps; ++ts) {
        const VECTOR3* src = sol->GetTimestepData(ts);
        std::memcpy(&out[static_cast<size_t>(ts) * n_vert_nodes],
                    src,
                    static_cast<size_t>(n_vert_nodes) * sizeof(mpaso_vec3));
    }
    return out;
}

} // namespace

namespace mpaso_gpu_host {

bool isAvailable()
{
    int count = 0;
    cudaError_t err = cudaGetDeviceCount(&count);
    if (err != cudaSuccess) {
        // Clear the sticky error so later cudaGetLastError calls aren't confused.
        (void)cudaGetLastError();
        return false;
    }
    return count > 0;
}


void TraceParticles(MPASOGrid*      grid,
                    Solution*       pSolution,
                    Solution*       vSolution,
                    const VECTOR3*  seeds,
                    const int*      seed_cell_id,
                    int             n_particles,
                    int             n_steps,
                    double          dt,
                    double          t_start,
                    bool            use_euler,
                    VECTOR3*        traces_out,
                    const int*      seed_max_steps,
                    const int*      seed_step_offset,
                    int*            steps_taken_out,
                    int*            final_cell_out,
                    int             save_interval,
                    int             max_saved_points,
                    int*            saved_counts_out)
{
    if (grid == nullptr || pSolution == nullptr || vSolution == nullptr) return;
    if (n_particles <= 0 || n_steps <= 0) return;

    MPASODeviceField field;
    field.n_cells            = grid->getNumCell();
    field.n_local_vertices   = grid->getNumLocalVert();
    field.n_max_edges        = grid->getNMaxEdges();
    field.n_vert_levels      = grid->getNVertLevels();
    field.n_vert_levels_p1   = field.n_vert_levels + 1;
    field.n_timesteps_loaded = grid->getNTimestepsLoaded();
    field.earth_radius       = grid->getEarthRadius();

    // ---- topology upload (one-shot) ----
    field.uploadTopology(
        reinterpret_cast<const mpaso_vec3*>(grid->getCellCoord()),
        reinterpret_cast<const mpaso_vec3*>(grid->getVertexCoord()),
        grid->getVerticesOnCell(),
        grid->getCellsOnCell(),
        grid->getNumVerticesOnCell(),
        grid->getMaxLevelCell());

    // ---- velocity / ztop window upload ----
    const int n_vert_nodes = field.n_local_vertices * field.n_vert_levels;
    const int n_ts         = field.n_timesteps_loaded;

    std::vector<mpaso_vec3> h_vel  = flattenSolution(pSolution, n_vert_nodes, n_ts);
    std::vector<mpaso_vec3> h_vvel = flattenSolution(vSolution, n_vert_nodes, n_ts);

    // zTop is already a flat [n_ts * n_vert_nodes] buffer owned by MPASOGrid.
    const double* h_ztop = grid->getZTopAll();

    // Timestamps: Solution / MPASOGrid carry real-time stamps internally; for
    // the stub roundtrip we just pass zeros — sample_velocity currently
    // ignores t_rel anyway.
    std::vector<double> h_timestamps(n_ts, 0.0);

    field.uploadVelocityWindow(h_vel.data(), h_vvel.data(), h_ztop, h_timestamps.data());

    // ---- launch ----
    mpaso_gpu::LaunchTracer(
        field,
        reinterpret_cast<const mpaso_vec3*>(seeds),
        seed_cell_id,
        n_particles,
        n_steps,
        dt,
        t_start,
        use_euler,
        reinterpret_cast<mpaso_vec3*>(traces_out),
        seed_max_steps,
        seed_step_offset,
        steps_taken_out,
        final_cell_out,
        save_interval,
        max_saved_points,
        saved_counts_out);

    field.release();
}

void TracePathlineBatch(MPASOGrid*      grid,
                        Solution*       pSolution,
                        Solution*       vSolution,
                        const VECTOR3*  seeds,
                        const int*      seed_cell_id,
                        const double*   seed_t_start,
                        const int*      seed_max_steps,
                        int             n_particles,
                        int             n_steps,
                        double          dt,
                        VECTOR3*        traces_out,
                        double*         final_time_out,
                        int*            steps_taken_out,
                        int*            final_cell_out)
{
    if (grid == nullptr || pSolution == nullptr || vSolution == nullptr) return;
    if (n_particles <= 0 || n_steps <= 0) return;

    MPASODeviceField field;
    field.n_cells            = grid->getNumCell();
    field.n_local_vertices   = grid->getNumLocalVert();
    field.n_max_edges        = grid->getNMaxEdges();
    field.n_vert_levels      = grid->getNVertLevels();
    field.n_vert_levels_p1   = field.n_vert_levels + 1;
    field.n_timesteps_loaded = grid->getNTimestepsLoaded();
    field.earth_radius       = grid->getEarthRadius();

    field.uploadTopology(
        reinterpret_cast<const mpaso_vec3*>(grid->getCellCoord()),
        reinterpret_cast<const mpaso_vec3*>(grid->getVertexCoord()),
        grid->getVerticesOnCell(),
        grid->getCellsOnCell(),
        grid->getNumVerticesOnCell(),
        grid->getMaxLevelCell());

    const int n_vert_nodes = field.n_local_vertices * field.n_vert_levels;
    const int n_ts         = field.n_timesteps_loaded;

    std::vector<mpaso_vec3> h_vel  = flattenSolution(pSolution, n_vert_nodes, n_ts);
    std::vector<mpaso_vec3> h_vvel = flattenSolution(vSolution, n_vert_nodes, n_ts);
    const double* h_ztop = grid->getZTopAll();

    // Real timestamps for time-varying blending. zTopTimestamps is set by
    // MPASOReader alongside the velocity window load.
    const std::vector<double>& ts = grid->getZTopTimestamps();
    std::vector<double> h_timestamps(n_ts, 0.0);
    for (int i = 0; i < n_ts && i < (int)ts.size(); ++i) h_timestamps[i] = ts[i];

    field.uploadVelocityWindow(h_vel.data(), h_vvel.data(), h_ztop, h_timestamps.data());

    mpaso_gpu::LaunchPathlineTracer(
        field,
        reinterpret_cast<const mpaso_vec3*>(seeds),
        seed_cell_id,
        seed_t_start,
        seed_max_steps,
        n_particles,
        n_steps,
        dt,
        reinterpret_cast<mpaso_vec3*>(traces_out),
        final_time_out,
        steps_taken_out,
        final_cell_out);

    field.release();
}

} // namespace mpaso_gpu_host
