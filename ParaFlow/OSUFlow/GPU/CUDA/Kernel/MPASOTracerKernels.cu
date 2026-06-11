#include "MPASOTracerKernels.h"
#include "MPASOVelocity.cuh"

#include <cuda_runtime.h>
#include <chrono>
#include <cstdio>
#include <cstdlib>

// =============================================================================
// RK4 particle tracer for MPAS-O, ported (trimmed) from MOPS
// src/GPU/CUDA/Kernel/MPASOVisualizerKernels.cu :: KernelStreamLine.
// =============================================================================

namespace {

constexpr int kMaxCellNeighbors = 21;  // matches MOPS; largest n_max_edges + 1

// ---- device vec3 helpers --------------------------------------------------

__host__ __device__ inline mpaso_vec3 v3_add(mpaso_vec3 a, mpaso_vec3 b) {
    return { a.x + b.x, a.y + b.y, a.z + b.z };
}
__host__ __device__ inline mpaso_vec3 v3_scale(mpaso_vec3 a, double s) {
    return { a.x * s, a.y * s, a.z * s };
}
__host__ __device__ inline double v3_len(mpaso_vec3 a) {
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}
__device__ inline bool v3_same(mpaso_vec3 a, mpaso_vec3 b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
__device__ inline void trace_save_point(mpaso_vec3* traces_out,
                                        int row,
                                        int max_saved_points,
                                        int* saved_count,
                                        mpaso_vec3 p)
{
    if (*saved_count < max_saved_points)
        traces_out[row + (*saved_count)++] = p;
}
__device__ inline void trace_save_final(mpaso_vec3* traces_out,
                                        int row,
                                        int max_saved_points,
                                        int* saved_count,
                                        mpaso_vec3 p)
{
    if (*saved_count <= 0 || !v3_same(traces_out[row + *saved_count - 1], p))
        trace_save_point(traces_out, row, max_saved_points, saved_count, p);
}
__device__ inline void trace_pad_points(mpaso_vec3* traces_out,
                                        int row,
                                        int max_saved_points,
                                        int saved_count,
                                        mpaso_vec3 p)
{
    for (int k = saved_count; k < max_saved_points; ++k)
        traces_out[row + k] = p;
}

// ---- velocity sampling -----------------------------------------------------
// Runs phys_to_cell_dev + sample_velocity_dev (ports of MPASOGrid::phys_to_cell
// and CVectorField::at_phys / MPASOGrid::interpolate). `cell_id` is updated in
// place so subsequent RK4 sub-steps can use it as the walk hint.
//
// Returns horizontal velocity; vertical scalar written to *out_vert. On
// failure (outside domain, bad cell, etc.) returns zero and sets *cell_id = -1.
__device__ mpaso_vec3 sample_velocity(mpaso_vec3 p,
                                      int* cell_id,
                                      double /*ts_blend*/,
                                      const MPASODeviceField& f,
                                      double* out_vert)
{
    *out_vert = 0.0;
    if (*cell_id < 0) return { 0.0, 0.0, 0.0 };

    int vert_ids[mpaso_gpu_dev::DEV_MAX_EDGES];
    double omegas[mpaso_gpu_dev::DEV_MAX_EDGES + 1];
    int new_cell, v_level, nV;

    bool ok = mpaso_gpu_dev::phys_to_cell_dev(
        f, p, *cell_id, /*ts_idx=*/0,
        &new_cell, &v_level, &nV, omegas, vert_ids);
    if (!ok) { *cell_id = -1; return { 0.0, 0.0, 0.0 }; }
    *cell_id = new_cell;

    mpaso_vec3 h_out;
    double     v_out;
    ok = mpaso_gpu_dev::sample_velocity_dev(
        f, /*ts_idx=*/0, vert_ids, nV, v_level, omegas, &h_out, &v_out);
    if (!ok) { *cell_id = -1; return { 0.0, 0.0, 0.0 }; }

    *out_vert = v_out;
    return h_out;
}

// ---- one RK4 step ---------------------------------------------------------
// Mirrors vtCFieldLine::MPASO_rk4 (FieldLine.C:247-316):
//   1. sample k1 at pt0, walk to pt1 = geodesic(pt0, k1, dt/2)
//   2. sample k2 at pt1, walk to pt2 = geodesic(pt0, k2, dt/2)
//   3. sample k3 at pt2, walk to pt3 = geodesic(pt0, k3, dt)
//   4. sample k4 at pt3
//   5. v_avg = (k1 + 2 k2 + 2 k3 + k4) / 6
//   6. pt_final = geodesic(pt0, v_avg, dt)
// Returns false (cell_id = -1) if any sample or step leaves the domain.
__device__ bool rk4_step(mpaso_vec3& p,
                         int* cell_id,
                         double dt,
                         double t_rel,
                         const MPASODeviceField& f)
{
    const double r0 = v3_len(p);
    if (r0 <= 0.0) { *cell_id = -1; return false; }
    const mpaso_vec3 pt0 = p;

    double k1_v, k2_v, k3_v, k4_v;
    mpaso_vec3 k1_h, k2_h, k3_h, k4_h;

    // Stage 1
    k1_h = sample_velocity(pt0, cell_id, t_rel, f, &k1_v);
    if (*cell_id < 0) return false;
    // cell0 = the cell freshly resolved at pt0 by stage 1 (matches CPU
    // cell0 = ci_tmp.inCell); capturing it before stage 1 feeds stages 2-4 a
    // one-step-stale hint and makes the single-pass walk fail early.
    const int cell0 = *cell_id;

    // Stage 2
    mpaso_vec3 pt1; double r1;
    if (!mpaso_gpu_dev::geodesic_step_dev(pt0, r0, k1_h, k1_v, 0.5 * dt, &pt1, &r1)) {
        *cell_id = -1; return false;
    }
    int cell_walk = cell0;
    k2_h = sample_velocity(pt1, &cell_walk, t_rel, f, &k2_v);
    if (cell_walk < 0) { *cell_id = -1; return false; }

    // Stage 3
    mpaso_vec3 pt2; double r2;
    if (!mpaso_gpu_dev::geodesic_step_dev(pt0, r0, k2_h, k2_v, 0.5 * dt, &pt2, &r2)) {
        *cell_id = -1; return false;
    }
    cell_walk = cell0;
    k3_h = sample_velocity(pt2, &cell_walk, t_rel, f, &k3_v);
    if (cell_walk < 0) { *cell_id = -1; return false; }

    // Stage 4
    mpaso_vec3 pt3; double r3;
    if (!mpaso_gpu_dev::geodesic_step_dev(pt0, r0, k3_h, k3_v, dt, &pt3, &r3)) {
        *cell_id = -1; return false;
    }
    cell_walk = cell0;
    k4_h = sample_velocity(pt3, &cell_walk, t_rel, f, &k4_v);
    if (cell_walk < 0) { *cell_id = -1; return false; }

    // RK4 weighted average
    mpaso_vec3 v_avg_h = {
        (k1_h.x + 2.0*k2_h.x + 2.0*k3_h.x + k4_h.x) / 6.0,
        (k1_h.y + 2.0*k2_h.y + 2.0*k3_h.y + k4_h.y) / 6.0,
        (k1_h.z + 2.0*k2_h.z + 2.0*k3_h.z + k4_h.z) / 6.0
    };
    double v_avg_v = (k1_v + 2.0*k2_v + 2.0*k3_v + k4_v) / 6.0;

    mpaso_vec3 pt_final; double r_final;
    if (!mpaso_gpu_dev::geodesic_step_dev(pt0, r0, v_avg_h, v_avg_v, dt, &pt_final, &r_final)) {
        *cell_id = -1; return false;
    }
    p = pt_final;
    // cell_id carries the stage-1 walk result, matching CPU's ci.fromCell = ci.inCell
    return true;
}

// ---- one Euler step -------------------------------------------------------
// Mirrors vtCFieldLine::MPASO_euler (FieldLine.C:201-239).
__device__ bool euler_step(mpaso_vec3& p,
                           int* cell_id,
                           double dt,
                           double t_rel,
                           const MPASODeviceField& f)
{
    const double r0 = v3_len(p);
    if (r0 <= 0.0) { *cell_id = -1; return false; }

    double v_vert;
    mpaso_vec3 h = sample_velocity(p, cell_id, t_rel, f, &v_vert);
    if (*cell_id < 0) return false;

    mpaso_vec3 pt_dst; double r_dst;
    if (!mpaso_gpu_dev::geodesic_step_dev(p, r0, h, v_vert, dt, &pt_dst, &r_dst)) {
        *cell_id = -1; return false;
    }
    p = pt_dst;
    return true;
}

// ---- main kernel -----------------------------------------------------------
__global__ void kernelTraceRK4(MPASODeviceField f,
                               const mpaso_vec3* seeds,
                               const int*        seed_cell_id,
                               const int*        seed_max_steps,
                               const int*        seed_step_offset,
                               int               n_particles,
                               int               n_steps,
                               double            dt,
                               double            t_start,
                               bool              use_euler,
                               mpaso_vec3*       traces_out,
                               int*              steps_taken_out,
                               int*              final_cell_out,
                               int               save_interval,
                               int               max_saved_points,
                               int*              saved_counts_out)
{
    const int gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= n_particles) return;

    const int row   = gid * max_saved_points;
    mpaso_vec3 p    = seeds[gid];
    int cell_id    = seed_cell_id[gid];
    int max_steps_for_seed = seed_max_steps ? seed_max_steps[gid] : n_steps;
    int step_offset = seed_step_offset ? seed_step_offset[gid] : 0;
    int saved_count = 0;
    trace_save_point(traces_out, row, max_saved_points, &saved_count, p);

    int steps_taken = 0;
    if (max_steps_for_seed <= 0) {
        trace_pad_points(traces_out, row, max_saved_points, saved_count, p);
        if (steps_taken_out) steps_taken_out[gid] = 0;
        if (final_cell_out)  final_cell_out[gid]  = cell_id;
        if (saved_counts_out) saved_counts_out[gid] = saved_count;
        return;
    }

    double t = t_start;
    for (int s = 0; s < n_steps && s < max_steps_for_seed; ++s) {
        // TODO: compute t_rel in [0,1] against f.d_timestamps window for
        // time-varying velocity blending. Streamline ignores this.
        double t_rel = 0.0;

        int prev_cell = cell_id;
        bool ok = use_euler ? euler_step(p, &cell_id, dt, t_rel, f)
                            : rk4_step(p, &cell_id, dt, t_rel, f);
        if (!ok || cell_id < 0) {
            trace_save_final(traces_out, row, max_saved_points, &saved_count, p);
            trace_pad_points(traces_out, row, max_saved_points, saved_count, p);
            if (steps_taken_out) steps_taken_out[gid] = steps_taken;
            if (final_cell_out)  final_cell_out[gid]  = prev_cell;
            if (saved_counts_out) saved_counts_out[gid] = saved_count;
            return;
        }
        t += dt;
        ++steps_taken;
        if (((step_offset + steps_taken) % save_interval) == 0)
            trace_save_point(traces_out, row, max_saved_points, &saved_count, p);
    }
    trace_save_final(traces_out, row, max_saved_points, &saved_count, p);
    trace_pad_points(traces_out, row, max_saved_points, saved_count, p);
    if (steps_taken_out) steps_taken_out[gid] = steps_taken;
    if (final_cell_out)  final_cell_out[gid]  = cell_id;
    if (saved_counts_out) saved_counts_out[gid] = saved_count;
}

inline void cudaCheck(cudaError_t err, const char* what) {
    if (err != cudaSuccess) {
        std::fprintf(stderr, "[mpaso_gpu] CUDA error in %s: %s\n",
                     what, cudaGetErrorString(err));
        std::abort();
    }
}

using HostClock = std::chrono::steady_clock;

inline HostClock::time_point host_timer_now()
{
    return HostClock::now();
}

inline float host_timer_ms_since(HostClock::time_point start)
{
    return std::chrono::duration<float, std::milli>(HostClock::now() - start).count();
}

} // namespace

namespace mpaso_gpu {

void LaunchTracer(const MPASODeviceField& field,
                  const mpaso_vec3* h_seeds,
                  const int*        h_seed_cell_id,
                  int               n_particles,
                  int               n_steps,
                  double            dt,
                  double            t_start,
                  bool              use_euler,
                  mpaso_vec3*       h_traces_out,
                  const int*        h_seed_max_steps,
                  const int*        h_seed_step_offset,
                  int*              h_steps_taken_out,
                  int*              h_final_cell_out,
                  int               save_interval,
                  int               max_saved_points,
                  int*              h_saved_counts_out,
                  float*            kernel_ms_out,
                  float*            alloc_ms_out,
                  float*            upload_particles_ms_out,
                  float*            download_results_ms_out,
                  float*            free_ms_out)
{
    if (n_particles <= 0 || n_steps <= 0) return;

    mpaso_vec3* d_seeds        = nullptr;
    int*        d_seed_cell_id = nullptr;
    int*        d_max_steps    = nullptr;
    int*        d_step_offset  = nullptr;
    mpaso_vec3* d_traces       = nullptr;
    int*        d_steps_taken  = nullptr;
    int*        d_final_cell   = nullptr;
    int*        d_saved_counts = nullptr;

    const size_t seed_bytes  = static_cast<size_t>(n_particles) * sizeof(mpaso_vec3);
    const size_t cellid_bytes= static_cast<size_t>(n_particles) * sizeof(int);
    int save_interval_eff = save_interval > 0 ? save_interval : 1;
    int max_saved_eff = max_saved_points;
    if (max_saved_eff <= 0) {
        max_saved_eff = (save_interval_eff == 1)
            ? n_steps + 1
            : n_steps / save_interval_eff + 2;
    }
    const size_t trace_bytes = static_cast<size_t>(n_particles)
                             * static_cast<size_t>(max_saved_eff)
                             * sizeof(mpaso_vec3);

    auto t_alloc = host_timer_now();
    cudaCheck(cudaMalloc(&d_seeds,        seed_bytes),   "alloc seeds");
    cudaCheck(cudaMalloc(&d_seed_cell_id, cellid_bytes), "alloc seed_cell_id");
    if (h_seed_max_steps)
        cudaCheck(cudaMalloc(&d_max_steps, cellid_bytes), "alloc seed_max_steps");
    if (h_seed_step_offset)
        cudaCheck(cudaMalloc(&d_step_offset, cellid_bytes), "alloc seed_step_offset");
    cudaCheck(cudaMalloc(&d_traces,       trace_bytes),  "alloc traces");
    if (h_steps_taken_out)
        cudaCheck(cudaMalloc(&d_steps_taken, cellid_bytes), "alloc steps_taken");
    if (h_final_cell_out)
        cudaCheck(cudaMalloc(&d_final_cell, cellid_bytes), "alloc final_cell");
    if (h_saved_counts_out)
        cudaCheck(cudaMalloc(&d_saved_counts, cellid_bytes), "alloc saved_counts");
    if (alloc_ms_out)
        *alloc_ms_out += host_timer_ms_since(t_alloc);

    auto t_upload_particles = host_timer_now();
    cudaCheck(cudaMemcpy(d_seeds,        h_seeds,        seed_bytes,   cudaMemcpyHostToDevice), "copy seeds");
    cudaCheck(cudaMemcpy(d_seed_cell_id, h_seed_cell_id, cellid_bytes, cudaMemcpyHostToDevice), "copy seed_cell_id");
    if (h_seed_max_steps)
        cudaCheck(cudaMemcpy(d_max_steps, h_seed_max_steps, cellid_bytes, cudaMemcpyHostToDevice), "copy seed_max_steps");
    if (h_seed_step_offset)
        cudaCheck(cudaMemcpy(d_step_offset, h_seed_step_offset, cellid_bytes, cudaMemcpyHostToDevice), "copy seed_step_offset");
    if (upload_particles_ms_out)
        *upload_particles_ms_out += host_timer_ms_since(t_upload_particles);

    const int block = 128;
    const int grid  = (n_particles + block - 1) / block;
    cudaEvent_t ev_start, ev_stop;
    cudaCheck(cudaEventCreate(&ev_start), "event create start");
    cudaCheck(cudaEventCreate(&ev_stop),  "event create stop");
    cudaCheck(cudaEventRecord(ev_start),  "event record start");
    kernelTraceRK4<<<grid, block>>>(field, d_seeds, d_seed_cell_id, d_max_steps, d_step_offset,
                                    n_particles, n_steps, dt, t_start,
                                    use_euler, d_traces,
                                    d_steps_taken, d_final_cell,
                                    save_interval_eff, max_saved_eff, d_saved_counts);
    cudaCheck(cudaGetLastError(),        "kernel launch");
    cudaCheck(cudaEventRecord(ev_stop),  "event record stop");
    cudaCheck(cudaEventSynchronize(ev_stop), "event sync");
    if (kernel_ms_out) {
        float ms = 0.0f;
        cudaCheck(cudaEventElapsedTime(&ms, ev_start, ev_stop), "event elapsed");
        *kernel_ms_out += ms;
    }
    cudaEventDestroy(ev_start);
    cudaEventDestroy(ev_stop);

    auto t_download_results = host_timer_now();
    cudaCheck(cudaMemcpy(h_traces_out, d_traces, trace_bytes, cudaMemcpyDeviceToHost), "copy traces");
    if (h_steps_taken_out)
        cudaCheck(cudaMemcpy(h_steps_taken_out, d_steps_taken, cellid_bytes, cudaMemcpyDeviceToHost), "copy steps_taken");
    if (h_final_cell_out)
        cudaCheck(cudaMemcpy(h_final_cell_out, d_final_cell, cellid_bytes, cudaMemcpyDeviceToHost), "copy final_cell");
    if (h_saved_counts_out)
        cudaCheck(cudaMemcpy(h_saved_counts_out, d_saved_counts, cellid_bytes, cudaMemcpyDeviceToHost), "copy saved_counts");
    if (download_results_ms_out)
        *download_results_ms_out += host_timer_ms_since(t_download_results);

    auto t_free = host_timer_now();
    cudaFree(d_seeds);
    cudaFree(d_seed_cell_id);
    cudaFree(d_max_steps);
    cudaFree(d_step_offset);
    cudaFree(d_traces);
    cudaFree(d_steps_taken);
    cudaFree(d_final_cell);
    cudaFree(d_saved_counts);
    if (free_ms_out)
        *free_ms_out += host_timer_ms_since(t_free);
}

} // namespace mpaso_gpu

// ============================================================================
// Pathline RK4 (time-varying velocity, temporal blend between two timesteps)
// ============================================================================
namespace {

// Sample velocity at particle time t_part. Resolves (lowT, highT, ratio) from
// f.d_timestamps and runs phys_to_cell_dev_blend + sample_velocity_dev_blend.
// `cell_id` acts as a walk hint; updated in place. Fails → sets *cell_id = -1.
__device__ mpaso_vec3 sample_velocity_t(mpaso_vec3 p,
                                        int* cell_id,
                                        double t_part,
                                        const MPASODeviceField& f,
                                        double* out_vert)
{
    *out_vert = 0.0;
    if (*cell_id < 0) return { 0.0, 0.0, 0.0 };

    int lowT, highT;
    double ratio;
    if (!mpaso_gpu_dev::resolve_time_window_dev(f, t_part, &lowT, &highT, &ratio)) {
        *cell_id = -1;
        return { 0.0, 0.0, 0.0 };
    }

    int vert_ids[mpaso_gpu_dev::DEV_MAX_EDGES];
    double omegas[mpaso_gpu_dev::DEV_MAX_EDGES + 1];
    int new_cell, v_level, nV;

    bool ok = mpaso_gpu_dev::phys_to_cell_dev_blend(
        f, p, *cell_id, lowT, highT, ratio,
        &new_cell, &v_level, &nV, omegas, vert_ids);
    if (!ok) { *cell_id = -1; return { 0.0, 0.0, 0.0 }; }
    *cell_id = new_cell;

    mpaso_vec3 h_out;
    double     v_out;
    ok = mpaso_gpu_dev::sample_velocity_dev_blend(
        f, lowT, highT, ratio, vert_ids, nV, v_level, omegas, &h_out, &v_out);
    if (!ok) { *cell_id = -1; return { 0.0, 0.0, 0.0 }; }

    *out_vert = v_out;
    return h_out;
}

// Pathline RK4: unlike streamline, each sub-stage evaluates velocity at a
// different time (t, t+dt/2, t+dt), matching vtCFieldLine::MPASO_rk4's
// UNSTEADY branch (FieldLine.C:273,292).
__device__ bool rk4_step_pathline(mpaso_vec3& p,
                                  int* cell_id,
                                  double dt,
                                  double t_cur,
                                  const MPASODeviceField& f)
{
    const double r0 = v3_len(p);
    if (r0 <= 0.0) { *cell_id = -1; return false; }
    const mpaso_vec3 pt0 = p;

    double k1_v, k2_v, k3_v, k4_v;
    mpaso_vec3 k1_h, k2_h, k3_h, k4_h;

    // Stage 1: (pt0, t_cur)
    k1_h = sample_velocity_t(pt0, cell_id, t_cur, f, &k1_v);
    if (*cell_id < 0) return false;
    // cell0 = the cell FRESHLY resolved at pt0 by stage 1 (matches the CPU's
    // cell0 = ci_tmp.inCell). Capturing it before stage 1 would feed stages 2-4
    // a one-step-stale hint, making the single-pass walk fail earlier than the
    // CPU and the pathline terminate prematurely.
    const int cell0 = *cell_id;

    // Stage 2: (pt1, t_cur + dt/2)
    mpaso_vec3 pt1; double r1;
    if (!mpaso_gpu_dev::geodesic_step_dev(pt0, r0, k1_h, k1_v, 0.5 * dt, &pt1, &r1)) {
        *cell_id = -1; return false;
    }
    int cell_walk = cell0;
    const double t_mid = t_cur + 0.5 * dt;
    k2_h = sample_velocity_t(pt1, &cell_walk, t_mid, f, &k2_v);
    // On a sub-stage sample failure the CPU advances to the sub-stage point
    // (FieldLine.C:276: ci.phyCoord = pt1) and stores it, so the GPU must too.
    if (cell_walk < 0) { p = pt1; *cell_id = -1; return false; }

    // Stage 3: (pt2, t_cur + dt/2)
    mpaso_vec3 pt2; double r2;
    if (!mpaso_gpu_dev::geodesic_step_dev(pt0, r0, k2_h, k2_v, 0.5 * dt, &pt2, &r2)) {
        *cell_id = -1; return false;
    }
    cell_walk = cell0;
    k3_h = sample_velocity_t(pt2, &cell_walk, t_mid, f, &k3_v);
    if (cell_walk < 0) { p = pt2; *cell_id = -1; return false; }  // FieldLine.C:285

    // Stage 4: (pt3, t_cur + dt)
    mpaso_vec3 pt3; double r3;
    if (!mpaso_gpu_dev::geodesic_step_dev(pt0, r0, k3_h, k3_v, dt, &pt3, &r3)) {
        *cell_id = -1; return false;
    }
    cell_walk = cell0;
    const double t_end = t_cur + dt;
    k4_h = sample_velocity_t(pt3, &cell_walk, t_end, f, &k4_v);
    if (cell_walk < 0) { p = pt3; *cell_id = -1; return false; }  // FieldLine.C:295

    mpaso_vec3 v_avg_h = {
        (k1_h.x + 2.0*k2_h.x + 2.0*k3_h.x + k4_h.x) / 6.0,
        (k1_h.y + 2.0*k2_h.y + 2.0*k3_h.y + k4_h.y) / 6.0,
        (k1_h.z + 2.0*k2_h.z + 2.0*k3_h.z + k4_h.z) / 6.0
    };
    double v_avg_v = (k1_v + 2.0*k2_v + 2.0*k3_v + k4_v) / 6.0;

    mpaso_vec3 pt_final; double r_final;
    if (!mpaso_gpu_dev::geodesic_step_dev(pt0, r0, v_avg_h, v_avg_v, dt, &pt_final, &r_final)) {
        *cell_id = -1; return false;
    }
    p = pt_final;
    return true;
}

// Main pathline kernel.
// Per-particle inputs: seed position, cell hint, start time.
// Per-particle outputs (row = gid * max_saved_points):
//   traces_out   — downsampled trajectory; one point every save_interval steps
//                  plus the final position; remainder padded with final point.
//   final_time_out[P]    — simulation time after the last successful step
//   steps_taken_out[P]   — RK4 steps completed before termination
//   final_cell_out[P]    — local cell id at termination (-1 = outside domain)
//   saved_counts_out[P]  — actual points written into traces_out row (optional)
// Termination: (a) n_steps reached, (b) time window exhausted, (c) sample fails.
__global__ void kernelTracePathline(MPASODeviceField f,
                                    const mpaso_vec3* seeds,
                                    const int*        seed_cell_id,
                                    const double*     seed_t_start,
                                    const int*        seed_max_steps,
                                    int               n_particles,
                                    int               n_steps,
                                    double            dt,
                                    mpaso_vec3*       traces_out,
                                    double*           final_time_out,
                                    int*              steps_taken_out,
                                    int*              final_cell_out,
                                    int               save_interval,
                                    int               max_saved_points,
                                    int*              saved_counts_out)
{
    const int gid = blockIdx.x * blockDim.x + threadIdx.x;
    if (gid >= n_particles) return;

    const int row   = gid * max_saved_points;
    mpaso_vec3 p    = seeds[gid];
    int cell_id     = seed_cell_id[gid];
    double t        = seed_t_start[gid];
    int max_steps_for_seed = seed_max_steps[gid];
    int saved_count = 0;
    const double kStationaryCutoff = 1.0e-5;  // vtCFieldLine::m_fStationaryCutoff

    // CPU seed-time check (TimeVaryingFieldLine.C:239-252): emit no trajectory if
    // the seed is out of the domain or stationary.
    {
        double seed_vv; int seed_cell = cell_id;
        mpaso_vec3 seed_vh = sample_velocity_t(p, &seed_cell, t, f, &seed_vv);
        if (seed_cell < 0 ||
            (fabs(seed_vh.x) < kStationaryCutoff &&
             fabs(seed_vh.y) < kStationaryCutoff &&
             fabs(seed_vh.z) < kStationaryCutoff)) {
            final_time_out[gid]  = t;
            steps_taken_out[gid] = 0;
            final_cell_out[gid]  = (seed_cell < 0) ? -1 : cell_id;
            if (saved_counts_out) saved_counts_out[gid] = 0;
            return;
        }
    }
    trace_save_point(traces_out, row, max_saved_points, &saved_count, p);

    int steps_taken = 0;
    if (max_steps_for_seed <= 0) {
        trace_pad_points(traces_out, row, max_saved_points, saved_count, p);
        final_time_out[gid]  = t;
        steps_taken_out[gid] = 0;
        final_cell_out[gid]  = cell_id;
        if (saved_counts_out) saved_counts_out[gid] = saved_count;
        return;
    }

    for (int s = 0; s < n_steps && s < max_steps_for_seed; ++s) {
        if (f.n_timesteps_loaded > 1 &&
            t + dt > f.d_timestamps[f.n_timesteps_loaded - 1] + 1.0e-9)
            break;  // out of loaded time window; driver reinjects into next window
        int prev_cell = cell_id;
        bool ok = rk4_step_pathline(p, &cell_id, dt, t, f);
        if (!ok || cell_id < 0) {
            trace_save_final(traces_out, row, max_saved_points, &saved_count, p);
            trace_pad_points(traces_out, row, max_saved_points, saved_count, p);
            final_time_out[gid]  = t;
            steps_taken_out[gid] = steps_taken;
            final_cell_out[gid]  = prev_cell;
            if (saved_counts_out) saved_counts_out[gid] = saved_count;
            return;
        }
        t += dt;
        ++steps_taken;
        if ((steps_taken % save_interval) == 0)
            trace_save_point(traces_out, row, max_saved_points, &saved_count, p);

        // CPU post-step check (TimeVaryingFieldLine.C:350-388): resample at the
        // new position/time, refresh the cell hint, and stop on out-of-bound or
        // a stationary (critical) point so the GPU terminates where the CPU does.
        double post_vv;
        int post_cell = cell_id;
        mpaso_vec3 post_vh = sample_velocity_t(p, &post_cell, t, f, &post_vv);
        if (post_cell < 0)
            break;
        cell_id = post_cell;
        if (fabs(post_vh.x) < kStationaryCutoff &&
            fabs(post_vh.y) < kStationaryCutoff &&
            fabs(post_vh.z) < kStationaryCutoff)
            break;
    }
    trace_save_final(traces_out, row, max_saved_points, &saved_count, p);
    trace_pad_points(traces_out, row, max_saved_points, saved_count, p);
    final_time_out[gid]  = t;
    steps_taken_out[gid] = steps_taken;
    final_cell_out[gid]  = cell_id;
    if (saved_counts_out) saved_counts_out[gid] = saved_count;
}

} // namespace

namespace mpaso_gpu {

void LaunchPathlineTracer(const MPASODeviceField& field,
                          const mpaso_vec3* h_seeds,
                          const int*        h_seed_cell_id,
                          const double*     h_seed_t_start,
                          const int*        h_seed_max_steps,
                          int               n_particles,
                          int               n_steps,
                          double            dt,
                          mpaso_vec3*       h_traces_out,
                          double*           h_final_time_out,
                          int*              h_steps_taken_out,
                          int*              h_final_cell_out,
                          int               save_interval,
                          int               max_saved_points,
                          int*              h_saved_counts_out,
                          float*            kernel_ms_out,
                          float*            alloc_ms_out,
                          float*            upload_particles_ms_out,
                          float*            download_results_ms_out,
                          float*            free_ms_out)
{
    if (n_particles <= 0 || n_steps <= 0) return;

    const int save_interval_eff = save_interval > 0 ? save_interval : 1;
    const int max_saved_eff = max_saved_points > 0
        ? max_saved_points
        : n_steps / save_interval_eff + 2;

    mpaso_vec3* d_seeds        = nullptr;
    int*        d_cell         = nullptr;
    int*        d_max_steps    = nullptr;
    double*     d_t_start      = nullptr;
    mpaso_vec3* d_traces       = nullptr;
    double*     d_final_time   = nullptr;
    int*        d_steps_taken  = nullptr;
    int*        d_final_cell   = nullptr;
    int*        d_saved_counts = nullptr;

    const size_t sb  = (size_t)n_particles * sizeof(mpaso_vec3);
    const size_t cb  = (size_t)n_particles * sizeof(int);
    const size_t db  = (size_t)n_particles * sizeof(double);
    const size_t tb  = (size_t)n_particles * (size_t)max_saved_eff * sizeof(mpaso_vec3);

    auto t_alloc = host_timer_now();
    cudaCheck(cudaMalloc(&d_seeds,       sb), "pl alloc seeds");
    cudaCheck(cudaMalloc(&d_cell,        cb), "pl alloc cell");
    cudaCheck(cudaMalloc(&d_max_steps,   cb), "pl alloc max_steps");
    cudaCheck(cudaMalloc(&d_t_start,     db), "pl alloc t_start");
    cudaCheck(cudaMalloc(&d_traces,      tb), "pl alloc traces");
    cudaCheck(cudaMalloc(&d_final_time,  db), "pl alloc final_time");
    cudaCheck(cudaMalloc(&d_steps_taken, cb), "pl alloc steps_taken");
    cudaCheck(cudaMalloc(&d_final_cell,  cb), "pl alloc final_cell");
    if (h_saved_counts_out)
        cudaCheck(cudaMalloc(&d_saved_counts, cb), "pl alloc saved_counts");
    if (alloc_ms_out)
        *alloc_ms_out += host_timer_ms_since(t_alloc);

    auto t_upload_particles = host_timer_now();
    cudaCheck(cudaMemcpy(d_seeds,     h_seeds,           sb, cudaMemcpyHostToDevice), "pl cpy seeds");
    cudaCheck(cudaMemcpy(d_cell,      h_seed_cell_id,    cb, cudaMemcpyHostToDevice), "pl cpy cell");
    cudaCheck(cudaMemcpy(d_max_steps, h_seed_max_steps,  cb, cudaMemcpyHostToDevice), "pl cpy max_steps");
    cudaCheck(cudaMemcpy(d_t_start,   h_seed_t_start,    db, cudaMemcpyHostToDevice), "pl cpy t_start");
    if (upload_particles_ms_out)
        *upload_particles_ms_out += host_timer_ms_since(t_upload_particles);

    const int block = 128;
    const int grid  = (n_particles + block - 1) / block;
    cudaEvent_t ev_start, ev_stop;
    cudaCheck(cudaEventCreate(&ev_start), "pl event create start");
    cudaCheck(cudaEventCreate(&ev_stop),  "pl event create stop");
    cudaCheck(cudaEventRecord(ev_start),  "pl event record start");
    kernelTracePathline<<<grid, block>>>(field, d_seeds, d_cell, d_t_start, d_max_steps,
                                         n_particles, n_steps, dt,
                                         d_traces, d_final_time, d_steps_taken, d_final_cell,
                                         save_interval_eff, max_saved_eff, d_saved_counts);
    cudaCheck(cudaGetLastError(),            "pl kernel launch");
    cudaCheck(cudaEventRecord(ev_stop),      "pl event record stop");
    cudaCheck(cudaEventSynchronize(ev_stop), "pl event sync");
    if (kernel_ms_out) {
        float ms = 0.0f;
        cudaCheck(cudaEventElapsedTime(&ms, ev_start, ev_stop), "pl event elapsed");
        *kernel_ms_out += ms;
    }
    cudaEventDestroy(ev_start);
    cudaEventDestroy(ev_stop);

    auto t_download_results = host_timer_now();
    cudaCheck(cudaMemcpy(h_traces_out,      d_traces,      tb, cudaMemcpyDeviceToHost), "pl cpy traces");
    cudaCheck(cudaMemcpy(h_final_time_out,  d_final_time,  db, cudaMemcpyDeviceToHost), "pl cpy final_time");
    cudaCheck(cudaMemcpy(h_steps_taken_out, d_steps_taken, cb, cudaMemcpyDeviceToHost), "pl cpy steps_taken");
    cudaCheck(cudaMemcpy(h_final_cell_out,  d_final_cell,  cb, cudaMemcpyDeviceToHost), "pl cpy final_cell");
    if (h_saved_counts_out)
        cudaCheck(cudaMemcpy(h_saved_counts_out, d_saved_counts, cb, cudaMemcpyDeviceToHost), "pl cpy saved_counts");
    if (download_results_ms_out)
        *download_results_ms_out += host_timer_ms_since(t_download_results);

    auto t_free = host_timer_now();
    cudaFree(d_seeds); cudaFree(d_cell); cudaFree(d_max_steps); cudaFree(d_t_start);
    cudaFree(d_traces); cudaFree(d_final_time); cudaFree(d_steps_taken); cudaFree(d_final_cell);
    cudaFree(d_saved_counts);
    if (free_ms_out)
        *free_ms_out += host_timer_ms_since(t_free);
}

} // namespace mpaso_gpu
