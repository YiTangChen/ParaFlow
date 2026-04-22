#ifndef MPASO_GPU_TRACER_H
#define MPASO_GPU_TRACER_H

// Host-side entry point for GPU particle tracing. This header is
// deliberately plain C++ (no CUDA types) so callers like OSUFlow.C or
// Streamline.C can include it without nvcc.
//
// Typical use:
//
//   std::vector<VECTOR3> seeds = ...;
//   std::vector<int>     seed_cells = ...;   // phys_to_cell for each seed
//   std::vector<VECTOR3> traces(seeds.size() * (n_steps + 1));
//
//   mpaso_gpu_host::TraceParticles(
//       grid, pSolution, vSolution,
//       seeds.data(), seed_cells.data(), seeds.size(),
//       n_steps, dt, /*t_start=*/0.0, /*use_euler=*/false,
//       traces.data());
//
// Semantics: launches one GPU thread per seed, advances `n_steps` RK4
// (or Euler) steps of size `dt` seconds. Output layout is row-major,
// [n_particles * (n_steps + 1)] with index 0 per row = seed position.

#include "VectorMatrix.h"

class MPASOGrid;
class Solution;

namespace mpaso_gpu_host {

// Returns true if at least one CUDA device is visible at runtime.
// Safe to call without a CUDA context; returns false on any cudaError.
bool isAvailable();


void TraceParticles(MPASOGrid*       grid,
                    Solution*        pSolution,       // horizontal vertex velocity
                    Solution*        vSolution,       // vertical   vertex velocity
                    const VECTOR3*   seeds,
                    const int*       seed_cell_id,
                    int              n_particles,
                    int              n_steps,
                    double           dt,
                    double           t_start,
                    bool             use_euler,
                    VECTOR3*         traces_out,
                    const int*       seed_max_steps = nullptr,
                    const int*       seed_step_offset = nullptr,
                    int*             steps_taken_out = nullptr,
                    int*             final_cell_out = nullptr,
                    int              save_interval = 1,
                    int              max_saved_points = 0,
                    int*             saved_counts_out = nullptr,
                    float*           kernel_ms_out = nullptr);  // CUDA Event kernel time (ms), accumulated

// Pathline batch tracer. Like TraceParticles, but per-seed start times and
// temporal velocity blending. The uploaded window's timestamps come from
// Solution::getTimestamps() / MPASOGrid::getZTopTimestamps (real seconds).
//
// save_interval:     record one point every N steps (default 1 = every step)
// max_saved_points:  capacity per particle row in traces_out; 0 = auto-compute
// saved_counts_out:  actual saved points per particle (optional)
// traces_out layout: [P * max_saved_points_eff] row-major, row 0 = seed
void TracePathlineBatch(MPASOGrid*       grid,
                        Solution*        pSolution,
                        Solution*        vSolution,
                        const VECTOR3*   seeds,
                        const int*       seed_cell_id,
                        const double*    seed_t_start,
                        const int*       seed_max_steps,
                        int              n_particles,
                        int              n_steps,
                        double           dt,
                        VECTOR3*         traces_out,
                        double*          final_time_out,    // [P]
                        int*             steps_taken_out,   // [P]
                        int*             final_cell_out,    // [P], local cell id or -1
                        int              save_interval = 1,
                        int              max_saved_points = 0,
                        int*             saved_counts_out = nullptr,
                        float*           kernel_ms_out = nullptr);  // CUDA Event kernel time (ms)

} // namespace mpaso_gpu_host

#endif // MPASO_GPU_TRACER_H
