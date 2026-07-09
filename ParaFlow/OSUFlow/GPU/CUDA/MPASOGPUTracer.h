#ifndef MPASO_GPU_TRACER_H
#define MPASO_GPU_TRACER_H

// Host-side entry point for GPU particle tracing. This header is
// deliberately plain C++ (no CUDA types) so callers like OSUFlow.C or
// Streamline.C can include it without nvcc.
//
// Semantics: launches one GPU thread per seed, advances `n_steps` RK4
// (or Euler) steps of size `dt` seconds. Output layout is row-major,
// [n_particles * (n_steps + 1)] with index 0 per row = seed position.

#include "VectorMatrix.h"

class MPASOGrid;
class Solution;

namespace mpaso_gpu_host {

struct GPUBlockContext;

struct GPUTimingBreakdown {
    float host_prepare_ms     = 0.0f;  // host metadata packing / Solution flattening
    float upload_topology_ms  = 0.0f;  // block-local static mesh upload
    float upload_velocity_ms  = 0.0f;  // velocity / zTop / timestamp window upload
    float alloc_ms            = 0.0f;  // per-launch CUDA allocations
    float upload_particles_ms = 0.0f;  // particle arrays host-to-device copies
    float kernel_ms           = 0.0f;  // CUDA Event kernel time
    float download_results_ms = 0.0f;  // device-to-host output copies
    float free_ms             = 0.0f;  // per-launch CUDA frees
    float field_release_ms    = 0.0f;  // release uploaded topology/velocity buffers
};

// Returns true if at least one CUDA device is visible at runtime.
// Safe to call without a CUDA context; returns false on any cudaError.
bool isAvailable();

GPUBlockContext* CreateGPUBlockContext(MPASOGrid* grid,
                                       GPUTimingBreakdown* timing_out = nullptr);

void UploadGPUVelocityWindow(GPUBlockContext* context,
                             MPASOGrid*       grid,
                             Solution*        pSolution,
                             Solution*        vSolution,
                             bool             use_real_timestamps,
                             GPUTimingBreakdown* timing_out = nullptr);

void TraceParticlesOnGPUContext(GPUBlockContext* context,
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
                                float*           kernel_ms_out = nullptr,
                                GPUTimingBreakdown* timing_out = nullptr);

void TracePathlineBatchOnGPUContext(GPUBlockContext* context,
                                    const VECTOR3*   seeds,
                                    const int*       seed_cell_id,
                                    const double*    seed_t_start,
                                    const int*       seed_max_steps,
                                    int              n_particles,
                                    int              n_steps,
                                    double           dt,
                                    VECTOR3*         traces_out,
                                    double*          final_time_out,
                                    int*             steps_taken_out,
                                    int*             final_cell_out,
                                    int              save_interval = 1,
                                    int              max_saved_points = 0,
                                    int*             saved_counts_out = nullptr,
                                    float*           kernel_ms_out = nullptr,
                                    GPUTimingBreakdown* timing_out = nullptr);

void DestroyGPUBlockContext(GPUBlockContext* context,
                            GPUTimingBreakdown* timing_out = nullptr);

} // namespace mpaso_gpu_host

#endif // MPASO_GPU_TRACER_H
