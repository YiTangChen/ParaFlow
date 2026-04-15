#ifndef MPASO_TRACER_KERNELS_H
#define MPASO_TRACER_KERNELS_H

#include "../MPASODeviceField.h"

// Host-visible entry point for the CUDA tracer. Launches one thread per
// particle and runs RK4 (or Euler when use_euler=true) for `n_steps` steps of
// size `dt` seconds, starting at `t_start`.
//
// Inputs / outputs live in host memory; this function handles the staging
// malloc/copy. A future optimized version can accept device pointers directly.
//
// Layout of traces_out: [n_particles * (n_steps + 1)] mpaso_vec3, row-major
// per particle, index 0 = seed position.
namespace mpaso_gpu {

void LaunchTracer(const MPASODeviceField& field,
                  const mpaso_vec3* h_seeds,
                  const int*        h_seed_cell_id,  // initial cell hint per particle
                  int               n_particles,
                  int               n_steps,
                  double            dt,
                  double            t_start,
                  bool              use_euler,
                  mpaso_vec3*       h_traces_out,
                  const int*        h_seed_max_steps = nullptr,
                  const int*        h_seed_step_offset = nullptr,
                  int*              h_steps_taken_out = nullptr,
                  int*              h_final_cell_out = nullptr,
                  int               save_interval = 1,
                  int               max_saved_points = 0,
                  int*              h_saved_counts_out = nullptr);

// Pathline launcher. Velocity is blended between two timesteps of the uploaded
// window based on each particle's current time. Each particle starts at its
// own t (from h_seed_t_start) and advances by dt per RK4 step.
// Outputs: per-particle full trajectory, final time, and steps actually taken
// before either completing n_steps or failing a sample.
void LaunchPathlineTracer(const MPASODeviceField& field,
                          const mpaso_vec3* h_seeds,
                          const int*        h_seed_cell_id,
                          const double*     h_seed_t_start,
                          const int*        h_seed_max_steps,
                          int               n_particles,
                          int               n_steps,
                          double            dt,
                          mpaso_vec3*       h_traces_out,        // [P * (n_steps + 1)]
                          double*           h_final_time_out,    // [P]
                          int*              h_steps_taken_out,   // [P]
                          int*              h_final_cell_out);   // [P], local cell id or -1
} // namespace mpaso_gpu

#endif // MPASO_TRACER_KERNELS_H
