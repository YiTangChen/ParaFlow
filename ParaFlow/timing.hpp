#ifndef _TIMING_HPP
#define _TIMING_HPP

#include <mpi.h>
#include <fstream>
#include <string>
#include <cstdio>

// Per-block timing accumulators (stored in each Block instance). All times are in seconds.
struct BlockTiming {
    // Initialization phase
    double t_block_load         = 0.0;  // set_data total (sum of the 4 sub-phases below + residual)
    double t_load_grid_build    = 0.0;  // MPASOGrid build + cell/vertex index mappings
    double t_load_read          = 0.0;  // NetCDF field reads (velocity/zTop/vertVel), all timesteps
    double t_load_convert       = 0.0;  // cell->vertex conversion, all timesteps
    double t_load_seed_filter   = 0.0;  // inBlock() seed filtering

    // Tracing phase (accumulated over all DIY iexchange loop iterations)
    double t_trace_dequeue       = 0.0;  // dequeue particles already received by DIY
    double t_trace_local_wall    = 0.0;  // local block tracing wall time, excluding dequeue
    double t_trace_prepare       = 0.0;  // prepare per-particle arrays, hints, and launch batches
    double t_trace_transfer      = 0.0;  // GPU host<->device transfer+setup (Solution flatten + H2D/D2H + alloc/free); 0 on CPU
    double t_trace_integrate_cpu = 0.0;  // pure CPU OSUFlow integration call time
    double t_trace_integrate_gpu = 0.0;  // pure CUDA integration kernel time
    double t_trace_postprocess   = 0.0;  // build segments, update cells/steps, classify exits
    double t_trace_enqueue       = 0.0;  // enqueue outgoing particles to DIY

    // Output phase
    double t_output_write       = 0.0;  // write_trajectory

    // Work counters (useful for normalizing time and detecting load imbalance)
    int    n_seeds_initial      = 0;    // seeds assigned to this block at startup
    long   n_steps_total        = 0;    // total RK4 integration steps completed
    int    n_particles_received = 0;    // particles received from neighbor blocks

    // GPU transfer drill-down: children of t_trace_transfer. They sum to it
    // (single axis, safe to add), and are all zero on CPU-only runs.
    double t_transfer_flatten         = 0.0;  // host Solution flatten (∝ local_cells × timesteps)
    double t_transfer_upload_topology = 0.0;  // static MPAS-O topology H2D (one-time)
    double t_transfer_upload_velocity = 0.0;  // velocity/zTop window H2D (per time-window)
    double t_transfer_misc            = 0.0;  // per-launch alloc + particle H2D + result D2H + free/release

    // Cell counts (set at init, used to derive Block-owned index array sizes)
    int    n_local_cells        = 0;    // number of cells owned by this block
    int    n_global_cells       = 0;    // total cells in the global mesh

    // Memory measurements (computed once at init, analytical formula)
    size_t mem_grid_bytes       = 0;    // MPASOGrid topology + coordinates + Block index arrays
    size_t mem_solution_bytes   = 0;    // Solution: velocity field data (scales with timesteps)

    // Memory measurements (OS-reported, from /proc/self/status)
    long   mem_vmrss_before_kb  = 0;    // VmRSS before set_data (kB)
    long   mem_vmrss_after_kb   = 0;    // VmRSS after  set_data (kB)
    long   mem_peak_vmhwm_kb    = 0;    // VmHWM peak after tracing (kB)
};

inline double pf_now(bool enabled) {
    return enabled ? MPI_Wtime() : 0.0;
}

inline void pf_accum(double& acc, double start, bool enabled) {
    if (enabled) acc += MPI_Wtime() - start;
}

inline long pf_read_vmrss_kb() {
    std::ifstream f("/proc/self/status");
    std::string line;
    while (std::getline(f, line)) {
        if (line.rfind("VmRSS:", 0) == 0) {
            long kb = 0;
            std::sscanf(line.c_str(), "VmRSS: %ld kB", &kb);
            return kb;
        }
    }
    return 0;
}

inline long pf_read_vmhwm_kb() {
    std::ifstream f("/proc/self/status");
    std::string line;
    while (std::getline(f, line)) {
        if (line.rfind("VmHWM:", 0) == 0) {
            long kb = 0;
            std::sscanf(line.c_str(), "VmHWM: %ld kB", &kb);
            return kb;
        }
    }
    return 0;
}

#endif  // _TIMING_HPP
