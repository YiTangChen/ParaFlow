#ifndef _TIMING_HPP
#define _TIMING_HPP

#include <mpi.h>
#include <fstream>
#include <string>
#include <cstdio>

// Per-block timing accumulators (stored in each Block instance). All times are in seconds.
struct BlockTiming {
    // Initialization phase
    double t_block_load         = 0.0;  // set_data: NetCDF read + seed filter

    // Tracing phase (accumulated over all iexchange loop iterations)
    double t_trace_compute      = 0.0;  // GenStreamLines / GenPathLines (RK4 integration)
    double t_trace_comm         = 0.0;  // deq_incoming_iexchange (receive particles)

    // Output phase
    double t_output_write       = 0.0;  // write_trajectory

    // Work counters (useful for normalizing time and detecting load imbalance)
    int    n_seeds_initial      = 0;    // seeds assigned to this block at startup
    long   n_steps_total        = 0;    // total RK4 integration steps completed
    int    n_particles_received = 0;    // particles received from neighbor blocks

    // GPU-specific timing (zero on CPU-only runs)
    double t_gpu_kernel_ms      = 0.0;  // pure CUDA kernel time (CUDA Events), accumulated across all launches

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
