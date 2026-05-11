#ifndef _TIMING_HPP
#define _TIMING_HPP

#include <mpi.h>
#include <time.h>
#include <fstream>
#include <string>
#include <cstdio>

// Per-block timing accumulators (stored in each Block instance). All times are in seconds.
struct BlockTiming {
    // ── Initialization phase ──────────────────────────────────────────────────
    // (block-load total is derived at parse time as t_netcdf_read + t_seed_filter)
    double t_netcdf_read        = 0.0;  // LoadMPASOData: NetCDF I/O + mesh topology build
    double t_seed_filter        = 0.0;  // inBlock() scan over all broadcast seeds

    // ── Tracing phase (accumulated over all iexchange rounds) ─────────────────
    double t_trace_compute      = 0.0;  // GenStreamLines / GenPathLines (RK4 integration)
    double t_trace_enqueue      = 0.0;  // Segment build/packing + cp.enqueue()
    double t_dequeue_local      = 0.0;  // deq_incoming_iexchange: local dequeue/receive work (no network)
    double t_fill_incoming      = 0.0;  // cp.fill_incoming(): exposed wait/progress/completion

    // ── Output phase ─────────────────────────────────────────────────────────
    double t_output_write       = 0.0;  // write_trajectory

    // ── Work counters ─────────────────────────────────────────────────────────
    int    n_iex_rounds         = 0;    // iexchange do-while iterations (for per-round averages)
    int    n_seeds_initial      = 0;    // seeds assigned to this block at startup
    long   n_steps_total        = 0;    // total RK4 integration steps completed
    int    n_particles_received = 0;    // particles received from neighbor blocks
    int    n_particles_sent     = 0;    // particles forwarded to neighbor blocks

    // ── Memory measurements (analytical formula) ──────────────────────────────
    size_t mem_grid_bytes       = 0;    // MPASOGrid: topology arrays + coordinates
    size_t mem_solution_bytes   = 0;    // Solution: velocity field data (scales with timesteps)

    // ── Memory measurements (OS-reported, from /proc/self/status) ────────────
    long   mem_vmrss_before_kb  = 0;    // VmRSS before set_data (kB)
    long   mem_vmrss_after_kb   = 0;    // VmRSS after  set_data (kB)
    long   mem_peak_vmhwm_kb    = 0;    // VmHWM peak after tracing (kB)
};

// Use clock_gettime(CLOCK_MONOTONIC) instead of MPI_Wtime().
// On Linux this is a vDSO call (~3-10 ns) vs MPI_Wtime (~100-200 ns).
// The difference matters when pf_now/pf_accum are called thousands of times
// per block in tight iexchange loops — especially with small seed counts where
// total compute time may be only a few milliseconds.
inline double pf_now(bool enabled) {
    if (!enabled) return 0.0;
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

inline void pf_accum(double& acc, double start, bool enabled) {
    if (!enabled) return;
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    acc += ((double)ts.tv_sec + (double)ts.tv_nsec * 1e-9) - start;
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
