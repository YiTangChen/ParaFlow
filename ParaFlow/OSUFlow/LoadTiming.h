#ifndef _LOAD_TIMING_H
#define _LOAD_TIMING_H

// Lightweight accumulator for decomposing block_load (set_data) into its
// sub-phases. The stages live deep in OSUFlow/MPASOReader, so rather than thread
// a timing pointer through every signature we accumulate into one process-global
// instance: DIY runs a single thread per rank and blocks are loaded sequentially,
// so there is no concurrency. ParaFlow resets it before each set_data, reads it
// after, and clears `enabled` so later pathline window reloads are not counted.
//
// std::chrono (not MPI_Wtime) keeps OSUFlow free of an MPI dependency; the two
// clocks agree to well within the resolution of these second-scale measurements.

#include <chrono>

namespace load_timing {

struct Breakdown {
    bool   enabled       = false;  // gate: zero accumulation when false
    double grid_build_s  = 0.0;    // MPASOGrid build + cell/vertex index mappings
    double read_s        = 0.0;    // NetCDF field reads per timestep (velocity/zTop/vertVel)
    double convert_s     = 0.0;    // cell->vertex conversion per timestep
    double seed_filter_s = 0.0;    // inBlock() seed filtering
    void reset() { grid_build_s = read_s = convert_s = seed_filter_s = 0.0; }
};

extern Breakdown g;   // defined once in MPASOReader.C

using clock = std::chrono::steady_clock;
inline clock::time_point tic() { return clock::now(); }
inline double toc(clock::time_point t0) {
    return std::chrono::duration<double>(clock::now() - t0).count();
}

}  // namespace load_timing

#endif  // _LOAD_TIMING_H
