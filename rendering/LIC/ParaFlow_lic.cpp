// ParaFlow_lic.cpp
// -----------------------------------------------------------------------------
// CPU streamline-LIC. A thin driver, just like ParaFlow_streamline.cpp /
// ParaFlow_pathline.cpp: construct a ParaFlow, RUN the streamlines with the core
// method ParaFlow::GenStreamLines, then hand them to the LIC method (lic_render.hpp,
// same folder) to fold into an image and write one PNG.
//
// The streamline "run" is still the core GenStreamLines (outside this folder); this
// driver + lic_render.hpp are the LIC feature. CPU tracer is selected by
// `useGPU: false` in the config (see ParaFlow_lic_gpu.cpp for the GPU counterpart).
//
// Both directions are traced so the LIC line runs THROUGH each seed (config key
// lic_direction: forward | backward | both). Each rank's segments are gathered
// straight to rank 0 over MPI (never written to disk), which reassembles whole
// streamlines (read_traces.py logic) and folds those.
//
// Build:  bash rendering/LIC/build_lic.sh
// Run:    mpirun -np <nproc> ./ParaFlow_lic conf/lic_streamline.yaml
// -----------------------------------------------------------------------------
#include <iostream>
#include <list>
#include <vector>
#include <mpi.h>
#include <ParaFlow.hpp>
#include "lic_render.hpp"

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cout << "Needs 1 argument!\n";
        return 0;
    }
    ParaFlow pf(argc, argv, argv[1]);

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    LicRun run = lic_begin(argv[1], lic_count_seeds(argv[1]), "cpu");

    // Generate (writeToDisk=false -- never touches disk), gather each direction's
    // segments straight to rank 0 over MPI, then reassemble those into whole
    // streamlines and fold. The in-memory point-list copy is unused here
    // (segBytesOut carries the pid/sid-tagged data lic_finish actually needs),
    // so drop it after each direction to keep memory down.
    std::list<std::vector<VECTOR3>> discard;
    std::vector<LicSeg> fwdSegs, bwdSegs;
    for (double sign : lic_directions(argv[1], rank)) {
        std::vector<char> localBytes;
        pf.GenStreamLines(discard, sign, false, &localBytes);   // run streamline (outside rendering)
        discard.clear();
        lic_gather_segments(localBytes, sign < 0.0 ? bwdSegs : fwdSegs);
    }
    lic_finish(run, argv[1], fwdSegs, bwdSegs);
    return 0;
}
