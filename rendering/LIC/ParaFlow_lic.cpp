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
// lic_direction: forward | backward | both).
//
// Two strategies, chosen by config key `lic_combine` (for A/B comparison):
//   false = fold each direction's per-block SEGMENTS locally, reduce the images
//   true  = rank 0 reads the <gid>.bin files, reassembles whole streamlines
//           (read_traces.py logic), and folds those  [use a single direction]
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

    if (lic_combine(argv[1])) {
        // METHOD 1 (combine): generate (writes per-block <gid>.bin), then rank 0
        // reassembles those into whole streamlines and folds. The in-memory copy is
        // unused here, so drop it after each direction to keep memory down.
        std::list<std::vector<VECTOR3>> discard;
        for (double sign : lic_directions(argv[1], rank)) {
            pf.GenStreamLines(discard, sign);        // run streamline (outside rendering)
            discard.clear();
        }
        lic_combine_finish(run, argv[1]);
    } else {
        // METHOD 2 (no combine): fold each direction's segments locally, reduce images.
        // streamline_storage: memory -> the traces live ONLY here, in RAM (writeBin=false
        // skips the <gid>.bin write entirely, see block.hpp write_trajectory): fold them
        // into the image and free them immediately, so peak memory is ONE direction's
        // traces (not both) and nothing ever touches disk.
        const bool writeBin = lic_write_bin_to_disk(argv[1]);
        for (double sign : lic_directions(argv[1], rank)) {
            std::list<std::vector<VECTOR3>> sls;
            pf.GenStreamLines(sls, sign, writeBin);  // run streamline (outside rendering)
            lic_fold(run, sls);                      // fold this direction ...
            sls.clear();                             // ... then free it before the next one
        }
        lic_finish(run, argv[1]);
    }
    return 0;
}
