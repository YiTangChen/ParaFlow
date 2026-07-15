// ParaFlow_lic_gpu.cpp
// -----------------------------------------------------------------------------
// GPU streamline-LIC. Identical in shape to ParaFlow_lic.cpp (and to
// ParaFlow_streamline.cpp): construct a ParaFlow, RUN the streamlines with the
// core method ParaFlow::GenStreamLines, then hand them to the SAME LIC method
// (lic_render.hpp, same folder). The only difference is that this binary
// integrates on the GPU -- it requires an OSUFlow built with -DOSUFLOW_ENABLE_CUDA
// and `useGPU: true` in the config, and aborts if no CUDA device is visible (so it
// never silently falls back to a CPU-traced image).
//
// The streamline "run" is still the core GenStreamLines; this driver + lic_render.hpp
// are the LIC feature.
//
// Build:  bash rendering/LIC/build_lic_gpu.sh
// Run:    srun -n <nproc> --gpus-per-task=1 ./ParaFlow_lic_gpu conf/lic_streamline_gpu.yaml
// -----------------------------------------------------------------------------
#include <iostream>
#include <list>
#include <vector>
#include <mpi.h>
#include <ParaFlow.hpp>
#include "GPU/CUDA/MPASOGPUTracer.h"          // mpaso_gpu_host::isAvailable
#include "lic_render.hpp"

#ifndef OSUFLOW_ENABLE_CUDA
#error "ParaFlow_lic_gpu.cpp requires an OSUFlow built with -DOSUFLOW_ENABLE_CUDA (see build_lic_gpu.sh)."
#endif

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cout << "Needs 1 argument!\n";
        return 0;
    }
    ParaFlow pf(argc, argv, argv[1]);

    int rank = 0, nranks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);

    // This binary integrates on the GPU by definition. If no CUDA device is
    // visible, GenStreamLines would silently fall back to the CPU tracer, so abort
    // loudly instead of quietly producing a CPU-traced image.
    if (rank == 0)
        std::cerr << "[lic] MPI world size = " << nranks
                  << "  (if this is 1 with -n>1, MPI didn't form a world -- check --mpi=pmi2)\n";
    if (!mpaso_gpu_host::isAvailable()) {
        std::cerr << "[lic] FATAL rank " << rank << ": no CUDA device visible. This binary "
                     "REQUIRES a GPU (and useGPU:true) -- run on a --gpus node. Aborting.\n";
        MPI_Abort(MPI_COMM_WORLD, 2);
    }

    LicRun run = lic_begin(argv[1], lic_count_seeds(argv[1]), "gpu");

    // Generate on the GPU (writeToDisk=false -- never touches disk), gather each
    // direction's segments straight to rank 0 over MPI, then reassemble those
    // into whole streamlines and fold. In-memory point-list copy unused here
    // (segBytesOut carries the pid/sid-tagged data lic_finish actually needs),
    // so drop it after each direction.
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
