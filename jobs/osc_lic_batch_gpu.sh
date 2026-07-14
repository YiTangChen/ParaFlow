#!/bin/bash
# osc_lic_batch_gpu.sh — run ParaFlow_lic_gpu once per config listed in a
# batch manifest (see conf/lic_batch_gpu.yaml), one after another.
#
# Usage (inside a salloc, same allocation lic_streamline_gpu.yaml expects):
#   bash jobs/osc_lic_batch_gpu.sh conf/lic_batch_gpu.yaml

set -euo pipefail

# ParaFlow_lic_gpu needs libcudart (from the CUDA module) and libnetcdf-cxx4
# (linked by absolute path with no rpath in rendering/LIC/build_lic_gpu.sh, so
# it won't resolve at runtime without this) -- same fix jobs/osc_run.sh applies
# for the ascend cluster's libnetcdf-cxx4.
module load cuda/12.4.1
export LD_LIBRARY_PATH=/apps/spack/0.21/cardinal/linux-rhel9-sapphirerapids/netcdf-cxx4/gcc/12.3.0/mvapich/3.0/4.3.1-tgr36hp/lib64:/apps/spack/0.21/cardinal/linux-rhel9-sapphirerapids/netcdf-c/gcc/12.3.0/mvapich/3.0/4.9.2-qd5jghn/lib:${LD_LIBRARY_PATH:-}

MANIFEST="${1:-conf/lic_batch_gpu.yaml}"

# Pull the "  - path" lines under the manifest's configs: list. Flat-list
# parsing only — no nested YAML support needed for this.
mapfile -t CONFIGS < <(grep -E '^\s*-\s' "$MANIFEST" | sed -E 's/^\s*-\s*//')

if [[ ${#CONFIGS[@]} -eq 0 ]]; then
    echo "No configs found in $MANIFEST" >&2
    exit 1
fi

for cfg in "${CONFIGS[@]}"; do
    echo "==> Running ParaFlow_lic_gpu $cfg"
    # --mpi=pmix (NOT pmi2): this cluster's mvapich links libpmix, not libpmi2 --
    # pmi2 leaves every rank in its own singleton MPI_COMM_WORLD (see
    # conf/lic_streamline_gpu.yaml's header comment for the failure mode).
    srun --export=ALL --mpi=pmix -n 16 --ntasks-per-node=4 --gpus-per-task=1 ./ParaFlow_lic_gpu "$cfg"
done
