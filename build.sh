#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
###### NERSC
CC=/opt/cray/pe/mpich/9.0.1/ofi/gnu/12.3/bin/mpicc
CXX=/opt/cray/pe/mpich/9.0.1/ofi/gnu/12.3/bin/mpicxx
CUDA_MODULE_NERSC=cudatoolkit
CUDA_ARCH_DEFAULT_NERSC=80

####### ascend
# CC=/apps/spack/0.21/ascend/linux-rhel9-zen2/mvapich/gcc/12.3.0/3.0-fndadii/bin/mpicc
# CXX=/apps/spack/0.21/ascend/linux-rhel9-zen2/mvapich/gcc/12.3.0/3.0-fndadii/bin/mpicxx
# CUDA_MODULE_ASCEND=cuda/12.4.1
# CUDA_ARCH_DEFAULT_ASCEND=80

####### cardinal
# CC=/apps/spack/0.21/cardinal/linux-rhel9-sapphirerapids/mvapich/gcc/12.3.0/3.0-mfh5vjl/bin/mpicc
# CXX=/apps/spack/0.21/cardinal/linux-rhel9-sapphirerapids/mvapich/gcc/12.3.0/3.0-mfh5vjl/bin/mpicxx
# CUDA_MODULE_CARDINAL=cuda/12.4.1
# CUDA_ARCH_DEFAULT_CARDINAL=90
JOBS=${JOBS:-$(nproc)}

# ---- optional CUDA backend ----
# Enable with: OSUFLOW_ENABLE_CUDA=1 ./build.sh
# Override arch (e.g. H100 = 90) with: CUDA_ARCH=90 ./build.sh
OSUFLOW_ENABLE_CUDA=${OSUFLOW_ENABLE_CUDA:-0}

# Pick up whichever per-cluster defaults are currently active (the uncommented
# block above). Only one of the *_NERSC / *_ASCEND / *_CARDINAL variables will
# actually be defined after sourcing switch.sh's edits.
CUDA_MODULE=${CUDA_MODULE:-${CUDA_MODULE_NERSC:-${CUDA_MODULE_ASCEND:-${CUDA_MODULE_CARDINAL:-}}}}
CUDA_ARCH=${CUDA_ARCH:-${CUDA_ARCH_DEFAULT_NERSC:-${CUDA_ARCH_DEFAULT_ASCEND:-${CUDA_ARCH_DEFAULT_CARDINAL:-80}}}}

CUDA_CMAKE_ARGS=()
if [[ "$OSUFLOW_ENABLE_CUDA" == "1" ]]; then
    # Try to load the cluster's CUDA module if `module` is available and a
    # module name was configured. Non-fatal: user may already have loaded it.
    if [[ -n "$CUDA_MODULE" ]] && command -v module >/dev/null 2>&1; then
        echo "==> module load $CUDA_MODULE"
        module load "$CUDA_MODULE" || echo "    (module load failed — continuing; ensure nvcc is in PATH)"
    fi
    if ! command -v nvcc >/dev/null 2>&1; then
        echo "ERROR: OSUFLOW_ENABLE_CUDA=1 but nvcc not found in PATH." >&2
        echo "       Either load a CUDA module manually or set CUDA_MODULE=<name>." >&2
        exit 1
    fi
    CUDA_CMAKE_ARGS+=(-DOSUFLOW_ENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES="$CUDA_ARCH")
    echo "==> CUDA backend enabled (arch=$CUDA_ARCH, nvcc=$(command -v nvcc))"
fi

####### NERSC
NETCDF_DIR=${NETCDF_DIR:-/opt/cray/pe/netcdf/4.9.2.1/gnu/12.3}
NETCDF_CXX_DIR=${NETCDF_CXX_DIR:-/opt/cray/pe/netcdf/4.9.2.1/gnu/12.3}

####### ascend
# NETCDF_DIR=${NETCDF_DIR:-/apps/spack/0.21/ascend/linux-rhel9-zen2/netcdf-c/gcc/12.3.0/mvapich/3.0/4.9.2-j2kquab/}
# NETCDF_CXX_DIR=${NETCDF_CXX_DIR:-/apps/spack/0.21/ascend/linux-rhel9-zen2/netcdf-cxx4/gcc/12.3.0/mvapich/3.0/4.3.1-aenx4lr/}

####### cardinal
# NETCDF_DIR=${NETCDF_DIR:-/apps/spack/0.21/cardinal/linux-rhel9-sapphirerapids/netcdf-c/gcc/12.3.0/mvapich/3.0/4.9.2-qd5jghn/}
# NETCDF_CXX_DIR=${NETCDF_CXX_DIR:-/apps/spack/0.21/cardinal/linux-rhel9-sapphirerapids/netcdf-cxx4/gcc/12.3.0/mvapich/3.0/4.3.1-tgr36hp/}

cmake_build() {
    local dir="$1"; shift
    echo "==> Building: $dir"
    cd "$dir"
    cmake -DCMAKE_C_COMPILER="$CC" -DCMAKE_CXX_COMPILER="$CXX" \
          -DNETCDF_DIR="$NETCDF_DIR" -DNETCDF_CXX_DIR="$NETCDF_CXX_DIR" "$@" .
    make clean
    make -j"$JOBS"
    cd "$SCRIPT_DIR"
}

cmake_build ParaFlow/OSUFlow "${CUDA_CMAKE_ARGS[@]}"

# yaml-cpp is an external library; skip make clean to avoid unnecessary rebuilds
echo "==> Building: ParaFlow/yaml-cpp"
cd ParaFlow/yaml-cpp
cmake -DCMAKE_C_COMPILER="$CC" -DCMAKE_CXX_COMPILER="$CXX" .
make -j"$JOBS"
cd "$SCRIPT_DIR"

# cmake_build ParaFlow/yaml-cpp -DYAML_BUILD_SHARED_LIBS=ON

cmake_build ParaFlow -DDIY_INCLUDE_DIRS=diy/include -Dyaml-cpp_DIR="$SCRIPT_DIR/ParaFlow/yaml-cpp" "${CUDA_CMAKE_ARGS[@]}"
cmake_build . -DDIY_INCLUDE_DIRS=ParaFlow/diy/include -Dyaml-cpp_DIR="$SCRIPT_DIR/ParaFlow/yaml-cpp" "${CUDA_CMAKE_ARGS[@]}"

echo "==> Done."
