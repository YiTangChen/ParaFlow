#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
###### NERSC
CC=/opt/cray/pe/mpich/9.0.1/ofi/gnu/12.3/bin/mpicc
CXX=/opt/cray/pe/mpich/9.0.1/ofi/gnu/12.3/bin/mpicxx

####### ascend
# CC=/apps/spack/0.21/ascend/linux-rhel9-zen2/mvapich/gcc/12.3.0/3.0-fndadii/bin/mpicc
# CXX=/apps/spack/0.21/ascend/linux-rhel9-zen2/mvapich/gcc/12.3.0/3.0-fndadii/bin/mpicxx

####### cardinal
# CC=/apps/spack/0.21/cardinal/linux-rhel9-sapphirerapids/mvapich/gcc/12.3.0/3.0-mfh5vjl/bin/mpicc
# CXX=/apps/spack/0.21/cardinal/linux-rhel9-sapphirerapids/mvapich/gcc/12.3.0/3.0-mfh5vjl/bin/mpicxx
JOBS=${JOBS:-$(nproc)}

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

cmake_build ParaFlow/OSUFlow

# yaml-cpp is an external library; skip make clean to avoid unnecessary rebuilds
echo "==> Building: ParaFlow/yaml-cpp"
cd ParaFlow/yaml-cpp
cmake -DCMAKE_C_COMPILER="$CC" -DCMAKE_CXX_COMPILER="$CXX" .
make -j"$JOBS"
cd "$SCRIPT_DIR"

# cmake_build ParaFlow/yaml-cpp -DYAML_BUILD_SHARED_LIBS=ON

cmake_build ParaFlow -DDIY_INCLUDE_DIRS=diy/include -Dyaml-cpp_DIR="$SCRIPT_DIR/ParaFlow/yaml-cpp"
cmake_build . -DDIY_INCLUDE_DIRS=ParaFlow/diy/include -Dyaml-cpp_DIR="$SCRIPT_DIR/ParaFlow/yaml-cpp"

echo "==> Done."
