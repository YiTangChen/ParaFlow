#!/usr/bin/env bash
# Build the GPU inline-LIC driver (ParaFlow_lic_gpu) against the CUDA-enabled
# OSUFlow in the GPUParaFlow tree. Unlike build_lic.sh (which links the CPU-only
# root OSUFlow), this links libOSUFlow.a built with -DOSUFLOW_ENABLE_CUDA plus
# the CUDA runtime, because the streamlines are integrated on the GPU.
#
# Prereqs (do these first, from THIS repo root — self-contained, no sibling tree needed):
#   1) bash switch.sh cardinal                                 # select the OSC Cardinal toolchain
#   2) OSUFLOW_ENABLE_CUDA=1 CUDA_ARCH=90 bash build.sh        # builds ParaFlow/OSUFlow/libOSUFlow.a (+yaml-cpp)
#   3) module load <mpi> <cudatoolkit>                         # so mpicxx + nvcc/cudart are visible
#      module load cuda/12.4.1   
# Then:
#   bash rendering/LIC/build_lic_gpu.sh                        # ROOT auto-resolves to this repo's own ParaFlow/
#
# Produces ./ParaFlow_lic_gpu at the repo root.
set -e

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# ROOT = the CUDA-enabled ParaFlow library dir (block.hpp, utils.hpp, OSUFlow/, yaml-cpp/, diy/).
# Self-contained layout: the core lib is THIS repo's own ParaFlow/ (i.e. $HERE/../../ParaFlow),
# searched FIRST so the clone builds using only its own tree. The GPUParaFlow / GPU_PARAFLOW
# spots (inside or beside the checkout) are kept only as fallbacks for the old layout.
# Override with ROOT=/path/to/ParaFlow if yours is elsewhere.
if [[ -z "${ROOT:-}" ]]; then
  for cand in \
      "$HERE/../../ParaFlow" \
      "$HERE/../../GPUParaFlow/ParaFlow"    "$HERE/../../GPU_PARAFLOW/ParaFlow" \
      "$HERE/../../../GPUParaFlow/ParaFlow" "$HERE/../../../GPU_PARAFLOW/ParaFlow"; do
    if [[ -f "$cand/OSUFlow/libOSUFlow.a" ]]; then ROOT="$(cd "$cand" && pwd)"; break; fi
  done
fi
if [[ -z "${ROOT:-}" || ! -f "$ROOT/OSUFlow/libOSUFlow.a" ]]; then
  echo "ERROR: could not find the CUDA-enabled ParaFlow tree (OSUFlow/libOSUFlow.a)." >&2
  echo "       Searched inside/beside ParaFlowClean for GPUParaFlow / GPU_PARAFLOW." >&2
  echo "       Set ROOT=/path/to/GPU_PARAFLOW/ParaFlow explicitly." >&2
  exit 1
fi
# Refuse to build a non-GPU binary: the archive MUST contain the CUDA tracer
# symbols, i.e. it was built with OSUFLOW_ENABLE_CUDA. Otherwise the link would
# fail confusingly (or, worse, mislead you into thinking it's a GPU build).
if command -v nm >/dev/null 2>&1 && \
   ! nm "$ROOT/OSUFlow/libOSUFlow.a" 2>/dev/null | grep -q CreateGPUBlockContext; then
  echo "ERROR: $ROOT/OSUFlow/libOSUFlow.a has NO GPU symbols (not built with CUDA)." >&2
  echo "       Rebuild it first:  cd \"$ROOT/..\" && OSUFLOW_ENABLE_CUDA=1 CUDA_ARCH=90 bash build.sh" >&2
  exit 1
fi

# Cardinal NetCDF (same trees build_lic.sh links); override via env if different.
NETCDF_DIR="${NETCDF_DIR:-/apps/spack/0.21/cardinal/linux-rhel9-sapphirerapids/netcdf-c/gcc/12.3.0/mvapich/3.0/4.9.2-qd5jghn}"
NETCDF_CXX_DIR="${NETCDF_CXX_DIR:-/apps/spack/0.21/cardinal/linux-rhel9-sapphirerapids/netcdf-cxx4/gcc/12.3.0/mvapich/3.0/4.3.1-tgr36hp}"

# CUDA runtime: locate libcudart via nvcc unless CUDA_HOME is set.
if [[ -z "$CUDA_HOME" ]] && command -v nvcc >/dev/null 2>&1; then
    CUDA_HOME="$(dirname "$(dirname "$(command -v nvcc)")")"
fi
if [[ -z "$CUDA_HOME" ]]; then
    echo "ERROR: CUDA not found. Set CUDA_HOME=<cuda toolkit root> or put nvcc in PATH." >&2
    exit 1
fi
CUDA_LIBDIR="$CUDA_HOME/lib64"
[[ -d "$CUDA_LIBDIR" ]] || CUDA_LIBDIR="$CUDA_HOME/lib"

REPO="$(cd "$HERE/../.." && pwd)"          # repo root (where the binary lands)
OUT="${OUT:-$REPO/ParaFlow_lic_gpu}"

echo "ROOT=$ROOT"
echo "REPO=$REPO"
echo "CUDA_HOME=$CUDA_HOME"
echo "building $OUT ..."

# The driver ParaFlow_lic_gpu.cpp lives HERE in rendering/LIC and #includes
# lic_render.hpp (the LIC method) from the same dir. ParaFlow.cpp is compiled in
# directly (not a prebuilt libParaFlow.a) so the driver calls the SAME
# ParaFlow::GenStreamLines used by ./ParaFlow_streamline and always picks up the
# current source. -DOSUFLOW_ENABLE_CUDA enables the GPU tracer path.
mpicxx -O2 -std=c++17 -DLIC_WITH_MPI -DOSUFLOW_ENABLE_CUDA \
  "$HERE/ParaFlow_lic_gpu.cpp" \
  "$ROOT/ParaFlow.cpp" \
  -I"$HERE" -I"$ROOT" -I"$ROOT/OSUFlow" -I"$ROOT/diy/include" \
  -I"$ROOT/yaml-cpp/include" -I"$NETCDF_DIR/include" -I"$CUDA_HOME/include" \
  "$ROOT/OSUFlow/libOSUFlow.a" \
  "$ROOT/yaml-cpp/libyaml-cpp.a" \
  "$NETCDF_DIR/lib/libnetcdf.so" \
  "$NETCDF_CXX_DIR/lib64/libnetcdf-cxx4.so" \
  -L"$CUDA_LIBDIR" -lcudart \
  -lpthread \
  -o "$OUT"

# Guard against something in the editor/sync layer stripping the execute bit
# after the build (observed happening between builds and runs on this setup).
chmod +x "$OUT"

echo "built $OUT"
