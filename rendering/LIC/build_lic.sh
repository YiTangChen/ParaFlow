#!/usr/bin/env bash
# Build the CPU streamline-LIC driver (ParaFlow_lic) against the already-built
# OSUFlow + yaml-cpp libraries. The driver delegates tracing to the SAME core
# method ParaFlow::GenStreamLines used by ./ParaFlow_streamline, so ParaFlow.cpp
# is compiled in directly (not a prebuilt libParaFlow.a) and always picks up the
# current source. No CUDA: this is the CPU-only OSUFlow path (see build_lic_gpu.sh
# for the GPU build of the SAME ParaFlow_lic.cpp source).
#
# Prereqs (do these first):
#   1) bash build.sh                              # builds libOSUFlow.a + libyaml-cpp.a
#   2) module reset && module load gcc/12.3.0 mvapich/3.0
# Then:
#   bash rendering/LIC/build_lic.sh               # (run from repo root, or anywhere)
#
# Produces ./ParaFlow_lic at the repo root, next to ParaFlow_streamline.
set -e

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# ROOT = the ParaFlow library dir (has ParaFlow.cpp, block.hpp, utils.hpp, OSUFlow/,
# yaml-cpp/, diy/). Self-contained layout: this repo's own ParaFlow/ (i.e.
# $HERE/../../ParaFlow), searched FIRST. Override with ROOT=/path/to/ParaFlow.
if [[ -z "${ROOT:-}" ]]; then
  for cand in "$HERE/../../ParaFlow" "$HERE/../.."; do
    if [[ -f "$cand/ParaFlow.cpp" ]]; then ROOT="$(cd "$cand" && pwd)"; break; fi
  done
fi
if [[ -z "${ROOT:-}" || ! -f "$ROOT/ParaFlow.cpp" ]]; then
  echo "ERROR: could not find the ParaFlow core (ParaFlow.cpp). Set ROOT=/path/to/ParaFlow." >&2
  exit 1
fi

# Cardinal NetCDF (same paths build.sh uses); override via env if different.
NETCDF_DIR="${NETCDF_DIR:-/apps/spack/0.21/cardinal/linux-rhel9-sapphirerapids/netcdf-c/gcc/12.3.0/mvapich/3.0/4.9.2-qd5jghn}"
NETCDF_CXX_DIR="${NETCDF_CXX_DIR:-/apps/spack/0.21/cardinal/linux-rhel9-sapphirerapids/netcdf-cxx4/gcc/12.3.0/mvapich/3.0/4.3.1-tgr36hp}"

REPO="$(cd "$HERE/../.." && pwd)"          # repo root (where the binary lands)
OUT="${OUT:-$REPO/ParaFlow_lic}"           # binary lands next to ParaFlow_streamline

echo "ROOT=$ROOT"
echo "REPO=$REPO"
echo "building $OUT ..."

# The driver ParaFlow_lic.cpp lives HERE in rendering/LIC and #includes lic_render.hpp
# (the LIC method) from the same dir. ParaFlow.cpp is compiled in directly (not a
# prebuilt libParaFlow.a) so the driver calls the SAME ParaFlow::GenStreamLines used
# by ./ParaFlow_streamline and always picks up the current source. No
# -DOSUFLOW_ENABLE_CUDA here: this is the CPU build.
mpicxx -O2 -std=c++17 -DLIC_WITH_MPI \
  "$HERE/ParaFlow_lic.cpp" \
  "$ROOT/ParaFlow.cpp" \
  -I"$HERE" -I"$ROOT" -I"$ROOT/OSUFlow" -I"$ROOT/diy/include" \
  -I"$ROOT/yaml-cpp/include" -I"$NETCDF_DIR/include" \
  "$ROOT/OSUFlow/libOSUFlow.a" \
  "$ROOT/yaml-cpp/libyaml-cpp.a" \
  "$NETCDF_DIR/lib/libnetcdf.so" \
  "$NETCDF_CXX_DIR/lib64/libnetcdf-cxx4.so" \
  -lpthread \
  -o "$OUT"

# Guard against something in the editor/sync layer stripping the execute bit
# after the build (observed happening between builds and runs on this setup).
chmod +x "$OUT"

echo "built $OUT"
