#!/usr/bin/env bash
# switch.sh – Switch build configuration between OSC ascend, OSC cardinal, and NERSC
#
# Usage: ./switch.sh {nersc|ascend|cardinal}
#
# Files modified:
#   build.sh
#   CMakeLists.txt
#   ParaFlow/CMakeLists.txt
#   ParaFlow/OSUFlow/CMakeLists.txt

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    echo "Usage: $0 {nersc|ascend|cardinal}"
    echo ""
    echo "  nersc    – NERSC Perlmutter (Cray MPI + Cray NetCDF)"
    echo "  ascend   – OSC Ascend       (MVAPICH + Spack NetCDF)"
    echo "  cardinal – OSC Cardinal     (MVAPICH + Spack NetCDF)"
    exit 1
}

[[ $# -eq 1 ]] || usage
TARGET=$(echo "$1" | tr '[:upper:]' '[:lower:]')
case "$TARGET" in nersc|ascend|cardinal) ;; *) usage ;; esac

BUILD="$SCRIPT_DIR/build.sh"
CMAKE_TOP="$SCRIPT_DIR/CMakeLists.txt"
CMAKE_PF="$SCRIPT_DIR/ParaFlow/CMakeLists.txt"
CMAKE_OS="$SCRIPT_DIR/ParaFlow/OSUFlow/CMakeLists.txt"

echo "==> Switching to: $TARGET"

# ── build.sh ─────────────────────────────────────────────────────────────────
# Step 1: Comment out any currently active CC/CXX/NETCDF_* assignment lines.
#         Lines already commented are left untouched (no double-commenting).
sed -i 's|^CC=|# CC=|'                 "$BUILD"
sed -i 's|^CXX=|# CXX=|'              "$BUILD"
sed -i 's|^NETCDF_DIR=|# NETCDF_DIR=|'         "$BUILD"
sed -i 's|^NETCDF_CXX_DIR=|# NETCDF_CXX_DIR=|' "$BUILD"

# Step 2: Uncomment the target machine's lines (matched by unique path prefix).
case "$TARGET" in
    nersc)
        sed -i 's|^# \(CC=/opt/cray/\)|\1|'    "$BUILD"
        sed -i 's|^# \(CXX=/opt/cray/\)|\1|'   "$BUILD"
        sed -i 's|^# \(NETCDF_DIR=${NETCDF_DIR:-/opt/cray/\)|\1|'         "$BUILD"
        sed -i 's|^# \(NETCDF_CXX_DIR=${NETCDF_CXX_DIR:-/opt/cray/\)|\1|' "$BUILD"
        ;;
    ascend)
        sed -i 's|^# \(CC=/apps/spack/0.21/ascend/\)|\1|'    "$BUILD"
        sed -i 's|^# \(CXX=/apps/spack/0.21/ascend/\)|\1|'   "$BUILD"
        sed -i 's|^# \(NETCDF_DIR=${NETCDF_DIR:-/apps/spack/0.21/ascend/\)|\1|'         "$BUILD"
        sed -i 's|^# \(NETCDF_CXX_DIR=${NETCDF_CXX_DIR:-/apps/spack/0.21/ascend/\)|\1|' "$BUILD"
        ;;
    cardinal)
        sed -i 's|^# \(CC=/apps/spack/0.21/cardinal/\)|\1|'    "$BUILD"
        sed -i 's|^# \(CXX=/apps/spack/0.21/cardinal/\)|\1|'   "$BUILD"
        sed -i 's|^# \(NETCDF_DIR=${NETCDF_DIR:-/apps/spack/0.21/cardinal/\)|\1|'         "$BUILD"
        sed -i 's|^# \(NETCDF_CXX_DIR=${NETCDF_CXX_DIR:-/apps/spack/0.21/cardinal/\)|\1|' "$BUILD"
        ;;
esac

# ── CMakeLists.txt files ──────────────────────────────────────────────────────
# Each file has two set(NetCDF_LIBRARIES ...) lines — one for NERSC, one for OSC —
# with one active and one commented out.
#
# Distinguishing strings:
#   NERSC: libnetcdf_c++4  (main + OSUFlow_MPASO_Parallel CMakeLists)
#          libnetcdf-c++4  (ParaFlow CMakeLists)
#   OSC:   libnetcdf-cxx4  (all three files)
#
# Note: sed BRE treats '+' as a literal character, so c++4 matches the
# literal string "c++4" without any escaping needed.

cmake_switch() {
    local file="$1"

    # Step 1: Comment out the currently active set(NetCDF_LIBRARIES ...) line.
    #         The address /^set[[:space:]]*(NetCDF_LIBRARIES/ only matches
    #         uncommented lines, so re-running is safe.
    sed -i '/^set[[:space:]]*(NetCDF_LIBRARIES/s/^/# /' "$file"

    # Step 2: Uncomment the target's set(NetCDF_LIBRARIES ...) line.
    if [[ "$TARGET" == "nersc" ]]; then
        # NERSC lines contain c++4 (either _c++4 or -c++4)
        sed -i 's|^# \(set.*NetCDF_LIBRARIES.*c++4\)|\1|' "$file"
    else
        # OSC (ascend/cardinal) lines contain cxx4
        sed -i 's|^# \(set.*NetCDF_LIBRARIES.*cxx4\)|\1|' "$file"
    fi
}

cmake_switch "$CMAKE_TOP"
cmake_switch "$CMAKE_PF"
cmake_switch "$CMAKE_OS"

echo "==> Done. Active configuration: $TARGET"
echo "    build.sh"
echo "    CMakeLists.txt"
echo "    ParaFlow/CMakeLists.txt"
echo "    ParaFlow/OSUFlow/CMakeLists.txt"
