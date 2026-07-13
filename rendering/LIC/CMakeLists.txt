# Standalone post-processing tool: streamline-footprint LIC renderer.
# Decoupled from the main MPI ParaFlow build (no MPI / NetCDF / OSUFlow deps).
#
#   cmake -S . -B build && cmake --build build
# or just:
#   g++ -O2 -std=c++17 streamline_lic.cpp -o streamline_lic

cmake_minimum_required(VERSION 3.13)
project(streamline_lic CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(streamline_lic streamline_lic.cpp)

if(NOT MSVC)
  target_compile_options(streamline_lic PRIVATE -O2)
endif()
