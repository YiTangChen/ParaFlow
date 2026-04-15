#include "MPASODeviceField.h"

#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>

namespace {

inline void cudaCheck(cudaError_t err, const char* what) {
    if (err != cudaSuccess) {
        std::fprintf(stderr, "[MPASODeviceField] CUDA error in %s: %s\n",
                     what, cudaGetErrorString(err));
        std::abort();
    }
}

template <typename T>
void uploadArray(T*& d_ptr, const T* h_src, size_t n, const char* tag) {
    if (d_ptr == nullptr) {
        cudaCheck(cudaMalloc(reinterpret_cast<void**>(&d_ptr), n * sizeof(T)), tag);
    }
    cudaCheck(cudaMemcpy(d_ptr, h_src, n * sizeof(T), cudaMemcpyHostToDevice), tag);
}

template <typename T>
void freeArray(T*& d_ptr) {
    if (d_ptr != nullptr) {
        cudaFree(d_ptr);
        d_ptr = nullptr;
    }
}

} // namespace

void MPASODeviceField::uploadTopology(const mpaso_vec3* h_cell_coord,
                                      const mpaso_vec3* h_vertex_coord,
                                      const int*        h_vertices_on_cell,
                                      const int*        h_cells_on_cell,
                                      const int*        h_num_vertices_on_cell,
                                      const int*        h_max_level_cell)
{
    const size_t n_cell_edge = static_cast<size_t>(n_cells) * n_max_edges;

    uploadArray(d_cell_coord,           h_cell_coord,           static_cast<size_t>(n_cells),          "cell_coord");
    uploadArray(d_vertex_coord,         h_vertex_coord,         static_cast<size_t>(n_local_vertices), "vertex_coord");
    uploadArray(d_vertices_on_cell,     h_vertices_on_cell,     n_cell_edge,                           "vertices_on_cell");
    uploadArray(d_cells_on_cell,        h_cells_on_cell,        n_cell_edge,                           "cells_on_cell");
    uploadArray(d_num_vertices_on_cell, h_num_vertices_on_cell, static_cast<size_t>(n_cells),          "num_vertices_on_cell");
    uploadArray(d_max_level_cell,       h_max_level_cell,       static_cast<size_t>(n_cells),          "max_level_cell");
}

void MPASODeviceField::uploadVelocityWindow(const mpaso_vec3* h_vertex_velocity,
                                            const mpaso_vec3* h_vertex_vert_velocity,
                                            const double*     h_vertex_ztop,
                                            const double*     h_timestamps)
{
    const size_t n_vert_nodes = static_cast<size_t>(n_local_vertices) * n_vert_levels;
    const size_t n_window     = n_vert_nodes * static_cast<size_t>(n_timesteps_loaded);

    uploadArray(d_vertex_velocity,      h_vertex_velocity,      n_window,                                     "vertex_velocity");
    uploadArray(d_vertex_vert_velocity, h_vertex_vert_velocity, n_window,                                     "vertex_vert_velocity");
    uploadArray(d_vertex_ztop,          h_vertex_ztop,          n_window,                                     "vertex_ztop");
    uploadArray(d_timestamps,           h_timestamps,           static_cast<size_t>(n_timesteps_loaded),      "timestamps");
}

void MPASODeviceField::release()
{
    freeArray(d_cell_coord);
    freeArray(d_vertex_coord);
    freeArray(d_vertices_on_cell);
    freeArray(d_cells_on_cell);
    freeArray(d_num_vertices_on_cell);
    freeArray(d_max_level_cell);

    freeArray(d_vertex_velocity);
    freeArray(d_vertex_vert_velocity);
    freeArray(d_vertex_ztop);
    freeArray(d_timestamps);
}
