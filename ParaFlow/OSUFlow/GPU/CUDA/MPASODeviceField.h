#ifndef MPASO_DEVICE_FIELD_H
#define MPASO_DEVICE_FIELD_H

#include <cstddef>

// Plain-old-data vec3 used on both host and device.
// Intentionally does NOT include VectorMatrix.h so that this header can be
// pulled into .cu translation units without dragging OSUFlow host code through
// nvcc.
struct mpaso_vec3 {
    double x;
    double y;
    double z;
};

// Device-resident copy of the MPAS-O grid + one velocity window.
//
// Lifetime model:
//   - Topology buffers (cell_coord, vertex_coord, vertices_on_cell, ...) are
//     uploaded once per block via uploadTopology() and stay on the device.
//   - Velocity / ztop / vertVel buffers are uploaded per time window via
//     uploadVelocityWindow(). Pathline re-uploads as the window slides; a
//     streamline calls it exactly once with n_timesteps == 1.
//
// All "vertex" arrays are indexed as [ts * n_vert_nodes + vertexIdx] where
// n_vert_nodes = n_local_vertices * n_vert_levels, matching the layout that
// MPASOReader::InitSolutions already produces on the host.
struct MPASODeviceField {
    // ----- dimensions -----
    int n_cells            = 0;
    int n_local_vertices   = 0;
    int n_max_edges        = 0;
    int n_vert_levels      = 0;
    int n_vert_levels_p1   = 0;
    int n_timesteps_loaded = 0;

    // ----- global scalars -----
    double earth_radius    = 0.0;  // mirrors MPASOGrid::earth_radius
    // lat/lon bounding box (radians), mirrors MPASOGrid::isInBBox. The domain is
    // not global in longitude, so the CPU terminates trajectories that leave
    // this box; the GPU must do the same to match.
    double lat_min         = 0.0;
    double lat_max         = 0.0;
    double lon_center      = 0.0;
    double lon_half_width  = 0.0;

    // ----- topology (time-independent, device pointers) -----
    mpaso_vec3* d_cell_coord           = nullptr;  // [n_cells]
    mpaso_vec3* d_vertex_coord         = nullptr;  // [n_local_vertices]
    int*        d_vertices_on_cell     = nullptr;  // [n_cells * n_max_edges]
    int*        d_cells_on_cell        = nullptr;  // [n_cells * n_max_edges]
    int*        d_num_vertices_on_cell = nullptr;  // [n_cells]
    int*        d_max_level_cell       = nullptr;  // [n_cells]

    // ----- velocity window (time-varying, device pointers) -----
    // Indexing: [ts * (n_local_vertices * n_vert_levels) + vertNodeIdx]
    mpaso_vec3* d_vertex_velocity      = nullptr;  // horizontal, per-layer
    mpaso_vec3* d_vertex_vert_velocity = nullptr;  // vertical,   per-layer
    double*     d_vertex_ztop          = nullptr;  // [ts * n_local_vertices * n_vert_levels]

    // Real-time timestamps for the loaded window (seconds since t=0).
    double* d_timestamps = nullptr;                // [n_timesteps_loaded]

    // ----- lifecycle -----
    // Upload topology from host pointers (MPASOGrid). Safe to call once per
    // block. Sizes must match the n_* fields which the caller should set first.
    void uploadTopology(const mpaso_vec3* h_cell_coord,
                        const mpaso_vec3* h_vertex_coord,
                        const int*        h_vertices_on_cell,
                        const int*        h_cells_on_cell,
                        const int*        h_num_vertices_on_cell,
                        const int*        h_max_level_cell);

    // Upload (or refresh) one velocity window. n_timesteps_loaded must be set
    // before the call. h_vertex_velocity / h_vertex_vert_velocity point to
    // [n_timesteps_loaded * n_local_vertices * n_vert_levels] vec3s.
    void uploadVelocityWindow(const mpaso_vec3* h_vertex_velocity,
                              const mpaso_vec3* h_vertex_vert_velocity,
                              const double*     h_vertex_ztop,
                              const double*     h_timestamps);

    // Free all device buffers. Safe to call multiple times.
    void release();
};

#endif // MPASO_DEVICE_FIELD_H
