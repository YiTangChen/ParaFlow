#ifndef MPASO_VELOCITY_CUH
#define MPASO_VELOCITY_CUH

// Device mirror of OSUFlow's CPU MPAS-O velocity sampling pipeline. The
// functions here are a line-by-line port of:
//
//   MPASOGrid::phys_to_cell()  in Grid.C
//   CVectorField::at_phys(...VECTOR4...) in Field.C
//   MPASOGrid::getCellVertices() (T5_CELL branch) in Grid.C
//   MPASOGrid::interpolate(VECTOR3*, int, double*) in Grid.C
//
// Correctness target is option (b): acceptably close, not bit-identical.
// Arithmetic order is kept close to CPU wherever practical so parity
// tests are meaningful.

#include "../MPASODeviceField.h"

namespace mpaso_gpu_dev {

// Match CPU's stack-array limits. See Grid.C:1315 and Grid.C:1341.
constexpr int DEV_MAX_EDGES      = 10;   // MPASO_MAX_EDGES
constexpr int DEV_MAX_VERT_LEV   = 120;  // MAX_VERT_LEVELS

// ---------- vec3 helpers ----------

__device__ inline mpaso_vec3 v_sub(mpaso_vec3 a, mpaso_vec3 b) {
    return { a.x - b.x, a.y - b.y, a.z - b.z };
}
__device__ inline mpaso_vec3 v_cross(mpaso_vec3 a, mpaso_vec3 b) {
    return { a.y * b.z - a.z * b.y,
             a.z * b.x - a.x * b.z,
             a.x * b.y - a.y * b.x };
}
__device__ inline double v_dot(mpaso_vec3 a, mpaso_vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
__device__ inline double v_mag(mpaso_vec3 a) {
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}
__device__ inline mpaso_vec3 v_scale(mpaso_vec3 a, double s) {
    return { a.x * s, a.y * s, a.z * s };
}
__device__ inline mpaso_vec3 v_normalize_scale(mpaso_vec3 a, double target_len) {
    double m = v_mag(a);
    if (m == 0.0) return { 0.0, 0.0, 0.0 };
    double s = target_len / m;
    return { a.x * s, a.y * s, a.z * s };
}

// Device mirror of MPASOGrid::isInBBox (Grid.C:1946) + xyz2latlon/wrapToPi
// (Grid.h:471,486). The MPAS domain is not global in longitude; the CPU
// rejects points outside this lat/lon box and terminates the trajectory there,
// so the GPU must too.
__device__ inline bool in_bbox_dev(const MPASODeviceField& f, mpaso_vec3 pos) {
    double lon = atan2(pos.y, pos.x);
    double lat = atan2(pos.z, hypot(pos.x, pos.y));
    lon = fmod(lon, 2.0 * M_PI);
    if (lon < 0.0) lon += 2.0 * M_PI;
    if (lat < f.lat_min || lat > f.lat_max) return false;
    double dlon = fmod((lon - f.lon_center) + M_PI, 2.0 * M_PI);  // wrapToPi
    if (dlon < 0.0) dlon += 2.0 * M_PI;
    dlon -= M_PI;
    if (fabs(dlon) > f.lon_half_width) return false;
    return true;
}

// Mirrors Grid.C:1140 triangle_area().
__device__ inline double triangle_area_dev(mpaso_vec3 v1, mpaso_vec3 v2, mpaso_vec3 v3) {
    mpaso_vec3 u = v_sub(v2, v1);
    mpaso_vec3 w = v_sub(v3, v1);
    return v_mag(v_cross(u, w)) * 0.5;
}

// ---------- geodesic step on device ----------
//
// Mirrors vtCFieldLine::geodesic_step (FieldLine.C:173-195). Forward-only
// (time_dir = +1); streamline/pathline tracing never flips direction on GPU.
// Uses Rodrigues' formula instead of building a 3x3 rotation matrix:
//   v_rot = v cos t + (k x v) sin t + k (k . v)(1 - cos t)
// Returns false if the particle dives below the planetary centre.
__device__ inline bool geodesic_step_dev(
    mpaso_vec3 pt_src, double r_src,
    mpaso_vec3 h_vel,  double v_vel,
    double fdt,
    mpaso_vec3* pt_dst, double* r_dst)
{
    const double vel_mag = v_mag(h_vel) * fdt;
    if (vel_mag == 0.0) {
        *pt_dst = pt_src;
        *r_dst  = r_src;
        return true;
    }
    const double r_new = r_src + v_vel * fdt;
    if (r_new <= 0.0) return false;

    const double omega  = vel_mag / r_src;
    mpaso_vec3  axis    = v_cross(pt_src, h_vel);
    const double amag   = v_mag(axis);
    if (amag == 0.0) {
        *pt_dst = pt_src;
        *r_dst  = r_new;
        // renormalise to new radius
        const double s = r_new / r_src;
        pt_dst->x = pt_src.x * s;
        pt_dst->y = pt_src.y * s;
        pt_dst->z = pt_src.z * s;
        return true;
    }
    const double inv_a = 1.0 / amag;
    axis.x *= inv_a; axis.y *= inv_a; axis.z *= inv_a;

    const double c = cos(omega);
    const double s = sin(omega);
    const double one_minus_c = 1.0 - c;
    const double k_dot_v = v_dot(axis, pt_src);
    mpaso_vec3 k_cross_v = v_cross(axis, pt_src);

    mpaso_vec3 rot = {
        pt_src.x * c + k_cross_v.x * s + axis.x * k_dot_v * one_minus_c,
        pt_src.y * c + k_cross_v.y * s + axis.y * k_dot_v * one_minus_c,
        pt_src.z * c + k_cross_v.z * s + axis.z * k_dot_v * one_minus_c
    };
    // normalise then scale to r_new (mirrors CPU's pt_dst.Normalize(); pt_dst.scale(r_new))
    *pt_dst = v_normalize_scale(rot, r_new);
    *r_dst  = r_new;
    return true;
}

// ---------- phys_to_cell on device ----------
//
// Exact mirror of MPASOGrid::phys_to_cell (Grid.C:1257-1411) specialised to
// the fast-path case (a valid `from_cell` hint is always available during
// RK4 substepping). Bounding-box check is omitted: outside points are caught
// by the "inside cell" sign test below.
//
// Returns true on success and writes:
//   out_cell_idx            : flat cell index in [0, n_cells)
//   out_v_level             : vertical layer in [0, n_vert_levels - 1)
//   out_n_vertices          : number of vertices for this cell
//   out_omegas[0..nVert-1]  : Wachspress horizontal weights
//   out_omegas[nVert]       : vertical lerp weight (upper layer weight)
//   out_vert_ids[0..nVert-1]: local vertex indices for this cell
__device__ inline bool phys_to_cell_dev(
    const MPASODeviceField& f,
    mpaso_vec3 phys,
    int from_cell_idx,
    int ts_idx,
    int* out_cell_idx,
    int* out_v_level,
    int* out_n_vertices,
    double* out_omegas /* size DEV_MAX_EDGES + 1 */,
    int* out_vert_ids  /* size DEV_MAX_EDGES */)
{
    const int n_cells      = f.n_cells;
    const int n_max_edges  = f.n_max_edges;
    const int n_vert_lev   = f.n_vert_levels;
    const int n_local_vert = f.n_local_vertices;

    // --- project onto sphere surface ---
    const double r_mag = v_mag(phys);
    if (r_mag > f.earth_radius || r_mag == 0.0) return false;
    if (!in_bbox_dev(f, phys)) return false;  // match CPU MPASOGrid::phys_to_cell
    const mpaso_vec3 phys_ontop = v_normalize_scale(phys, f.earth_radius);

    // --- nearest-cell walk via cellsOnCell neighbours ---
    if (from_cell_idx < 0 || from_cell_idx >= n_cells) return false;

    // Single 1-ring pass, matching the CPU MPASOGrid::phys_to_cell single-pass
    // neighbour walk: terminate the line when a step lands more than one ring
    // away, so the GPU stops where the CPU does.
    int nearest = from_cell_idx;
    mpaso_vec3 d0 = v_sub(phys_ontop, f.d_cell_coord[from_cell_idx]);
    double min_dist2 = v_dot(d0, d0);
    const int from_nvoc = f.d_num_vertices_on_cell[from_cell_idx];
    const int from_off  = from_cell_idx * n_max_edges;
    for (int i = 0; i < from_nvoc; ++i) {
        int nb = f.d_cells_on_cell[from_off + i];
        if (nb < 0) continue;
        mpaso_vec3 d = v_sub(phys_ontop, f.d_cell_coord[nb]);
        double dd = v_dot(d, d);
        if (dd < min_dist2) { min_dist2 = dd; nearest = nb; }
    }
    if (nearest < 0 || nearest >= n_cells) return false;

    const int curr_nv  = f.d_num_vertices_on_cell[nearest];
    const int nvoc_off = nearest * n_max_edges;
    if (curr_nv <= 0 || curr_nv > DEV_MAX_EDGES) return false;

    // --- fetch cell vertex coords and verify phys is inside ---
    mpaso_vec3 vertices[DEV_MAX_EDGES];
    int vertexidx[DEV_MAX_EDGES];
    for (int i = 0; i < curr_nv; ++i) {
        int vi = f.d_vertices_on_cell[nvoc_off + i];
        if (vi < 0 || vi >= n_local_vert) return false;
        vertexidx[i] = vi;
        vertices[i] = f.d_vertex_coord[vi];
    }
    bool dot_sign_pos = false;
    for (int i = 0; i < curr_nv; ++i) {
        mpaso_vec3 normal = v_cross(vertices[i], vertices[(i + 1) % curr_nv]);
        double d = v_dot(normal, phys);
        bool pos = (d > 0.0);
        if (i == 0) dot_sign_pos = pos;
        else if (pos != dot_sign_pos) return false;
    }

    // --- Wachspress omegas (Grid.C:1345-1361) ---
    double Aj[DEV_MAX_EDGES];
    double Bk[DEV_MAX_EDGES];
    for (int j = 0; j < curr_nv; ++j)
        Aj[j] = triangle_area_dev(vertices[j], vertices[(j + 1) % curr_nv], phys_ontop);
    for (int k = 0; k < curr_nv; ++k)
        Bk[k] = triangle_area_dev(vertices[k], vertices[(k + 1) % curr_nv], vertices[(k + 2) % curr_nv]);

    double omegas[DEV_MAX_EDGES];
    double omega_sum = 0.0;
    for (int k = 0; k < curr_nv; ++k) {
        double w = Bk[k];
        for (int i = 2; i < curr_nv; ++i)
            w *= Aj[(k + i) % curr_nv];
        omegas[(k + 1) % curr_nv] = w;
        omega_sum += w;
    }
    if (!(omega_sum > 1e-30) && !(omega_sum < -1e-30)) return false;
    for (int i = 0; i < curr_nv; ++i)
        omegas[i] /= omega_sum;

    // --- blend vertex zTop columns into a per-cell currztop[vLevel] ---
    // Single-timestep path for now: use ts_idx directly (no time blend).
    if (n_vert_lev > DEV_MAX_VERT_LEV) return false;
    double currztop[DEV_MAX_VERT_LEV];
    for (int vl = 0; vl < n_vert_lev; ++vl) currztop[vl] = 0.0;

    const double* ztop_ts = f.d_vertex_ztop
        + static_cast<size_t>(ts_idx) * n_local_vert * n_vert_lev;
    for (int i = 0; i < curr_nv; ++i) {
        int lv_off = vertexidx[i] * n_vert_lev;
        for (int vl = 0; vl < n_vert_lev; ++vl)
            currztop[vl] += ztop_ts[lv_off + vl] * omegas[i];
    }

    // --- vertical layer lookup ---
    // Mirror CPU's upper_bound(greater<>) walk on currztop.
    const double depth = r_mag - v_mag(f.d_cell_coord[nearest]);
    int v_level = -1;
    for (int vl = 0; vl < n_vert_lev; ++vl) {
        if (currztop[vl] < depth) { v_level = (vl > 0) ? vl - 1 : 0; break; }
    }
    if (v_level == -1) return false;
    if (v_level >= n_vert_lev - 1) return false;

    const double dz = currztop[v_level + 1] - currztop[v_level];
    const double vert_w = (dz != 0.0) ? (currztop[v_level + 1] - depth) / dz : 0.5;

    // --- write outputs ---
    *out_cell_idx   = nearest;
    *out_v_level    = v_level;
    *out_n_vertices = curr_nv;
    for (int i = 0; i < curr_nv; ++i) {
        out_omegas[i]   = omegas[i];
        out_vert_ids[i] = vertexidx[i];
    }
    out_omegas[curr_nv] = vert_w;
    return true;
}

// ---------- horizontal + vertical velocity sampling ----------
//
// Mirrors CVectorField::at_phys(VECTOR4) + MPASOGrid::interpolate().
// The caller has already run phys_to_cell_dev to fill vert_ids, omegas,
// nV and v_level.
//
// Output: vec3 whose (x,y,z) is the horizontal velocity in xyz and ...
// (On CPU, VECTOR4[3] is the vertical scalar. We return it via out_vert.)
__device__ inline bool sample_velocity_dev(
    const MPASODeviceField& f,
    int ts_idx,
    const int* vert_ids,
    int nV,
    int v_level,
    const double* omegas,       // size nV + 1; [nV] is vertical weight
    mpaso_vec3* out_horizontal,
    double* out_vertical)
{
    const int n_vert_nodes = f.n_local_vertices * f.n_vert_levels;
    const size_t ts_off    = static_cast<size_t>(ts_idx) * n_vert_nodes;
    const mpaso_vec3* h_base = f.d_vertex_velocity      + ts_off;
    const mpaso_vec3* v_base = f.d_vertex_vert_velocity + ts_off;

    // Upper-layer weighted sum and lower-layer weighted sum, mirroring
    // MPASOGrid::interpolate(VECTOR3*, int, double*).
    mpaso_vec3 upper_h = {0.0, 0.0, 0.0};
    mpaso_vec3 lower_h = {0.0, 0.0, 0.0};
    double     upper_v = 0.0;
    double     lower_v = 0.0;

    for (int i = 0; i < nV; ++i) {
        int vi         = vert_ids[i];
        int upper_node = vi * f.n_vert_levels + v_level;
        int lower_node = upper_node + 1;
        double w = omegas[i];

        mpaso_vec3 uh = h_base[upper_node];
        mpaso_vec3 lh = h_base[lower_node];
        upper_h.x += w * uh.x; upper_h.y += w * uh.y; upper_h.z += w * uh.z;
        lower_h.x += w * lh.x; lower_h.y += w * lh.y; lower_h.z += w * lh.z;

        // CPU stores vertical velocity as VECTOR3 with only [0] populated.
        upper_v += w * v_base[upper_node].x;
        lower_v += w * v_base[lower_node].x;
    }

    const double vw = omegas[nV];          // upper-layer weight
    const double lw = 1.0 - vw;            // lower-layer weight

    out_horizontal->x = vw * upper_h.x + lw * lower_h.x;
    out_horizontal->y = vw * upper_h.y + lw * lower_h.y;
    out_horizontal->z = vw * upper_h.z + lw * lower_h.z;
    *out_vertical     = vw * upper_v    + lw * lower_v;
    return true;
}

// ---------- time-window resolution ----------
//
// Given a particle time t_part and the uploaded window's real-time stamps,
// find (lowT, highT, ratio) such that:
//   timestamps[lowT] <= t_part <= timestamps[highT]
//   ratio = (t_part - timestamps[lowT]) / (timestamps[highT] - timestamps[lowT])
// Mirrors Solution::GetValue(UNSTEADY)'s window resolution + clamping.
// Returns false if t_part is outside the window or there is < 2 loaded ts.
__device__ inline bool resolve_time_window_dev(
    const MPASODeviceField& f, double t_part,
    int* out_lowT, int* out_highT, double* out_ratio)
{
    const int n_ts = f.n_timesteps_loaded;
    if (n_ts <= 0) return false;
    if (n_ts == 1) {
        *out_lowT  = 0; *out_highT = 0; *out_ratio = 0.0;
        return true;
    }
    const double t0 = f.d_timestamps[0];
    const double tN = f.d_timestamps[n_ts - 1];
    if (t_part < t0 || t_part > tN) return false;

    int lowT = 0;
    for (int i = 0; i + 1 < n_ts; ++i) {
        if (t_part >= f.d_timestamps[i]) lowT = i;
        else break;
    }
    int highT = lowT + 1;
    if (highT >= n_ts) { highT = lowT; *out_ratio = 0.0; }
    else {
        double span = f.d_timestamps[highT] - f.d_timestamps[lowT];
        *out_ratio = (span > 0.0) ? (t_part - f.d_timestamps[lowT]) / span : 0.0;
    }
    *out_lowT  = lowT;
    *out_highT = highT;
    return true;
}

// ---------- phys_to_cell with ztop temporal blend ----------
//
// Mirrors phys_to_cell_dev but blends ztop between two timesteps using
// (lowT, highT, ratio). Wachspress omegas are time-independent.
__device__ inline bool phys_to_cell_dev_blend(
    const MPASODeviceField& f,
    mpaso_vec3 phys,
    int from_cell_idx,
    int lowT, int highT, double ratio,
    int* out_cell_idx,
    int* out_v_level,
    int* out_n_vertices,
    double* out_omegas,
    int* out_vert_ids)
{
    const int n_cells      = f.n_cells;
    const int n_max_edges  = f.n_max_edges;
    const int n_vert_lev   = f.n_vert_levels;
    const int n_local_vert = f.n_local_vertices;

    const double r_mag = v_mag(phys);
    if (r_mag > f.earth_radius || r_mag == 0.0) return false;
    if (!in_bbox_dev(f, phys)) return false;  // match CPU MPASOGrid::phys_to_cell
    const mpaso_vec3 phys_ontop = v_normalize_scale(phys, f.earth_radius);

    if (from_cell_idx < 0 || from_cell_idx >= n_cells) return false;

    // Single 1-ring pass, matching the CPU MPASOGrid::phys_to_cell single-pass
    // neighbour walk (terminate when a step lands more than one ring away).
    int nearest = from_cell_idx;
    mpaso_vec3 d0 = v_sub(phys_ontop, f.d_cell_coord[from_cell_idx]);
    double min_dist2 = v_dot(d0, d0);

    const int from_nvoc = f.d_num_vertices_on_cell[from_cell_idx];
    const int from_off  = from_cell_idx * n_max_edges;
    for (int i = 0; i < from_nvoc; ++i) {
        int nb = f.d_cells_on_cell[from_off + i];
        if (nb < 0) continue;
        mpaso_vec3 d = v_sub(phys_ontop, f.d_cell_coord[nb]);
        double dd = v_dot(d, d);
        if (dd < min_dist2) { min_dist2 = dd; nearest = nb; }
    }
    if (nearest < 0 || nearest >= n_cells) return false;

    const int curr_nv  = f.d_num_vertices_on_cell[nearest];
    const int nvoc_off = nearest * n_max_edges;
    if (curr_nv <= 0 || curr_nv > DEV_MAX_EDGES) return false;

    mpaso_vec3 vertices[DEV_MAX_EDGES];
    int vertexidx[DEV_MAX_EDGES];
    for (int i = 0; i < curr_nv; ++i) {
        int vi = f.d_vertices_on_cell[nvoc_off + i];
        if (vi < 0 || vi >= n_local_vert) return false;
        vertexidx[i] = vi;
        vertices[i] = f.d_vertex_coord[vi];
    }
    bool dot_sign_pos = false;
    for (int i = 0; i < curr_nv; ++i) {
        mpaso_vec3 normal = v_cross(vertices[i], vertices[(i + 1) % curr_nv]);
        double d = v_dot(normal, phys);
        bool pos = (d > 0.0);
        if (i == 0) dot_sign_pos = pos;
        else if (pos != dot_sign_pos) return false;
    }

    double Aj[DEV_MAX_EDGES];
    double Bk[DEV_MAX_EDGES];
    for (int j = 0; j < curr_nv; ++j)
        Aj[j] = triangle_area_dev(vertices[j], vertices[(j + 1) % curr_nv], phys_ontop);
    for (int k = 0; k < curr_nv; ++k)
        Bk[k] = triangle_area_dev(vertices[k], vertices[(k + 1) % curr_nv], vertices[(k + 2) % curr_nv]);

    double omegas[DEV_MAX_EDGES];
    double omega_sum = 0.0;
    for (int k = 0; k < curr_nv; ++k) {
        double w = Bk[k];
        for (int i = 2; i < curr_nv; ++i)
            w *= Aj[(k + i) % curr_nv];
        omegas[(k + 1) % curr_nv] = w;
        omega_sum += w;
    }
    if (!(omega_sum > 1e-30) && !(omega_sum < -1e-30)) return false;
    for (int i = 0; i < curr_nv; ++i) omegas[i] /= omega_sum;

    if (n_vert_lev > DEV_MAX_VERT_LEV) return false;
    double currztop[DEV_MAX_VERT_LEV];
    for (int vl = 0; vl < n_vert_lev; ++vl) currztop[vl] = 0.0;

    const size_t ts_stride = (size_t)n_local_vert * n_vert_lev;
    const double* ztop_lo = f.d_vertex_ztop + (size_t)lowT  * ts_stride;
    const double* ztop_hi = f.d_vertex_ztop + (size_t)highT * ts_stride;
    const bool    blend   = (highT != lowT);
    for (int i = 0; i < curr_nv; ++i) {
        int lv_off = vertexidx[i] * n_vert_lev;
        for (int vl = 0; vl < n_vert_lev; ++vl) {
            double z = ztop_lo[lv_off + vl];
            if (blend) z += (ztop_hi[lv_off + vl] - z) * ratio;
            currztop[vl] += z * omegas[i];
        }
    }

    const double depth = r_mag - v_mag(f.d_cell_coord[nearest]);
    int v_level = -1;
    for (int vl = 0; vl < n_vert_lev; ++vl) {
        if (currztop[vl] < depth) { v_level = (vl > 0) ? vl - 1 : 0; break; }
    }
    if (v_level == -1) return false;
    if (v_level >= n_vert_lev - 1) return false;

    const double dz = currztop[v_level + 1] - currztop[v_level];
    const double vert_w = (dz != 0.0) ? (currztop[v_level + 1] - depth) / dz : 0.5;

    *out_cell_idx   = nearest;
    *out_v_level    = v_level;
    *out_n_vertices = curr_nv;
    for (int i = 0; i < curr_nv; ++i) {
        out_omegas[i]   = omegas[i];
        out_vert_ids[i] = vertexidx[i];
    }
    out_omegas[curr_nv] = vert_w;
    return true;
}

// ---------- sample velocity with temporal blend ----------
//
// Mirrors sample_velocity_dev but blends each vertex's (h,v) between lowT/highT
// before the Wachspress sum, matching CPU's Solution::GetValue per-vertex lerp.
__device__ inline bool sample_velocity_dev_blend(
    const MPASODeviceField& f,
    int lowT, int highT, double ratio,
    const int* vert_ids,
    int nV,
    int v_level,
    const double* omegas,
    mpaso_vec3* out_horizontal,
    double* out_vertical)
{
    const int n_vert_nodes = f.n_local_vertices * f.n_vert_levels;
    const size_t ts_off_lo = (size_t)lowT  * n_vert_nodes;
    const size_t ts_off_hi = (size_t)highT * n_vert_nodes;
    const mpaso_vec3* h_lo = f.d_vertex_velocity      + ts_off_lo;
    const mpaso_vec3* h_hi = f.d_vertex_velocity      + ts_off_hi;
    const mpaso_vec3* v_lo = f.d_vertex_vert_velocity + ts_off_lo;
    const mpaso_vec3* v_hi = f.d_vertex_vert_velocity + ts_off_hi;
    const bool blend = (highT != lowT);

    mpaso_vec3 upper_h = {0.0, 0.0, 0.0};
    mpaso_vec3 lower_h = {0.0, 0.0, 0.0};
    double     upper_v = 0.0;
    double     lower_v = 0.0;

    for (int i = 0; i < nV; ++i) {
        int vi         = vert_ids[i];
        int upper_node = vi * f.n_vert_levels + v_level;
        int lower_node = upper_node + 1;
        double w = omegas[i];

        mpaso_vec3 uh = h_lo[upper_node];
        mpaso_vec3 lh = h_lo[lower_node];
        if (blend) {
            mpaso_vec3 uh_hi = h_hi[upper_node];
            mpaso_vec3 lh_hi = h_hi[lower_node];
            uh.x += (uh_hi.x - uh.x) * ratio;
            uh.y += (uh_hi.y - uh.y) * ratio;
            uh.z += (uh_hi.z - uh.z) * ratio;
            lh.x += (lh_hi.x - lh.x) * ratio;
            lh.y += (lh_hi.y - lh.y) * ratio;
            lh.z += (lh_hi.z - lh.z) * ratio;
        }
        upper_h.x += w * uh.x; upper_h.y += w * uh.y; upper_h.z += w * uh.z;
        lower_h.x += w * lh.x; lower_h.y += w * lh.y; lower_h.z += w * lh.z;

        double uv = v_lo[upper_node].x;
        double lv = v_lo[lower_node].x;
        if (blend) {
            uv += (v_hi[upper_node].x - uv) * ratio;
            lv += (v_hi[lower_node].x - lv) * ratio;
        }
        upper_v += w * uv;
        lower_v += w * lv;
    }

    const double vw = omegas[nV];
    const double lw = 1.0 - vw;

    out_horizontal->x = vw * upper_h.x + lw * lower_h.x;
    out_horizontal->y = vw * upper_h.y + lw * lower_h.y;
    out_horizontal->z = vw * upper_h.z + lw * lower_h.z;
    *out_vertical     = vw * upper_v    + lw * lower_v;
    return true;
}

} // namespace mpaso_gpu_dev

#endif // MPASO_VELOCITY_CUH
