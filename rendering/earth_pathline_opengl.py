#!/usr/bin/env python3
"""
3D animated pathline visualization using PyVista.

Renders MPAS-O ocean pathlines as animated head+tail particles on a 3D globe.
The globe is a texture-mapped sphere: ocean cells are coloured by bottomDepth
(light blue = shallow → dark blue = deep) and land is filled with an earth tone.
Particle heads and tails are coloured by their instantaneous depth below the
ocean surface so depth changes are visible throughout the animation.

Prerequisites:
    pip install pyvista --user   (netCDF4 is already on NERSC)

    If pathlines.bin does not yet exist, run first:
        python scripts/plot/read_traces.py 256 pathlines/nersc_highres_256_cpu_10k

Usage:
    # Interactive preview — rotate/zoom with mouse, close to print camera pos:
    python rendering/render_3d_animation.py --interactive

    # Static mid-point PNG (fast sanity check):
    python rendering/render_3d_animation.py ... rendering/preview.png

    # Full animation:
    python rendering/render_3d_animation.py ... rendering/pathlines_3d.mp4

All positional args are optional; see defaults below.
"""

import struct
import os
import argparse

import numpy as np
from datetime import datetime, timedelta
import netCDF4 as nc_mod

# ---------------------------------------------------------------------------
# Tunable constants
# ---------------------------------------------------------------------------

MESH_SUBSAMPLE     = 1       # keep 1-in-N particles for animation speed tests
PARTICLE_SUBSAMPLE = 1       # (ocean texture always uses ALL cells)
TAIL_LEN           = 500     # trailing timesteps shown per particle
FRAME_SKIP         = 20      # timesteps per animation frame  (→ ~53 s @ 24 fps)
FPS                = 24
WINDOW_SIZE        = (1920, 1080)

HEAD_POINT_SIZE    = 4       # sphere point size for particle heads
TAIL_LINE_WIDTH    = 1.5
OCEAN_DEPTH_CLIM   = 5000    # depth (m) at which ocean colour saturates to darkest blue
PARTICLE_DEPTH_CLIM = None   # depth (m) for particle colourmap upper bound;
                             # None = auto-compute from the 95th-percentile depth

# CAMERA_POSITION = [(x,y,z), (fx,fy,fz), (ux,uy,uz)]
# Pacific view: lon≈220°, equatorial, pulled back
# CAMERA_POSITION = [(-22_000_000.0, -18_000_000.0, 3_000_000.0), (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)]
# Atlantic view: lon≈310° (50°W), shows mid-Atlantic + American east coast on right
CAMERA_POSITION = [(18_000_000.0, -21_000_000.0, 2_000_000.0), (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)]

EARTH_RADIUS = 6_371_229.0   # metres (MPAS-O reference)

def read_sim_time_info(nc_path, paraflow_yaml=None):
    """
    Read simulation start time from MPAS NetCDF (xtime.orig) and the
    particle output interval from the ParaFlow yaml config.

    The MPAS config_dt is the ocean model integration step (e.g. 5 min) —
    NOT the particle output interval.  The correct interval is:
        dt (seconds) × record_interval  (from the ParaFlow yaml)
    e.g. dt=60 s, record_interval=60 → 1 hour per stored particle step.
    """
    with nc_mod.Dataset(nc_path) as ds:
        xtime_chars = ds.variables["xtime.orig"][0].data
        xtime_str   = b"".join(xtime_chars).decode().strip("\x00")
    sim_start = datetime.strptime(xtime_str, "%Y-%m-%d_%H:%M:%S")

    dt_minutes = None
    if paraflow_yaml and os.path.exists(paraflow_yaml):
        import re
        with open(paraflow_yaml) as f:
            txt = f.read()
        m_dt  = re.search(r'^\s*dt\s*:\s*([\d.]+)', txt, re.MULTILINE)
        m_rec = re.search(r'^\s*record_interval\s*:\s*(\d+)', txt, re.MULTILINE)
        if m_dt and m_rec:
            dt_sec     = float(m_dt.group(1))
            rec_int    = int(m_rec.group(1))
            dt_minutes = dt_sec * rec_int / 60.0

    if dt_minutes is None:
        # Fallback: guess from nPts — if user says ~3 years default to 60 min
        dt_minutes = 60.0
        print("  Warning: ParaFlow yaml not found; defaulting dt_minutes=60 (1 hr/step).")

    return sim_start, dt_minutes, xtime_str


def fmt_mpas_time(sim_start, dt_minutes, t):
    """Format a timestep as a MPAS-style timestamp (YYYY-MM-DD_HH:MM:SS)."""
    dt = sim_start + timedelta(minutes=t * dt_minutes)
    return f"{dt.year:04d}-{dt.month:02d}-{dt.day:02d}_{dt.hour:02d}:{dt.minute:02d}:{dt.second:02d}"

# CHANGE THESE to your own files (or pass them as command-line args) before
# running -- the defaults point at specific paths that won't exist for you.
DEFAULT_MESH_NC = (
    "/global/cfs/projectdirs/m4259/milena/forParticleTracing/"
    "mpaso.RRSwISC6to18E3r5.rstFromG-chrysalis.20240603.nc"
)
DEFAULT_PATHLINES = "pathlines/nersc_highres_256_cpu_10k/pathlines.bin"


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def read_pathlines(path):
    """Read pathlines.bin → list of (pid, Nx3 float64 array)."""
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"{path} not found.\n"
            "Run:  python scripts/plot/read_traces.py 256 pathlines/nersc_highres_256_cpu_10k"
        )
    print(f"Reading pathlines: {path}")
    with open(path, "rb") as f:
        data = f.read()
    offset = 0
    n_pl = struct.unpack_from("<i", data, offset)[0]
    offset += 4
    pathlines = []
    for _ in range(n_pl):
        pid, nPts = struct.unpack_from("<ii", data, offset)
        offset += 8
        coords = (
            np.frombuffer(data, dtype="<f8", count=nPts * 3, offset=offset)
            .reshape(nPts, 3)
            .copy()
        )
        offset += nPts * 24
        pathlines.append((pid, coords))
    print(f"  {len(pathlines):,} pathlines loaded")
    return pathlines


def build_ocean_globe(nc_path, depth_clim, earth_radius):
    """
    Rasterise MPAS bottomDepth onto a 2400×1200 equirectangular texture and
    map it onto a sphere.

    Ocean cells → Blues colormap by depth (light = shallow, dark = deep).
    Land (unvisited pixels) → muted olive-green earth tone.

    Returns (sphere_mesh, texture).
    """
    import netCDF4 as nc_mod
    import pyvista as pv
    from matplotlib import colormaps

    LAT_RES, LON_RES = 1200, 2400

    print(f"Building ocean texture from {nc_path} ...")
    with nc_mod.Dataset(nc_path) as ds:
        lat_rad = ds.variables["latCell"][:].data
        lon_rad = ds.variables["lonCell"][:].data
        depth   = ds.variables["bottomDepth"][:].data.astype(np.float32)

    lat_deg = np.degrees(lat_rad)        # -90..90
    lon_deg = np.degrees(lon_rad) % 360  #   0..360

    # Image layout we will use:
    #   row 0 = south pole (lat=-90),  row increases northward
    #   col 0 = lon=0°,                col increases eastward
    row = np.clip(((lat_deg + 90.0) / 180.0 * LAT_RES).astype(int), 0, LAT_RES - 1)
    col = np.clip((lon_deg / 360.0 * LON_RES).astype(int), 0, LON_RES - 1)

    depth_grid = np.zeros((LAT_RES, LON_RES), dtype=np.float32)
    depth_grid[row, col] = depth

    # Gap-fill: MPAS cells don't cover every pixel at this resolution.
    # Expand ocean into neighboring empty pixels with 8 passes of 4-neighbour
    # max-fill (~8 × 17 km ≈ 135 km spread, enough for all coastal gaps).
    depth_fill = depth_grid.copy()
    for _ in range(8):
        nb = np.maximum.reduce([
            np.roll(depth_fill,  1, axis=0),   # from south
            np.roll(depth_fill, -1, axis=0),   # from north
            np.roll(depth_fill,  1, axis=1),   # from west
            np.roll(depth_fill, -1, axis=1),   # from east
        ])
        depth_fill = np.where(depth_fill > 0, depth_fill, nb)

    # Build RGB texture image
    ocean_mask = depth_fill > 0
    cmap       = colormaps["Blues"]
    land_rgb   = np.array([110, 130, 80], dtype=np.uint8)
    img        = np.full((LAT_RES, LON_RES, 3), land_rgb, dtype=np.uint8)
    depth_norm = np.clip(depth_fill[ocean_mask] / depth_clim, 0, 1)
    img[ocean_mask] = (cmap(depth_norm)[:, :3] * 255).astype(np.uint8)

    pre  = (depth_grid > 0).mean() * 100
    post = ocean_mask.mean() * 100
    print(f"  Ocean coverage: {pre:.1f}% raw → {post:.1f}% after gap-fill")

    # Build sphere and assign UV texture coordinates directly from each
    # vertex's geographic position — no dependence on VTK auto-mapping.
    #
    # Coordinate convention matching the image layout above:
    #   u = lon / 360          (u=0 at lon=0°,  u=1 at lon=360°)
    #   v = (lat + 90) / 180   (v=0 at south pole, v=1 at north pole)
    #
    # PyVista passes the numpy array through PIL before handing to VTK.
    # PIL stores row 0 at the top of the image, but VTK/OpenGL expects row 0
    # at the bottom (v=0).  PyVista flips the array during that conversion,
    # so our row-0=south-pole image ends up at OpenGL bottom = v=0 ✓.
    sphere = pv.Sphere(
        radius=earth_radius * 0.997,
        theta_resolution=360,
        phi_resolution=180,
    )
    sph_pts = sphere.points                                              # (N, 3)
    r_pts   = np.sqrt(np.einsum("ij,ij->i", sph_pts, sph_pts))
    lat_sph = np.degrees(np.arcsin(np.clip(sph_pts[:, 2] / r_pts, -1, 1)))
    lon_sph = np.degrees(np.arctan2(sph_pts[:, 1], sph_pts[:, 0])) % 360

    u = (lon_sph / 360.0).astype(np.float32)
    v = ((lat_sph + 90.0) / 180.0).astype(np.float32)
    sphere.active_texture_coordinates = np.column_stack([u, v])

    # Seam fix: triangles bridging u≈0 and u≈1 would otherwise interpolate UVs
    # across the entire texture (visible blue streak at lon=0°).  For each such
    # triangle, replace its low-U vertex with a duplicate whose U is shifted
    # up by 1.0 — so interpolation goes from u≈0.997 → u≈1.001 (a tiny step).
    faces      = sphere.faces.reshape(-1, 4).copy()
    tri_u      = u[faces[:, 1:]]
    seam_tris  = np.where((tri_u.max(axis=1) - tri_u.min(axis=1)) > 0.5)[0]
    if seam_tris.size:
        pts_list = list(sphere.points)
        u_list   = list(u)
        v_list   = list(v)
        remap    = {}
        for ti in seam_tris:
            for ci in (1, 2, 3):
                vi = int(faces[ti, ci])
                if u[vi] < 0.5:
                    if vi not in remap:
                        remap[vi] = len(pts_list)
                        pts_list.append(sphere.points[vi])
                        u_list.append(u[vi] + 1.0)
                        v_list.append(v[vi])
                    faces[ti, ci] = remap[vi]
        sphere.points = np.array(pts_list)
        sphere.faces  = faces.ravel()
        sphere.active_texture_coordinates = np.column_stack([u_list, v_list]).astype(np.float32)
        print(f"  Seam fix: duplicated {len(remap)} vertices across {len(seam_tris)} triangles")

    # Pre-flip: PyVista._from_array internally flips the row axis, so we flip
    # first so the two flips cancel and V=0 correctly lands at the south pole.
    texture = pv.Texture(img[::-1])
    # Enable U-direction tiling so duplicated vertices with u>1 sample the
    # correct wrapped texel instead of clamping to the rightmost column.
    texture.SetRepeat(True)
    return sphere, texture


# ---------------------------------------------------------------------------
# Particle depth helper
# ---------------------------------------------------------------------------

def compute_depth(pts):
    """Depth below Earth surface for Nx3 Cartesian points (metres)."""
    return EARTH_RADIUS - np.sqrt(np.einsum("ij,ij->i", pts, pts))


def visible_mask(pts, camera_pos, margin=0.05):
    """
    Return a boolean mask: True if the point is on the camera-facing hemisphere.

    Z-fighting at globe scale makes the depth buffer unreliable for back-side
    particles (sphere back-face and particle are only ~20 km apart in depth
    space out of 24 M km camera distance).  Explicit hemisphere culling is
    more robust: a point is visible when its outward normal has a positive
    dot-product with the camera direction.

    margin > 0 hides points near the horizon to avoid edge artefacts.
    """
    cam_dir   = np.asarray(camera_pos, dtype=np.float64)
    cam_dir   = cam_dir / np.linalg.norm(cam_dir)          # unit vector toward camera
    pt_norms  = np.sqrt(np.einsum("ij,ij->i", pts, pts))
    pt_norms  = np.maximum(pt_norms, 1.0)
    pt_unit   = pts / pt_norms[:, None]
    return (pt_unit @ cam_dir) > margin


# ---------------------------------------------------------------------------
# Animation buffer helpers
# ---------------------------------------------------------------------------

def build_tail_mesh(n_parts, tail_len):
    """
    Pre-allocate a PyVista PolyData for all particle tails.

    Points layout: particle k → rows [k*tail_len : (k+1)*tail_len].
    Per-point 'particle_depth' scalar is updated each frame so the tail
    colour shows actual depth, not trajectory age.

    Returns (mesh, pts_buffer, depth_buffer).
    """
    import pyvista as pv

    n_pts = n_parts * tail_len
    pts   = np.zeros((n_pts, 3), dtype=np.float64)

    counts  = np.full((n_parts, 1), tail_len, dtype=np.int32)
    offsets = (np.arange(n_parts, dtype=np.int32) * tail_len)[:, None]
    indices = offsets + np.arange(tail_len, dtype=np.int32)[None, :]
    conn    = np.hstack([counts, indices]).ravel()

    depth_buf = np.zeros(n_pts, dtype=np.float32)

    mesh              = pv.PolyData()
    mesh.points       = pts
    mesh.lines        = conn
    mesh["particle_depth"] = depth_buf
    return mesh, pts, depth_buf


def fill_tail_buffer(pts_buf, depth_buf, pathlines, t, tail_len, vis_mask):
    """Fill tail positions/depths; collapse back-side tails to origin."""
    window = np.arange(t - tail_len + 1, t + 1, dtype=np.int64)
    for k, (_, coords) in enumerate(pathlines):
        npts = len(coords)
        if not vis_mask[k]:
            # Particle on back hemisphere — collapse inside the globe
            pts_buf  [k * tail_len : (k + 1) * tail_len] = 0.0
            depth_buf[k * tail_len : (k + 1) * tail_len] = 0.0
            continue
        idxs = np.clip(window, 0, npts - 1)
        seg  = coords[idxs]
        pts_buf  [k * tail_len : (k + 1) * tail_len] = seg
        depth_buf[k * tail_len : (k + 1) * tail_len] = compute_depth(seg)


def fill_head_buffer(pts_buf, depth_buf, pathlines, t, vis_mask):
    """Fill head positions/depths; collapse back-side heads to origin."""
    for k, (_, coords) in enumerate(pathlines):
        if not vis_mask[k]:
            pts_buf[k]   = 0.0
            depth_buf[k] = 0.0
            continue
        pos          = coords[min(t, len(coords) - 1)]
        pts_buf[k]   = pos
        depth_buf[k] = EARTH_RADIUS - float(np.sqrt(pos @ pos))


# ---------------------------------------------------------------------------
# ffmpeg stitching (libx264 → libvpx-vp9 fallback)
# ---------------------------------------------------------------------------

def _ffmpeg_stitch(frames_dir, output_path, fps):
    """Stitch PNG frames into a video. Falls back to VP9/WebM if libx264 is absent."""
    import subprocess

    pattern = os.path.join(frames_dir, "frame_%06d.png")
    codecs  = [
        (["-c:v", "libx264", "-preset", "fast", "-crf", "18", "-pix_fmt", "yuv420p"],
         output_path),
        (["-c:v", "libvpx-vp9", "-crf", "30", "-b:v", "0"],
         os.path.splitext(output_path)[0] + ".webm"),
    ]
    for codec_args, out in codecs:
        cmd    = ["ffmpeg", "-y", "-framerate", str(fps), "-i", pattern] + codec_args + [out]
        result = subprocess.run(cmd, capture_output=True, text=True)
        print(f"Stitching ({codec_args[1]}): {out}")
        if result.returncode == 0:
            return out
        print(f"  failed: {result.stderr.splitlines()[-1] if result.stderr else 'unknown'}")
    raise RuntimeError("ffmpeg failed with all codecs. Run: ffmpeg -encoders | grep -E 'libx264|vp9'")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="3D pathline animation — PyVista")
    parser.add_argument("pathlines_bin", nargs="?", default=DEFAULT_PATHLINES)
    parser.add_argument("mesh_nc",       nargs="?", default=DEFAULT_MESH_NC)
    parser.add_argument("output",        nargs="?", default="rendering/pathlines_3d.mp4")
    parser.add_argument("--interactive", action="store_true",
                        help="Open interactive window (requires display / X11 forwarding). "
                             "Rotate/zoom with mouse, close to print camera position.")
    parser.add_argument("--final", action="store_true",
                        help="Render a single summary image with all complete trajectories "
                             "(every particle's full path from start to end).")
    parser.add_argument("--config", default="conf/nersc_highres_cpu.yaml",
                        help="ParaFlow yaml config — used to read dt and record_interval "
                             "for correct real-world timestamps.")
    args = parser.parse_args()

    is_png = args.output.lower().endswith(".png")

    # --- Simulation time ---
    sim_start, dt_minutes, xtime_str = read_sim_time_info(args.mesh_nc, args.config)
    print(f"  Sim start: {xtime_str}  dt={dt_minutes} min")

    # --- Load data ---
    pathlines = read_pathlines(args.pathlines_bin)
    pathlines = pathlines[::PARTICLE_SUBSAMPLE]
    n_parts   = len(pathlines)
    max_npts  = max(len(c) for _, c in pathlines)
    print(f"  {n_parts} particles, {max_npts} max timesteps")

    # Auto-compute particle depth range from all trajectory points
    global PARTICLE_DEPTH_CLIM
    if PARTICLE_DEPTH_CLIM is None:
        all_depths = np.concatenate([
            compute_depth(coords) for _, coords in pathlines
        ])
        PARTICLE_DEPTH_CLIM = float(np.percentile(all_depths, 95))
        print(f"  Particle depth range: {all_depths.min():.1f}–{all_depths.max():.1f} m  "
              f"(clim set to 95th-pct = {PARTICLE_DEPTH_CLIM:.1f} m)")

    globe, texture = build_ocean_globe(args.mesh_nc, OCEAN_DEPTH_CLIM, EARTH_RADIUS)

    # --- Plotter setup ---
    import pyvista as pv

    if not args.interactive:
        try:
            pv.start_xvfb()
        except Exception:
            pass

    plotter = pv.Plotter(off_screen=not args.interactive, window_size=WINDOW_SIZE)
    plotter.set_background("black")

    # Textured globe (ocean depth + land fill)
    plotter.add_mesh(globe, texture=texture, smooth_shading=True, lighting=True)

    # Scalar bar for ocean depth (static reference — invisible ghost mesh drives the legend)
    dummy = pv.PolyData(np.zeros((2, 3)))
    dummy["depth"] = np.array([0.0, float(OCEAN_DEPTH_CLIM)])
    plotter.add_mesh(dummy, scalars="depth", cmap="Blues",
                     clim=[0, OCEAN_DEPTH_CLIM], opacity=0,
                     show_scalar_bar=True,
                     scalar_bar_args={"title": "Ocean\nDepth (m)", "color": "white",
                                      "fmt": "%.0f", "vertical": True,
                                      "width": 0.035, "height": 0.22,
                                      "position_x": 0.09, "position_y": 0.05})

    # Particle tails
    tail_mesh, tail_pts_buf, tail_depth_buf = build_tail_mesh(n_parts, TAIL_LEN)

    # Camera — resolve now so visibility culling can use the camera position
    if CAMERA_POSITION is not None:
        plotter.camera_position = CAMERA_POSITION
    else:
        cam_dist = EARTH_RADIUS * 2.8
        plotter.camera_position = [
            (cam_dist * 0.7, cam_dist * 0.5, cam_dist * 0.55),
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 1.0),
        ]
        plotter.camera.zoom(1.1)

    cam_pos = np.array(plotter.camera_position[0], dtype=np.float64)

    # Pre-compute per-particle visibility (fixed camera → mask is constant)
    t0_heads = np.array([coords[0] for _, coords in pathlines])
    vis = visible_mask(t0_heads, cam_pos)
    print(f"  Visible particles (front hemisphere): {vis.sum():,} / {n_parts:,}")

    if not (is_png and args.final):
        fill_tail_buffer(tail_pts_buf, tail_depth_buf, pathlines, 0, TAIL_LEN, vis)
        tail_mesh.points            = tail_pts_buf
        tail_mesh["particle_depth"] = tail_depth_buf
        plotter.add_mesh(tail_mesh, scalars="particle_depth", cmap="plasma_r",
                         clim=[0, PARTICLE_DEPTH_CLIM],
                         line_width=TAIL_LINE_WIDTH, lighting=False, show_scalar_bar=True,
                         scalar_bar_args={"title": "Particle\nDepth (m)", "color": "white",
                                          "fmt": "%.0f", "vertical": True,
                                          "width": 0.035, "height": 0.22,
                                          "position_x": 0.02, "position_y": 0.05})

        # Particle heads
        head_pts_buf   = np.zeros((n_parts, 3), dtype=np.float64)
        head_depth_buf = np.zeros(n_parts, dtype=np.float32)
        fill_head_buffer(head_pts_buf, head_depth_buf, pathlines, 0, vis)
        head_mesh = pv.PolyData(head_pts_buf.copy())
        head_mesh["particle_depth"] = head_depth_buf
        plotter.add_mesh(head_mesh, scalars="particle_depth", cmap="plasma_r",
                         clim=[0, PARTICLE_DEPTH_CLIM],
                         point_size=HEAD_POINT_SIZE, render_points_as_spheres=True,
                         lighting=False, show_scalar_bar=False)

    # --- Interactive mode ---
    if args.interactive:
        print("\nInteractive window open — use mouse to rotate/zoom/pan.")
        print("Close the window to print the camera position.\n")
        plotter.show()
        pos = list(plotter.camera_position)
        print("\n--- Copy this into the script as CAMERA_POSITION ---")
        print(f"CAMERA_POSITION = {pos}")
        return

    # --- Final summary image: all complete trajectories ---
    if is_png and args.final:
        sub_pathlines = pathlines[::5]
        print(f"  --final subsample: 1-in-5 particles → {len(sub_pathlines)}/{len(pathlines)}")

        end_heads = np.array([coords[-1] for _, coords in sub_pathlines])
        vis_final = visible_mask(end_heads, cam_pos)
        print(f"  Visible particles at end: {vis_final.sum():,} / {len(sub_pathlines):,}")

        # Build one PolyData with every particle's full path (front hemisphere only)
        all_pts, all_depths = [], []
        counts, offsets = [], []
        cumulative = 0
        for k, (_, coords) in enumerate(sub_pathlines):
            if not vis_final[k]:
                continue
            depths = compute_depth(coords)
            all_pts.append(coords)
            all_depths.append(depths)
            n = len(coords)
            counts.append(n)
            offsets.append(cumulative)
            cumulative += n

        pts_cat   = np.concatenate(all_pts,   axis=0).astype(np.float64)
        depth_cat = np.concatenate(all_depths, axis=0).astype(np.float32)

        lines_list = []
        for n, off in zip(counts, offsets):
            lines_list.append(n)
            lines_list.extend(range(off, off + n))
        lines_arr = np.array(lines_list, dtype=np.int32)

        full_mesh = pv.PolyData()
        full_mesh.points = pts_cat
        full_mesh.lines  = lines_arr
        full_mesh["particle_depth"] = depth_cat

        plotter.add_mesh(full_mesh, scalars="particle_depth", cmap="plasma_r",
                         clim=[0, PARTICLE_DEPTH_CLIM],
                         line_width=TAIL_LINE_WIDTH, lighting=False, show_scalar_bar=True,
                         scalar_bar_args={"title": "Particle\nDepth (m)", "color": "white",
                                          "fmt": "%.0f", "vertical": True,
                                          "width": 0.035, "height": 0.22,
                                          "position_x": 0.02, "position_y": 0.05})

        end_label = fmt_mpas_time(sim_start, dt_minutes, max_npts - 1)
        plotter.add_text(f"All pathlines  {xtime_str} → {end_label}",
                         position="upper_left", font_size=12, color="white")
        os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
        plotter.screenshot(args.output, transparent_background=False)
        print(f"Final summary image saved: {args.output}")
        plotter.close()
        return

    # --- Static PNG preview (mid-point snapshot) ---
    if is_png:
        mid     = max_npts // 2
        vis_mid = visible_mask(
            np.array([coords[min(mid, len(coords)-1)] for _, coords in pathlines]),
            cam_pos)
        fill_tail_buffer(tail_pts_buf, tail_depth_buf, pathlines, mid, TAIL_LEN, vis_mid)
        tail_mesh.points            = tail_pts_buf
        tail_mesh["particle_depth"] = tail_depth_buf
        fill_head_buffer(head_pts_buf, head_depth_buf, pathlines, mid, vis_mid)
        head_mesh.points            = head_pts_buf.copy()
        head_mesh["particle_depth"] = head_depth_buf
        plotter.add_text(fmt_mpas_time(sim_start, dt_minutes, mid),
                         position="upper_left", font_size=14, color="white", name="ts_label")
        os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
        plotter.screenshot(args.output, transparent_background=False)
        print(f"Preview saved: {args.output}")
        plotter.close()
        return

    # --- MP4 animation ---
    frames     = list(range(0, max_npts, FRAME_SKIP))
    out_dir    = os.path.dirname(args.output) or "."
    stem       = os.path.splitext(os.path.basename(args.output))[0]
    frames_dir = os.path.join(out_dir, f"{stem}_frames")
    os.makedirs(frames_dir, exist_ok=True)

    print(f"Rendering {len(frames)} frames at {FPS} fps → {args.output}")
    print(f"  Estimated duration: {len(frames) / FPS:.1f} s")

    for fi, t in enumerate(frames):
        cur_heads = np.array([coords[min(t, len(coords)-1)] for _, coords in pathlines])
        vis_t     = visible_mask(cur_heads, cam_pos)

        fill_tail_buffer(tail_pts_buf, tail_depth_buf, pathlines, t, TAIL_LEN, vis_t)
        tail_mesh.points            = tail_pts_buf
        tail_mesh["particle_depth"] = tail_depth_buf

        fill_head_buffer(head_pts_buf, head_depth_buf, pathlines, t, vis_t)
        head_mesh.points            = head_pts_buf.copy()
        head_mesh["particle_depth"] = head_depth_buf

        plotter.add_text(fmt_mpas_time(sim_start, dt_minutes, t),
                         position="upper_left", font_size=14, color="white", name="ts_label")

        plotter.screenshot(
            os.path.join(frames_dir, f"frame_{fi:06d}.png"),
            transparent_background=False,
        )

        if fi % 50 == 0:
            print(f"  Frame {fi:5d}/{len(frames)}  t={t}")

    plotter.close()

    out_path = _ffmpeg_stitch(frames_dir, args.output, FPS)
    print(f"Done → {out_path}")
    print(f"  (PNG frames kept in {frames_dir}/)")


if __name__ == "__main__":
    main()
