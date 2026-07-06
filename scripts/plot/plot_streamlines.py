"""
Plot reassembled streamlines on top of the DrawSubdomain ocean mask.

Reads streamlines.bin (output of read_traces.py) and draws each streamline
as a lat/lon polyline over the ocean background map.

Binary format of streamlines.bin:
  [int32 nStreamlines]
  For each streamline:
    [int32 pid][int32 nPts][float64 x0][float64 y0][float64 z0] ...

Usage:
  python plot_streamlines.py [streamlines_bin] [draw_dir] [output_png]
  python plot_streamlines.py trace_outputs/streamlines.bin drawSubdomain streamlines_map.png
"""

import struct
import math
import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize

EARTH_RADIUS = 6371229.0  # metres (MPAS-O reference)

LAT_RES = 1200
LON_RES = 2400


def read_ocean_mask(draw_dir):
    path = os.path.join(draw_dir, f"{LAT_RES}_{LON_RES}.bin")
    data = np.fromfile(path, dtype=np.float64).reshape((LAT_RES, LON_RES))
    return data


def read_streamlines_bin(path):
    """
    Read streamlines.bin produced by read_traces.py.
    Returns list of dicts: {'pid': int, 'pts': [(lat, lon, depth_m), ...]}
    depth_m = EARTH_RADIUS - r  (positive = below surface)
    """
    streamlines = []
    with open(path, 'rb') as f:
        n_sl, = struct.unpack('i', f.read(4))
        for _ in range(n_sl):
            pid, nPts = struct.unpack('ii', f.read(8))
            raw = f.read(nPts * 24)
            if len(raw) < nPts * 24:
                print(f'Warning: truncated data for pid={pid}')
                break
            pts = []
            for x, y, z in struct.iter_unpack('ddd', raw):
                lon   = math.degrees(math.atan2(y, x))
                lat   = math.degrees(math.atan2(z, math.hypot(x, y)))
                depth = EARTH_RADIUS - math.sqrt(x*x + y*y + z*z)
                pts.append((lat, lon, depth))
            streamlines.append({'pid': pid, 'pts': pts})
    return streamlines


def split_at_dateline(lons, lats):
    """
    Split a polyline into segments wherever a longitude jump > 180°
    occurs (dateline crossing), to avoid spurious horizontal lines.
    Returns list of (lons_seg, lats_seg).
    """
    segments = []
    seg_lons, seg_lats = [lons[0]], [lats[0]]
    for i in range(1, len(lons)):
        if abs(lons[i] - lons[i - 1]) > 180.0:
            if len(seg_lons) > 1:
                segments.append((seg_lons, seg_lats))
            seg_lons, seg_lats = [lons[i]], [lats[i]]
        else:
            seg_lons.append(lons[i])
            seg_lats.append(lats[i])
    if len(seg_lons) > 1:
        segments.append((seg_lons, seg_lats))
    return segments


def main():
    sl_file  = sys.argv[1] if len(sys.argv) > 1 else 'streamlines/streamlines.bin'
    draw_dir = sys.argv[2] if len(sys.argv) > 2 else 'drawSubdomain'
    out_png  = sys.argv[3] if len(sys.argv) > 3 else 'streamlines_map.png'

    print(f"Ocean mask   : {draw_dir}/{LAT_RES}_{LON_RES}.bin")
    print(f"Streamlines  : {sl_file}")

    ocean_mask  = read_ocean_mask(draw_dir)
    streamlines = read_streamlines_bin(sl_file)
    print(f"Loaded {len(streamlines)} streamlines")
    npts_list = [len(sl['pts']) for sl in streamlines]
    if npts_list:
        print(f"Points per streamline — min:{min(npts_list)}  max:{max(npts_list)}  mean:{sum(npts_list)/len(npts_list):.1f}")
        print(f"Streamlines with >1 pt: {sum(1 for n in npts_list if n > 1)} / {len(npts_list)}")

    # --- plot ---
    fig, ax = plt.subplots(figsize=(16, 8))

    # Background: ocean mask, rolled so column 0 = -180°
    rolled  = np.roll(ocean_mask, LON_RES // 2, axis=1)
    display = np.where(rolled != 0, 0.5, 0)
    ax.imshow(display, cmap='ocean', origin='lower', vmin=0, vmax=1,
              extent=[-180, 180, -90, 90], aspect='auto')

    # Color each segment by depth
    lw    = 0.5
    alpha = 0.8
    depth_cmap = matplotlib.colormaps['plasma']

    # Global depth range for a consistent colormap across all streamlines
    all_depths = [p[2] for sl in streamlines for p in sl['pts']]
    depth_min  = min(all_depths) if all_depths else 0.0
    depth_max  = max(all_depths) if all_depths else 1.0
    norm = Normalize(vmin=depth_min, vmax=depth_max*0.1)

    for sl in streamlines:
        pts = sl['pts']
        if len(pts) < 2:
            continue
        segments   = []
        seg_depths = []
        for i in range(len(pts) - 1):
            lat0, lon0, d0 = pts[i]
            lat1, lon1, d1 = pts[i + 1]
            if abs(lon1 - lon0) > 180.0:   # skip dateline crossings
                continue
            segments.append([(lon0, lat0), (lon1, lat1)])
            seg_depths.append((d0 + d1) * 0.5)
        if not segments:
            continue
        lc = LineCollection(segments, cmap=depth_cmap, norm=norm,
                            linewidths=lw, alpha=alpha, zorder=3)
        lc.set_array(np.array(seg_depths))
        ax.add_collection(lc)

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=depth_cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm, ax=ax, label='Depth (m)', shrink=0.6, pad=0.02)

    # Mark start points
    start_lons = [sl['pts'][0][1] for sl in streamlines]
    start_lats = [sl['pts'][0][0] for sl in streamlines]
    # ax.scatter(start_lons, start_lats, s=6, c='yellow', zorder=5,
    #            label='Seed (start)')

    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel('Longitude (degrees)')
    ax.set_ylabel('Latitude (degrees)')
    ax.set_title(f'Streamlines on ocean mask  (n={len(streamlines)})')
    ax.legend(loc='lower left', fontsize=8)
    ax.grid(True, linewidth=0.1, alpha=0.4)

    plt.tight_layout()
    plt.savefig(out_png, dpi=150)
    print(f"Saved: {out_png}")


if __name__ == '__main__':
    main()