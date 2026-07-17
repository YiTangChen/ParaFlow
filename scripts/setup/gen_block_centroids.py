#!/usr/bin/env python3
"""
gen_block_centroids.py — per-block geographic centroids from a METIS partition
file + an MPAS-O mesh. Feeds the L5 spatial-work heatmap in plot_timing.py.

A partition file (`<prefix>.info.part.<N>`, produced by gen_partition.py) holds one
integer per global cell: the block/part id (0..N-1), in mesh `nCells` order. Each
block's centroid is the mean of its cells' coordinates — circular mean for longitude
so the dateline doesn't distort it.

Usage:
    python3 gen_block_centroids.py <mesh.nc> <partition.info.part.N> <out.csv> [--append]
        [--cells-out=block_cells.csv] [--cell-stride=25]

Output CSV columns: n_blocks, gid, lon, lat, n_cells
  gid = block id (matches ParaFlow gid and block_detail.csv gid)
  lon = circular-mean longitude (degrees, -180..180) ; lat = mean latitude (degrees)

--cells-out also writes a subsampled cell->block table (n_blocks, gid, lon, lat) so
plots can draw each block's actual geographic FOOTPRINT, not just its centroid. This
matters: on a variable-resolution mesh, block area varies several-fold, and centroid
dots hide that entirely.

Then point plot_timing.py at it (it looks for <results_dir>/block_centroids.csv):
    cp out.csv results/<run>/block_centroids.csv
Use --append to accumulate several partition sizes into one file (L5 filters by
n_blocks), e.g. run once per .part.16 / .part.32 / ... into the same out.csv.
"""

import sys
import csv
import numpy as np
import netCDF4 as nc


def read_partition(part_file):
    """One int per global cell → int array (fast whole-file read)."""
    with open(part_file) as f:
        return np.array(f.read().split(), dtype=np.int64)


def cell_lonlat_deg(mesh_file):
    """Return (lon_rad, lat_deg) for every cell. Prefers latCell/lonCell (radians);
    falls back to Cartesian xCell/yCell/zCell if the angular vars are absent."""
    ds = nc.Dataset(mesh_file, "r")
    try:
        if "latCell" in ds.variables and "lonCell" in ds.variables:
            lat = np.ma.filled(ds.variables["latCell"][:], np.nan).astype(float)
            lon = np.ma.filled(ds.variables["lonCell"][:], np.nan).astype(float)
            return lon, np.degrees(lat)
        x = np.ma.filled(ds.variables["xCell"][:], np.nan).astype(float)
        y = np.ma.filled(ds.variables["yCell"][:], np.nan).astype(float)
        z = np.ma.filled(ds.variables["zCell"][:], np.nan).astype(float)
        r = np.sqrt(x * x + y * y + z * z)
        return np.arctan2(y, x), np.degrees(np.arcsin(np.clip(z / r, -1, 1)))
    finally:
        ds.close()


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    opts = {a.split("=")[0]: a.split("=")[1]
            for a in sys.argv[1:] if a.startswith("--") and "=" in a}
    append = "--append" in sys.argv
    cells_out = opts.get("--cells-out")
    cell_stride = int(opts.get("--cell-stride", 25))
    if len(args) != 3:
        print(__doc__)
        sys.exit(1)
    mesh_file, part_file, out_csv = args

    part = read_partition(part_file)
    n_blocks = int(part.max()) + 1
    lon_rad, lat_deg = cell_lonlat_deg(mesh_file)

    if len(part) != len(lat_deg):
        sys.exit(f"ERROR: cell-count mismatch — partition {len(part)} vs mesh {len(lat_deg)}. "
                 "Is this the same mesh the partition was built from?")

    # Vectorised per-block reduction. Longitude via circular mean (mean of unit vectors).
    cnt = np.bincount(part, minlength=n_blocks).astype(float)
    lat = np.bincount(part, weights=lat_deg, minlength=n_blocks) / np.where(cnt > 0, cnt, np.nan)
    cx  = np.bincount(part, weights=np.cos(lon_rad), minlength=n_blocks)
    cy  = np.bincount(part, weights=np.sin(lon_rad), minlength=n_blocks)
    lon = np.degrees(np.arctan2(cy, cx))

    mode = "a" if (append and __import__("os").path.exists(out_csv)) else "w"
    with open(out_csv, mode, newline="") as f:
        w = csv.writer(f)
        if mode == "w":
            w.writerow(["n_blocks", "gid", "lon", "lat", "n_cells"])
        for gid in range(n_blocks):
            if cnt[gid] > 0:
                w.writerow([n_blocks, gid, round(float(lon[gid]), 6),
                            round(float(lat[gid]), 6), int(cnt[gid])])

    print(f"{'Appended to' if mode == 'a' else 'Wrote'} {out_csv}: "
          f"{int((cnt > 0).sum())} block centroids (n_blocks={n_blocks}, {len(part)} cells)")

    # Optional: subsampled cell -> block footprint table, so plots can draw the
    # actual extent of each block instead of a single centroid dot.
    if cells_out:
        import os
        lon_deg = np.degrees(np.arctan2(np.sin(lon_rad), np.cos(lon_rad)))   # wrap to -180..180
        sel = np.arange(0, len(part), cell_stride)
        cmode = "a" if (append and os.path.exists(cells_out)) else "w"
        # `cell` is the GLOBAL cell index — needed to join against a per-cell work
        # map (cell_workmap.npy), which is indexed by global cell id.
        rows = np.column_stack([np.full(len(sel), n_blocks), sel, part[sel],
                                lon_deg[sel].round(4), lat_deg[sel].round(4)])
        with open(cells_out, cmode, newline="") as f:
            if cmode == "w":
                f.write("n_blocks,cell,gid,lon,lat\n")
            np.savetxt(f, rows, fmt="%d,%d,%d,%.4f,%.4f")
        print(f"  {'appended' if cmode == 'a' else 'wrote'} cell footprints -> {cells_out} "
              f"(stride={cell_stride}, {len(sel):,} cells)")


if __name__ == "__main__":
    main()
