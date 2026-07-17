#!/usr/bin/env python3
"""
gen_cell_workmap.py — per-cell work map from already-written trajectory output.

Phase 0 of the load-imbalance work: we need a *cell-level* weight to re-partition
by work instead of by cell count. That weight can be recovered offline from the
trajectories ParaFlow already writes — no new instrumentation, no re-run, and DIY
is untouched.

How it works: each `<gid>.bin` is a stream of segments
    [int32 pid][int32 gid][int32 sid][int32 nPts] [float64 x,y,z] * nPts
Points are Cartesian on the sphere, the same frame as the mesh's xCell/yCell/zCell.
MPAS cells are Voronoi, so the *nearest cell centre is the containing cell* — a
KD-tree nearest-neighbour query maps every trajectory point to its cell exactly.
Binning the points per cell gives visit counts ∝ RK4 step density = the work map.
(Points are saved every `record_interval` steps, so counts are a fixed fraction of
the true step count — fine for relative weights.)

Usage:
    python3 gen_cell_workmap.py <mesh.nc> <trace_dir> <out.npy> [--stride N] [--limit N]

  --stride N   keep every Nth point (uniform subsampling preserves relative
               density; use it to trade accuracy for speed on huge runs)
  --limit N    only read the first N block files (quick smoke test)

Output: <out.npy>, an int64 array of length nCells = visits per cell.
Feed it to pymetis as vertex weights to build a work-balanced partition.
"""

import sys
import struct
import numpy as np
import netCDF4 as nc
from pathlib import Path
from scipy.spatial import cKDTree


def iter_block_points(path, stride=1):
    """Yield the (N,3) point array of one <gid>.bin. Offsets stay 8-byte aligned
    (header 16B, coords 24B/point), so frombuffer views are zero-copy."""
    data = Path(path).read_bytes()
    n, off, chunks = len(data), 0, []
    while off + 16 <= n:
        pid, gid, sid, npts = struct.unpack_from("iiii", data, off)
        off += 16
        nbytes = npts * 24
        if npts < 0 or off + nbytes > n:
            print(f"  warning: truncated segment in {Path(path).name} (pid={pid} sid={sid})")
            break
        if npts:
            pts = np.frombuffer(data, dtype=np.float64, count=npts * 3, offset=off).reshape(-1, 3)
            chunks.append(pts[::stride] if stride > 1 else pts)
        off += nbytes
    return np.concatenate(chunks) if chunks else np.empty((0, 3))


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    opts = {a.split("=")[0]: a.split("=")[1] for a in sys.argv[1:] if a.startswith("--") and "=" in a}
    if len(args) != 3:
        print(__doc__)
        sys.exit(1)
    mesh_file, trace_dir, out_npy = args
    stride = int(opts.get("--stride", 1))
    limit = int(opts.get("--limit", 0))

    ds = nc.Dataset(mesh_file, "r")
    xyz = np.column_stack([
        np.ma.filled(ds.variables[v][:], np.nan).astype(np.float64) for v in ("xCell", "yCell", "zCell")
    ])
    ds.close()
    n_cells = len(xyz)
    print(f"mesh: {n_cells} cells — building KD-tree...")
    tree = cKDTree(xyz)

    files = sorted(Path(trace_dir).glob("*.bin"), key=lambda p: int(p.stem) if p.stem.isdigit() else 1 << 30)
    files = [f for f in files if f.stem.isdigit()]
    if limit:
        files = files[:limit]
    if not files:
        sys.exit(f"ERROR: no <gid>.bin files in {trace_dir}")

    work = np.zeros(n_cells, dtype=np.int64)
    total = 0
    for i, f in enumerate(files):
        pts = iter_block_points(f, stride)
        if len(pts) == 0:
            continue
        _, idx = tree.query(pts, k=1, workers=-1)
        work += np.bincount(idx, minlength=n_cells)
        total += len(pts)
        print(f"  [{i+1}/{len(files)}] {f.name}: {len(pts):>9,} pts  (cum {total:>11,})")

    np.save(out_npy, work)
    hit = int((work > 0).sum())
    nz = work[work > 0]
    print(f"\nWrote {out_npy}")
    print(f"  points binned : {total:,}  (stride={stride})")
    print(f"  cells visited : {hit:,} / {n_cells:,}  ({hit/n_cells*100:.1f}%)")
    print(f"  work per cell : max={work.max():,}  mean(visited)={nz.mean():.1f}  "
          f"total={work.sum():,}")


if __name__ == "__main__":
    main()
