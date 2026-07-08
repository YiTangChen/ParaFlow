"""
All seeds are released *simultaneously* at t=0, but a particle's recorded
trajectory can end early when it:
  * beaches on a land cell (no ocean velocity available),
  * exits the model domain,
  * hits NaN / undefined velocity (e.g. under-ice cells), or
  * is pruned by a stuck-particle filter.
ParaFlow stops appending positions in those cases, so each particle stores
its own `nPts` — the renderer reads that value per particle.

This script walks the binary, collects every nPts, and reports the
distribution so you can tell whether the dataset is uniform (all particles
ran the full sim) or attrition-heavy.
"""

import os
import struct

import numpy as np


# Real-world minutes per stored step.
#   = ParaFlow yaml's dt × record_interval / 60
# For nersc_highres_256_cpu_10k.yaml: dt=60s, record_interval=60 → 60 min/step.
DT_MINUTES = 60.0

BIN_PATH = os.path.join(os.path.dirname(__file__), "pathlines_10k_active.bin")


def read_nPts(path):
    """
    Walk the binary file and return an int64 array of every particle's nPts.

    Layout per particle:  [int32 pid][int32 nPts][float64 × 3 × nPts coords]
    We read the header but skip past the coord block to keep memory low.
    """
    with open(path, "rb") as f:
        data = f.read()
    n_particles = struct.unpack_from("<i", data, 0)[0]
    offset = 4
    nps    = np.empty(n_particles, dtype=np.int64)
    for i in range(n_particles):
        _pid, nPts = struct.unpack_from("<ii", data, offset)
        nps[i]     = nPts
        offset    += 8 + nPts * 24       # 8-byte header + nPts × 3×float64
    return nps


def fmt_steps(steps):
    """Render a step count as 'N steps → days (years)'."""
    days  = steps * DT_MINUTES / 1440.0
    years = days / 365.25
    return f"{steps:>6,} steps  →  {days:>7,.1f} days  ({years:>4.2f} yr)"


def main():
    print(f"Reading {BIN_PATH}")
    nps = read_nPts(BIN_PATH)
    print(f"  {len(nps):,} particles\n")

    nmax = int(nps.max())

    print("Trajectory-length statistics:")
    print(f"  min     : {fmt_steps(int(nps.min()))}")
    print(f"  max     : {fmt_steps(nmax)}")
    print(f"  mean    : {fmt_steps(int(nps.mean()))}")
    print(f"  median  : {fmt_steps(int(np.median(nps)))}\n")

    print("Percentiles:")
    for p in (10, 25, 50, 75, 90, 99):
        print(f"  p{p:02d}    : {fmt_steps(int(np.percentile(nps, p)))}")
    print()

    full   = int((nps == nmax).sum())
    short  = int((nps < nmax * 0.5).sum())
    print(f"Particles reaching max length ({nmax:,} steps):")
    print(f"  {full:,} / {len(nps):,}  ({full / len(nps):.1%})")
    print(f"Particles ending before half of max length:")
    print(f"  {short:,} / {len(nps):,}  ({short / len(nps):.1%})\n")

    print("Histogram (10 bins across [0, max]):")
    bins        = np.linspace(0, nmax, 11)
    hist, _     = np.histogram(nps, bins=bins)
    bar_unit    = max(hist.max() // 50, 1)
    for i, count in enumerate(hist):
        lo, hi = int(bins[i]), int(bins[i + 1])
        bar    = "█" * int(count / bar_unit) if count > 0 else ""
        print(f"  [{lo:>6,} – {hi:>6,}]  {count:>5,}  {bar}")


if __name__ == "__main__":
    main()
