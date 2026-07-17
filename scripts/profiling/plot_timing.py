#!/usr/bin/env python3
"""
Plot ParaFlow timing/memory results from CSV files.  Outputs: <results>/plots/*.png

Usage:  python3 plot_timing.py <results_dir>

The analysis is ordered load-imbalance-first, since that is the current focus:

Section 1 — Per-run load-imbalance diagnostics (one set per run):
  L1  Work distribution + Lorenz curve + Gini        (how bad is the skew?)
  L2  Work vs trace-time scatter (throughput CoV)     (work-bound or hardware?)
  L3  Partition weight (cells) vs actual work         (did we balance the wrong thing?)
  L4  Rebalance ceiling: current vs work-balanced      (how much is recoverable?)
  L5  Spatial work map: block footprints (optional)    (where are the hot blocks?)
  L6  Per-cell work density (needs cell_workmap.npy)   (work *inside* the blocks)
  L7  Phase Gantt per rank                             (timeline)

Section 2 — Scaling (varies n_blocks). S1/S2/S3/S6 all use the end-to-end wall
(t_run_max = block_load + exchange + write), falling back to exchange-only on old CSVs:
  S1  Strong scaling, end-to-end wall  S2  Parallel efficiency    S3  Speedup S(P)
  S4  Imbalance vs scale: TIME (barrier-flattened) vs WORK (the real skew)
  S5  Idle fraction vs scale          S6  Weak scaling (needs seeds ∝ blocks)
  S7  Phase breakdown vs scale        S8  Memory per block vs scale

Section 3 — CPU vs GPU:
  G1  End-to-end wall time   G2  Aligned phase breakdown   G3  Speedup by phase
  G4  Peak memory            G5  GPU seed scalability
"""

import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
from matplotlib.colors import LogNorm
import numpy as np
from pathlib import Path

# ── config ────────────────────────────────────────────────────────────────────
RESULTS_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("results/gpu_analysis")
PLOTS_DIR   = RESULTS_DIR / "plots"
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

SUMMARY_CSV   = RESULTS_DIR / "run_summary.csv"
BLOCK_CSV     = RESULTS_DIR / "block_detail.csv"
CENTROIDS_CSV = RESULTS_DIR / "block_centroids.csv"   # optional: n_blocks,gid,lon,lat
CELLS_CSV     = RESULTS_DIR / "block_cells.csv"       # optional: n_blocks,cell,gid,lon,lat
WORKMAP_NPY   = RESULTS_DIR / "cell_workmap.npy"      # optional: per-cell work, indexed by cell id

C_COMPUTE = "#2196F3"
C_IDLE    = "#FF5722"
C_LOAD    = "#4CAF50"
C_WRITE   = "#9C27B0"
C_GRID    = "#00BCD4"
C_SOL     = "#F44336"
C_IDEAL   = "#9E9E9E"
C_GPU     = "#F57C00"
C_CPU     = "#1565C0"
C_KERNEL  = "#FF9800"
C_WORK    = "#6A1B9A"

SEED_COLORS = ["#1565C0", "#2E7D32", "#6A1B9A"]
KEYS = ["n_blocks", "n_seeds", "run_id", "device", "mesh"]


def save(fig, name):
    path = PLOTS_DIR / name
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved: {path}")


# ── helpers ───────────────────────────────────────────────────────────────────
def series(df, name, default=0.0):
    """Column as float ndarray, or a filled default if absent (old-CSV safe)."""
    if name in df.columns:
        return pd.to_numeric(df[name], errors="coerce").fillna(default).values.astype(float)
    return np.full(len(df), float(default))


def _cov(v):
    v = np.asarray(v, float)
    m = v.mean() if v.size else 0.0
    return (v.std() / m) if m > 0 else np.nan


def _moa(v):
    v = np.asarray(v, float)
    m = v.mean() if v.size else 0.0
    return (v.max() / m) if m > 0 else np.nan


def _gini(v):
    v = np.sort(np.asarray(v, float))
    n, s = len(v), v.sum()
    if n == 0 or s == 0:
        return 0.0
    idx = np.arange(1, n + 1)
    return (2.0 * np.sum(idx * v)) / (n * s) - (n + 1.0) / n


def _lorenz(v):
    v = np.sort(np.asarray(v, float))
    n, s = len(v), v.sum()
    if n == 0 or s == 0:
        return np.array([0.0, 1.0]), np.array([0.0, 1.0])
    return np.arange(0, n + 1) / n, np.concatenate([[0.0], np.cumsum(v) / s])


# ── load data ─────────────────────────────────────────────────────────────────
summary = pd.read_csv(SUMMARY_CSV)
blocks  = pd.read_csv(BLOCK_CSV)
for df in (summary, blocks):
    if "device" not in df.columns: df["device"] = "cpu"
    if "mesh"   not in df.columns: df["mesh"]   = "lowres"

centroids = pd.read_csv(CENTROIDS_CSV) if CENTROIDS_CSV.exists() else None
cells_map = pd.read_csv(CELLS_CSV) if CELLS_CSV.exists() else None
cell_work = np.load(WORKMAP_NPY) if WORKMAP_NPY.exists() else None

# Scaling metric = end-to-end wall (run_streamline/run_pathline total: block_load +
# exchange + write), which is what a strong-scaling study should report — the phases
# that do not scale are exactly what the measurement must expose. Falls back to the
# exchange-only wall for older CSVs that predate t_run_max.
summary["wall"] = 0.0
for c in ("t_run_max", "t_iex_max", "t_dist_trace_total_max"):
    if c in summary.columns:
        summary["wall"] = summary["wall"].where(summary["wall"] > 0, series(summary, c))

# Idle within the exchange = exchange_wall − local compute (both avg-per-rank).
ex   = series(summary, "t_iex_avg")
comp = series(summary, "t_compute_avg")   # compat alias == t_local_wall_avg
summary["idle_fraction"] = np.clip((ex - comp) / np.where(ex > 0, ex, np.nan), 0, 1)

# Per-run imbalance, computed from block_detail so it works on any recent CSV.
def _run_imbalance(g):
    work  = series(g, "n_steps_total")
    wall  = series(g, "t_trace_local_wall")
    cells = series(g, "n_local_cells")
    return pd.Series({
        "work_moa":  _moa(work),   "work_cov": _cov(work),
        "time_moa":  _moa(wall),   "work_gini": _gini(work),
        "cells_cov": _cov(cells),
        "ceiling":   (wall.max() / wall.mean()) if wall.mean() > 0 else np.nan,
    })

imb = (blocks.groupby(KEYS, group_keys=False).apply(_run_imbalance).reset_index())
summary = summary.merge(imb, on=KEYS, how="left")

valid = summary[summary["wall"] > 0].copy()

print(f"Loaded {len(summary)} runs ({len(valid)} with wall>0), {len(blocks)} block rows.")
print(f"Devices:  {sorted(valid.device.unique())}")
print(f"Meshes:   {sorted(valid.mesh.unique())}")
print(f"n_blocks: {sorted(valid.n_blocks.unique())}")
print(f"n_seeds:  {sorted(valid.n_seeds.unique())}")
print(f"Centroids for spatial plot: {'yes' if centroids is not None else 'no (skipping L5)'}")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — per-run load-imbalance diagnostics
# ═══════════════════════════════════════════════════════════════════════════════

def plot_per_run(df_blocks, n_blocks, n_seeds, run_id, device, mesh):
    tag   = f"{device}_{mesh}_b{n_blocks:03d}_s{n_seeds:07d}_r{run_id}"
    title = f"{device.upper()} {mesh} · {n_blocks} blocks · {n_seeds:,} seeds"
    df    = df_blocks.sort_values("gid").reset_index(drop=True)
    n     = len(df)
    w     = max(8, n * 0.5 + 2)
    work  = series(df, "n_steps_total")
    wall  = series(df, "t_trace_local_wall")
    cells = series(df, "n_local_cells")

    # ── L1: Work distribution + Lorenz + Gini ─────────────────────────────────
    if work.sum() > 0:
        avg, mx = work.mean(), work.max()
        moa = mx / avg if avg > 0 else 0.0
        gini = _gini(work)
        lx, ly = _lorenz(work)
        fig, (axL, axR) = plt.subplots(1, 2, figsize=(min(15, w + 5), 4.5))
        axL.bar(range(n), np.sort(work)[::-1] / 1e6, color=C_WORK, zorder=3)
        axL.axhline(avg / 1e6, color="red", ls="--", lw=1.4, label=f"avg = {avg/1e6:.2f} M")
        axL.set_xlabel("Block (sorted by work)"); axL.set_ylabel("RK4 steps (millions)")
        axL.set_title(f"Work per block — max/avg = {moa:.2f}×")
        axL.legend(); axL.grid(axis="y", alpha=0.3, zorder=0)

        axR.plot([0, 1], [0, 1], color=C_IDEAL, ls="--", lw=1.2, label="perfect balance")
        axR.plot(lx, ly, color=C_WORK, lw=2, label="work")
        axR.fill_between(lx, ly, lx, alpha=0.15, color=C_WORK)
        axR.set_xlabel("Cumulative fraction of blocks")
        axR.set_ylabel("Cumulative fraction of work")
        axR.set_title(f"Lorenz curve — Gini = {gini:.3f}")
        axR.legend(fontsize=8, loc="upper left"); axR.grid(alpha=0.3)
        fig.suptitle(f"Work Imbalance — {title}")
        plt.tight_layout()
        save(fig, f"{tag}_L1_work_lorenz.png")

    # ── L2: Work vs trace-time (throughput CoV) ───────────────────────────────
    ok = wall > 0
    if ok.sum() > 0 and work.sum() > 0:
        tput = np.zeros_like(work); tput[ok] = work[ok] / wall[ok] / 1e6
        tmean = tput[ok].mean(); tcov = tput[ok].std() / tmean if tmean > 0 else 0.0
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.scatter(work / 1e6, wall, s=40, color=C_COMPUTE, zorder=3)
        if tmean > 0:
            xs = np.array([0.0, work.max() / 1e6])
            ax.plot(xs, xs / tmean, color=C_IDEAL, ls="--", lw=1.3,
                    label=f"uniform throughput ({tmean:.1f} M steps/s)")
        verdict = ("work-bound → re-partition by work" if tcov < 0.15
                   else "throughput varies → NUMA / GPU contention")
        ax.set_xlabel("Work (M steps)"); ax.set_ylabel("trace_local_wall (s)")
        ax.set_title(f"Work vs time — throughput CoV = {tcov:.2f}\n{verdict}")
        ax.legend(fontsize=8); ax.grid(alpha=0.3)
        save(fig, f"{tag}_L2_work_vs_throughput.png")

    # ── L3: Partition weight (cells) vs actual work ───────────────────────────
    if cells.sum() > 0 and work.sum() > 0:
        ccov, wcov = _cov(cells), _cov(work)
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.scatter(cells, work / 1e6, s=40, color=C_LOAD, zorder=3)
        note = ("  → cells balanced, work is not: re-weight the partition by work"
                if (np.isfinite(ccov) and np.isfinite(wcov) and wcov > max(ccov, 0.02) * 1.5) else "")
        ax.set_xlabel("n_local_cells  (quantity the graph partition balances)")
        ax.set_ylabel("Work — RK4 steps (millions)")
        ax.set_title(f"Partition weight vs actual work — {title}\n"
                     f"cells CoV = {ccov:.2f}   work CoV = {wcov:.2f}{note}")
        ax.grid(alpha=0.3)
        save(fig, f"{tag}_L3_partition_vs_work.png")

    # ── L4: Rebalance ceiling (current slowest vs work-balanced ideal) ────────
    if wall.sum() > 0:
        current, ideal = wall.max(), wall.mean()
        ceil = current / ideal if ideal > 0 else 0.0
        fig, ax = plt.subplots(figsize=(5, 5))
        bars = ax.bar(["current\n(slowest block sets wall)", "work-balanced\n(ideal)"],
                      [current, ideal], color=[C_IDLE, C_COMPUTE], width=0.5, zorder=3)
        for b, v in zip(bars, [current, ideal]):
            ax.text(b.get_x() + b.get_width() / 2, v, f"{v:.1f} s",
                    ha="center", va="bottom", fontsize=10, fontweight="bold")
        ax.set_ylabel("Local trace wall (s)")
        ax.set_title(f"Rebalance ceiling — up to {ceil:.2f}× faster\n{title}")
        ax.grid(axis="y", alpha=0.3, zorder=0)
        save(fig, f"{tag}_L4_imbalance_ceiling.png")

    # ── L5: Spatial work map — each block's FOOTPRINT coloured by its work ────
    # Centroid dots alone are misleading: every dot renders the same size, hiding
    # that block area varies several-fold on a variable-resolution mesh (equatorial
    # blocks are geographically huge, polar ones tiny). Drawing the real cell
    # footprint shows the area difference that actually drives the seed/work skew.
    if centroids is not None:
        cen = centroids[centroids.n_blocks == n_blocks] if "n_blocks" in centroids.columns else centroids
        if len(cen) and {"gid", "lon", "lat"}.issubset(cen.columns):
            workby = dict(zip(df.gid.astype(int).values, work))
            cen = cen[cen.gid.isin(workby.keys())]
            fig, ax = plt.subplots(figsize=(11, 5.2))
            sc = None
            if cells_map is not None:
                cm = cells_map[cells_map.n_blocks == n_blocks] if "n_blocks" in cells_map.columns else cells_map
                cm = cm[cm.gid.isin(workby.keys())]
                if len(cm):
                    sc = ax.scatter(cm.lon, cm.lat, c=cm.gid.map(workby).astype(float) / 1e6,
                                    cmap="inferno", s=1.5, marker=".", linewidths=0,
                                    rasterized=True)
            if sc is None:      # no footprint file — fall back to centroid dots
                sc = ax.scatter(cen.lon, cen.lat, c=cen.gid.map(workby).astype(float) / 1e6,
                                cmap="inferno", s=70, zorder=3)
                note = "centroids only — run gen_block_centroids.py --cells-out for footprints"
            else:               # footprints drawn: overlay the centroid inside each block
                ax.scatter(cen.lon, cen.lat, s=22, facecolor="none", edgecolor="#00E5FF",
                           linewidths=1.0, zorder=4)
                if n_blocks <= 32:
                    for _, r in cen.iterrows():
                        ax.annotate(str(int(r.gid)), (r.lon, r.lat), fontsize=6,
                                    color="#00E5FF", ha="center", va="center", zorder=5)
                note = "block footprint coloured by that block's work · ring = centroid"
            fig.colorbar(sc, ax=ax, label="Block work (M steps)")
            ax.set_xlim(-180, 180); ax.set_ylim(-90, 90)
            ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
            ax.set_title(f"Spatial work map — {title}\n{note}")
            save(fig, f"{tag}_L5_spatial_work.png")

    # ── L6: Per-cell work map — the work density INSIDE the blocks ────────────
    # L5 gives each block one flat colour (its total work). This shows where the
    # work actually sits within a block, and how much of the mesh does nothing at
    # all. The map is partition-independent (a function of flow field + seeds), so
    # only the overlaid block centroids change between runs.
    if cells_map is not None and cell_work is not None and "cell" in cells_map.columns:
        cm = cells_map[cells_map.n_blocks == n_blocks] if "n_blocks" in cells_map.columns else cells_map
        idx = cm.cell.values.astype(np.int64)
        keep = idx < len(cell_work)
        idx, cm = idx[keep], cm[keep]
        if len(idx):
            wv = cell_work[idx].astype(float)
            lon, lat = cm.lon.values, cm.lat.values
            dead = wv <= 0
            fig, ax = plt.subplots(figsize=(11, 5.2))
            if dead.any():
                ax.scatter(lon[dead], lat[dead], c="#D9D9D9", s=1.0, marker=".",
                           linewidths=0, rasterized=True)
            live = ~dead
            sc = ax.scatter(lon[live], lat[live], c=wv[live], cmap="inferno",
                            norm=LogNorm(vmin=1, vmax=max(wv.max(), 2.0)),
                            s=1.5, marker=".", linewidths=0, rasterized=True)
            fig.colorbar(sc, ax=ax, label="Work per cell — trajectory points (log)")
            if centroids is not None:
                cen = centroids[centroids.n_blocks == n_blocks] if "n_blocks" in centroids.columns else centroids
                ax.scatter(cen.lon, cen.lat, s=22, facecolor="none", edgecolor="#00E5FF",
                           linewidths=1.0, zorder=4)
            ax.set_xlim(-180, 180); ax.set_ylim(-90, 90)
            ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
            ax.set_title(f"Per-cell work density — {title}\n"
                         f"grey = never visited ({dead.mean()*100:.0f}% of cells) · "
                         f"ring = block centroid")
            save(fig, f"{tag}_L6_cell_work_map.png")

    # ── L7: Phase Gantt per rank ──────────────────────────────────────────────
    phases = [("t_block_load", "Block load", C_LOAD),
              ("t_trace_compute", "Compute", C_COMPUTE),
              ("t_output_write", "Write", C_WRITE),
              ("t_trace_idle", "Idle", C_IDLE)]
    if any(p in df.columns for p, _, _ in phases):
        fig, ax = plt.subplots(figsize=(10, max(4, n * 0.28 + 1.5)))
        for _, row in df.iterrows():
            start = 0.0
            for c, _l, color in phases:
                dur = float(row[c]) if c in df.columns and pd.notna(row[c]) else 0.0
                if dur > 1e-4:
                    ax.barh(row.gid, dur, left=start, height=0.6, color=color, alpha=0.9)
                start += dur
        ax.set_yticks(df.gid.values)
        ax.set_yticklabels([f"rank {int(r)}" for r in df.get("rank", df.gid)], fontsize=8)
        ax.set_xlabel("Approximate elapsed time (s)")
        ax.set_title(f"Phase Timeline per Rank — {title}")
        ax.legend(handles=[mpatches.Patch(color=c, label=l) for _, l, c in phases],
                  loc="lower right", fontsize=8)
        ax.grid(axis="x", alpha=0.3)
        save(fig, f"{tag}_L7_gantt_timeline.png")

    gv = _gini(work) if work.sum() > 0 else 0.0
    cv = (wall.max() / wall.mean()) if wall.sum() > 0 and wall.mean() > 0 else 0.0
    print(f"  → work Gini = {gv:.3f}   rebalance ceiling = {cv:.2f}×")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — scaling
# ═══════════════════════════════════════════════════════════════════════════════

def plot_scaling(valid, device_filter=None):
    if device_filter:
        valid = valid[valid.device == device_filter]
    suffix = f"_{device_filter}" if device_filter else ""
    dev_tag = f" — {device_filter.upper()}" if device_filter else ""

    seed_counts  = sorted(valid.n_seeds.unique())
    block_counts = sorted(valid.n_blocks.unique())
    if len(block_counts) < 2:
        print(f"  [S{suffix}] only 1 n_blocks value — scaling skipped")
        return
    sc = (SEED_COLORS * 4)[:max(1, len(seed_counts))]

    # ── S1: Strong scaling ────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7.5, 5))
    for i, ns in enumerate(seed_counts):
        df = valid[valid.n_seeds == ns].sort_values("n_blocks")
        if len(df) < 2: continue
        wmin = df.wall.values / 60.0                       # seconds -> minutes
        p0, t0 = df.iloc[0].n_blocks, wmin[0]
        ax.plot(df.n_blocks, wmin, "o-", color=sc[i], label=f"{ns:,} seeds")
        ax.plot(df.n_blocks, t0 * p0 / df.n_blocks.values, "--", alpha=0.3, color=sc[i])
        for xv, yv in zip(df.n_blocks.values, wmin):       # label every point
            ax.annotate(f"{yv:.0f}", xy=(xv, yv), xytext=(5, 5), textcoords="offset points",
                        fontsize=7.5, color=sc[i], fontweight="bold")
    ax.set_xscale("log", base=2); ax.set_yscale("log")
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter()); ax.set_xticks(block_counts)
    # Denser, plainly-readable y ticks (log axis would otherwise only label 10/100/1000).
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10, subs=(1.0,), numticks=12))
    ax.yaxis.set_minor_locator(ticker.LogLocator(base=10, subs=(2, 3, 4, 5, 6, 8), numticks=12))
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.yaxis.set_minor_formatter(ticker.ScalarFormatter())
    ax.tick_params(axis="y", which="minor", labelsize=7)
    ax.set_xlabel("Number of blocks (MPI ranks)"); ax.set_ylabel("End to end wall (minutes)")
    ax.set_title(f"Strong Scaling end to end{dev_tag}\ndashed = ideal")
    ax.legend(); ax.grid(True, alpha=0.3, which="both")
    save(fig, f"S1_strong_scaling{suffix}.png")

    # ── S2: Parallel efficiency ───────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 5))
    for i, ns in enumerate(seed_counts):
        df = valid[valid.n_seeds == ns].sort_values("n_blocks")
        if len(df) < 2: continue
        p0, t0 = df.iloc[0].n_blocks, df.iloc[0].wall
        ax.plot(df.n_blocks, (t0 * p0) / (df.wall * df.n_blocks), "o-", color=sc[i],
                label=f"{ns:,} seeds")
    ax.axhline(1.0, color="black", ls="--", lw=1, label="ideal (100%)")
    ax.set_xscale("log", base=2); ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.set_xticks(block_counts); ax.set_ylim(0, 1.15)
    ax.set_xlabel("Number of blocks (MPI ranks)"); ax.set_ylabel("Parallel efficiency")
    ax.set_title(f"Parallel Efficiency end to end{dev_tag}")
    ax.legend(); ax.grid(True, alpha=0.3)
    save(fig, f"S2_parallel_efficiency{suffix}.png")

    # ── S3: Speedup S(P) ──────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 5))
    for i, ns in enumerate(seed_counts):
        df = valid[valid.n_seeds == ns].sort_values("n_blocks")
        if len(df) < 2: continue
        ax.plot(df.n_blocks, df.iloc[0].wall / df.wall.values, "o-", color=sc[i],
                label=f"{ns:,} seeds")
    p_all = np.array(block_counts)
    ax.plot(p_all, p_all / p_all[0], "--", color=C_IDEAL, lw=1.5, label="ideal")
    ax.set_xscale("log", base=2); ax.set_yscale("log", base=2)
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter()); ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.set_xticks(block_counts)
    ax.set_xlabel("Number of blocks (MPI ranks)"); ax.set_ylabel("Speedup S(P) = T(P₀)/T(P)")
    ax.set_title(f"Speedup Curve{dev_tag}")
    ax.legend(); ax.grid(True, alpha=0.3, which="both")
    save(fig, f"S3_speedup{suffix}.png")

    # ── S4: Imbalance vs scale — TIME (barrier-flattened) vs WORK (real) ──────
    # The headline load-imbalance plot: exchange-wall max/min stays ~1 while the
    # underlying work skew grows with scale. Data from block_detail (work_moa).
    fig, ax = plt.subplots(figsize=(7.5, 5))
    for i, ns in enumerate(seed_counts):
        df = valid[valid.n_seeds == ns].sort_values("n_blocks")
        if len(df) < 2: continue
        if "work_moa" in df.columns:
            ax.plot(df.n_blocks, df.work_moa, "o-", color=sc[i], lw=2,
                    label=f"work skew (max/avg) — {ns:,} seeds")
        if "time_moa" in df.columns:
            ax.plot(df.n_blocks, df.time_moa, "s--", color=sc[i], alpha=0.6,
                    label=f"time skew — {ns:,} seeds")
    ax.axhline(1.0, color="black", ls=":", lw=1, label="perfect balance")
    ax.set_xscale("log", base=2); ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.set_xticks(block_counts)
    ax.set_xlabel("Number of blocks"); ax.set_ylabel("Imbalance  (max / avg)")
    ax.set_title(f"Load Imbalance vs Scale{dev_tag}\n"
                 "solid = work skew (real) · dashed = time skew (hidden by the barrier)")
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3)
    save(fig, f"S4_imbalance_vs_scale{suffix}.png")

    # ── S5: Idle fraction vs scale ────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 5))
    for i, ns in enumerate(seed_counts):
        df = valid[valid.n_seeds == ns].sort_values("n_blocks")
        if len(df) < 2: continue
        ax.plot(df.n_blocks, df.idle_fraction * 100, "o-", color=sc[i], label=f"{ns:,} seeds")
    ax.set_xscale("log", base=2); ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.set_xticks(block_counts); ax.set_ylim(0, 105)
    ax.set_xlabel("Number of blocks"); ax.set_ylabel("Idle fraction (% of exchange wall)")
    ax.set_title(f"Idle Fraction vs Scale{dev_tag}")
    ax.legend(); ax.grid(True, alpha=0.3)
    save(fig, f"S5_idle_fraction{suffix}.png")

    # ── S6: Weak scaling (needs seeds ∝ blocks; error bars from repetitions) ──
    v = valid.copy()
    v["spb"] = (v.n_seeds / v.n_blocks).round(3)
    drew = False
    fig, ax = plt.subplots(figsize=(7, 5))
    for i, spb in enumerate(sorted(v.spb.unique())):
        grp = v[v.spb == spb]
        if grp.n_blocks.nunique() < 2:
            continue
        agg = grp.groupby("n_blocks").wall.agg(["mean", "std", "count"]).sort_index()
        base = agg["mean"].iloc[0]
        eff = base / agg["mean"]                       # weak-scaling: ideal is flat 1.0
        yerr = (agg["std"] / agg["mean"] * eff).fillna(0)
        ax.errorbar(agg.index, eff, yerr=yerr, fmt="o-", capsize=3,
                    color=sc[i % len(sc)], label=f"{int(spb)} seeds/block")
        drew = True
    if drew:
        ax.axhline(1.0, color="black", ls="--", lw=1, label="ideal")
        ax.set_xscale("log", base=2); ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
        ax.set_xticks(block_counts); ax.set_ylim(0, 1.15)
        ax.set_xlabel("Number of blocks"); ax.set_ylabel("Weak-scaling efficiency")
        ax.set_title(f"Weak Scaling{dev_tag}  (error bars = run-to-run std)")
        ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
        save(fig, f"S6_weak_scaling{suffix}.png")
    else:
        plt.close(fig)
        print(f"  [S6{suffix}] no weak-scaling data (need seeds ∝ blocks) — skipped")

    # ── S7: Phase breakdown vs scale ──────────────────────────────────────────
    v = valid.copy()
    v["t_idle_avg"] = (series(v, "t_iex_avg") - series(v, "t_compute_avg")).clip(min=0)
    for ns in seed_counts:
        df = v[v.n_seeds == ns].sort_values("n_blocks")
        if len(df) < 2: continue
        xpos = np.arange(len(df)); bottom = np.zeros(len(df))
        fig, ax = plt.subplots(figsize=(8, 4))
        for col, label, color in [("t_blockload_avg", "Block load", C_LOAD),
                                   ("t_compute_avg", "Compute", C_COMPUTE),
                                   ("t_write_avg", "Write", C_WRITE),
                                   ("t_idle_avg", "Idle", C_IDLE)]:
            vals = series(df, col)
            ax.bar(xpos, vals, bottom=bottom, label=label, color=color,
                   alpha=0.88 if col == "t_idle_avg" else 1.0)
            bottom = bottom + vals
        ax.set_xticks(xpos); ax.set_xticklabels([str(int(b)) for b in df.n_blocks])
        ax.set_xlabel("Number of blocks"); ax.set_ylabel("Avg time per rank (s)")
        ax.set_title(f"Phase Breakdown vs Scale — {ns:,} seeds{dev_tag}")
        ax.legend(loc="upper right", fontsize=8); ax.grid(axis="y", alpha=0.3, zorder=0)
        save(fig, f"S7_phase_breakdown_s{ns:07d}{suffix}.png")

    # ── S8: Memory per block vs scale ─────────────────────────────────────────
    for ns in seed_counts:
        df = valid[valid.n_seeds == ns].sort_values("n_blocks")
        if len(df) < 2 or "mem_grid_bytes_avg" not in df.columns: continue
        grid_mb = series(df, "mem_grid_bytes_avg") / 1e6
        sol_mb  = series(df, "mem_solution_bytes_avg") / 1e6
        peak_mb = series(df, "mem_peak_vmhwm_kb_avg") / 1024
        xpos = np.arange(len(df))
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.bar(xpos, grid_mb, label="Grid topology", color=C_GRID)
        ax.bar(xpos, sol_mb, bottom=grid_mb, label="Velocity solution", color=C_SOL)
        ax2 = ax.twinx()
        ax2.plot(xpos, peak_mb, "k^--", markersize=6, lw=1.2, label="Peak RSS (VmHWM)")
        ax2.set_ylabel("Peak RSS per rank (MB)")
        ax.set_xticks(xpos); ax.set_xticklabels([str(int(b)) for b in df.n_blocks])
        ax.set_xlabel("Number of blocks"); ax.set_ylabel("Analytical data memory / block (MB)")
        ax.set_title(f"Memory per Block vs Scale — {ns:,} seeds{dev_tag}")
        l1, la1 = ax.get_legend_handles_labels(); l2, la2 = ax2.get_legend_handles_labels()
        ax.legend(l1 + l2, la1 + la2, loc="upper right", fontsize=8)
        ax.grid(axis="y", alpha=0.3, zorder=0)
        save(fig, f"S8_memory_vs_scale_s{ns:07d}{suffix}.png")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — CPU vs GPU
# ═══════════════════════════════════════════════════════════════════════════════

def plot_gpu_vs_cpu(valid, blocks):
    print("\n── Section 3: CPU vs GPU ──")
    # Auto-pick a (mesh, blocks, seeds) config present for BOTH devices; prefer 16 / 10k.
    combos = valid.groupby(["mesh", "n_blocks", "n_seeds"]).device.nunique()
    both = [k for k, v in combos.items() if v >= 2]
    if not both:
        print("  no config has both CPU and GPU — Section 3 skipped")
        return
    both.sort(key=lambda k: (k[1] != 16, k[2] != 10000, k[1]))
    mesh_r, nb_r, ns_r = both[0]
    reftag = f"{mesh_r} · {nb_r} blocks · {ns_r:,} seeds"
    print(f"  reference config: {reftag}")
    ref = valid[(valid.mesh == mesh_r) & (valid.n_blocks == nb_r) & (valid.n_seeds == ns_r)]
    cpu_ref, gpu_ref = ref[ref.device == "cpu"], ref[ref.device == "gpu"]
    have_ref = len(cpu_ref) and len(gpu_ref)

    def g(row, col):
        try: return float(row.get(col, 0.0) or 0.0)
        except Exception: return 0.0

    def e2e(row):
        r = float(row.get("t_run_max", 0) or 0)
        return r if r > 0 else float(row.get("wall", 0) or 0)

    # ── G1: End-to-end wall time ──────────────────────────────────────────────
    if have_ref:
        cw, gw = e2e(cpu_ref.iloc[0]), e2e(gpu_ref.iloc[0])
        sp = cw / gw if gw > 0 else 0.0
        fig, ax = plt.subplots(figsize=(6, 5))
        bars = ax.bar(["CPU (16 ranks)", "GPU (16 ranks)"], [cw / 3600, gw / 3600],
                      color=[C_CPU, C_GPU], width=0.45)
        for b, v in zip(bars, [cw, gw]):
            ax.text(b.get_x() + b.get_width() / 2, v / 3600, f"{v/3600:.2f} h",
                    ha="center", va="bottom", fontsize=11, fontweight="bold")
        ax.set_ylabel("End-to-end wall time (hours)")
        ax.set_title(f"CPU vs GPU End-to-End Wall\n({reftag} · {sp:.1f}× speedup)")
        ax.grid(axis="y", alpha=0.3)
        save(fig, "G1_wall_time.png")

    # ── G2: Aligned phase breakdown ───────────────────────────────────────────
    if have_ref:
        cr, gr = cpu_ref.iloc[0], gpu_ref.iloc[0]
        trace = [("t_prepare_avg", "Prepare", "#8E24AA"),
                 ("t_transfer_avg", "Transfer", C_KERNEL),
                 ("t_integrate_avg", "Integrate", C_COMPUTE),
                 ("t_postprocess_avg", "Postprocess", "#00897B"),
                 ("t_enqueue_avg", "Enqueue", "#546E7A")]
        e2e_ph = [("t_blockload_avg", "Block load", C_LOAD)] + trace + [("t_output_write_avg", "Write", C_WRITE)]
        fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5))
        for ax, phs, ttl in [(axL, e2e_ph, "End-to-end"), (axR, trace, "Trace buckets (aligned)")]:
            bc = bg = 0.0
            for col, label, color in phs:
                cv, gv = g(cr, col), g(gr, col)
                ax.bar([0], [cv], bottom=[bc], color=color, label=label, width=0.45)
                ax.bar([1], [gv], bottom=[bg], color=color, width=0.45, alpha=0.85)
                bc += cv; bg += gv
            ax.set_xticks([0, 1]); ax.set_xticklabels(["CPU", "GPU"])
            ax.set_ylabel("Avg time per rank (s)"); ax.set_title(ttl)
            ax.legend(loc="upper right", fontsize=8); ax.grid(axis="y", alpha=0.3)
        fig.suptitle(f"Aligned Phase Breakdown: CPU vs GPU ({reftag})")
        plt.tight_layout()
        save(fig, "G2_phase_breakdown.png")

    # ── G3: Speedup by phase ──────────────────────────────────────────────────
    if have_ref:
        cr, gr = cpu_ref.iloc[0], gpu_ref.iloc[0]
        def ratio(col): d = g(gr, col); return (g(cr, col) / d) if d > 0 else 0.0
        metrics = ["Wall", "Block load", "Integrate", "Write"]
        vals = [e2e(cr) / e2e(gr) if e2e(gr) > 0 else 0,
                ratio("t_blockload_avg"), ratio("t_integrate_avg"), ratio("t_output_write_avg")]
        fig, ax = plt.subplots(figsize=(7, 5))
        bars = ax.bar(metrics, vals, color=[C_GPU, C_LOAD, C_COMPUTE, C_WRITE], width=0.5)
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v, f"{v:.1f}×",
                    ha="center", va="bottom", fontsize=11, fontweight="bold")
        ax.axhline(1.0, color="grey", ls="--", lw=1, label="1× (no speedup)")
        ax.set_ylabel("GPU speedup over CPU"); ax.set_title("GPU Speedup by Phase")
        ax.legend(fontsize=8); ax.grid(axis="y", alpha=0.3)
        save(fig, "G3_speedup_by_phase.png")

    # ── G4: Peak memory CPU vs GPU ────────────────────────────────────────────
    rb = blocks[(blocks.mesh == mesh_r) & (blocks.n_blocks == nb_r) & (blocks.n_seeds == ns_r)]
    cb, gb = rb[rb.device == "cpu"].sort_values("gid"), rb[rb.device == "gpu"].sort_values("gid")
    if len(cb) and len(gb) and "mem_peak_vmhwm_kb" in rb.columns:
        gids = sorted(set(cb.gid) & set(gb.gid))
        cr = [float(cb[cb.gid == i]["mem_peak_vmhwm_kb"].iloc[0]) / 1024 for i in gids]
        gr = [float(gb[gb.gid == i]["mem_peak_vmhwm_kb"].iloc[0]) / 1024 for i in gids]
        x = np.arange(len(gids)); wbar = 0.35
        fig, ax = plt.subplots(figsize=(max(8, len(gids) * 0.6 + 2), 4))
        ax.bar(x - wbar / 2, cr, width=wbar, label="CPU", color=C_CPU, alpha=0.85)
        ax.bar(x + wbar / 2, gr, width=wbar, label="GPU", color=C_GPU, alpha=0.85)
        ax.set_xticks(x); ax.set_xticklabels([f"gid {i}" for i in gids], rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("Peak RSS per rank (MB)"); ax.set_title("Peak Memory per Block: CPU vs GPU")
        ax.legend(); ax.grid(axis="y", alpha=0.3)
        save(fig, "G4_memory.png")

    # ── G5: GPU seed scalability ──────────────────────────────────────────────
    g16 = valid[(valid.device == "gpu") & (valid.mesh == mesh_r) & (valid.n_blocks == nb_r)].sort_values("n_seeds")
    if len(g16) >= 2 and "total_steps" in g16.columns:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))
        axes[0].plot(g16.n_seeds, g16.wall / 3600, "o-", color=C_GPU, lw=2, markersize=8)
        axes[0].set_xscale("log"); axes[0].set_xlabel("Number of seeds")
        axes[0].set_ylabel("Wall time (hours)"); axes[0].set_title("GPU Wall vs Seeds")
        axes[0].grid(True, alpha=0.3)
        tput = series(g16, "total_steps") / np.where(g16.wall > 0, g16.wall, np.nan) / 1e6
        axes[1].plot(g16.n_seeds, tput, "s-", color=C_GPU, lw=2, markersize=8)
        axes[1].set_xscale("log"); axes[1].set_xlabel("Number of seeds")
        axes[1].set_ylabel("Throughput (M steps/s)"); axes[1].set_title("GPU Throughput vs Seeds")
        axes[1].grid(True, alpha=0.3)
        plt.tight_layout()
        save(fig, "G5_gpu_seed_scalability.png")


# ── run ────────────────────────────────────────────────────────────────────────
print("\n── Section 1: per-run load-imbalance diagnostics ──")
for (nb, ns, rid, dev, mesh), grp in blocks.groupby(KEYS):
    print(f"  run: {dev} {mesh} blocks={nb} seeds={ns} run={rid}")
    plot_per_run(grp, nb, ns, rid, dev, mesh)

print("\n── Section 2: scaling (CPU) ──");  plot_scaling(valid, "cpu")
print("\n── Section 2: scaling (GPU) ──");  plot_scaling(valid, "gpu")

plot_gpu_vs_cpu(valid, blocks)

print(f"\nDone. All plots saved to {PLOTS_DIR}/")
