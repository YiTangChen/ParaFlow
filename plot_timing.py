#!/usr/bin/env python3
"""
Plot ParaFlow timing and memory results from CSV files.
Outputs: results/plots/*.png

Section A — single-run diagnostics (one (n_blocks, n_seeds) config):
  A1  Compute vs Exposed Wait per block  — where observed time goes per block
  A2  Load imbalance sorted              — which ranks are hot and by how much
  A3  RK4 steps per block                — explains WHY imbalance exists
  A4  Phase Gantt (approximate timeline) — visual timeline per rank

Section B — scaling analysis (varies n_blocks):
  B1  Strong scaling wall-clock time     — overall runtime vs ideal
  B2  Parallel efficiency                — quantifies scaling quality
  B3  Load imbalance ratio vs scale      — imbalance grows with P
  B4  Iexchange breakdown vs scale       — compute / local message work / exposed wait
  B5  Memory per block vs scale          — data memory ∝ 1/n (positive result)
  B6  Speedup curve S(P)                 — T_base / T(P) vs ideal
  B7  Compute vs wall-clock speedup      — gap = exposed async overhead
  B8  Exposed wait fraction vs scale     — unhidden wait/progress grows with P
"""

import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path

# ── config ────────────────────────────────────────────────────────────────────
RESULTS_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("results")
PLOTS_DIR   = RESULTS_DIR / "plots"
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

SUMMARY_CSV = RESULTS_DIR / "run_summary.csv"
BLOCK_CSV   = RESULTS_DIR / "block_detail.csv"

# Colors — consistent across all plots.
# Exposed wait is orange-red: it is time not hidden by asynchronous overlap.
C_COMPUTE = "#2196F3"   # blue
C_IDLE    = "#FF5722"   # deep orange — wasted time, make it visible
C_LOAD    = "#4CAF50"   # green
C_WRITE   = "#9C27B0"   # purple
C_GRID    = "#00BCD4"   # teal
C_SOL     = "#F44336"   # red
C_IDEAL   = "#9E9E9E"   # grey — reference/ideal lines

SEED_COLORS = ["#1565C0", "#2E7D32", "#6A1B9A"]   # blue / green / purple per seed count


def save(fig, name):
    path = PLOTS_DIR / name
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved: {path}")


# ── load data ─────────────────────────────────────────────────────────────────
summary = pd.read_csv(SUMMARY_CSV)
blocks  = pd.read_csv(BLOCK_CSV)

for col in ("t_enqueue_avg", "t_dequeue_local_avg", "t_fill_incoming_avg",
            "t_iex_barrier_avg", "t_residual_avg"):
    if col not in summary:
        summary[col] = 0.0
# Accept legacy column name `t_comm_avg` from CSVs produced before the
# t_trace_comm → t_dequeue_local rename.
if "t_comm_avg" in summary and (summary["t_dequeue_local_avg"] == 0).all():
    summary["t_dequeue_local_avg"] = summary["t_comm_avg"]
for col in ("t_trace_enqueue", "t_dequeue_local", "t_fill_incoming",
            "t_trace_idle"):
    if col not in blocks:
        blocks[col] = 0.0
if "t_trace_comm" in blocks and (blocks["t_dequeue_local"] == 0).all():
    blocks["t_dequeue_local"] = blocks["t_trace_comm"]

# Backward-compatible derived columns.  DIY iexchange is asynchronous, so these
# are exposed-time categories, not a full physical network-transfer timeline.
# Iexchange timing excludes block load, output write, and post-iexchange barrier.
if "t_iex_with_barrier_avg" not in summary:
    summary["t_iex_with_barrier_avg"] = summary["t_iex_avg"] + summary["t_iex_barrier_avg"]
if "t_msg_overhead_avg" not in summary:
    summary["t_msg_overhead_avg"] = summary["t_enqueue_avg"] + summary["t_dequeue_local_avg"]
if "t_wait_sync_avg" not in summary:
    summary["t_wait_sync_avg"] = summary["t_fill_incoming_avg"]
if (summary["t_residual_avg"] == 0).all():
    summary["t_residual_avg"] = (
        summary["t_iex_avg"]
        - summary["t_compute_avg"]
        - summary["t_enqueue_avg"]
        - summary["t_dequeue_local_avg"]
        - summary["t_fill_incoming_avg"]
    ).clip(lower=0)

summary["t_idle_avg"] = summary["t_wait_sync_avg"] + summary["t_residual_avg"]

summary["idle_fraction"] = (
    summary["t_idle_avg"] / summary["t_iex_avg"].replace(0, np.nan)
)

if "t_msg_overhead" not in blocks:
    blocks["t_msg_overhead"] = blocks["t_trace_enqueue"] + blocks["t_dequeue_local"]
if "t_wait_sync" not in blocks:
    blocks["t_wait_sync"] = blocks["t_fill_incoming"]
if (blocks["t_wait_sync"] == 0).all() and (blocks["t_trace_idle"] > 0).any():
    # Old CSVs used t_trace_idle for all unaccounted/exposed overhead.
    blocks["t_wait_sync"] = blocks["t_trace_idle"]

# Keep only runs with actual data (t_iex_avg > 0)
valid = summary[summary["t_iex_avg"] > 0].copy()

print(f"Loaded {len(summary)} rows ({len(valid)} valid runs), {len(blocks)} block rows.")
print(f"n_blocks: {sorted(valid.n_blocks.unique())}")
print(f"n_seeds:  {sorted(valid.n_seeds.unique())}")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION A — per-run single-case plots
# ═══════════════════════════════════════════════════════════════════════════════

def plot_per_run(df_blocks, n_blocks, n_seeds, run_id):
    tag   = f"b{n_blocks:03d}_s{n_seeds:07d}"
    title = f"{n_blocks} blocks · {n_seeds:,} seeds"
    df    = df_blocks.sort_values("gid").reset_index(drop=True)
    x     = df.gid.values
    n     = len(df)
    w     = max(8, n * 0.55 + 2)   # figure width scales with number of blocks

    # ── A1: Compute vs exposed wait/progress per block ────────────────────────
    # Compute = useful work.  Exposed wait/progress is mostly fill_incoming().
    fig, ax = plt.subplots(figsize=(w, 4))
    ax.bar(x, df.t_trace_compute, label="Compute (RK4 integration)", color=C_COMPUTE, zorder=3)
    ax.bar(x, df.t_wait_sync,
           bottom=df.t_trace_compute,
           label="Exposed wait/progress", color=C_IDLE, alpha=0.88, zorder=3)
    residual_bottom = df.t_trace_compute + df.t_wait_sync
    ax.bar(x, df.t_trace_idle,
           bottom=residual_bottom,
           label="Residual uninstrumented", color="#795548", alpha=0.70, zorder=3)

    hot_idx = df.t_trace_compute.idxmax()
    hot     = df.loc[hot_idx]
    ax.annotate(
        f"gid {int(hot.gid)}\n{hot.t_trace_compute:.0f} s",
        xy=(hot.gid, hot.t_trace_compute),
        xytext=(0, 10), textcoords="offset points",
        ha="center", fontsize=8, color="darkblue",
        arrowprops=dict(arrowstyle="->", color="darkblue", lw=0.8),
    )
    ax.set_xlabel("Block gid")
    ax.set_ylabel("Time (s)")
    ax.set_title(f"Compute vs Exposed Wait per Block — {title}")
    ax.legend(loc="upper right")
    ax.set_xticks(x)
    ax.grid(axis="y", alpha=0.3, zorder=0)
    save(fig, f"{tag}_A1_compute_vs_wait.png")

    # ── A2: Load imbalance — sorted by compute time ───────────────────────────
    # Purpose: quantify how unequal the workload distribution is.
    # Sorted descending so the imbalance shape is immediately obvious.
    df_s    = df.sort_values("t_trace_compute", ascending=False).reset_index(drop=True)
    avg     = df.t_trace_compute.mean()
    mx      = df.t_trace_compute.max()
    mn      = df.t_trace_compute.min()
    imb_ma  = mx / avg
    imb_mm  = mx / mn if mn > 0 else float("inf")

    fig, ax = plt.subplots(figsize=(w, 4))
    bar_colors = [C_COMPUTE if t >= avg else "#90CAF9" for t in df_s.t_trace_compute]
    ax.bar(range(n), df_s.t_trace_compute, color=bar_colors, zorder=3)
    ax.axhline(avg, color="red",   linestyle="--", linewidth=1.5, label=f"avg = {avg:.0f} s")
    ax.axhline(mn,  color="green", linestyle=":",  linewidth=1.2, label=f"min = {mn:.0f} s")

    ax.set_xticks(range(n))
    ax.set_xticklabels([f"gid {int(g)}" for g in df_s.gid],
                       rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("Compute time (s)")
    ax.set_title(
        f"Load Imbalance (sorted by compute time) — {title}\n"
        f"max/avg = {imb_ma:.2f}×   ·   max/min = {imb_mm:.2f}×"
    )
    ax.legend()
    ax.grid(axis="y", alpha=0.3, zorder=0)
    save(fig, f"{tag}_A2_load_imbalance.png")

    # ── A3: RK4 steps per block ───────────────────────────────────────────────
    # Purpose: explain WHY imbalance exists.
    # More steps = more ocean particles integrated = harder domain geometry.
    fig, ax = plt.subplots(figsize=(w, 4))
    steps_M = df.n_steps_total / 1e6
    ax.bar(x, steps_M, color=C_COMPUTE, alpha=0.85, zorder=3)
    ax.set_xlabel("Block gid")
    ax.set_ylabel("RK4 steps (millions)")
    ax.set_title(f"Total Integration Steps per Block — {title}")
                #  f"(unequal steps → unequal compute time → load imbalance)")
    ax.set_xticks(x)
    ymax = steps_M.max()
    for _, row in df.iterrows():
        ax.text(row.gid, row.n_steps_total / 1e6 + ymax * 0.012,
                f"{int(row.n_seeds_initial)}s", ha="center", fontsize=6, color="#444")
    ax.annotate("label = initial seeds assigned to block",
                xy=(0.98, 0.97), xycoords="axes fraction",
                ha="right", fontsize=8, color="grey")
    ax.grid(axis="y", alpha=0.3, zorder=0)
    save(fig, f"{tag}_A3_steps_per_block.png")

    # ── A4: Phase Gantt — approximate timeline per rank ───────────────────────
    # Purpose: visually show that ranks finish at different times.
    # Phases stacked sequentially as approximation (no absolute timestamps).
    phases = [
        ("t_block_load",    "Block load",  C_LOAD),
        ("t_trace_compute", "Compute",     C_COMPUTE),
        ("t_msg_overhead",  "Local message work", "#607D8B"),
        ("t_wait_sync",     "Exposed wait/progress", C_IDLE),
        ("t_trace_idle",    "Residual",    "#795548"),
        ("t_output_write",  "Write",       C_WRITE),
    ]
    fig, ax = plt.subplots(figsize=(10, max(4, n * 0.30 + 1.5)))
    for _, row in df.iterrows():
        start = 0.0
        for col, _label, color in phases:
            dur = row[col]
            if dur > 1e-4:
                ax.barh(row.gid, dur, left=start, height=0.6, color=color, alpha=0.90)
            start += dur
    ax.set_yticks(df.gid.values)
    ax.set_yticklabels([f"rank {int(r)}" for r in df["rank"].values], fontsize=8)
    ax.set_xlabel("Approximate elapsed time (s)")
    ax.set_title(
        f"Phase Timeline per Rank — {title}"
        # f"(sequential phase approximation; comm < 1 ms omitted)"
    )
    legend_patches = [mpatches.Patch(color=c, label=l) for _, l, c in phases]
    ax.legend(handles=legend_patches, loc="lower right", fontsize=8)
    ax.grid(axis="x", alpha=0.3)
    save(fig, f"{tag}_A4_gantt_timeline.png")

    print(f"  → imbalance  max/avg = {imb_ma:.2f}×   max/min = {imb_mm:.2f}×  "
          f"  hottest = gid {int(hot.gid)}")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION B — scaling plots
# ═══════════════════════════════════════════════════════════════════════════════

def plot_scaling(valid):
    seed_counts  = sorted(valid.n_seeds.unique())
    block_counts = sorted(valid.n_blocks.unique())

    if len(block_counts) < 2:
        print("  [B] Only 1 n_blocks value — scaling plots skipped (need ≥2 runs)")
        return

    sc = SEED_COLORS[:len(seed_counts)]   # one colour per seed count

    # ── B1: Strong scaling — wall-clock time ──────────────────────────────────
    # Wall time = t_iex_max (max across ranks = the slowest rank determines runtime).
    # Each seed count gets its own ideal-scaling dashed line.
    fig, ax = plt.subplots(figsize=(7, 5))
    for i, ns in enumerate(seed_counts):
        df = valid[valid.n_seeds == ns].sort_values("n_blocks")
        if len(df) < 2:
            continue
        p0, t0 = df.iloc[0]["n_blocks"], df.iloc[0]["t_iex_max"]
        ax.plot(df.n_blocks, df.t_iex_max, "o-", color=sc[i], label=f"{ns:,} seeds")
        ax.plot(df.n_blocks, t0 * p0 / df.n_blocks.values,
                "--", alpha=0.25, color=sc[i])
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.set_xticks(block_counts)
    ax.set_xlabel("Number of blocks (MPI ranks)")
    ax.set_ylabel("Wall-clock max time (s)")
    ax.set_title("Strong Scaling\n(dashed = ideal linear scaling)")
    ax.legend()
    ax.grid(True, alpha=0.3, which="both")
    save(fig, "B1_strong_scaling.png")

    # ── B2: Parallel efficiency ───────────────────────────────────────────────
    # E(P) = T_base · P_base / (T(P) · P).  Ideal = 1.0.
    # Efficiency drop tells how much of the added hardware is wasted.
    fig, ax = plt.subplots(figsize=(7, 5))
    for i, ns in enumerate(seed_counts):
        df = valid[valid.n_seeds == ns].sort_values("n_blocks")
        if len(df) < 2:
            continue
        p0, t0 = df.iloc[0]["n_blocks"], df.iloc[0]["t_iex_max"]
        eff = (t0 * p0) / (df.t_iex_max * df.n_blocks)
        ax.plot(df.n_blocks, eff, "o-", color=sc[i], label=f"{ns:,} seeds")
    ax.axhline(1.0, color="black", linestyle="--", linewidth=1,   label="Ideal (100%)")
    ax.set_xscale("log", base=2)
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.set_xticks(block_counts)
    ax.set_ylim(0, 1.15)
    ax.set_xlabel("Number of blocks (MPI ranks)")
    ax.set_ylabel("Parallel efficiency")
    ax.set_title("Parallel Efficiency")
    ax.legend()
    ax.grid(True, alpha=0.3)
    save(fig, "B2_parallel_efficiency.png")

    # ── B3: Load imbalance ratio vs scale ─────────────────────────────────────
    # max/avg compute time across ranks.  Perfect balance = 1.0.
    # Key finding: ratio grows as more blocks → smaller, more unequal sub-domains.
    fig, ax = plt.subplots(figsize=(7, 5))
    for i, ns in enumerate(seed_counts):
        df = valid[valid.n_seeds == ns].sort_values("n_blocks")
        if len(df) < 2:
            continue
        imb = df.t_compute_max / df.t_compute_avg
        ax.plot(df.n_blocks, imb, "o-", color=sc[i], label=f"{ns:,} seeds")
        last = df.iloc[-1]
        ax.annotate(
            f"{last.t_compute_max / last.t_compute_avg:.1f}×",
            xy=(last.n_blocks, last.t_compute_max / last.t_compute_avg),
            xytext=(5, 0), textcoords="offset points",
            fontsize=8, color=sc[i],
        )
    ax.axhline(1.0, color="black", linestyle="--", linewidth=1, label="Perfect balance")
    ax.set_xscale("log", base=2)
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.set_xticks(block_counts)
    ax.set_xlabel("Number of blocks")
    ax.set_ylabel("max / avg compute time  (imbalance ratio)")
    ax.set_title("Load Imbalance Grows with Scale")
                #  "(larger ratio = hottest rank takes much longer than average)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    save(fig, "B3_load_imbalance_vs_scale.png")

    # ── B4: Iexchange exposed-time breakdown vs scale ─────────────────────────
    # Because DIY communication is asynchronous, hidden network transfer may
    # overlap with compute.  This plot shows exposed local work and waits only.
    for i, ns in enumerate(seed_counts):
        df = valid[valid.n_seeds == ns].sort_values("n_blocks")
        if len(df) < 2:
            continue
        fig, ax = plt.subplots(figsize=(8, 4))
        xpos   = np.arange(len(df))
        labels = [str(int(n)) for n in df.n_blocks]

        b0 = np.zeros(len(df))
        ax.bar(xpos, df.t_compute_avg.values,   bottom=b0,
               label="Compute",     color=C_COMPUTE)
        b1 = b0 + df.t_compute_avg.values
        ax.bar(xpos, df.t_msg_overhead_avg.values, bottom=b1,
               label="Local message work", color="#607D8B")
        b2 = b1 + df.t_msg_overhead_avg.values
        ax.bar(xpos, df.t_wait_sync_avg.values, bottom=b2,
               label="Exposed wait/progress", color=C_IDLE, alpha=0.88)
        b3 = b2 + df.t_wait_sync_avg.values
        ax.bar(xpos, df.t_residual_avg.values,  bottom=b3,
               label="Residual uninstrumented", color="#795548", alpha=0.70)

        ax.set_xticks(xpos)
        ax.set_xticklabels(labels)
        ax.set_xlabel("Number of blocks (MPI ranks)")
        ax.set_ylabel("Avg time per rank (s)")
        ax.set_title(
            f"Iexchange Exposed-Time Breakdown — {ns:,} seeds\n"
            f"(asynchronous transfer may overlap with compute)"
        )
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(axis="y", alpha=0.3, zorder=0)
        save(fig, f"B4_phase_breakdown_s{ns:07d}.png")

    # ── B5: Memory per block vs scale ─────────────────────────────────────────
    # Analytical data memory (grid + solution) scales exactly ∝ 1/n_blocks.
    # This is a positive result: no memory overhead from parallelism.
    # Peak RSS (VmHWM) includes OS + MPI overhead and decreases more slowly.
    for i, ns in enumerate(seed_counts):
        df = valid[valid.n_seeds == ns].sort_values("n_blocks")
        if len(df) < 2:
            continue
        grid_mb = df.mem_grid_bytes_avg    / 1e6
        sol_mb  = df.mem_solution_bytes_avg / 1e6
        peak_mb = df.mem_peak_vmhwm_kb_avg  / 1024   # KB → MB
        xpos    = np.arange(len(df))

        fig, ax = plt.subplots(figsize=(7, 4))
        ax.bar(xpos, grid_mb.values, label="Grid topology (analytical)", color=C_GRID)
        ax.bar(xpos, sol_mb.values,  bottom=grid_mb.values,
               label="Velocity solution (analytical)", color=C_SOL)
        ax2 = ax.twinx()
        ax2.plot(xpos, peak_mb.values, "k^--", markersize=6, linewidth=1.2,
                 label="Peak process RSS (VmHWM)")
        ax2.set_ylabel("Peak process RSS per rank (MB)")
        ax.set_xticks(xpos)
        ax.set_xticklabels([str(int(n)) for n in df.n_blocks])
        ax.set_xlabel("Number of blocks")
        ax.set_ylabel("Analytical data memory per block (MB)")
        ax.set_title(
            f"Memory per Block vs Scale — {ns:,} seeds\n"
            f"(left: data memory ∝ 1/n as expected;  right: OS+MPI overhead)"
        )
        l1, la1 = ax.get_legend_handles_labels()
        l2, la2 = ax2.get_legend_handles_labels()
        ax.legend(l1 + l2, la1 + la2, loc="upper right", fontsize=8)
        ax.grid(axis="y", alpha=0.3, zorder=0)
        save(fig, f"B5_memory_vs_scale_s{ns:07d}.png")

    # ── B6: Speedup curve  S(P) = T(P_base) / T(P) ───────────────────────────
    # The primary scaling result.  Both axes log-scale so ideal = straight line.
    fig, ax = plt.subplots(figsize=(7, 5))
    for i, ns in enumerate(seed_counts):
        df = valid[valid.n_seeds == ns].sort_values("n_blocks")
        if len(df) < 2:
            continue
        p0  = df.iloc[0]["n_blocks"]
        t0  = df.iloc[0]["t_iex_max"]
        spd = t0 / df.t_iex_max.values
        ax.plot(df.n_blocks, spd, "o-", color=sc[i], label=f"{ns:,} seeds")

    # ideal speedup relative to smallest block count in valid data
    p_all = np.array(sorted(valid.n_blocks.unique()))
    p_min = p_all[0]
    ax.plot(p_all, p_all / p_min, "--", color=C_IDEAL, linewidth=1.5, label="Ideal")

    ax.set_xscale("log", base=2)
    ax.set_yscale("log", base=2)
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.set_xticks(block_counts)
    ax.set_xlabel("Number of blocks (MPI ranks)")
    ax.set_ylabel("Speedup  S(P) = T(P₀) / T(P)")
    ax.set_title(f"Speedup Curve  (baseline = {p_min} blocks)")
    ax.legend()
    ax.grid(True, alpha=0.3, which="both")
    save(fig, "B6_speedup.png")

    # ── B7: Compute-only vs wall-clock speedup (one plot per seed count) ────────
    # This is the KEY finding plot.
    # Compute speedup ≈ ideal → the algorithm itself scales perfectly.
    # Wall-clock speedup << ideal → exposed async overhead is the bottleneck.
    # The shaded gap between the two lines is non-compute overhead.
    #
    # Generated per seed count so each has its own baseline (critical for 100k
    # seeds which only has data from n=32 onward — a shared baseline would mix
    # incompatible reference points).
    #
    # Additional insight across plots: larger problems show a smaller gap because
    # each rank has more compute work to hide exchange/progress overhead.
    for i, ns in enumerate(seed_counts):
        df = valid[(valid.n_seeds == ns) & (valid.t_compute_avg > 0)].sort_values("n_blocks")
        if len(df) < 2:
            continue
        p0      = df.iloc[0]["n_blocks"]   # baseline = smallest valid block count
        tc0     = df.iloc[0]["t_compute_avg"]
        tw0     = df.iloc[0]["t_iex_max"]
        sp_comp = tc0 / df.t_compute_avg.values
        sp_wall = tw0 / df.t_iex_max.values
        ideal   = df.n_blocks.values / p0

        fig, ax = plt.subplots(figsize=(7, 5))
        ax.plot(df.n_blocks, ideal,   "--", color=C_IDEAL,   linewidth=1.5, label="Ideal")
        ax.plot(df.n_blocks, sp_comp, "o-", color=C_COMPUTE, linewidth=2,
                label="Compute-only speedup")
        ax.plot(df.n_blocks, sp_wall, "s-", color=C_IDLE,    linewidth=2,
                label="Wall-clock speedup")
        ax.fill_between(df.n_blocks, sp_wall, sp_comp,
                        alpha=0.18, color=C_IDLE, label="Gap = exposed overhead")

        ax.set_xscale("log", base=2)
        ax.set_yscale("log", base=2)
        ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
        ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
        ax.set_xticks([n for n in block_counts if n >= p0])
        ax.set_xlabel("Number of blocks (MPI ranks)")
        ax.set_ylabel("Speedup")
        ax.set_title(
            f"Compute vs Wall-Clock Speedup — {ns:,} seeds  (baseline = {p0} blocks)\n"
            f"(gap = exposed wait/progress and residual overhead)"
        )
        ax.legend()
        ax.grid(True, alpha=0.3, which="both")
        print(f"  [B7] n_seeds={ns:,}: baseline={p0} blocks, "
              f"compute S={sp_comp[-1]:.1f}×  wall S={sp_wall[-1]:.1f}×  "
              f"gap={sp_comp[-1]-sp_wall[-1]:.1f}×")
        save(fig, f"B7_compute_vs_total_speedup_s{ns:07d}.png")

    # ── B8: Exposed wait fraction vs scale ────────────────────────────────────
    # Shows the fraction of pure iexchange time not hidden by async overlap.
    # Directly links B3 (imbalance) to B2 (efficiency loss).
    fig, ax = plt.subplots(figsize=(7, 5))
    for i, ns in enumerate(seed_counts):
        df = valid[valid.n_seeds == ns].sort_values("n_blocks")
        if len(df) < 2:
            continue
        wait_pct = df.idle_fraction * 100
        ax.plot(df.n_blocks, wait_pct, "o-", color=sc[i], label=f"{ns:,} seeds")
        last = df.iloc[-1]
        ax.annotate(
            f"{last.idle_fraction * 100:.0f}%",
            xy=(last.n_blocks, last.idle_fraction * 100),
            xytext=(5, 0), textcoords="offset points",
            fontsize=8, color=sc[i],
        )
    ax.set_xscale("log", base=2)
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.set_xticks(block_counts)
    ax.set_ylim(0, 105)
    ax.set_xlabel("Number of blocks (MPI ranks)")
    ax.set_ylabel("Exposed overhead fraction (% of pure iexchange)")
    ax.set_title("Exposed Wait/Progress Fraction vs Scale\n"
                 "(fill_incoming + residual; barrier reported separately)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    save(fig, "B8_wait_fraction.png")


# ── run ────────────────────────────────────────────────────────────────────────
print("\n── Section A: per-block plots ──")
for (nb, ns, rid), grp in blocks.groupby(["n_blocks", "n_seeds", "run_id"]):
    print(f"  run: n_blocks={nb}  n_seeds={ns}  run_id={rid}")
    plot_per_run(grp, nb, ns, rid)

print("\n── Section B: scaling plots ──")
plot_scaling(valid)

print(f"\nDone. All plots saved to {PLOTS_DIR}/")
