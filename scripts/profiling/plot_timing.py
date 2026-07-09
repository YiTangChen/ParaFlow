#!/usr/bin/env python3
"""
Plot ParaFlow timing and memory results from CSV files.
Outputs: results/plots/*.png

Section A — single-run diagnostics (one (n_blocks, n_seeds) config):
  A1  Compute vs Idle per block
  A2  Load imbalance sorted
  A3  RK4 steps per block
  A4  Phase Gantt (approximate timeline)
  A5  Imbalance diagnosis: work distribution + work-vs-throughput scatter

Section B — scaling analysis (varies n_blocks, same device):
  B1  Strong scaling wall-clock time
  B2  Parallel efficiency
  B3  Load imbalance ratio vs scale
  B4  Phase breakdown vs scale
  B5  Memory per block vs scale
  B6  Speedup curve S(P)
  B7  Compute vs wall-clock speedup
  B8  Idle fraction vs scale

Section C — GPU vs CPU comparison:
  C1  Wall-clock time: CPU vs GPU (same 16 blocks, 10k seeds, lowres)
  C2  Phase breakdown comparison: CPU vs GPU stacked bars
  C3  GPU speedup: overall and compute-only
  C4  GPU seed scalability (1k / 10k / 100k)
  C5  GPU kernel fraction of compute time
  C6  Memory: CPU vs GPU peak RSS per block
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
RESULTS_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("results/gpu_analysis")
PLOTS_DIR   = RESULTS_DIR / "plots"
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

SUMMARY_CSV = RESULTS_DIR / "run_summary.csv"
BLOCK_CSV   = RESULTS_DIR / "block_detail.csv"

C_COMPUTE = "#2196F3"
C_IDLE    = "#FF5722"
C_LOAD    = "#4CAF50"
C_WRITE   = "#9C27B0"
C_GRID    = "#00BCD4"
C_SOL     = "#F44336"
C_IDEAL   = "#9E9E9E"
C_GPU     = "#F57C00"   # orange for GPU
C_CPU     = "#1565C0"   # dark blue for CPU
C_KERNEL  = "#FF9800"   # amber for GPU kernel overhead

SEED_COLORS = ["#1565C0", "#2E7D32", "#6A1B9A"]


def save(fig, name):
    path = PLOTS_DIR / name
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved: {path}")


# ── load data ─────────────────────────────────────────────────────────────────
summary = pd.read_csv(SUMMARY_CSV)
blocks  = pd.read_csv(BLOCK_CSV)

# Ensure device/mesh columns exist (backwards compat with old CSVs)
if "device" not in summary.columns:
    summary["device"] = "cpu"
if "mesh" not in summary.columns:
    summary["mesh"] = "lowres"
if "device" not in blocks.columns:
    blocks["device"] = "cpu"
if "mesh" not in blocks.columns:
    blocks["mesh"] = "lowres"

summary["t_idle_avg"] = (
    summary["t_iex_avg"]
    - summary["t_blockload_avg"]
    - summary["t_compute_avg"]
    - summary["t_comm_avg"]
    - summary["t_write_avg"]
).clip(lower=0)

summary["idle_fraction"] = (
    summary["t_idle_avg"] / summary["t_iex_avg"].replace(0, np.nan)
)

valid = summary[summary["t_iex_avg"] > 0].copy()

print(f"Loaded {len(summary)} rows ({len(valid)} valid runs), {len(blocks)} block rows.")
print(f"Devices:  {sorted(valid.device.unique())}")
print(f"Meshes:   {sorted(valid.mesh.unique())}")
print(f"n_blocks: {sorted(valid.n_blocks.unique())}")
print(f"n_seeds:  {sorted(valid.n_seeds.unique())}")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION A — per-run single-case plots
# ═══════════════════════════════════════════════════════════════════════════════

def plot_per_run(df_blocks, n_blocks, n_seeds, run_id, device, mesh):
    tag   = f"{device}_{mesh}_b{n_blocks:03d}_s{n_seeds:07d}"
    title = f"{device.upper()} {mesh} · {n_blocks} blocks · {n_seeds:,} seeds"
    df    = df_blocks.sort_values("gid").reset_index(drop=True)
    x     = df.gid.values
    n     = len(df)
    w     = max(8, n * 0.55 + 2)

    # ── A1: Compute vs Idle per block ─────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(w, 4))
    ax.bar(x, df.t_trace_compute, label="Compute", color=C_COMPUTE, zorder=3)
    ax.bar(x, df.t_trace_idle,
           bottom=df.t_trace_compute,
           label="Idle (waiting for slowest rank)", color=C_IDLE, alpha=0.88, zorder=3)

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
    ax.set_title(f"Compute vs Idle per Block — {title}")
    ax.legend(loc="upper right")
    ax.set_xticks(x)
    ax.grid(axis="y", alpha=0.3, zorder=0)
    save(fig, f"{tag}_A1_compute_vs_idle.png")

    # ── A2: Load imbalance ────────────────────────────────────────────────────
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
        f"Load Imbalance (sorted) — {title}\n"
        f"max/avg = {imb_ma:.2f}×   ·   max/min = {imb_mm:.2f}×"
    )
    ax.legend()
    ax.grid(axis="y", alpha=0.3, zorder=0)
    save(fig, f"{tag}_A2_load_imbalance.png")

    # ── A3: RK4 steps per block ───────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(w, 4))
    steps_M = df.n_steps_total / 1e6
    ax.bar(x, steps_M, color=C_COMPUTE, alpha=0.85, zorder=3)
    ax.set_xlabel("Block gid")
    ax.set_ylabel("RK4 steps (millions)")
    ax.set_title(f"Total Integration Steps per Block — {title}")
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

    # ── A4: Phase Gantt ───────────────────────────────────────────────────────
    phases = [
        ("t_block_load",    "Block load",  C_LOAD),
        ("t_trace_compute", "Compute",     C_COMPUTE),
        ("t_output_write",  "Write",       C_WRITE),
        ("t_trace_idle",    "Idle",        C_IDLE),
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
    ax.set_title(f"Phase Timeline per Rank — {title}")
    legend_patches = [mpatches.Patch(color=c, label=l) for _, l, c in phases]
    ax.legend(handles=legend_patches, loc="lower right", fontsize=8)
    ax.grid(axis="x", alpha=0.3)
    save(fig, f"{tag}_A4_gantt_timeline.png")

    # ── A5: Imbalance diagnosis — work distribution + work-vs-throughput ───────
    # The scatter is the key: a tight line ⇒ throughput is uniform, so imbalance is
    # pure work distribution (fixable by weighting the partition by predicted steps).
    # A cloud ⇒ throughput varies ⇒ hardware (NUMA / GPU contention), a different fix.
    steps = df.n_steps_total.astype(float).values
    walls = df.t_trace_local_wall.astype(float).values
    ok    = walls > 0
    tput  = np.zeros_like(steps)
    tput[ok] = steps[ok] / walls[ok] / 1e6                 # M steps/s
    work_avg = steps.mean() if len(steps) else 0.0
    work_moa = (steps.max() / work_avg) if work_avg > 0 else 0.0
    tmean = tput[ok].mean() if ok.any() else 0.0
    tcov  = (tput[ok].std() / tmean) if tmean > 0 else 0.0

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 4.5))
    order = np.argsort(-steps)
    axL.bar(range(n), steps[order] / 1e6, color=C_COMPUTE, zorder=3)
    axL.axhline(work_avg / 1e6, color="red", ls="--", lw=1.4, label=f"avg = {work_avg/1e6:.2f} M")
    axL.set_xlabel("Block (sorted by work)")
    axL.set_ylabel("RK4 steps (millions)")
    axL.set_title(f"Work distribution — max/avg = {work_moa:.2f}×")
    axL.legend(); axL.grid(axis="y", alpha=0.3, zorder=0)

    axR.scatter(steps / 1e6, walls, s=40, color=C_COMPUTE, zorder=3)
    if tmean > 0:
        xs = np.array([0.0, steps.max() / 1e6])
        axR.plot(xs, xs / tmean, color=C_IDEAL, ls="--", lw=1.3,
                 label=f"uniform throughput ({tmean:.1f} M steps/s)")
    axR.set_xlabel("Work (M steps)")
    axR.set_ylabel("trace_local_wall (s)")
    verdict = ("work-bound → re-partition by predicted steps"
               if tcov < 0.15 else "throughput varies → check NUMA / GPU contention")
    axR.set_title(f"Work vs time — throughput CoV = {tcov:.2f}\n{verdict}")
    axR.legend(fontsize=8); axR.grid(alpha=0.3, zorder=0)

    fig.suptitle(f"Imbalance Diagnosis — {title}")
    plt.tight_layout()
    save(fig, f"{tag}_A5_imbalance_diagnosis.png")

    print(f"  → imbalance  max/avg = {imb_ma:.2f}×   max/min = {imb_mm:.2f}×  "
          f"  hottest = gid {int(hot.gid)}  throughput_CoV = {tcov:.2f}")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION B — scaling plots (per device)
# ═══════════════════════════════════════════════════════════════════════════════

def plot_scaling(valid, device_filter=None):
    if device_filter:
        valid = valid[valid.device == device_filter]
        suffix = f"_{device_filter}"
    else:
        suffix = ""

    seed_counts  = sorted(valid.n_seeds.unique())
    block_counts = sorted(valid.n_blocks.unique())

    if len(block_counts) < 2:
        print(f"  [B{suffix}] Only 1 n_blocks value — scaling plots skipped (need ≥2 runs)")
        return

    sc = SEED_COLORS[:len(seed_counts)]

    # ── B1: Strong scaling ────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 5))
    for i, ns in enumerate(seed_counts):
        df = valid[valid.n_seeds == ns].sort_values("n_blocks")
        if len(df) < 2:
            continue
        p0, t0 = df.iloc[0]["n_blocks"], df.iloc[0]["t_iex_max"]
        ax.plot(df.n_blocks, df.t_iex_max, "o-", color=sc[i], label=f"{ns:,} seeds")
        ax.plot(df.n_blocks, t0 * p0 / df.n_blocks.values, "--", alpha=0.25, color=sc[i])
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.set_xticks(block_counts)
    ax.set_xlabel("Number of blocks (MPI ranks)")
    ax.set_ylabel("Wall-clock max time (s)")
    ax.set_title(f"Strong Scaling{' — ' + device_filter.upper() if device_filter else ''}\n"
                 "(dashed = ideal)")
    ax.legend()
    ax.grid(True, alpha=0.3, which="both")
    save(fig, f"B1_strong_scaling{suffix}.png")

    # ── B2: Parallel efficiency ───────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 5))
    for i, ns in enumerate(seed_counts):
        df = valid[valid.n_seeds == ns].sort_values("n_blocks")
        if len(df) < 2:
            continue
        p0, t0 = df.iloc[0]["n_blocks"], df.iloc[0]["t_iex_max"]
        eff = (t0 * p0) / (df.t_iex_max * df.n_blocks)
        ax.plot(df.n_blocks, eff, "o-", color=sc[i], label=f"{ns:,} seeds")
    ax.axhline(1.0, color="black", linestyle="--", linewidth=1, label="Ideal (100%)")
    ax.axhline(0.7, color="red",   linestyle=":",  linewidth=1.2, label="70% threshold")
    ax.set_xscale("log", base=2)
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.set_xticks(block_counts)
    ax.set_ylim(0, 1.15)
    ax.set_xlabel("Number of blocks (MPI ranks)")
    ax.set_ylabel("Parallel efficiency")
    ax.set_title(f"Parallel Efficiency{' — ' + device_filter.upper() if device_filter else ''}")
    ax.legend()
    ax.grid(True, alpha=0.3)
    save(fig, f"B2_parallel_efficiency{suffix}.png")

    # ── B3: Load imbalance vs scale ───────────────────────────────────────────
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
    ax.set_title(f"Load Imbalance vs Scale{' — ' + device_filter.upper() if device_filter else ''}")
    ax.legend()
    ax.grid(True, alpha=0.3)
    save(fig, f"B3_load_imbalance_vs_scale{suffix}.png")

    # ── B4: Phase breakdown ───────────────────────────────────────────────────
    for i, ns in enumerate(seed_counts):
        df = valid[valid.n_seeds == ns].sort_values("n_blocks")
        if len(df) < 2:
            continue
        fig, ax = plt.subplots(figsize=(8, 4))
        xpos   = np.arange(len(df))
        labels = [str(int(n)) for n in df.n_blocks]

        b0 = np.zeros(len(df))
        ax.bar(xpos, df.t_blockload_avg.values, bottom=b0, label="Block load", color=C_LOAD)
        b1 = b0 + df.t_blockload_avg.values
        ax.bar(xpos, df.t_compute_avg.values,   bottom=b1, label="Compute",    color=C_COMPUTE)
        b2 = b1 + df.t_compute_avg.values
        ax.bar(xpos, df.t_write_avg.values,     bottom=b2, label="Write",      color=C_WRITE)
        b3 = b2 + df.t_write_avg.values
        ax.bar(xpos, df.t_idle_avg.values,      bottom=b3,
               label="Idle (imbalance overhead)", color=C_IDLE, alpha=0.88)

        ax.set_xticks(xpos)
        ax.set_xticklabels(labels)
        ax.set_xlabel("Number of blocks (MPI ranks)")
        ax.set_ylabel("Avg time per rank (s)")
        ax.set_title(
            f"Phase Breakdown vs Scale — {ns:,} seeds"
            + (f" — {device_filter.upper()}" if device_filter else "")
        )
        ax.legend(loc="upper right", fontsize=8)
        ax.grid(axis="y", alpha=0.3, zorder=0)
        save(fig, f"B4_phase_breakdown_s{ns:07d}{suffix}.png")

    # ── B5: Memory per block vs scale ─────────────────────────────────────────
    for i, ns in enumerate(seed_counts):
        df = valid[valid.n_seeds == ns].sort_values("n_blocks")
        if len(df) < 2:
            continue
        grid_mb = df.mem_grid_bytes_avg    / 1e6
        sol_mb  = df.mem_solution_bytes_avg / 1e6
        peak_mb = df.mem_peak_vmhwm_kb_avg  / 1024
        xpos    = np.arange(len(df))

        fig, ax = plt.subplots(figsize=(7, 4))
        ax.bar(xpos, grid_mb.values, label="Grid topology", color=C_GRID)
        ax.bar(xpos, sol_mb.values, bottom=grid_mb.values,
               label="Velocity solution", color=C_SOL)
        ax2 = ax.twinx()
        ax2.plot(xpos, peak_mb.values, "k^--", markersize=6, linewidth=1.2,
                 label="Peak RSS (VmHWM)")
        ax2.set_ylabel("Peak RSS per rank (MB)")
        ax.set_xticks(xpos)
        ax.set_xticklabels([str(int(n)) for n in df.n_blocks])
        ax.set_xlabel("Number of blocks")
        ax.set_ylabel("Analytical data memory per block (MB)")
        ax.set_title(f"Memory per Block vs Scale — {ns:,} seeds"
                     + (f" — {device_filter.upper()}" if device_filter else ""))
        l1, la1 = ax.get_legend_handles_labels()
        l2, la2 = ax2.get_legend_handles_labels()
        ax.legend(l1 + l2, la1 + la2, loc="upper right", fontsize=8)
        ax.grid(axis="y", alpha=0.3, zorder=0)
        save(fig, f"B5_memory_vs_scale_s{ns:07d}{suffix}.png")

    # ── B6: Speedup ───────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 5))
    for i, ns in enumerate(seed_counts):
        df = valid[valid.n_seeds == ns].sort_values("n_blocks")
        if len(df) < 2:
            continue
        p0  = df.iloc[0]["n_blocks"]
        t0  = df.iloc[0]["t_iex_max"]
        spd = t0 / df.t_iex_max.values
        ax.plot(df.n_blocks, spd, "o-", color=sc[i], label=f"{ns:,} seeds")

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
    ax.set_title(f"Speedup Curve{' — ' + device_filter.upper() if device_filter else ''}")
    ax.legend()
    ax.grid(True, alpha=0.3, which="both")
    save(fig, f"B6_speedup{suffix}.png")

    # ── B7: Compute vs wall speedup ───────────────────────────────────────────
    for i, ns in enumerate(seed_counts):
        df = valid[(valid.n_seeds == ns) & (valid.t_compute_avg > 0)].sort_values("n_blocks")
        if len(df) < 2:
            continue
        p0      = df.iloc[0]["n_blocks"]
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
                        alpha=0.18, color=C_IDLE, label="Gap = idle overhead")

        ax.set_xscale("log", base=2)
        ax.set_yscale("log", base=2)
        ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
        ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
        ax.set_xticks([n for n in block_counts if n >= p0])
        ax.set_xlabel("Number of blocks (MPI ranks)")
        ax.set_ylabel("Speedup")
        ax.set_title(
            f"Compute vs Wall-Clock Speedup — {ns:,} seeds  (baseline = {p0} blocks)"
            + (f" — {device_filter.upper()}" if device_filter else "")
        )
        ax.legend()
        ax.grid(True, alpha=0.3, which="both")
        save(fig, f"B7_compute_vs_total_speedup_s{ns:07d}{suffix}.png")

    # ── B8: Idle fraction ─────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 5))
    for i, ns in enumerate(seed_counts):
        df = valid[valid.n_seeds == ns].sort_values("n_blocks")
        if len(df) < 2:
            continue
        idle_pct = df.idle_fraction * 100
        ax.plot(df.n_blocks, idle_pct, "o-", color=sc[i], label=f"{ns:,} seeds")
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
    ax.set_ylabel("Idle fraction (% of wall time)")
    ax.set_title(f"Idle Time Fraction vs Scale"
                 + (f" — {device_filter.upper()}" if device_filter else ""))
    ax.legend()
    ax.grid(True, alpha=0.3)
    save(fig, f"B8_idle_fraction{suffix}.png")


# ═══════════════════════════════════════════════════════════════════════════════
# SECTION C — GPU vs CPU comparison
# ═══════════════════════════════════════════════════════════════════════════════

def plot_gpu_vs_cpu(valid, blocks):
    """Direct GPU vs CPU comparison plots."""
    print("\n── Section C: GPU vs CPU comparison plots ──")

    # Shared reference: lowres, 16 blocks, 10k seeds
    ref = valid[(valid.mesh == "lowres") & (valid.n_blocks == 16) & (valid.n_seeds == 10000)]
    cpu_ref = ref[ref.device == "cpu"]
    gpu_ref = ref[ref.device == "gpu"]

    have_ref = len(cpu_ref) > 0 and len(gpu_ref) > 0

    # End-to-end driver wall (block_load + exchange + write) if available,
    # else fall back to the exchange-only wall from older CSVs.
    def _wall(row):
        r = float(row.get("t_run_max", 0) or 0)
        return r if r > 0 else float(row["t_iex_max"])

    # ── C1: Wall-clock time bar chart ─────────────────────────────────────────
    # Compare all runs grouped by device, at 16 blocks / 10k seeds
    if have_ref:
        cpu_wall = _wall(cpu_ref.iloc[0])
        gpu_wall = _wall(gpu_ref.iloc[0])
        speedup  = cpu_wall / gpu_wall

        fig, ax = plt.subplots(figsize=(6, 5))
        bars = ax.bar(["CPU (16 ranks)", "GPU (16 ranks)"],
                      [cpu_wall / 3600, gpu_wall / 3600],
                      color=[C_CPU, C_GPU], width=0.45)
        for bar, val in zip(bars, [cpu_wall, gpu_wall]):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.2,
                    f"{val/3600:.1f} h", ha="center", va="bottom", fontsize=11, fontweight="bold")
        ax.set_ylabel("End-to-end wall time (hours)")
        ax.set_title(
            f"End-to-End Wall Time: CPU vs GPU\n"
            f"(lowres · 16 blocks · 10,000 seeds · GPU speedup = {speedup:.1f}×)"
        )
        ax.grid(axis="y", alpha=0.3)
        save(fig, "C1_wall_time_cpu_vs_gpu.png")
        print(f"  [C1] CPU={cpu_wall/3600:.2f} h  GPU={gpu_wall/3600:.2f} h  speedup={speedup:.1f}×")

    # ── C2: Aligned phase breakdown, CPU vs GPU ───────────────────────────────
    # Both devices decompose into the SAME buckets so the bars are directly
    # comparable. transfer is 0 on CPU by construction; integrate is the unified
    # cpu+gpu bucket. Left panel = end-to-end (block_load dominates); right panel
    # = trace-only buckets zoomed in, to expose where the GPU win is spent.
    if have_ref:
        cpu_row = cpu_ref.iloc[0]
        gpu_row = gpu_ref.iloc[0]

        def g(row, col):
            try:
                return float(row.get(col, 0.0) or 0.0)
            except Exception:
                return 0.0

        trace_phases = [
            ("t_prepare_avg",     "Prepare (host setup)", "#8E24AA"),
            ("t_transfer_avg",    "Transfer (H2D/D2H)",   C_KERNEL),
            ("t_integrate_avg",   "Integrate",            C_COMPUTE),
            ("t_postprocess_avg", "Postprocess",          "#00897B"),
            ("t_enqueue_avg",     "Enqueue (MPI)",        "#546E7A"),
        ]
        e2e_phases = [("t_blockload_avg", "Block load", C_LOAD)] + trace_phases \
                     + [("t_output_write_avg", "Write", C_WRITE)]

        fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5))

        # Left: end-to-end (seconds)
        bc, bg = 0.0, 0.0
        for col, label, color in e2e_phases:
            cv, gv = g(cpu_row, col), g(gpu_row, col)
            axL.bar([0], [cv], bottom=[bc], color=color, label=label, width=0.45)
            axL.bar([1], [gv], bottom=[bg], color=color, width=0.45, alpha=0.85)
            bc += cv; bg += gv
        axL.set_xticks([0, 1]); axL.set_xticklabels(["CPU", "GPU"])
        axL.set_ylabel("Avg time per rank (s)")
        axL.set_title("End-to-end (block_load dominates)")
        axL.legend(loc="upper right", fontsize=8)
        axL.grid(axis="y", alpha=0.3)

        # Right: trace-only buckets, zoomed
        bc, bg = 0.0, 0.0
        for col, label, color in trace_phases:
            cv, gv = g(cpu_row, col), g(gpu_row, col)
            axR.bar([0], [cv], bottom=[bc], color=color, label=label, width=0.45)
            axR.bar([1], [gv], bottom=[bg], color=color, width=0.45, alpha=0.85)
            bc += cv; bg += gv
        axR.set_xticks([0, 1]); axR.set_xticklabels(["CPU", "GPU"])
        axR.set_ylabel("Avg time per rank (s)")
        axR.set_title("Trace buckets only (aligned)")
        axR.legend(loc="upper right", fontsize=8)
        axR.grid(axis="y", alpha=0.3)

        fig.suptitle("Aligned Phase Breakdown: CPU vs GPU  (lowres · 16 blocks · 10,000 seeds)")
        plt.tight_layout()
        save(fig, "C2_phase_breakdown_cpu_vs_gpu.png")

    # ── C3: GPU speedup bar (overall + compute-only) ──────────────────────────
    if have_ref:
        cpu_row   = cpu_ref.iloc[0]
        gpu_row   = gpu_ref.iloc[0]
        sp_wall   = float(cpu_row["t_iex_max"])   / float(gpu_row["t_iex_max"])
        sp_load   = float(cpu_row["t_blockload_avg"]) / float(gpu_row["t_blockload_avg"]) \
                    if float(gpu_row["t_blockload_avg"]) > 0 else 0
        sp_comp   = float(cpu_row["t_compute_avg"])   / float(gpu_row["t_compute_avg"]) \
                    if float(gpu_row["t_compute_avg"])   > 0 else 0
        sp_write  = float(cpu_row["t_write_avg"])     / float(gpu_row["t_write_avg"]) \
                    if float(gpu_row["t_write_avg"])     > 0 else 0

        metrics = ["Wall-clock", "Block load", "Compute", "Write"]
        values  = [sp_wall, sp_load, sp_comp, sp_write]
        bar_colors = [C_GPU, C_LOAD, C_COMPUTE, C_WRITE]

        fig, ax = plt.subplots(figsize=(7, 5))
        bars = ax.bar(metrics, values, color=bar_colors, width=0.5)
        for bar, v in zip(bars, values):
            ax.text(bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.1,
                    f"{v:.1f}×", ha="center", va="bottom", fontsize=11, fontweight="bold")
        ax.axhline(1.0, color="grey", linestyle="--", linewidth=1, label="1× (no speedup)")
        ax.set_ylabel("GPU speedup over CPU")
        ax.set_title("GPU Speedup by Phase\n(lowres · 16 blocks · 10,000 seeds)")
        ax.legend(fontsize=8)
        ax.grid(axis="y", alpha=0.3)
        save(fig, "C3_gpu_speedup_by_phase.png")
        print(f"  [C3] Speedups — wall={sp_wall:.1f}× load={sp_load:.1f}× "
              f"compute={sp_comp:.1f}× write={sp_write:.1f}×")

    # ── C4: GPU seed scalability (lowres, 16 blocks) ──────────────────────────
    gpu_16 = valid[(valid.device == "gpu") & (valid.mesh == "lowres") & (valid.n_blocks == 16)
                  ].sort_values("n_seeds")
    if len(gpu_16) >= 2:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4))

        # Wall time vs seeds
        ax = axes[0]
        ax.plot(gpu_16.n_seeds, gpu_16.t_iex_max / 3600, "o-", color=C_GPU, linewidth=2,
                markersize=8, label="GPU wall time")
        ax.set_xlabel("Number of seeds")
        ax.set_ylabel("Wall-clock time (hours)")
        ax.set_xscale("log")
        ax.set_title("GPU Wall Time vs Seed Count\n(lowres · 16 blocks)")
        ax.grid(True, alpha=0.3)
        for _, row in gpu_16.iterrows():
            ax.annotate(f"{row.t_iex_max/3600:.2f} h",
                        xy=(row.n_seeds, row.t_iex_max / 3600),
                        xytext=(5, 5), textcoords="offset points", fontsize=8)

        # Throughput: million steps / second
        ax = axes[1]
        # Total steps / wall time
        throughput = gpu_16.total_steps / gpu_16.t_iex_max / 1e6  # M steps/s
        ax.plot(gpu_16.n_seeds, throughput, "s-", color=C_GPU, linewidth=2,
                markersize=8, label="Throughput")
        ax.set_xlabel("Number of seeds")
        ax.set_ylabel("Integration throughput (M steps/s)")
        ax.set_xscale("log")
        ax.set_title("GPU Throughput vs Seed Count\n(lowres · 16 blocks)")
        ax.grid(True, alpha=0.3)
        for _, row in gpu_16.iterrows():
            tp = row.total_steps / row.t_iex_max / 1e6
            ax.annotate(f"{tp:.1f}",
                        xy=(row.n_seeds, tp),
                        xytext=(5, 5), textcoords="offset points", fontsize=8)

        plt.tight_layout()
        save(fig, "C4_gpu_seed_scalability.png")
        print(f"  [C4] GPU seed scalability: {list(gpu_16.n_seeds)} seeds "
              f"→ {list(gpu_16.t_iex_max.round(0).astype(int))} s wall")

    # ── C5: GPU kernel fraction of compute time ───────────────────────────────
    gpu_blocks = blocks[blocks.device == "gpu"].copy()
    if len(gpu_blocks) > 0 and "t_gpu_kernel" in gpu_blocks.columns:
        gpu_blocks = gpu_blocks[gpu_blocks.t_gpu_kernel > 0].copy()
    if len(gpu_blocks) > 0 and "t_gpu_kernel" in gpu_blocks.columns:
        gpu_blocks["gpu_frac"] = gpu_blocks["t_gpu_kernel"] / gpu_blocks["t_trace_compute"].replace(0, np.nan)

        seed_vals = sorted(gpu_blocks.n_seeds.unique())
        fig, ax = plt.subplots(figsize=(8, 4))
        for i, ns in enumerate(seed_vals):
            df = gpu_blocks[(gpu_blocks.n_seeds == ns) & (gpu_blocks.mesh == "lowres")
                           ].sort_values("gid")
            if len(df) == 0:
                continue
            ax.plot(df.gid, df.gpu_frac * 100, "o-",
                    label=f"{ns:,} seeds", color=SEED_COLORS[i % len(SEED_COLORS)],
                    markersize=5)

        ax.set_xlabel("Block gid")
        ax.set_ylabel("GPU kernel time / compute time (%)")
        ax.set_ylim(90, 101)
        ax.set_title("GPU Kernel Fraction of Compute Time\n"
                     "(≈100% means nearly all compute is on GPU)")
        ax.legend(fontsize=8)
        ax.grid(axis="y", alpha=0.3)
        save(fig, "C5_gpu_kernel_fraction.png")

    # ── C6: Memory: CPU vs GPU peak RSS ───────────────────────────────────────
    # Compare peak RSS for same domain (lowres, 16 blocks, 10k seeds)
    ref_blocks = blocks[
        (blocks.mesh == "lowres") & (blocks.n_blocks == 16) & (blocks.n_seeds == 10000)
    ]
    cpu_blk = ref_blocks[ref_blocks.device == "cpu"].sort_values("gid")
    gpu_blk = ref_blocks[ref_blocks.device == "gpu"].sort_values("gid")

    if len(cpu_blk) > 0 and len(gpu_blk) > 0:
        gids = sorted(set(cpu_blk.gid) & set(gpu_blk.gid))
        cpu_rss = [float(cpu_blk[cpu_blk.gid == g]["mem_peak_vmhwm_kb"].iloc[0]) / 1024
                   for g in gids]
        gpu_rss = [float(gpu_blk[gpu_blk.gid == g]["mem_peak_vmhwm_kb"].iloc[0]) / 1024
                   for g in gids]

        x = np.arange(len(gids))
        w = 0.35
        fig, ax = plt.subplots(figsize=(max(8, len(gids) * 0.6 + 2), 4))
        ax.bar(x - w/2, cpu_rss, width=w, label="CPU", color=C_CPU, alpha=0.85)
        ax.bar(x + w/2, gpu_rss, width=w, label="GPU", color=C_GPU, alpha=0.85)
        ax.set_xticks(x)
        ax.set_xticklabels([f"gid {g}" for g in gids], rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("Peak RSS per rank (MB)")
        ax.set_title("Peak Memory per Block: CPU vs GPU\n(lowres · 16 blocks · 10,000 seeds)")
        ax.legend()
        ax.grid(axis="y", alpha=0.3)

        avg_cpu = np.mean(cpu_rss)
        avg_gpu = np.mean(gpu_rss)
        ax.axhline(avg_cpu, color=C_CPU, linestyle="--", linewidth=1,
                   label=f"CPU avg {avg_cpu:.0f} MB")
        ax.axhline(avg_gpu, color=C_GPU, linestyle="--", linewidth=1,
                   label=f"GPU avg {avg_gpu:.0f} MB")
        ax.legend(fontsize=8)
        save(fig, "C6_memory_cpu_vs_gpu.png")
        print(f"  [C6] Peak RSS — CPU avg={avg_cpu:.0f} MB  GPU avg={avg_gpu:.0f} MB  "
              f"ratio={avg_gpu/avg_cpu:.2f}×")

    # ── C7: highres vs lowres CPU comparison ──────────────────────────────────
    cpu_runs = valid[(valid.device == "cpu") & (valid.n_blocks == 16) & (valid.n_seeds == 10000)]
    if len(cpu_runs) >= 2:
        fig, ax = plt.subplots(figsize=(6, 5))
        labels = [f"{row['mesh'].capitalize()}\nCPU 16b 10k seeds"
                  for _, row in cpu_runs.iterrows()]
        walls  = [row["t_iex_max"] / 3600 for _, row in cpu_runs.iterrows()]
        ax.bar(labels, walls, color=[C_CPU, "#5C8FE8"], width=0.4)
        for i, (label, wall) in enumerate(zip(labels, walls)):
            ax.text(i, wall + 0.05, f"{wall:.2f} h", ha="center", fontsize=11,
                    fontweight="bold")
        ax.set_ylabel("Wall-clock time (hours)")
        ax.set_title("CPU: Lowres vs Highres Mesh\n(16 blocks · 10,000 seeds)")
        ax.grid(axis="y", alpha=0.3)
        save(fig, "C7_cpu_lowres_vs_highres.png")


# ── run ────────────────────────────────────────────────────────────────────────
print("\n── Section A: per-block plots ──")
for (nb, ns, rid, dev, mesh), grp in blocks.groupby(
        ["n_blocks", "n_seeds", "run_id", "device", "mesh"]):
    print(f"  run: n_blocks={nb}  n_seeds={ns}  run_id={rid}  device={dev}  mesh={mesh}")
    plot_per_run(grp, nb, ns, rid, dev, mesh)

print("\n── Section B: scaling plots (CPU) ──")
plot_scaling(valid, device_filter="cpu")

print("\n── Section B: scaling plots (GPU) ──")
plot_scaling(valid, device_filter="gpu")

print("\n── Section B: scaling plots (all devices combined) ──")
plot_scaling(valid)

print("\n── Section C: GPU vs CPU comparison ──")
plot_gpu_vs_cpu(valid, blocks)

print(f"\nDone. All plots saved to {PLOTS_DIR}/")
