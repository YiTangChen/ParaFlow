#!/usr/bin/env python3
"""
Parse ParaFlow timing + memory logs into CSV files for plotting.

Each log file is one run (one n_blocks / n_seeds combination).
The filename encodes the configuration using the convention:
    logs/b<nblocks>_s<nseeds>_run<N>.stderr
    e.g.  logs/b016_s010000_run1.stderr

Outputs (appended to, so you can run on many files):
    results/run_summary.csv   — one row per run   (scaling plots)
    results/rank_detail.csv   — one row per rank  (iexchange timing decomposition)
    results/block_detail.csv  — one row per block  (load imbalance plots)

Usage:
    # parse one file
    python3 parse_timing.py logs/b016_s002000_run1.stderr

    # parse all logs at once
    python3 parse_timing.py logs/*.stderr

    # reset output files first (fresh experiment)
    python3 parse_timing.py --reset logs/*.stderr

    # write to a separate result directory
    python3 parse_timing.py --results-dir=results/small_async --reset logs/*.stderr
    
    # parse files with non-standard names (override n_blocks/n_seeds/run_id)
    python3 parse_timing.py --blocks=256 --seeds=10000 --run=1 slurm-51151129.out
"""

import re
import sys
import csv
import os
from pathlib import Path
from collections import defaultdict


# ─── default output file paths ────────────────────────────────────────────────
RESULTS_DIR    = Path("results/slurm1129")
SUMMARY_CSV    = RESULTS_DIR / "run_summary.csv"
RANK_CSV       = RESULTS_DIR / "rank_detail.csv"
BLOCK_CSV      = RESULTS_DIR / "block_detail.csv"

SUMMARY_FIELDS = [
    "n_blocks", "n_seeds", "run_id",
    "t_seed_read",
    "t_seed_bcast_min",        "t_seed_bcast_max",        "t_seed_bcast_avg",
    "t_blockload_min",         "t_blockload_max",         "t_blockload_avg",
    "t_netcdf_min",            "t_netcdf_max",            "t_netcdf_avg",
    "t_filter_min",            "t_filter_max",            "t_filter_avg",
    "t_compute_min",           "t_compute_max",           "t_compute_avg",
    "t_enqueue_min",           "t_enqueue_max",           "t_enqueue_avg",
    "t_dequeue_local_min",     "t_dequeue_local_max",     "t_dequeue_local_avg",
    "t_msg_overhead_min",      "t_msg_overhead_max",      "t_msg_overhead_avg",
    "t_fill_incoming_min",     "t_fill_incoming_max",     "t_fill_incoming_avg",
    "t_wait_sync_min",         "t_wait_sync_max",         "t_wait_sync_avg",
    "t_residual_min",          "t_residual_max",          "t_residual_avg",
    "t_write_min",             "t_write_max",             "t_write_avg",
    "t_iex_min",               "t_iex_max",               "t_iex_avg",
    "t_iex_barrier_min",       "t_iex_barrier_max",       "t_iex_barrier_avg",
    "t_iex_with_barrier_min",  "t_iex_with_barrier_max",  "t_iex_with_barrier_avg",
    "iex_imbalance",
    "mem_grid_bytes_min",      "mem_grid_bytes_max",      "mem_grid_bytes_avg",
    "mem_solution_bytes_min",  "mem_solution_bytes_max",  "mem_solution_bytes_avg",
    "mem_delta_kb_min",        "mem_delta_kb_max",        "mem_delta_kb_avg",
    "mem_peak_vmhwm_kb_min",   "mem_peak_vmhwm_kb_max",  "mem_peak_vmhwm_kb_avg",
    "total_seeds_assigned",    "total_steps",
    "total_particles_received", "total_particles_sent",
    "total_iex_rounds",
]

RANK_FIELDS = [
    "n_blocks", "n_seeds", "run_id", "rank",
    "t_iex", "t_iex_barrier", "t_iex_with_barrier",
    "t_block_load", "t_netcdf_read", "t_seed_filter",
    "t_trace_compute", "t_trace_enqueue", "t_dequeue_local", "t_fill_incoming",
    "t_msg_overhead", "t_wait_sync", "t_residual",
    "t_output_write",
    "n_local_blocks", "n_iex_rounds", "n_seeds_initial", "n_steps_total",
    "n_particles_received", "n_particles_sent",
]

BLOCK_FIELDS = [
    "n_blocks", "n_seeds", "run_id", "gid", "rank",
    "t_block_load", "t_netcdf_read", "t_seed_filter",
    "t_trace_compute", "t_trace_enqueue", "t_dequeue_local", "t_fill_incoming",
    "t_msg_overhead", "t_wait_sync",
    "t_output_write",
    "t_trace_idle",          # residual: rank iexchange - explicit timed phases
    "n_iex_rounds", "n_seeds_initial", "n_steps_total",
    "n_particles_received", "n_particles_sent",
    "mem_grid_bytes", "mem_solution_bytes", "mem_total_bytes",
    "mem_delta_kb", "mem_peak_vmhwm_kb",
]


# ─── line parsers ─────────────────────────────────────────────────────────────
def parse_kv(line: str) -> dict:
    """Parse 'key=value key=value ...' tokens into a dict (strings)."""
    return dict(m.groups() for m in re.finditer(r'(\w+)=([\w.+\-]+)', line))


def parse_log(path: str):
    """
    Read one log file and return:
        global_timing  — dict of scalar values from non-block lines
        ranks          — dict keyed by rank, each a dict of rank-level values
        blocks         — dict keyed by gid, each a dict of accumulated values
    """
    global_timing = {}   # seed_read, seed_bcast list, iexchange_global
    ranks = defaultdict(dict)
    blocks = defaultdict(dict)  # gid -> field -> value

    seed_bcasts = []     # one entry per rank
    iex_totals  = {}     # rank -> iexchange time, excluding post-iexchange barrier
    iex_barriers = {}    # rank -> post-iexchange barrier time

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not (line.startswith("TIMING") or line.startswith("MEM_")):
                continue
            kv = parse_kv(line)
            phase = kv.get("phase", "")

            # ── global / rank-level lines ──────────────────────────────────
            if phase == "seed_read":
                global_timing["t_seed_read"]  = float(kv.get("t", 0))
                global_timing["n_seeds_total"] = int(kv.get("n_seeds", 0))

            elif phase == "seed_bcast":
                seed_bcasts.append(float(kv.get("t", 0)))

            elif phase == "iexchange_global":
                global_timing["iex_min"]       = float(kv.get("iex_min", 0))
                global_timing["iex_max"]       = float(kv.get("iex_max", 0))
                global_timing["iex_avg"]       = float(kv.get("iex_avg", 0))
                global_timing["iex_barrier_min"] = float(kv.get("iex_barrier_min", 0))
                global_timing["iex_barrier_max"] = float(kv.get("iex_barrier_max", 0))
                global_timing["iex_barrier_avg"] = float(kv.get("iex_barrier_avg", 0))
                global_timing["iex_with_barrier_min"] = float(kv.get("iex_with_barrier_min", 0))
                global_timing["iex_with_barrier_max"] = float(kv.get("iex_with_barrier_max", 0))
                global_timing["iex_with_barrier_avg"] = float(kv.get("iex_with_barrier_avg", 0))
                global_timing["iex_imbalance"] = float(kv.get("iex_imbalance", 0))

            elif phase == "iexchange_total":
                rank = int(kv.get("rank", 0))
                iex_totals[rank] = float(kv.get("t", 0))
                iex_barriers[rank] = float(kv.get("t_barrier", 0))

            elif phase == "iexchange_barrier":
                rank = int(kv.get("rank", 0))
                iex_barriers[rank] = float(kv.get("t", 0))

            # ── per-block lines ────────────────────────────────────────────
            # `trace_comm` is an accepted alias for `dequeue_local` so that
            # log files written before the rename still parse correctly.
            elif phase in ("block_load", "trace_compute", "trace_comm",
                           "dequeue_local", "trace_enqueue", "fill_incoming",
                           "output_write"):
                gid = int(kv.get("gid", -1))
                b   = blocks[gid]
                b["gid"]  = gid
                b["rank"] = int(kv.get("rank", 0))
                if phase == "block_load":
                    b["n_seeds_initial"] = int(kv.get("nseeds", 0))
                    b["t_netcdf_read"]   = float(kv.get("t_netcdf", 0))
                    b["t_seed_filter"]   = float(kv.get("t_filter", 0))
                    # Old logs carried an explicit `t=` (= netcdf + filter).
                    # New logs omit it; we derive the block-load total here so
                    # the t_block_load CSV column stays populated for plotting.
                    if "t" in kv:
                        b["t_block_load"] = float(kv["t"])
                    else:
                        b["t_block_load"] = b["t_netcdf_read"] + b["t_seed_filter"]
                elif phase == "trace_compute":
                    b["t_trace_compute"]      = float(kv.get("t", 0))
                    b["n_steps_total"]        = int(kv.get("nsteps", 0))
                    b["n_particles_received"] = int(kv.get("nrecv", 0))
                    b["n_particles_sent"]     = int(kv.get("nsent", 0))
                    b["n_iex_rounds"]         = int(kv.get("rounds", 0))
                elif phase in ("dequeue_local", "trace_comm"):
                    b["t_dequeue_local"] = float(kv.get("t", 0))
                elif phase == "trace_enqueue":
                    b["t_trace_enqueue"] = float(kv.get("t", 0))
                elif phase == "fill_incoming":
                    b["t_fill_incoming"] = float(kv.get("t", 0))
                elif phase == "output_write":
                    b["t_output_write"]  = float(kv.get("t", 0))

            # ── memory lines ───────────────────────────────────────────────
            elif line.startswith("MEM_ANALYTICAL"):
                gid = int(kv.get("gid", -1))
                blocks[gid]["mem_grid_bytes"]     = int(kv.get("grid_bytes", 0))
                blocks[gid]["mem_solution_bytes"] = int(kv.get("solution_bytes", 0))

            elif line.startswith("MEM_DELTA"):
                gid = int(kv.get("gid", -1))
                blocks[gid]["mem_delta_kb"] = int(kv.get("delta_kb", 0))

            elif line.startswith("MEM_PEAK"):
                gid = int(kv.get("gid", -1))
                # Support both key=value format (vmhwm_kb=N) and
                # /proc/status format (VmHWM:\t<N> kB)
                vmhwm = kv.get("vmhwm_kb")
                if vmhwm is None:
                    m = re.search(r'VmHWM:\s*(\d+)', line)
                    vmhwm = m.group(1) if m else 0
                blocks[gid]["mem_peak_vmhwm_kb"] = int(vmhwm)

    # Store per-rank iexchange total on each block (same rank = same iex time).
    # t_iex excludes the post-iexchange barrier.  t_iex_barrier is tracked
    # separately and should not be folded into the asynchronous iexchange
    # breakdown.
    for gid, b in blocks.items():
        rank = b.get("rank", 0)
        b["t_iex_rank"] = iex_totals.get(rank, 0.0)
        b["t_iex_barrier_rank"] = iex_barriers.get(rank, 0.0)

    global_timing["seed_bcasts"] = seed_bcasts
    for rank in set(iex_totals) | set(iex_barriers):
        ranks[rank]["rank"] = rank
        ranks[rank]["t_iex"] = iex_totals.get(rank, 0.0)
        ranks[rank]["t_iex_barrier"] = iex_barriers.get(rank, 0.0)
        ranks[rank]["t_iex_with_barrier"] = (
            ranks[rank]["t_iex"] + ranks[rank]["t_iex_barrier"]
        )
    return global_timing, ranks, blocks


def _stats(values):
    """Return (min, max, avg) or (0,0,0) for empty list."""
    if not values:
        return 0.0, 0.0, 0.0
    return min(values), max(values), sum(values) / len(values)


def extract_n_blocks_n_seeds(path: str):
    """
    Extract n_blocks and n_seeds from filename convention:
        b<nblocks>_s<nseeds>_run<N>.stderr
    Returns (n_blocks, n_seeds, run_id) as ints, or (None, None, None).
    """
    name = Path(path).stem   # e.g. b016_s010000_run1
    m = re.match(r'b(\d+)_s(\d+)_run(\d+)', name)
    if m:
        return int(m.group(1)), int(m.group(2)), int(m.group(3))
    return None, None, None


def build_rank_rows(n_blocks, n_seeds, run_id, ranks, blocks):
    """Aggregate per-block timings to rank-level rows."""
    by_rank = defaultdict(list)
    for b in blocks.values():
        by_rank[b.get("rank", 0)].append(b)

    rows = []
    for rank in sorted(set(by_rank) | set(ranks)):
        bdata = by_rank.get(rank, [])

        def total(field):
            return sum(b.get(field, 0) for b in bdata)

        t_iex = ranks.get(rank, {}).get("t_iex", 0.0)
        t_barrier = ranks.get(rank, {}).get("t_iex_barrier", 0.0)
        t_comp = total("t_trace_compute")
        t_enq = total("t_trace_enqueue")
        t_deq = total("t_dequeue_local")
        t_fill = total("t_fill_incoming")
        t_msg = t_enq + t_deq
        # DIY iexchange is asynchronous; communication may overlap with compute.
        # t_wait_sync is therefore only the exposed wait/progress inside
        # fill_incoming(), not total network communication and not barrier time.
        t_wait = t_fill
        t_residual = max(0.0, t_iex - t_comp - t_enq - t_deq - t_fill)

        rows.append({
            "n_blocks": n_blocks, "n_seeds": n_seeds, "run_id": run_id,
            "rank": rank,
            "t_iex": t_iex,
            "t_iex_barrier": t_barrier,
            "t_iex_with_barrier": t_iex + t_barrier,
            "t_block_load": total("t_block_load"),
            "t_netcdf_read": total("t_netcdf_read"),
            "t_seed_filter": total("t_seed_filter"),
            "t_trace_compute": t_comp,
            "t_trace_enqueue": t_enq,
            "t_dequeue_local": t_deq,
            "t_fill_incoming": t_fill,
            "t_msg_overhead": t_msg,
            "t_wait_sync": t_wait,
            "t_residual": t_residual,
            "t_output_write": total("t_output_write"),
            "n_local_blocks": len(bdata),
            "n_iex_rounds": sum(b.get("n_iex_rounds", 0) for b in bdata),
            "n_seeds_initial": sum(b.get("n_seeds_initial", 0) for b in bdata),
            "n_steps_total": sum(b.get("n_steps_total", 0) for b in bdata),
            "n_particles_received": sum(b.get("n_particles_received", 0) for b in bdata),
            "n_particles_sent": sum(b.get("n_particles_sent", 0) for b in bdata),
        })
    return rows


def build_summary_row(n_blocks, n_seeds, run_id, global_timing, blocks, rank_rows):
    bdata = list(blocks.values())
    rdata = rank_rows

    def col(field):
        return [b.get(field, 0) for b in bdata]

    def rcol(field):
        return [r.get(field, 0) for r in rdata]

    bcasts = global_timing.get("seed_bcasts", [])
    bc_min, bc_max, bc_avg = _stats(bcasts)
    iex_min, iex_max, iex_avg = _stats(rcol("t_iex"))
    barrier_min, barrier_max, barrier_avg = _stats(rcol("t_iex_barrier"))
    with_barrier_min, with_barrier_max, with_barrier_avg = _stats(rcol("t_iex_with_barrier"))

    return {
        "n_blocks": n_blocks, "n_seeds": n_seeds, "run_id": run_id,
        "t_seed_read": global_timing.get("t_seed_read", 0),
        "t_seed_bcast_min": bc_min, "t_seed_bcast_max": bc_max, "t_seed_bcast_avg": bc_avg,
        **dict(zip(["t_blockload_min",     "t_blockload_max",     "t_blockload_avg"],
                   _stats(rcol("t_block_load")))),
        **dict(zip(["t_netcdf_min",        "t_netcdf_max",        "t_netcdf_avg"],
                   _stats(rcol("t_netcdf_read")))),
        **dict(zip(["t_filter_min",        "t_filter_max",        "t_filter_avg"],
                   _stats(rcol("t_seed_filter")))),
        **dict(zip(["t_compute_min",       "t_compute_max",       "t_compute_avg"],
                   _stats(rcol("t_trace_compute")))),
        **dict(zip(["t_enqueue_min",       "t_enqueue_max",       "t_enqueue_avg"],
                   _stats(rcol("t_trace_enqueue")))),
        **dict(zip(["t_dequeue_local_min", "t_dequeue_local_max", "t_dequeue_local_avg"],
                   _stats(rcol("t_dequeue_local")))),
        **dict(zip(["t_msg_overhead_min",  "t_msg_overhead_max",  "t_msg_overhead_avg"],
                   _stats(rcol("t_msg_overhead")))),
        **dict(zip(["t_fill_incoming_min", "t_fill_incoming_max", "t_fill_incoming_avg"],
                   _stats(rcol("t_fill_incoming")))),
        **dict(zip(["t_wait_sync_min",     "t_wait_sync_max",     "t_wait_sync_avg"],
                   _stats(rcol("t_wait_sync")))),
        **dict(zip(["t_residual_min",      "t_residual_max",      "t_residual_avg"],
                   _stats(rcol("t_residual")))),
        **dict(zip(["t_write_min",         "t_write_max",         "t_write_avg"],
                   _stats(rcol("t_output_write")))),
        "t_iex_min":     global_timing.get("iex_min", iex_min),
        "t_iex_max":     global_timing.get("iex_max", iex_max),
        "t_iex_avg":     global_timing.get("iex_avg", iex_avg),
        "t_iex_barrier_min":      global_timing.get("iex_barrier_min", barrier_min),
        "t_iex_barrier_max":      global_timing.get("iex_barrier_max", barrier_max),
        "t_iex_barrier_avg":      global_timing.get("iex_barrier_avg", barrier_avg),
        "t_iex_with_barrier_min": global_timing.get("iex_with_barrier_min", with_barrier_min),
        "t_iex_with_barrier_max": global_timing.get("iex_with_barrier_max", with_barrier_max),
        "t_iex_with_barrier_avg": global_timing.get("iex_with_barrier_avg", with_barrier_avg),
        "iex_imbalance": global_timing.get(
            "iex_imbalance",
            (iex_max / iex_min) if iex_min > 0.0 else 0.0,
        ),
        **dict(zip(
            ["mem_grid_bytes_min",     "mem_grid_bytes_max",     "mem_grid_bytes_avg"],
            _stats(col("mem_grid_bytes")))),
        **dict(zip(
            ["mem_solution_bytes_min", "mem_solution_bytes_max", "mem_solution_bytes_avg"],
            _stats(col("mem_solution_bytes")))),
        **dict(zip(
            ["mem_delta_kb_min",       "mem_delta_kb_max",       "mem_delta_kb_avg"],
            _stats(col("mem_delta_kb")))),
        **dict(zip(
            ["mem_peak_vmhwm_kb_min",  "mem_peak_vmhwm_kb_max",  "mem_peak_vmhwm_kb_avg"],
            _stats(col("mem_peak_vmhwm_kb")))),
        "total_seeds_assigned":      sum(col("n_seeds_initial")),
        "total_steps":               sum(col("n_steps_total")),
        "total_particles_received":  sum(col("n_particles_received")),
        "total_particles_sent":      sum(col("n_particles_sent")),
        "total_iex_rounds":          sum(col("n_iex_rounds")),
    }


def build_block_rows(n_blocks, n_seeds, run_id, blocks):
    rows = []
    rank_counts = defaultdict(int)
    for b in blocks.values():
        rank_counts[b.get("rank", 0)] += 1
    for gid, b in sorted(blocks.items()):
        t_comp  = b.get("t_trace_compute",  0.0)
        t_enq   = b.get("t_trace_enqueue",  0.0)
        t_deq   = b.get("t_dequeue_local",  0.0)
        t_fill  = b.get("t_fill_incoming",  0.0)
        t_iex   = b.get("t_iex_rank",       0.0)
        t_msg   = t_enq + t_deq
        # Block-level wait excludes the rank-level post-iexchange barrier.
        # Accurate barrier/residual accounting is reported in rank_detail.csv.
        t_wait  = t_fill
        t_idle  = (
            max(0.0, t_iex - t_comp - t_enq - t_deq - t_fill)
            if rank_counts[b.get("rank", 0)] == 1 else 0.0
        )
        rows.append({
            "n_blocks": n_blocks, "n_seeds": n_seeds, "run_id": run_id,
            "gid":  gid,
            "rank": b.get("rank", 0),
            "t_block_load":     b.get("t_block_load",    0.0),
            "t_netcdf_read":    b.get("t_netcdf_read",   0.0),
            "t_seed_filter":    b.get("t_seed_filter",   0.0),
            "t_trace_compute":  t_comp,
            "t_trace_enqueue":  t_enq,
            "t_dequeue_local":  t_deq,
            "t_fill_incoming":  t_fill,
            "t_msg_overhead":   t_msg,
            "t_wait_sync":      t_wait,
            "t_output_write":   b.get("t_output_write",  0.0),
            "t_trace_idle":     t_idle,
            "n_iex_rounds":          b.get("n_iex_rounds",          0),
            "n_seeds_initial":       b.get("n_seeds_initial",       0),
            "n_steps_total":         b.get("n_steps_total",         0),
            "n_particles_received":  b.get("n_particles_received",  0),
            "n_particles_sent":      b.get("n_particles_sent",      0),
            "mem_grid_bytes":     b.get("mem_grid_bytes",     0),
            "mem_solution_bytes": b.get("mem_solution_bytes", 0),
            "mem_total_bytes":    b.get("mem_grid_bytes", 0) + b.get("mem_solution_bytes", 0),
            "mem_delta_kb":       b.get("mem_delta_kb",     0),
            "mem_peak_vmhwm_kb":  b.get("mem_peak_vmhwm_kb", 0),
        })
    return rows


def write_csv(path, fields, rows, append=True):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not path.exists() or not append
    mode = "a" if append else "w"
    with open(path, mode, newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        if write_header:
            writer.writeheader()
        writer.writerows(rows)


# ─── main ─────────────────────────────────────────────────────────────────────
def main():
    global RESULTS_DIR, SUMMARY_CSV, RANK_CSV, BLOCK_CSV

    args = sys.argv[1:]
    reset = False
    override_blocks = None
    override_seeds  = None
    override_run    = None
    results_dir = RESULTS_DIR

    filtered = []
    for a in args:
        if a == "--reset":
            reset = True
        elif a.startswith("--results-dir="):
            results_dir = Path(a.split("=", 1)[1])
        elif a.startswith("--blocks="):
            override_blocks = int(a.split("=", 1)[1])
        elif a.startswith("--seeds="):
            override_seeds = int(a.split("=", 1)[1])
        elif a.startswith("--run="):
            override_run = int(a.split("=", 1)[1])
        else:
            filtered.append(a)
    args = filtered

    if not args:
        print(__doc__)
        sys.exit(1)

    RESULTS_DIR = results_dir
    SUMMARY_CSV = RESULTS_DIR / "run_summary.csv"
    RANK_CSV    = RESULTS_DIR / "rank_detail.csv"
    BLOCK_CSV   = RESULTS_DIR / "block_detail.csv"

    if reset:
        for p in (SUMMARY_CSV, RANK_CSV, BLOCK_CSV):
            if p.exists():
                p.unlink()
                print(f"Removed {p}")

    for path in args:
        n_blocks, n_seeds, run_id = extract_n_blocks_n_seeds(path)
        if n_blocks is None:
            if override_blocks is not None and override_seeds is not None:
                n_blocks = override_blocks
                n_seeds  = override_seeds
                run_id   = override_run if override_run is not None else 1
            else:
                print(f"WARN: cannot parse n_blocks/n_seeds from '{path}' — expected b<N>_s<N>_run<N>.stderr")
                print(f"      Use --blocks=N --seeds=N [--run=N] to override.")
                continue

        global_timing, ranks, blocks = parse_log(path)

        rank_rows    = build_rank_rows(n_blocks, n_seeds, run_id, ranks, blocks)
        summary_row  = build_summary_row(n_blocks, n_seeds, run_id, global_timing, blocks, rank_rows)
        block_rows   = build_block_rows(n_blocks, n_seeds, run_id, blocks)

        write_csv(SUMMARY_CSV, SUMMARY_FIELDS, [summary_row], append=not reset)
        write_csv(RANK_CSV,    RANK_FIELDS,    rank_rows,     append=not reset)
        write_csv(BLOCK_CSV,   BLOCK_FIELDS,   block_rows,    append=not reset)
        reset = False   # only reset on first file

        print(f"[{Path(path).name}] n_blocks={n_blocks} n_seeds={n_seeds} "
              f"blocks_parsed={len(blocks)} "
              f"iex_imbalance={summary_row['iex_imbalance']:.3f}")

    print(f"\nOutputs:")
    print(f"  {SUMMARY_CSV}  ({SUMMARY_CSV.stat().st_size if SUMMARY_CSV.exists() else 0} bytes)")
    print(f"  {RANK_CSV}     ({RANK_CSV.stat().st_size    if RANK_CSV.exists()    else 0} bytes)")
    print(f"  {BLOCK_CSV}    ({BLOCK_CSV.stat().st_size   if BLOCK_CSV.exists()   else 0} bytes)")


if __name__ == "__main__":
    main()
