#!/usr/bin/env python3
"""
Parse ParaFlow timing + memory logs into CSV files for plotting.

Each log file is one run (one n_blocks / n_seeds combination).
The filename encodes the configuration using the convention:
    logs/b<nblocks>_s<nseeds>_run<N>.stderr
    e.g.  logs/b016_s010000_run1.stderr

Outputs (appended to, so you can run on many files):
    results/run_summary.csv   — one row per run   (scaling plots)
    results/block_detail.csv  — one row per block  (load imbalance plots)

Usage:
    # parse one file
    python3 parse_timing.py logs/b016_s002000_run1.stderr

    # parse all logs at once
    python3 parse_timing.py logs/*.stderr

    # reset output files first (fresh experiment)
    python3 parse_timing.py --reset logs/*.stderr

    # parse files with non-standard names (override n_blocks/n_seeds/run_id)
    python3 parse_timing.py --blocks=256 --seeds=10000 --run=1 slurm-51151129.out

    # specify device and mesh type
    python3 parse_timing.py --blocks=16 --seeds=10000 --device=gpu --mesh=lowres run.err
"""

import re
import sys
import csv
import os
from pathlib import Path
from collections import defaultdict


# ─── output file paths ────────────────────────────────────────────────────────
RESULTS_DIR    = Path("results/gpu_analysis")
SUMMARY_CSV    = RESULTS_DIR / "run_summary.csv"
BLOCK_CSV      = RESULTS_DIR / "block_detail.csv"

SUMMARY_FIELDS = [
    "n_blocks", "n_seeds", "run_id", "device", "mesh",
    "t_seed_read",
    "t_seed_bcast_min", "t_seed_bcast_max", "t_seed_bcast_avg",
    "t_blockload_min",  "t_blockload_max",  "t_blockload_avg",
    "t_compute_min",    "t_compute_max",    "t_compute_avg",
    "t_gpukernel_min",  "t_gpukernel_max",  "t_gpukernel_avg",
    "t_comm_min",       "t_comm_max",       "t_comm_avg",
    "t_write_min",      "t_write_max",      "t_write_avg",
    "t_iex_min",        "t_iex_max",        "t_iex_avg",
    "iex_imbalance",
    "mem_grid_bytes_min",      "mem_grid_bytes_max",      "mem_grid_bytes_avg",
    "mem_solution_bytes_min",  "mem_solution_bytes_max",  "mem_solution_bytes_avg",
    "mem_delta_kb_min",        "mem_delta_kb_max",        "mem_delta_kb_avg",
    "mem_peak_vmhwm_kb_min",   "mem_peak_vmhwm_kb_max",  "mem_peak_vmhwm_kb_avg",
    "total_seeds_assigned",    "total_steps",             "total_particles_received",
]

BLOCK_FIELDS = [
    "n_blocks", "n_seeds", "run_id", "device", "mesh", "gid", "rank",
    "t_block_load", "t_trace_compute", "t_gpu_kernel", "t_trace_comm", "t_output_write",
    "t_trace_idle",
    "n_seeds_initial", "n_steps_total", "n_particles_received",
    "n_local_cells", "n_global_cells",
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
        blocks         — dict keyed by gid, each a dict of accumulated values
        has_gpu        — True if gpu_kernel lines were found
    """
    global_timing = {}
    blocks = defaultdict(dict)

    seed_bcasts = []
    iex_totals  = []
    has_gpu     = False

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
                global_timing["iex_imbalance"] = float(kv.get("iex_imbalance", 0))

            elif phase == "iexchange_total":
                iex_totals.append((int(kv.get("rank", 0)), float(kv.get("t", 0))))

            # ── per-block lines ────────────────────────────────────────────
            elif phase in ("block_load", "trace_compute", "gpu_kernel",
                           "trace_comm", "output_write"):
                gid = int(kv.get("gid", -1))
                b   = blocks[gid]
                b["gid"]  = gid
                b["rank"] = int(kv.get("rank", 0))
                if phase == "block_load":
                    b["t_block_load"]    = float(kv.get("t", 0))
                    b["n_seeds_initial"] = int(kv.get("nseeds", 0))
                elif phase == "trace_compute":
                    b["t_trace_compute"]      = float(kv.get("t", 0))
                    b["n_steps_total"]        = int(kv.get("nsteps", 0))
                    b["n_particles_received"] = int(kv.get("nrecv", 0))
                elif phase == "gpu_kernel":
                    b["t_gpu_kernel"] = float(kv.get("t", 0))
                    has_gpu = True
                elif phase == "trace_comm":
                    b["t_trace_comm"] = float(kv.get("t", 0))
                elif phase == "output_write":
                    b["t_output_write"] = float(kv.get("t", 0))

            # ── memory lines ───────────────────────────────────────────────
            elif line.startswith("MEM_ANALYTICAL"):
                gid = int(kv.get("gid", -1))
                blocks[gid]["mem_grid_bytes"]     = int(kv.get("grid_bytes", 0))
                blocks[gid]["mem_solution_bytes"] = int(kv.get("solution_bytes", 0))

            elif line.startswith("MEM_DELTA"):
                gid = int(kv.get("gid", -1))
                blocks[gid]["mem_delta_kb"] = int(kv.get("delta_kb", 0))

            elif line.startswith("MEM_CELLCOUNT"):
                gid = int(kv.get("gid", -1))
                blocks[gid]["n_local_cells"]  = int(kv.get("n_local_cells", 0))
                blocks[gid]["n_global_cells"] = int(kv.get("n_global_cells", 0))

            elif line.startswith("MEM_PEAK"):
                gid = int(kv.get("gid", -1))
                vmhwm = kv.get("vmhwm_kb")
                if vmhwm is None:
                    m = re.search(r'VmHWM:\s*(\d+)', line)
                    vmhwm = m.group(1) if m else 0
                blocks[gid]["mem_peak_vmhwm_kb"] = int(vmhwm)

    rank_to_iex = {rank: t for rank, t in iex_totals}
    for gid, b in blocks.items():
        b["t_iex_rank"] = rank_to_iex.get(b.get("rank", 0), 0.0)

    global_timing["seed_bcasts"] = seed_bcasts
    return global_timing, blocks, has_gpu


def _stats(values):
    """Return (min, max, avg) or (0,0,0) for empty list."""
    if not values:
        return 0.0, 0.0, 0.0
    return min(values), max(values), sum(values) / len(values)


def extract_n_blocks_n_seeds(path: str):
    name = Path(path).stem
    m = re.match(r'b(\d+)_s(\d+)_run(\d+)', name)
    if m:
        return int(m.group(1)), int(m.group(2)), int(m.group(3))
    return None, None, None


def build_summary_row(n_blocks, n_seeds, run_id, device, mesh, global_timing, blocks):
    bdata = list(blocks.values())

    def col(field):
        return [b.get(field, 0) for b in bdata]

    bcasts = global_timing.get("seed_bcasts", [])
    bc_min, bc_max, bc_avg = _stats(bcasts)
    bl_min, bl_max, bl_avg = _stats(col("t_block_load"))
    co_min, co_max, co_avg = _stats(col("t_trace_compute"))
    gk_min, gk_max, gk_avg = _stats([v for v in col("t_gpu_kernel") if v > 0])
    cm_min, cm_max, cm_avg = _stats(col("t_trace_comm"))
    wr_min, wr_max, wr_avg = _stats(col("t_output_write"))

    return {
        "n_blocks": n_blocks, "n_seeds": n_seeds, "run_id": run_id,
        "device": device, "mesh": mesh,
        "t_seed_read":        global_timing.get("t_seed_read", 0),
        "t_seed_bcast_min": bc_min, "t_seed_bcast_max": bc_max, "t_seed_bcast_avg": bc_avg,
        "t_blockload_min":  bl_min, "t_blockload_max":  bl_max, "t_blockload_avg":  bl_avg,
        "t_compute_min":    co_min, "t_compute_max":    co_max, "t_compute_avg":    co_avg,
        "t_gpukernel_min":  gk_min, "t_gpukernel_max":  gk_max, "t_gpukernel_avg":  gk_avg,
        "t_comm_min":       cm_min, "t_comm_max":       cm_max, "t_comm_avg":       cm_avg,
        "t_write_min":      wr_min, "t_write_max":      wr_max, "t_write_avg":      wr_avg,
        "t_iex_min":     global_timing.get("iex_min", 0),
        "t_iex_max":     global_timing.get("iex_max", 0),
        "t_iex_avg":     global_timing.get("iex_avg", 0),
        "iex_imbalance": global_timing.get("iex_imbalance", 0),
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
        "total_seeds_assigned":     sum(col("n_seeds_initial")),
        "total_steps":              sum(col("n_steps_total")),
        "total_particles_received": sum(col("n_particles_received")),
    }


def build_block_rows(n_blocks, n_seeds, run_id, device, mesh, blocks):
    rows = []
    for gid, b in sorted(blocks.items()):
        t_comp = b.get("t_trace_compute", 0.0)
        t_comm = b.get("t_trace_comm",    0.0)
        t_iex  = b.get("t_iex_rank",      0.0)
        t_idle = max(0.0, t_iex - t_comp - t_comm)
        rows.append({
            "n_blocks": n_blocks, "n_seeds": n_seeds, "run_id": run_id,
            "device": device, "mesh": mesh,
            "gid":  gid,
            "rank": b.get("rank", 0),
            "t_block_load":    b.get("t_block_load",   0.0),
            "t_trace_compute": t_comp,
            "t_gpu_kernel":    b.get("t_gpu_kernel",   0.0),
            "t_trace_comm":    t_comm,
            "t_output_write":  b.get("t_output_write", 0.0),
            "t_trace_idle":    t_idle,
            "n_seeds_initial":      b.get("n_seeds_initial",      0),
            "n_steps_total":        b.get("n_steps_total",        0),
            "n_particles_received": b.get("n_particles_received", 0),
            "n_local_cells":      b.get("n_local_cells",      0),
            "n_global_cells":     b.get("n_global_cells",     0),
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
    args = sys.argv[1:]
    reset = False
    override_blocks = None
    override_seeds  = None
    override_run    = None
    override_device = None
    override_mesh   = None

    filtered = []
    for a in args:
        if a == "--reset":
            reset = True
        elif a.startswith("--blocks="):
            override_blocks = int(a.split("=", 1)[1])
        elif a.startswith("--seeds="):
            override_seeds = int(a.split("=", 1)[1])
        elif a.startswith("--run="):
            override_run = int(a.split("=", 1)[1])
        elif a.startswith("--device="):
            override_device = a.split("=", 1)[1]
        elif a.startswith("--mesh="):
            override_mesh = a.split("=", 1)[1]
        else:
            filtered.append(a)
    args = filtered

    if not args:
        print(__doc__)
        sys.exit(1)

    if reset:
        for p in (SUMMARY_CSV, BLOCK_CSV):
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
                print(f"WARN: cannot parse n_blocks/n_seeds from '{path}'")
                print(f"      Use --blocks=N --seeds=N [--run=N] to override.")
                continue

        global_timing, blocks, has_gpu = parse_log(path)

        device = override_device if override_device else ("gpu" if has_gpu else "cpu")
        mesh   = override_mesh   if override_mesh   else "lowres"

        summary_row = build_summary_row(n_blocks, n_seeds, run_id, device, mesh,
                                        global_timing, blocks)
        block_rows  = build_block_rows(n_blocks, n_seeds, run_id, device, mesh, blocks)

        write_csv(SUMMARY_CSV, SUMMARY_FIELDS, [summary_row], append=not reset)
        write_csv(BLOCK_CSV,   BLOCK_FIELDS,   block_rows,    append=not reset)
        reset = False

        print(f"[{Path(path).name}] device={device} mesh={mesh} "
              f"n_blocks={n_blocks} n_seeds={n_seeds} "
              f"blocks_parsed={len(blocks)} "
              f"iex_imbalance={summary_row['iex_imbalance']:.3f} "
              f"wall={summary_row['t_iex_max']:.1f}s")

    print(f"\nOutputs:")
    print(f"  {SUMMARY_CSV}  ({SUMMARY_CSV.stat().st_size if SUMMARY_CSV.exists() else 0} bytes)")
    print(f"  {BLOCK_CSV}    ({BLOCK_CSV.stat().st_size   if BLOCK_CSV.exists()   else 0} bytes)")


if __name__ == "__main__":
    main()
