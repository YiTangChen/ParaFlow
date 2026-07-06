#!/usr/bin/env python3
"""
Parse ParaFlow timing + memory logs into CSV files for plotting.

This parser follows the newer timing names emitted by ParaFlow, including:
  run_seed_read, run_seed_bcast,
  block_load, trace_local_wall, trace_dequeue, trace_prepare,
  trace_integrate_cpu, trace_integrate_gpu, trace_postprocess, trace_enqueue,
  gpu_pipeline_wall, gpu_host_prepare, gpu_upload_topology, gpu_upload_velocity,
  gpu_alloc, gpu_upload_particles, gpu_download_results, gpu_free, gpu_field_release,
  output_write,
  dist_trace_window, window_load, window_reinject,
  dist_trace_total, dist_trace_global.

It also writes compatibility columns used by older plotting scripts
(e.g., t_compute_*, t_iex_*, iex_imbalance).

Filename convention for automatic metadata extraction:
  logs/b<nblocks>_s<nseeds>_run<N>.stderr

Usage:
  python3 parse_timing.py logs/b016_s0010000_run1.stderr
  python3 parse_timing.py --reset logs/*.stderr
  python3 parse_timing.py --blocks=1 --seeds=1000 --run=1 my.log
  python3 parse_timing.py --outdir=results/cpu_analysis logs/*.stderr
"""

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


DEFAULT_RESULTS_DIR = Path("results/gpu_analysis")
RESULTS_DIR = DEFAULT_RESULTS_DIR
SUMMARY_CSV = RESULTS_DIR / "run_summary.csv"
BLOCK_CSV = RESULTS_DIR / "block_detail.csv"


# Per-block time phases in current ParaFlow logs.
BLOCK_TIME_PHASES = {
    "block_load": "t_block_load",
    "trace_local_wall": "t_trace_local_wall",
    "trace_dequeue": "t_trace_dequeue",
    "trace_prepare": "t_trace_prepare",
    "trace_integrate_cpu": "t_trace_integrate_cpu",
    "trace_integrate_gpu": "t_trace_integrate_gpu",
    "trace_postprocess": "t_trace_postprocess",
    "trace_enqueue": "t_trace_enqueue",
    "output_write": "t_output_write",
    "gpu_pipeline_wall": "t_gpu_pipeline_wall",
    "gpu_host_prepare": "t_gpu_host_prepare",
    "gpu_upload_topology": "t_gpu_upload_topology",
    "gpu_upload_velocity": "t_gpu_upload_velocity",
    "gpu_alloc": "t_gpu_alloc",
    "gpu_upload_particles": "t_gpu_upload_particles",
    "gpu_download_results": "t_gpu_download_results",
    "gpu_free": "t_gpu_free",
    "gpu_field_release": "t_gpu_field_release",
}

GPU_BLOCK_KEYS = [
    "t_trace_integrate_gpu",
    "t_gpu_pipeline_wall",
    "t_gpu_host_prepare",
    "t_gpu_upload_topology",
    "t_gpu_upload_velocity",
    "t_gpu_alloc",
    "t_gpu_upload_particles",
    "t_gpu_download_results",
    "t_gpu_free",
    "t_gpu_field_release",
]


SUMMARY_FIELDS = [
    "n_blocks", "n_seeds", "run_id", "device", "mesh",

    # startup
    "t_run_seed_read",
    "t_run_seed_bcast_min", "t_run_seed_bcast_max", "t_run_seed_bcast_avg",

    # distributed totals (rank-level)
    "t_dist_trace_total_min", "t_dist_trace_total_max", "t_dist_trace_total_avg",
    "dist_trace_imbalance",

    # pathline window diagnostics (sum per rank)
    "t_dist_trace_window_sum_min", "t_dist_trace_window_sum_max", "t_dist_trace_window_sum_avg",
    "t_window_load_sum_min", "t_window_load_sum_max", "t_window_load_sum_avg",
    "t_window_reinject_sum_min", "t_window_reinject_sum_max", "t_window_reinject_sum_avg",

    # per-block phase stats
    "t_blockload_min", "t_blockload_max", "t_blockload_avg",
    "t_local_wall_min", "t_local_wall_max", "t_local_wall_avg",
    "t_dequeue_min", "t_dequeue_max", "t_dequeue_avg",
    "t_prepare_min", "t_prepare_max", "t_prepare_avg",
    "t_integrate_cpu_min", "t_integrate_cpu_max", "t_integrate_cpu_avg",
    "t_integrate_gpu_min", "t_integrate_gpu_max", "t_integrate_gpu_avg",
    "t_postprocess_min", "t_postprocess_max", "t_postprocess_avg",
    "t_enqueue_min", "t_enqueue_max", "t_enqueue_avg",
    "t_output_write_min", "t_output_write_max", "t_output_write_avg",

    # gpu pipeline stats
    "t_gpu_pipeline_wall_min", "t_gpu_pipeline_wall_max", "t_gpu_pipeline_wall_avg",
    "t_gpu_host_prepare_min", "t_gpu_host_prepare_max", "t_gpu_host_prepare_avg",
    "t_gpu_upload_topology_min", "t_gpu_upload_topology_max", "t_gpu_upload_topology_avg",
    "t_gpu_upload_velocity_min", "t_gpu_upload_velocity_max", "t_gpu_upload_velocity_avg",
    "t_gpu_alloc_min", "t_gpu_alloc_max", "t_gpu_alloc_avg",
    "t_gpu_upload_particles_min", "t_gpu_upload_particles_max", "t_gpu_upload_particles_avg",
    "t_gpu_download_results_min", "t_gpu_download_results_max", "t_gpu_download_results_avg",
    "t_gpu_free_min", "t_gpu_free_max", "t_gpu_free_avg",
    "t_gpu_field_release_min", "t_gpu_field_release_max", "t_gpu_field_release_avg",

    # memory
    "mem_grid_bytes_min", "mem_grid_bytes_max", "mem_grid_bytes_avg",
    "mem_solution_bytes_min", "mem_solution_bytes_max", "mem_solution_bytes_avg",
    "mem_delta_kb_min", "mem_delta_kb_max", "mem_delta_kb_avg",
    "mem_peak_vmhwm_kb_min", "mem_peak_vmhwm_kb_max", "mem_peak_vmhwm_kb_avg",

    # work counters
    "total_seeds_assigned", "total_steps", "total_particles_received",

    # compatibility aliases for older analysis scripts
    "t_seed_read",
    "t_seed_bcast_min", "t_seed_bcast_max", "t_seed_bcast_avg",
    "t_compute_min", "t_compute_max", "t_compute_avg",
    "t_gpukernel_min", "t_gpukernel_max", "t_gpukernel_avg",
    "t_comm_min", "t_comm_max", "t_comm_avg",
    "t_write_min", "t_write_max", "t_write_avg",
    "t_iex_min", "t_iex_max", "t_iex_avg", "iex_imbalance",
]


BLOCK_FIELDS = [
    "n_blocks", "n_seeds", "run_id", "device", "mesh", "gid", "rank",

    "t_block_load",
    "t_trace_local_wall",
    "t_trace_dequeue",
    "t_trace_prepare",
    "t_trace_integrate_cpu",
    "t_trace_integrate_gpu",
    "t_trace_postprocess",
    "t_trace_enqueue",
    "t_output_write",

    "t_gpu_pipeline_wall",
    "t_gpu_host_prepare",
    "t_gpu_upload_topology",
    "t_gpu_upload_velocity",
    "t_gpu_alloc",
    "t_gpu_upload_particles",
    "t_gpu_download_results",
    "t_gpu_free",
    "t_gpu_field_release",

    "t_dist_trace_total_rank",
    "t_dist_trace_window_sum_rank",
    "t_window_load_sum_rank",
    "t_window_reinject_sum_rank",
    "t_trace_idle",

    "n_seeds_initial", "n_steps_total", "n_particles_received",
    "n_local_cells", "n_global_cells",
    "mem_grid_bytes", "mem_solution_bytes", "mem_total_bytes",
    "mem_delta_kb", "mem_peak_vmhwm_kb",

    # compatibility aliases
    "t_trace_compute", "t_gpu_kernel", "t_trace_comm",
]


def parse_kv(line: str) -> dict:
    # Supports values like 1, -1, 3.14, 1e-5, inf.
    return dict(m.groups() for m in re.finditer(r"(\w+)=([^\s]+)", line))


def _safe_float(v, default=0.0):
    try:
        return float(v)
    except Exception:
        return default


def _safe_int(v, default=0):
    try:
        return int(v)
    except Exception:
        return default


def _stats(values):
    if not values:
        return 0.0, 0.0, 0.0
    return min(values), max(values), sum(values) / len(values)


def extract_n_blocks_n_seeds(path: str):
    name = Path(path).stem
    m = re.match(r"b(\d+)_s(\d+)_run(\d+)", name)
    if m:
        return int(m.group(1)), int(m.group(2)), int(m.group(3))
    return None, None, None


def parse_log(path: str):
    """
    Returns:
      global_data: dict
      blocks: dict[gid] -> dict
      has_gpu: bool
    """
    global_data = {
        "run_seed_read": [],
        "run_seed_bcast": [],
        "dist_trace_global": {},
        "dist_trace_total_per_rank": defaultdict(list),
        "dist_trace_window_per_rank": defaultdict(list),
        "window_load_per_rank": defaultdict(list),
        "window_reinject_per_rank": defaultdict(list),
    }
    blocks = defaultdict(dict)
    has_gpu = False

    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not (line.startswith("TIMING") or line.startswith("MEM_")):
                continue

            kv = parse_kv(line)
            phase = kv.get("phase", "")

            # ----- timing lines -----
            if phase == "run_seed_read":
                global_data["run_seed_read"].append(_safe_float(kv.get("t")))
                continue
            if phase == "run_seed_bcast":
                global_data["run_seed_bcast"].append(_safe_float(kv.get("t")))
                continue

            if phase == "dist_trace_total":
                rank = _safe_int(kv.get("rank", 0))
                global_data["dist_trace_total_per_rank"][rank].append(_safe_float(kv.get("t")))
                continue

            if phase == "dist_trace_window":
                rank = _safe_int(kv.get("rank", 0))
                global_data["dist_trace_window_per_rank"][rank].append(_safe_float(kv.get("t")))
                continue

            if phase == "window_load":
                rank = _safe_int(kv.get("rank", 0))
                global_data["window_load_per_rank"][rank].append(_safe_float(kv.get("t")))
                continue

            if phase == "window_reinject":
                rank = _safe_int(kv.get("rank", 0))
                global_data["window_reinject_per_rank"][rank].append(_safe_float(kv.get("t")))
                continue

            if phase == "dist_trace_global":
                # format: min=... max=... avg=... imbalance=...
                g = global_data["dist_trace_global"]
                g["min"] = _safe_float(kv.get("min"))
                g["max"] = _safe_float(kv.get("max"))
                g["avg"] = _safe_float(kv.get("avg"))
                g["imbalance"] = _safe_float(kv.get("imbalance"))
                continue

            if phase in BLOCK_TIME_PHASES:
                gid = _safe_int(kv.get("gid", -1))
                if gid < 0:
                    continue
                b = blocks[gid]
                b["gid"] = gid
                b["rank"] = _safe_int(kv.get("rank", b.get("rank", 0)))

                dst = BLOCK_TIME_PHASES[phase]
                b[dst] = b.get(dst, 0.0) + _safe_float(kv.get("t", 0.0))

                if phase == "block_load":
                    b["n_seeds_initial"] = _safe_int(kv.get("nseeds", b.get("n_seeds_initial", 0)))
                elif phase == "trace_local_wall":
                    b["n_steps_total"] = _safe_int(kv.get("nsteps", b.get("n_steps_total", 0)))
                    b["n_particles_received"] = _safe_int(kv.get("nrecv", b.get("n_particles_received", 0)))

                if dst in GPU_BLOCK_KEYS and b[dst] > 0.0:
                    has_gpu = True

                continue

            # ----- memory lines -----
            if line.startswith("MEM_ANALYTICAL"):
                gid = _safe_int(kv.get("gid", -1))
                if gid >= 0:
                    b = blocks[gid]
                    b["gid"] = gid
                    b["rank"] = _safe_int(kv.get("rank", b.get("rank", 0)))
                    b["mem_grid_bytes"] = _safe_int(kv.get("grid_bytes", 0))
                    b["mem_solution_bytes"] = _safe_int(kv.get("solution_bytes", 0))
                continue

            if line.startswith("MEM_DELTA"):
                gid = _safe_int(kv.get("gid", -1))
                if gid >= 0:
                    b = blocks[gid]
                    b["gid"] = gid
                    b["rank"] = _safe_int(kv.get("rank", b.get("rank", 0)))
                    b["mem_delta_kb"] = _safe_int(kv.get("delta_kb", 0))
                continue

            if line.startswith("MEM_CELLCOUNT"):
                gid = _safe_int(kv.get("gid", -1))
                if gid >= 0:
                    b = blocks[gid]
                    b["gid"] = gid
                    b["rank"] = _safe_int(kv.get("rank", b.get("rank", 0)))
                    b["n_local_cells"] = _safe_int(kv.get("n_local_cells", 0))
                    b["n_global_cells"] = _safe_int(kv.get("n_global_cells", 0))
                continue

            if line.startswith("MEM_PEAK"):
                gid = _safe_int(kv.get("gid", -1))
                if gid >= 0:
                    b = blocks[gid]
                    b["gid"] = gid
                    b["rank"] = _safe_int(kv.get("rank", b.get("rank", 0)))
                    vmhwm = kv.get("vmhwm_kb")
                    if vmhwm is None:
                        m = re.search(r"VmHWM:\s*(\d+)", line)
                        vmhwm = m.group(1) if m else "0"
                    b["mem_peak_vmhwm_kb"] = _safe_int(vmhwm, 0)
                continue

    # per-rank sums
    rank_dist_total = {
        rank: sum(vals) for rank, vals in global_data["dist_trace_total_per_rank"].items()
    }
    rank_dist_window = {
        rank: sum(vals) for rank, vals in global_data["dist_trace_window_per_rank"].items()
    }
    rank_window_load = {
        rank: sum(vals) for rank, vals in global_data["window_load_per_rank"].items()
    }
    rank_window_reinject = {
        rank: sum(vals) for rank, vals in global_data["window_reinject_per_rank"].items()
    }

    for gid, b in blocks.items():
        rank = b.get("rank", 0)
        b["t_dist_trace_total_rank"] = rank_dist_total.get(rank, 0.0)
        b["t_dist_trace_window_sum_rank"] = rank_dist_window.get(rank, 0.0)
        b["t_window_load_sum_rank"] = rank_window_load.get(rank, 0.0)
        b["t_window_reinject_sum_rank"] = rank_window_reinject.get(rank, 0.0)

        t_local = b.get("t_trace_local_wall", 0.0)
        t_total = b.get("t_dist_trace_total_rank", 0.0)
        b["t_trace_idle"] = max(0.0, t_total - t_local)

    # explicit dist_trace_global if present; otherwise compute from per-rank totals
    if not global_data["dist_trace_global"]:
        vals = list(rank_dist_total.values())
        if vals:
            mn, mx, av = _stats(vals)
            global_data["dist_trace_global"] = {
                "min": mn,
                "max": mx,
                "avg": av,
                "imbalance": (mx / mn) if mn > 0 else 0.0,
            }
        else:
            global_data["dist_trace_global"] = {
                "min": 0.0,
                "max": 0.0,
                "avg": 0.0,
                "imbalance": 0.0,
            }

    global_data["rank_dist_totals"] = list(rank_dist_total.values())
    global_data["rank_dist_window_sums"] = list(rank_dist_window.values())
    global_data["rank_window_load_sums"] = list(rank_window_load.values())
    global_data["rank_window_reinject_sums"] = list(rank_window_reinject.values())

    return global_data, blocks, has_gpu


def build_summary_row(n_blocks, n_seeds, run_id, device, mesh, global_data, blocks):
    bdata = list(blocks.values())

    def col(field):
        return [b.get(field, 0.0) for b in bdata]

    # startup
    seed_read = max(global_data.get("run_seed_read", [0.0])) if global_data.get("run_seed_read") else 0.0
    bcast_min, bcast_max, bcast_avg = _stats(global_data.get("run_seed_bcast", []))

    # rank-level distributed totals
    dist_rank_vals = global_data.get("rank_dist_totals", [])
    dist_total_min, dist_total_max, dist_total_avg = _stats(dist_rank_vals)

    # pathline window sums
    wdist_min, wdist_max, wdist_avg = _stats(global_data.get("rank_dist_window_sums", []))
    wload_min, wload_max, wload_avg = _stats(global_data.get("rank_window_load_sums", []))
    wreinj_min, wreinj_max, wreinj_avg = _stats(global_data.get("rank_window_reinject_sums", []))

    # per-block phases
    bl_min, bl_max, bl_avg = _stats(col("t_block_load"))
    lw_min, lw_max, lw_avg = _stats(col("t_trace_local_wall"))
    dq_min, dq_max, dq_avg = _stats(col("t_trace_dequeue"))
    pr_min, pr_max, pr_avg = _stats(col("t_trace_prepare"))
    ic_min, ic_max, ic_avg = _stats(col("t_trace_integrate_cpu"))
    ig_min, ig_max, ig_avg = _stats(col("t_trace_integrate_gpu"))
    po_min, po_max, po_avg = _stats(col("t_trace_postprocess"))
    en_min, en_max, en_avg = _stats(col("t_trace_enqueue"))
    ow_min, ow_max, ow_avg = _stats(col("t_output_write"))

    gpw_min, gpw_max, gpw_avg = _stats(col("t_gpu_pipeline_wall"))
    ghp_min, ghp_max, ghp_avg = _stats(col("t_gpu_host_prepare"))
    gut_min, gut_max, gut_avg = _stats(col("t_gpu_upload_topology"))
    guv_min, guv_max, guv_avg = _stats(col("t_gpu_upload_velocity"))
    gal_min, gal_max, gal_avg = _stats(col("t_gpu_alloc"))
    gup_min, gup_max, gup_avg = _stats(col("t_gpu_upload_particles"))
    gdr_min, gdr_max, gdr_avg = _stats(col("t_gpu_download_results"))
    gfr_min, gfr_max, gfr_avg = _stats(col("t_gpu_free"))
    grel_min, grel_max, grel_avg = _stats(col("t_gpu_field_release"))

    mg_min, mg_max, mg_avg = _stats(col("mem_grid_bytes"))
    ms_min, ms_max, ms_avg = _stats(col("mem_solution_bytes"))
    md_min, md_max, md_avg = _stats(col("mem_delta_kb"))
    mp_min, mp_max, mp_avg = _stats(col("mem_peak_vmhwm_kb"))

    gdist = global_data.get("dist_trace_global", {})
    imbalance = _safe_float(gdist.get("imbalance", 0.0))

    row = {
        "n_blocks": n_blocks,
        "n_seeds": n_seeds,
        "run_id": run_id,
        "device": device,
        "mesh": mesh,

        "t_run_seed_read": seed_read,
        "t_run_seed_bcast_min": bcast_min,
        "t_run_seed_bcast_max": bcast_max,
        "t_run_seed_bcast_avg": bcast_avg,

        "t_dist_trace_total_min": dist_total_min,
        "t_dist_trace_total_max": dist_total_max,
        "t_dist_trace_total_avg": dist_total_avg,
        "dist_trace_imbalance": imbalance,

        "t_dist_trace_window_sum_min": wdist_min,
        "t_dist_trace_window_sum_max": wdist_max,
        "t_dist_trace_window_sum_avg": wdist_avg,
        "t_window_load_sum_min": wload_min,
        "t_window_load_sum_max": wload_max,
        "t_window_load_sum_avg": wload_avg,
        "t_window_reinject_sum_min": wreinj_min,
        "t_window_reinject_sum_max": wreinj_max,
        "t_window_reinject_sum_avg": wreinj_avg,

        "t_blockload_min": bl_min,
        "t_blockload_max": bl_max,
        "t_blockload_avg": bl_avg,
        "t_local_wall_min": lw_min,
        "t_local_wall_max": lw_max,
        "t_local_wall_avg": lw_avg,
        "t_dequeue_min": dq_min,
        "t_dequeue_max": dq_max,
        "t_dequeue_avg": dq_avg,
        "t_prepare_min": pr_min,
        "t_prepare_max": pr_max,
        "t_prepare_avg": pr_avg,
        "t_integrate_cpu_min": ic_min,
        "t_integrate_cpu_max": ic_max,
        "t_integrate_cpu_avg": ic_avg,
        "t_integrate_gpu_min": ig_min,
        "t_integrate_gpu_max": ig_max,
        "t_integrate_gpu_avg": ig_avg,
        "t_postprocess_min": po_min,
        "t_postprocess_max": po_max,
        "t_postprocess_avg": po_avg,
        "t_enqueue_min": en_min,
        "t_enqueue_max": en_max,
        "t_enqueue_avg": en_avg,
        "t_output_write_min": ow_min,
        "t_output_write_max": ow_max,
        "t_output_write_avg": ow_avg,

        "t_gpu_pipeline_wall_min": gpw_min,
        "t_gpu_pipeline_wall_max": gpw_max,
        "t_gpu_pipeline_wall_avg": gpw_avg,
        "t_gpu_host_prepare_min": ghp_min,
        "t_gpu_host_prepare_max": ghp_max,
        "t_gpu_host_prepare_avg": ghp_avg,
        "t_gpu_upload_topology_min": gut_min,
        "t_gpu_upload_topology_max": gut_max,
        "t_gpu_upload_topology_avg": gut_avg,
        "t_gpu_upload_velocity_min": guv_min,
        "t_gpu_upload_velocity_max": guv_max,
        "t_gpu_upload_velocity_avg": guv_avg,
        "t_gpu_alloc_min": gal_min,
        "t_gpu_alloc_max": gal_max,
        "t_gpu_alloc_avg": gal_avg,
        "t_gpu_upload_particles_min": gup_min,
        "t_gpu_upload_particles_max": gup_max,
        "t_gpu_upload_particles_avg": gup_avg,
        "t_gpu_download_results_min": gdr_min,
        "t_gpu_download_results_max": gdr_max,
        "t_gpu_download_results_avg": gdr_avg,
        "t_gpu_free_min": gfr_min,
        "t_gpu_free_max": gfr_max,
        "t_gpu_free_avg": gfr_avg,
        "t_gpu_field_release_min": grel_min,
        "t_gpu_field_release_max": grel_max,
        "t_gpu_field_release_avg": grel_avg,

        "mem_grid_bytes_min": mg_min,
        "mem_grid_bytes_max": mg_max,
        "mem_grid_bytes_avg": mg_avg,
        "mem_solution_bytes_min": ms_min,
        "mem_solution_bytes_max": ms_max,
        "mem_solution_bytes_avg": ms_avg,
        "mem_delta_kb_min": md_min,
        "mem_delta_kb_max": md_max,
        "mem_delta_kb_avg": md_avg,
        "mem_peak_vmhwm_kb_min": mp_min,
        "mem_peak_vmhwm_kb_max": mp_max,
        "mem_peak_vmhwm_kb_avg": mp_avg,

        "total_seeds_assigned": int(sum(col("n_seeds_initial"))),
        "total_steps": int(sum(col("n_steps_total"))),
        "total_particles_received": int(sum(col("n_particles_received"))),

        # compatibility aliases
        "t_seed_read": seed_read,
        "t_seed_bcast_min": bcast_min,
        "t_seed_bcast_max": bcast_max,
        "t_seed_bcast_avg": bcast_avg,
        "t_compute_min": lw_min,
        "t_compute_max": lw_max,
        "t_compute_avg": lw_avg,
        "t_gpukernel_min": ig_min,
        "t_gpukernel_max": ig_max,
        "t_gpukernel_avg": ig_avg,
        "t_comm_min": en_min,
        "t_comm_max": en_max,
        "t_comm_avg": en_avg,
        "t_write_min": ow_min,
        "t_write_max": ow_max,
        "t_write_avg": ow_avg,
        "t_iex_min": dist_total_min,
        "t_iex_max": dist_total_max,
        "t_iex_avg": dist_total_avg,
        "iex_imbalance": imbalance,
    }

    return row


def build_block_rows(n_blocks, n_seeds, run_id, device, mesh, blocks):
    rows = []
    for gid, b in sorted(blocks.items()):
        row = {
            "n_blocks": n_blocks,
            "n_seeds": n_seeds,
            "run_id": run_id,
            "device": device,
            "mesh": mesh,
            "gid": gid,
            "rank": b.get("rank", 0),

            "t_block_load": b.get("t_block_load", 0.0),
            "t_trace_local_wall": b.get("t_trace_local_wall", 0.0),
            "t_trace_dequeue": b.get("t_trace_dequeue", 0.0),
            "t_trace_prepare": b.get("t_trace_prepare", 0.0),
            "t_trace_integrate_cpu": b.get("t_trace_integrate_cpu", 0.0),
            "t_trace_integrate_gpu": b.get("t_trace_integrate_gpu", 0.0),
            "t_trace_postprocess": b.get("t_trace_postprocess", 0.0),
            "t_trace_enqueue": b.get("t_trace_enqueue", 0.0),
            "t_output_write": b.get("t_output_write", 0.0),

            "t_gpu_pipeline_wall": b.get("t_gpu_pipeline_wall", 0.0),
            "t_gpu_host_prepare": b.get("t_gpu_host_prepare", 0.0),
            "t_gpu_upload_topology": b.get("t_gpu_upload_topology", 0.0),
            "t_gpu_upload_velocity": b.get("t_gpu_upload_velocity", 0.0),
            "t_gpu_alloc": b.get("t_gpu_alloc", 0.0),
            "t_gpu_upload_particles": b.get("t_gpu_upload_particles", 0.0),
            "t_gpu_download_results": b.get("t_gpu_download_results", 0.0),
            "t_gpu_free": b.get("t_gpu_free", 0.0),
            "t_gpu_field_release": b.get("t_gpu_field_release", 0.0),

            "t_dist_trace_total_rank": b.get("t_dist_trace_total_rank", 0.0),
            "t_dist_trace_window_sum_rank": b.get("t_dist_trace_window_sum_rank", 0.0),
            "t_window_load_sum_rank": b.get("t_window_load_sum_rank", 0.0),
            "t_window_reinject_sum_rank": b.get("t_window_reinject_sum_rank", 0.0),
            "t_trace_idle": b.get("t_trace_idle", 0.0),

            "n_seeds_initial": b.get("n_seeds_initial", 0),
            "n_steps_total": b.get("n_steps_total", 0),
            "n_particles_received": b.get("n_particles_received", 0),
            "n_local_cells": b.get("n_local_cells", 0),
            "n_global_cells": b.get("n_global_cells", 0),
            "mem_grid_bytes": b.get("mem_grid_bytes", 0),
            "mem_solution_bytes": b.get("mem_solution_bytes", 0),
            "mem_total_bytes": b.get("mem_grid_bytes", 0) + b.get("mem_solution_bytes", 0),
            "mem_delta_kb": b.get("mem_delta_kb", 0),
            "mem_peak_vmhwm_kb": b.get("mem_peak_vmhwm_kb", 0),

            # compatibility aliases
            "t_trace_compute": b.get("t_trace_local_wall", 0.0),
            "t_gpu_kernel": b.get("t_trace_integrate_gpu", 0.0),
            "t_trace_comm": b.get("t_trace_enqueue", 0.0),
        }
        rows.append(row)

    return rows


def write_csv(path, fields, rows, append=True):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not path.exists() or not append
    mode = "a" if append else "w"
    with open(path, mode, newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        if write_header:
            w.writeheader()
        w.writerows(rows)


def main():
    global RESULTS_DIR, SUMMARY_CSV, BLOCK_CSV

    args = sys.argv[1:]
    reset = False
    override_blocks = None
    override_seeds = None
    override_run = None
    override_device = None
    override_mesh = None

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
        elif a.startswith("--outdir="):
            RESULTS_DIR = Path(a.split("=", 1)[1])
            SUMMARY_CSV = RESULTS_DIR / "run_summary.csv"
            BLOCK_CSV = RESULTS_DIR / "block_detail.csv"
        else:
            filtered.append(a)
    args = filtered

    if not args:
        print(__doc__)
        raise SystemExit(1)

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
                n_seeds = override_seeds
                run_id = override_run if override_run is not None else 1
            else:
                print(f"WARN: cannot parse n_blocks/n_seeds from '{path}'")
                print("      Use --blocks=N --seeds=N [--run=N] to override.")
                continue

        global_data, blocks, has_gpu = parse_log(path)

        device = override_device if override_device else ("gpu" if has_gpu else "cpu")
        mesh = override_mesh if override_mesh else "lowres"

        summary_row = build_summary_row(n_blocks, n_seeds, run_id, device, mesh, global_data, blocks)
        block_rows = build_block_rows(n_blocks, n_seeds, run_id, device, mesh, blocks)

        write_csv(SUMMARY_CSV, SUMMARY_FIELDS, [summary_row], append=not reset)
        write_csv(BLOCK_CSV, BLOCK_FIELDS, block_rows, append=not reset)
        reset = False

        print(
            f"[{Path(path).name}] device={device} mesh={mesh} "
            f"n_blocks={n_blocks} n_seeds={n_seeds} blocks_parsed={len(blocks)} "
            f"dist_imbalance={summary_row['dist_trace_imbalance']:.3f} "
            f"dist_max={summary_row['t_dist_trace_total_max']:.3f}s"
        )

    print("\nOutputs:")
    print(f"  {SUMMARY_CSV}  ({SUMMARY_CSV.stat().st_size if SUMMARY_CSV.exists() else 0} bytes)")
    print(f"  {BLOCK_CSV}    ({BLOCK_CSV.stat().st_size if BLOCK_CSV.exists() else 0} bytes)")


if __name__ == "__main__":
    main()
