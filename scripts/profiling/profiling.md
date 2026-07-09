# ParaFlow Profiling — Design & Phase Reference

Reference for ParaFlow's built-in timing/memory profiling: how it is switched on,
what each phase measures and why, and how the CPU and GPU paths are decomposed into
the **same** buckets so they compare directly.

> This is a *design* document. It contains no measured results — profiling numbers
> live in the generated CSVs/plots, not here.

---

## 1. The switch: `enable_timing`

```yaml
enable_timing: true    # emit TIMING/MEM lines to stderr; false = zero overhead
```

- `false` (default): every timer is a no-op via `pf_now`/`pf_accum`
  (`ParaFlow/timing.hpp`) — nothing printed, no measurable overhead.
- `true`: each rank writes `TIMING phase=… t=…` / `MEM_… ` to **stderr** (kept off
  stdout, where the binary trajectories go). Every emission is gated by this flag.

**Where the code lives**

| File | Role |
|---|---|
| `ParaFlow/timing.hpp` | `BlockTiming` accumulators, `pf_now`/`pf_accum`, `/proc` readers |
| `ParaFlow/ParaFlow.cpp` | the only file that prints `TIMING`/`MEM_` |
| `ParaFlow/OSUFlow/LoadTiming.h` + `MPASOReader.C` / `OSUFlow.C` | `block_load` sub-phase timers |
| `ParaFlow/OSUFlow/GPU/CUDA/MPASOGPUTracer.cpp` | GPU per-stage transfer timings (`GPUTimingBreakdown`) |
| `scripts/profiling/parse_timing.py` | logs → `run_summary.csv` + `block_detail.csv` |
| `scripts/profiling/plot_timing.py` | CSVs → plots |
| `scripts/profiling/parse_mem.py` | `MEM_*` → memory report |

---

## 2. Naming convention

- `…_rank` — one line per rank, that rank's own value.
- `…_summary` — one line (rank 0) with `min max avg imbalance` across ranks
  (`imbalance = max/min`). A summary, **not** a sum.

Per-block lines carry `rank=` and `gid=`. Production runs use `nproc == nblocks`
(1 block/rank), so `rank ↔ gid` is 1:1.

---

## 3. Phase catalog

Nested levels; times in seconds.

### Level 0 — whole-run & distributed

| Phase | Scope | Measures |
|---|---|---|
| `run_streamline_{rank,summary}` | rank / reduce | End-to-end streamline driver = load + exchange + write (the headline) |
| `run_pathline_{rank,summary}` | rank / reduce | Same for pathlines (load + all windows + write) |
| `trace_exchange_{rank,summary}` | rank / reduce | Wall inside the DIY `iexchange` loop (excludes load & write) |
| `block_load` / `block_load_summary` | gid / reduce | Initial data load (`set_data`); carries `nseeds=`; decomposed below |
| `output_write` | gid | `write_trajectory` |
| `run_seed_read` / `run_seed_bcast` | startup | Seed file read / broadcast |

### `block_load` decomposition (children sum to `block_load`)

Instrumented in OSUFlow, not DIY (§7). `load_read` vs `load_convert` answers
*is loading I/O-bound or compute-bound*.

| Phase | Measures | Scales with |
|---|---|---|
| `load_grid_build` | MPASOGrid + cell/vertex mappings | mesh size |
| `load_read` | NetCDF reads (velocity/zTop/vertVel), all timesteps | I/O |
| `load_convert` | cell→vertex conversion, all timesteps | cells × timesteps |
| `load_seed_filter` | `inBlock()` seed filtering | seed count |

Residual (`block_load − Σ children`) = `readAllTimestamps` + allocations + neighbor setup.

### Pathline windowing (pathline driver only, per window)

| Phase | Measures |
|---|---|
| `dist_trace_window` | Trace wall for one time window |
| `window_load` | Load next window's velocity (`advanceTimestep`) |
| `window_reinject` | Re-inject particles crossing the temporal boundary |

### Level 1 — per-block trace (the aligned buckets)

Accumulate over all `iexchange` iterations; **close to `trace_local_wall`** (§4).
CPU and GPU write the same field names.

| Phase | CPU | GPU |
|---|---|---|
| `trace_dequeue` | receive particles (outside `trace_local_wall`) | same |
| `trace_local_wall` | all local work per block (carries `nsteps=`, `nrecv=`) | same |
| `trace_prepare` | step-budget array (≈0) | `at_phys` point-location + launch arrays |
| `trace_transfer` | 0 (by construction) | Solution flatten + H2D/D2H + alloc/free |
| `trace_integrate_cpu` / `_gpu` | OSUFlow RK4 | CUDA kernel |
| `trace_postprocess` | build segments, classify exits | same |
| `trace_enqueue` | send particles to neighbors | same |

The parser merges `integrate_cpu + integrate_gpu` → one `t_trace_integrate` column
(the device suffix is kept in the log).

### Level 2 — `trace_transfer` drill-down (GPU only; children sum to parent)

Single axis, safe to add. All sub-millisecond stages are folded into `transfer_misc`.

| Phase | Measures | Distinguishes |
|---|---|---|
| `transfer_flatten` | host `flattenSolution` | flatten vs PCIe |
| `transfer_upload_topology` | topology H2D (one-time) | one-time vs recurring |
| `transfer_upload_velocity` | velocity window H2D (per window) | " |
| `transfer_misc` | alloc + particle H2D + result D2H + free | — |

`trace_transfer = flatten + upload_topology + upload_velocity + misc`.

### Memory lines

| Line | Measures |
|---|---|
| `MEM_ANALYTICAL` | analytical `grid_bytes` + `solution_bytes` |
| `MEM_DELTA` | VmRSS before/after `set_data` |
| `MEM_PEAK` | peak VmHWM |
| `MEM_CELLCOUNT` | `n_local_cells` / `n_global_cells` |

---

## 4. CPU / GPU granularity alignment (the core)

**Why they differ.** The CPU integrator locates each particle's cell *during*
integration, so prepare is trivial and there is no device transfer. The GPU must
flatten+upload the field and pre-locate each seed's cell (`at_phys`) before the
kernel — host-side phases with **no CPU analog** — while the kernel is a separately
isolated cost.

**Two adjustments make them comparable:**

1. **Unified integrate** — parser sums `integrate_cpu + integrate_gpu` (one is 0).
2. **`trace_transfer`** — GPU flatten+upload+download, summed from its Level-2
   children and defined **0 on CPU**. The kernel is excluded (it lives in
   `integrate`), so there is no double-counting.

**One taxonomy, both devices** — closes to `trace_local_wall`:

```
trace_local_wall ≈ trace_prepare + trace_transfer + trace_integrate
                 + trace_postprocess + trace_enqueue
```

- CPU: `transfer = 0` → `prepare(≈0) + integrate_cpu + postprocess + enqueue`.
- GPU: `prepare + transfer + kernel + postprocess + enqueue`.
- Any leftover → an explicit `trace_other` bucket at analysis time.

Identical buckets on both sides make a CPU-vs-GPU stacked bar apples-to-apples: you
see exactly which buckets move between devices.

---

## 5. CSV outputs

`parse_timing.py` → `run_summary.csv` (one row/run) + `block_detail.csv`
(one row/block); default dir `results/gpu_analysis/`, override `--outdir=`. Summary
column groups:

- Identity: `n_blocks, n_seeds, run_id, device, mesh, driver`
- Headline / exchange: `t_run_*`, `run_imbalance`, `t_trace_exchange_*`
- Aligned buckets: `t_prepare_*`, `t_transfer_*`, `t_integrate_*`, `t_postprocess_*`,
  `t_enqueue_*`, `t_blockload_*`, `t_output_write_*`
- `block_load` decomposition: `t_load_{grid_build,read,convert,seed_filter}_avg`
- Transfer drill-down: `t_transfer_{flatten,upload_topology,upload_velocity,misc}_*`
- Memory: `mem_*` · Work: `total_{seeds_assigned,steps,particles_received}`
- Imbalance: `work_max_over_avg`, `work_cov`, `throughput_cov`, `idle_frac_avg`,
  `straggler_gid`, and `block_load_*` (`block_load_cov`, `read_per_cell_cov`, …)

The parser accepts current and legacy phase names.

### Why these imbalance metrics

Exchange-wall `max/min` is misleading — `iexchange` terminates globally, so it reads
~1.0 even under heavy skew. Measure the underlying distributions instead:

- **`work_max_over_avg`** — Amdahl ceiling (slowest block's work ÷ ideal).
- **`throughput_cov`** (per-block steps/s) — uniform ⇒ imbalance is *work
  distribution* (re-partition by predicted steps); varying ⇒ *processing-speed*
  (NUMA / GPU sharing / cache), which re-partitioning won't fix.
- **`read_per_cell_cov`** (`load_read / n_local_cells`) — uniform ⇒ each rank reads
  only its share (a cell partition balances I/O); varying ⇒ ranks read the whole
  field, so partitioning can't cut I/O (need collective/hyperslab reads).

Plot `A5_imbalance_diagnosis` shows the work-vs-`local_wall` scatter behind these.

---

## 6. Workflow

```bash
bash NERSC_job/submit_highres_16_pair.sh                                   # → logs/*.stderr
python3 scripts/profiling/parse_timing.py --reset --outdir=results/cpu_gpu logs/*.stderr
python3 scripts/profiling/plot_timing.py results/cpu_gpu                    # needs pandas+matplotlib
python3 scripts/profiling/parse_mem.py logs/cpu.stderr logs/gpu.stderr      # optional
```

Name logs `bNNN_sNNNNNNN_runN.stderr` for auto metadata, or pass
`--blocks/--seeds/--run/--device/--mesh`. Run each config ≥3× and take the median to
average out filesystem I/O variance.

---

## 7. Design notes

- **One nested tree, sum freely.** `trace_*` buckets close to `trace_local_wall`;
  `transfer_*` children close to `trace_transfer`. No overlapping second axis.
- **`_summary` ≠ total.** `_rank` is one rank's value; `_summary` is
  min/max/avg/imbalance across ranks. `trace_exchange_*` covers the `iexchange` loop
  only — not the whole run (use `run_streamline`/`run_pathline`).
- **DIY is never instrumented.** The external DIY submodule is used only for
  `master.iexchange`; we time either side of it (`block_load`, `output_write`) or
  inside our own callback (`trace_*`), which we bracket. Time blocked *inside* DIY is
  not observable — it surfaces only as `idle = trace_exchange − trace_local_wall` and
  cannot be split into "waiting for particles" vs "finished, waiting for termination".
