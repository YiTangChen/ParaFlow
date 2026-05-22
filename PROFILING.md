# ParaFlow Timing and Profiling Guide

This file defines the current `TIMING phase=...` log names and how to use
them for CPU, GPU, and CPU+GPU profiling. All reported times are seconds.

## Timing Model

ParaFlow has two timing levels:

1. End-to-end distributed cost: what the MPI + DIY application actually pays.
2. Integration engine cost: the pure particle tracing compute time.

Do not compare GPU pipeline time directly against CPU RK4 time and call that
kernel speedup. The fair pure-compute comparison is:

```text
CPU pure integration: trace_integrate_cpu
GPU pure integration: trace_integrate_gpu
```

The fair application-level comparison is:

```text
CPU distributed run: dist_trace_total
GPU distributed run: dist_trace_total
```

## Run Startup Phases

```text
run_seed_read
run_seed_bcast
```

- `run_seed_read`: rank 0 reads the seed file.
- `run_seed_bcast`: MPI broadcast of seeds from rank 0 to all ranks.

These are startup costs. They should usually be excluded from pure tracing
speedup, but included in end-to-end workflow plots if seed I/O is part of the
claimed runtime.

## Per-Block Setup

```text
block_load
```

- `block_load`: block-local data setup, including NetCDF read, MPAS-O grid and
  solution construction, and seed filtering for the block.

Memory lines:

```text
MEM_ANALYTICAL
MEM_DELTA
MEM_PEAK
MEM_CELLCOUNT
```

- `MEM_ANALYTICAL`: estimated block grid and solution bytes.
- `MEM_DELTA`: process RSS before/after block load. This is process-level, not
  strictly per-block if one MPI rank owns multiple blocks.
- `MEM_PEAK`: process peak RSS after tracing.
- `MEM_CELLCOUNT`: block-local and global cell counts for normalization.

## Distributed Tracing Phases

```text
dist_trace_window
dist_trace_total
dist_trace_global
```

- `dist_trace_window`: one pathline temporal-window iteration. It includes local
  compute, DIY particle movement, waiting, and imbalance for that window.
- `dist_trace_total`: total distributed tracing wall time for the rank. This is
  the old `iexchange_total`, renamed because it is not pure communication.
- `dist_trace_global`: min, max, average, and imbalance ratio across ranks.

Use `dist_trace_total` for application speedup and strong/weak scaling.

## Per-Block Local Tracing Phases

```text
trace_dequeue
trace_local_wall
trace_prepare
trace_integrate_cpu
trace_integrate_gpu
trace_postprocess
trace_enqueue
```

- `trace_dequeue`: time spent dequeuing particles that DIY has already received.
  This is not full network communication time.
- `trace_local_wall`: total local block tracing time after dequeue. It includes
  integration, GPU pipeline overhead if enabled, postprocess, and enqueue.
- `trace_prepare`: CPU-side setup before integration, including remaining-step
  calculation, time arrays, seed cell hints, output arrays, and GPU launch-batch
  packing.
- `trace_integrate_cpu`: pure CPU OSUFlow integration call time.
- `trace_integrate_gpu`: pure CUDA kernel integration time measured with CUDA
  events. This is the GPU equivalent of CPU RK4 integration time.
- `trace_postprocess`: CPU work after integration: segment construction, step
  updates, final-cell handling, temporal-boundary classification, and handoff
  preparation.
- `trace_enqueue`: time spent enqueueing outgoing particles into DIY.

For CPU-only runs, expect `trace_integrate_cpu` to be nonzero and
`trace_integrate_gpu` to be zero. For GPU runs, expect the opposite.

## GPU Pipeline Phases

```text
gpu_pipeline_wall
gpu_host_prepare
gpu_upload_topology
gpu_upload_velocity
gpu_alloc
gpu_upload_particles
gpu_download_results
gpu_free
gpu_field_release
```

- `gpu_pipeline_wall`: full GPU wrapper wall time, excluding ParaFlow
  postprocess. This is the GPU-side cost ParaFlow pays before it can rebuild
  distributed segments.
- `gpu_host_prepare`: CPU-side preparation for GPU input, mainly metadata setup,
  velocity flattening, and timestamp array construction.
- `gpu_upload_topology`: copy static MPAS-O mesh/topology arrays to the GPU.
- `gpu_upload_velocity`: copy velocity, vertical velocity, zTop, and timestamps.
- `gpu_alloc`: per-launch CUDA allocation cost.
- `gpu_upload_particles`: host-to-device copy of seeds, cell hints, and
  per-particle control arrays.
- `trace_integrate_gpu`: the pure CUDA RK4/pathline kernel time. It is not
  repeated as `gpu_kernel` to avoid duplicate timing fields.
- `gpu_download_results`: device-to-host copy of traces, final cells, times, and
  step counts.
- `gpu_free`: per-launch CUDA free cost.
- `gpu_field_release`: release topology and velocity buffers owned by the GPU
  field object.

If `gpu_upload_topology` or `gpu_upload_velocity` dominates, the next
optimization should be persistent GPU field buffers across launches or windows.
If `gpu_alloc` and `gpu_free` dominate, use persistent per-block particle
buffers or a workspace allocator.

## Pathline Temporal Window Phases

```text
window_load
window_reinject
```

- `window_load`: cost of advancing to the next time window, including field
  loading performed by `advanceTimestep()`.
- `window_reinject`: cost of moving temporal-boundary particles back into the
  next tracing window.

These are important for high-resolution time-varying data, where I/O and
temporal-window movement can hide the benefit of a faster integrator.

## Recommended Plots

### Plot 1: End-to-End Distributed Runtime

Use stacked bars per run configuration:

```text
block_load | dist_trace_total | output_write
```

For pathlines, add temporal I/O:

```text
block_load | window_load | dist_trace_total | output_write
```

This plot answers: does the whole ParaFlow workflow get faster?

### Plot 2: Local Block Runtime Anatomy

Use stacked bars aggregated across blocks:

```text
trace_dequeue | trace_prepare | trace_integrate_cpu/gpu |
trace_postprocess | trace_enqueue
```

For GPU runs, use `trace_integrate_gpu` for the compute slice and show GPU
overheads separately in Plot 3. Include `trace_prepare` when comparing the
whole local block cost because GPU cell-hint preparation can be nontrivial.

This plot answers: where does local block time go?

### Plot 3: GPU Pipeline Anatomy

Use stacked bars:

```text
gpu_host_prepare | gpu_upload_topology | gpu_upload_velocity |
gpu_alloc | gpu_upload_particles | trace_integrate_gpu |
gpu_download_results | gpu_free | gpu_field_release | trace_postprocess
```

This plot answers: is the GPU kernel fast but hidden by data movement or setup?

### Plot 4: Pure Integration Speedup

Use a simple bar or line plot:

```text
speedup = trace_integrate_cpu / trace_integrate_gpu
```

Normalize by total RK4 steps when the particle count or termination behavior is
not identical:

```text
seconds_per_step = trace_integrate_* / nsteps
```

This plot answers: how much faster is the integration engine itself?

### Plot 5: Rank Imbalance

Use `dist_trace_global`:

```text
min | avg | max | imbalance = max / min
```

This plot answers: are some ranks idle while others are still tracing?

## Strong Scaling Design

Strong scaling keeps the global problem fixed and increases resources.

Recommended setup:

```text
Fixed:
  dataset
  total seeds
  max_steps
  time windows
  output interval
  total blocks, unless testing block-count sensitivity

Vary:
  MPI ranks / compute nodes
  GPUs used
```

Suggested runs:

```text
CPU-only:  1, 2, 4, 8, 16, 32 ranks
GPU:       1, 2, 4, 8, 16 GPUs, with one MPI rank per GPU
CPU+GPU:   same rank/GPU layout, with GPU enabled
```

Report:

```text
speedup(P) = T(1) / T(P)
efficiency(P) = speedup(P) / P
```

Use `dist_trace_total` for application scaling and
`trace_integrate_cpu/gpu` for pure integration scaling. If strong scaling
flattens, check `dist_trace_global` imbalance and the GPU pipeline overheads.

## Weak Scaling Design

Weak scaling increases resources and problem size together.

Recommended setup:

```text
Fixed per rank or per GPU:
  seeds per rank/GPU
  blocks per rank/GPU
  max_steps
  data resolution per block

Vary:
  ranks/GPUs and total seeds/blocks proportionally
```

Suggested GPU layout:

```text
1 MPI rank per physical GPU
1 to N DIY blocks per rank
```

Report:

```text
weak_efficiency(P) = T(1) / T(P)
```

For ideal weak scaling, `dist_trace_total` should stay roughly flat as the
problem grows. If it rises, inspect:

- `trace_local_wall`: local work increased unexpectedly.
- `trace_dequeue` and `trace_enqueue`: more particle exchange.
- `dist_trace_global`: rank imbalance.
- `window_load`: high-resolution temporal data loading.
- `gpu_upload_topology` and `gpu_upload_velocity`: repeated GPU data movement.

## Presentation Interpretation

A good GPU result has two different claims:

1. The GPU kernel accelerates pure RK4/pathline integration:
   `trace_integrate_cpu` vs `trace_integrate_gpu`.
2. The distributed CPU+GPU system is efficient only if upload, allocation,
   download, postprocess, and DIY imbalance do not dominate:
   `dist_trace_total` and the GPU pipeline breakdown.

This distinction is important. It lets you say the GPU algorithm is fast while
also showing honestly whether the full ParaFlow distributed workflow benefits.
