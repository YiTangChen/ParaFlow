# ParaFlow Timing and Profiling Guide

This file defines the current `TIMING phase=...` log names and how to use
them for CPU-only and CPU+GPU profiling. All reported times are seconds.

## Timing Model

ParaFlow has two application modes:

1. CPU-only: MPI + DIY distributed blocks, with OSUFlow integration on CPU.
2. CPU+GPU: MPI + DIY distributed blocks, with local integration offloaded to
   CUDA while CPU still handles setup, downloads, segments, and DIY handoff.

ParaFlow also has two timing levels:

1. End-to-end distributed cost: what the full MPI + DIY application pays.
2. Integration engine cost: the particle tracing compute kernel/call only.

The pure-compute comparison is:

```text
CPU integration call: trace_integrate_cpu
GPU integration kernel inside CPU+GPU mode: trace_integrate_gpu
```

The application-level comparison is:

```text
CPU-only distributed run: dist_trace_total
CPU+GPU distributed run: dist_trace_total
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
`trace_integrate_gpu` to be zero. For CPU+GPU runs, expect the opposite.

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
  In the production ParaFlow GPU path this should be nonzero only the first
  time a block creates its persistent GPU context.
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
  field object. Persistent block contexts release this at block destruction, so
  it may not appear in the per-block log output.

If `gpu_upload_topology` or `gpu_upload_velocity` dominates, the next
optimization should be persistent velocity/window buffers or fewer window
uploads. Topology is already persistent in the production GPU path.
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

## Results Analysis Plots

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

For CPU+GPU runs, use `trace_integrate_gpu` for the compute slice and show GPU
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

### Plot 4: Integration Speedup

Use a simple bar or line plot:

```text
speedup = trace_integrate_cpu / trace_integrate_gpu
```

Normalize by total RK4 steps when the particle count or termination behavior is
not identical:

```text
seconds_per_step = trace_integrate_* / nsteps
```

This plot answers: how much faster is the GPU integration kernel inside
CPU+GPU mode than the CPU OSUFlow integration call?

### Plot 5: Rank Imbalance

Use `dist_trace_global`:

```text
min | avg | max | imbalance = max / min
```

This plot answers: are some ranks idle while others are still tracing?

## Strong Scaling Design

Strong scaling keeps the global science problem fixed and increases resources.
For the current high-resolution MPAS-O pathline data, do not use toy 1-block or
2-block tests. Existing NERSC scripts show that the high-res case needs at
least the 16-block layout to be practical:

```text
dataset: highres MPAS-O, 3.7M cells, 80 vertical levels
seed file: seeds/random_seeds_10k.bin
time-varying pathline data: 36 monthly files
loadNTimeSteps: 16
dt: 60.0
record_interval: 60
minimum practical highres layout: 16 blocks / 16 ranks
```

There are only two application modes:

```text
CPU-only: useGPU: false
CPU+GPU:  useGPU: true, one MPI rank per physical GPU
```

For strong scaling, report two related but different stories:

```text
application scaling: dist_trace_total
integration scaling: trace_integrate_cpu vs trace_integrate_gpu
```

Application speedup:

```text
speedup(P) = T_16 / T_P
efficiency(P) = speedup(P) / (P / 16)
```

Here `T_16` is the 16-rank or 16-GPU baseline because smaller high-res layouts
are not realistic for this dataset.

### Recommended High-Res CPU-Only Strong Scaling

Use the existing CPU memory-sweep scripts as the CPU-only strong-scaling series:

```text
16 ranks, 16 blocks, 4 CPU nodes
32 ranks, 32 blocks, 4 CPU nodes
64 ranks, 64 blocks, 4 CPU nodes
128 ranks, 128 blocks, 8 CPU nodes
256 ranks, 256 blocks, 16 CPU nodes
```

Existing files:

```text
NERSC_job/highres_16_cpu_10k.sh
NERSC_job/highres_32_cpu_10k.sh
NERSC_job/highres_64_cpu_10k.sh
NERSC_job/highres_128_cpu_10k.sh
NERSC_job/highres_256_cpu_10k.sh
```

These runs change both rank count and partition count. That is acceptable for
this code because each rank owns block-local data; strong scaling at high-res
must reduce block size enough to fit memory.

### Recommended High-Res CPU+GPU Strong Scaling

Create matching GPU scripts/configs for:

```text
16 ranks, 16 blocks, 16 GPUs, 4 GPU nodes
32 ranks, 32 blocks, 32 GPUs, 8 GPU nodes
64 ranks, 64 blocks, 64 GPUs, 16 GPU nodes
128 ranks, 128 blocks, 128 GPUs, 32 GPU nodes
256 ranks, 256 blocks, 256 GPUs, 64 GPU nodes
```

GPU launch rule:

```text
srun -n P -c 16 --gpus-per-task=1 --gpu-bind=map_gpu:0,1,2,3
```

Each GPU YAML should match the CPU YAML except:

```text
outputDir: pathlines/nersc_highres_${P}_gpu_10k
useGPU: true
nproc: P
nblocks: P
graph_partition_indices: partition/highres_6y.info.part.${P}
```

The 16-GPU version already exists:

```text
NERSC_job/highres_16_gpu_10k.sh
NERSC_job/nersc_highres_16_gpu_10k.yaml
```

For 32/64/128/256 GPU runs, copy the corresponding CPU YAML and script, change
`-C cpu` to `-C gpu`, set `useGPU: true`, use 4 ranks per GPU node, and keep
one GPU per rank.

### Strong Scaling Plots

Plot these first:

```text
1. dist_trace_total vs rank/GPU count
2. speedup relative to 16-rank or 16-GPU baseline
3. efficiency relative to 16-rank or 16-GPU baseline
4. dist_trace_global min/avg/max and imbalance
```

Then diagnose with:

```text
CPU-only: trace_integrate_cpu, trace_prepare, trace_postprocess, trace_dequeue, trace_enqueue
CPU+GPU: trace_integrate_gpu, gpu_pipeline_wall, gpu_upload_velocity,
         gpu_alloc, gpu_download_results, trace_prepare, trace_postprocess
```

Expected risks:

```text
16 blocks: enough work per block, but high memory pressure
256 blocks: much smaller blocks, but more DIY boundary exchange and possible GPU underutilization
GPU: many small blocks/seeds can make upload and CPU postprocess dominate kernel time
```

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
- `gpu_upload_velocity`: repeated GPU field-window movement.
- `gpu_upload_topology`: should remain one-time per block. If it scales with
  windows, the persistent GPU context is not being reused.

## Interpretation

A good GPU result has two different claims:

1. The GPU kernel accelerates pure RK4/pathline integration:
   `trace_integrate_cpu` vs `trace_integrate_gpu`.
2. The distributed CPU+GPU system is efficient only if upload, allocation,
   download, postprocess, and DIY imbalance do not dominate:
   `dist_trace_total` and the GPU pipeline breakdown.

This distinction is important. It lets you say the GPU algorithm is fast while
also showing honestly whether the full ParaFlow distributed workflow benefits.

## GPU Testing Plan

Use this order before collecting final experiment numbers. The point is to
avoid mixing correctness debugging, GPU pipeline debugging, and scaling results.

### Stage 0: Build and Environment Sanity

Goal: prove the binary is the CUDA build you think it is.

Run:

```bash
OSUFLOW_ENABLE_CUDA=1 JOBS=8 ./build.sh
```

Check:

- The build reports `OSUFlow: CUDA backend ENABLED`.
- `ParaFlow_streamline` and `ParaFlow_pathline` are rebuilt.
- The job script requests one MPI rank per physical GPU for GPU experiments.
- Each GPU run prints `GPU_RANK_MAP` lines with rank, block gid,
  `CUDA_VISIBLE_DEVICES`, selected CUDA device, visible device count, and GPU
  name.

### Stage 1: High-Res Smoke and Parity Correctness

Goal: verify the GPU kernel matches the CPU reference before running expensive
scaling jobs. For the high-resolution MPAS-O data, do not use 1-block or
2-block toy tests. Use the smallest practical high-res layout: 16 blocks,
16 ranks, 16 GPUs.

Use the built-in parity switches:

```bash
OSUFLOW_GPU_SMOKE=1
OSUFLOW_GPU_PATHLINE_SMOKE=1
```

Recommended correctness settings:

```text
dataset: highres
blocks: 16
MPI ranks: 16
GPUs: 16
seed file: seeds/random_seeds_10k.bin
enable_timing: true
smoke seeds: controlled by OSUFLOW_GPU_SMOKE_SEEDS if needed
smoke steps: controlled by OSUFLOW_GPU_SMOKE_STEPS if needed
```

Pass criteria:

- The smoke test reports a nonzero comparison count.
- Maximum relative error is below the configured tolerance.
- CPU and GPU step counts match for pathlines.

The smoke tests now use the same persistent context API as production:

```text
CreateGPUBlockContext
UploadGPUVelocityWindow
Trace*OnGPUContext
DestroyGPUBlockContext
```

There is no one-shot `TraceParticles` or `TracePathlineBatch` wrapper path left.

### Stage 2: Persistent Topology Verification

Goal: prove topology upload happens once per block, not once per batch/window.

Run:

```text
NERSC_job/highres_16_gpu_10k.sh
enable_timing: true
```

Check log behavior:

- `gpu_upload_topology` should be nonzero once per GPU block.
- `gpu_upload_velocity` may appear every block/window because velocity changes.
- `trace_integrate_gpu` should accumulate across launches.
- `gpu_pipeline_wall` should be larger than `trace_integrate_gpu`, but topology
  should no longer scale with the number of windows or chunks.

Recommended diagnostic plot:

```text
gpu_upload_topology vs number_of_windows
```

Expected result:

```text
flat per block, not growing with windows
```

### Stage 3: CPU-Only vs CPU+GPU Integration

Goal: compare the integration engines while remembering that the only full
application modes are CPU-only and CPU+GPU.

Run matched CPU and CPU+GPU cases:

```text
same dataset
same blocks
same seeds
same max_steps
same interval
```

Use the existing first pair:

```text
NERSC_job/highres_16_cpu_10k.sh
NERSC_job/highres_16_gpu_10k.sh
```

Compare:

```text
trace_integrate_cpu / nsteps
trace_integrate_gpu / nsteps
```

Report the kernel/call comparison:

```text
integration speedup = trace_integrate_cpu / trace_integrate_gpu
```

This answers: how much faster is the GPU RK4/pathline kernel than the CPU
OSUFlow integration call? It is not an application speedup number.

### Stage 4: GPU Pipeline Anatomy

Goal: find the overhead limiting full application speedup.

Use a high-resolution dataset and plot:

```text
gpu_host_prepare
gpu_upload_topology
gpu_upload_velocity
gpu_alloc
gpu_upload_particles
trace_integrate_gpu
gpu_download_results
gpu_free
trace_postprocess
```

Interpretation:

- Large `gpu_upload_velocity`: time-varying field movement dominates.
- Large `gpu_alloc` / `gpu_free`: persistent particle workspace should be next.
- Large `trace_postprocess`: CPU segment construction or handoff dominates.
- Large `trace_prepare`: CPU cell-hint setup is expensive.

### Stage 5: End-to-End CPU vs CPU+GPU

Goal: measure real ParaFlow benefit.

Compare CPU-only and CPU+GPU runs using:

```text
dist_trace_total
output_write
block_load
window_load
```

Report:

```text
application speedup = CPU-only dist_trace_total / CPU+GPU dist_trace_total
```

Do not use `trace_integrate_gpu` alone for application speedup. That only
measures the kernel.

### Stage 6: Oversubscription Check

Goal: ensure each MPI rank maps cleanly to GPU hardware.

For CPU+GPU runs, prefer:

```text
1 MPI rank per physical GPU
1 OpenMP thread per rank
1 to N DIY blocks per rank
```

Avoid:

```text
more MPI ranks than physical GPUs unless deliberately testing oversubscription
```

Initial high-res GPU mapping check:

```text
4 nodes, 16 GPUs, 16 ranks
```

Then use the strong-scaling matrix in Stage 7.

Check:

- `dist_trace_global imbalance` should not explode.
- `gpu_pipeline_wall` should not increase unexpectedly at fixed work per GPU.
- GPU memory should fit for the selected blocks per rank.

### Stage 7: High-Res Strong Scaling

Goal: fixed high-res pathline problem, increasing ranks/blocks/GPUs from the
practical 16-block baseline.

CPU-only runs:

```text
16, 32, 64, 128, 256 ranks
```

CPU+GPU runs:

```text
16, 32, 64, 128, 256 ranks/GPUs
```

Plot:

```text
dist_trace_total vs ranks/GPUs
speedup relative to 16
parallel efficiency relative to 16
dist_trace_global min/avg/max
GPU pipeline anatomy for each GPU scale
```

If scaling flattens, diagnose with:

```text
trace_prepare
trace_postprocess
trace_enqueue
trace_dequeue
gpu_upload_velocity
gpu_alloc
gpu_download_results
dist_trace_global imbalance
```

### Stage 8: Weak Scaling

Goal: fixed work per GPU, increasing total problem size.

Keep fixed per GPU:

```text
seeds per GPU
blocks per GPU
max_steps
time windows per block
```

Vary:

```text
GPUs/ranks: 16, 32, 64, 128, 256
total seeds scale with GPUs
total blocks scale with GPUs
```

Plot:

```text
dist_trace_total vs GPUs
weak scaling efficiency
gpu_pipeline_wall per block
trace_integrate_gpu per step
```

Expected ideal behavior:

```text
dist_trace_total stays roughly flat
```

### Stage 9: Stress and Failure Tests

Goal: find memory limits and bad scheduling choices before production runs.

Stress dimensions:

```text
blocks per GPU: 1, 2, 4, 8, ...
seeds per block: small, medium, large
save interval: sparse vs dense output
pathline windows: short vs long temporal runs
```

Watch:

- GPU out-of-memory errors.
- `gpu_upload_velocity` growing with high-resolution temporal windows.
- `gpu_download_results` growing when `save_interval` is too small.
- `trace_postprocess` growing when many trajectory points are saved.

Final choice:

```text
Use the largest blocks-per-GPU setting that does not hurt dist_trace_total,
does not cause GPU memory pressure, and keeps imbalance reasonable.
```
