# Streamline LIC (rendering/LIC)

Renders a line-integral-convolution (LIC) image from ParaFlow streamlines:
traces streamlines the same way `ParaFlow_streamline` does, then folds them
into a single grayscale PNG showing the flow pattern.

## Build

```bash
# GPU driver (what you'll normally use)
OSUFLOW_ENABLE_CUDA=1 CUDA_ARCH=90 bash build.sh
module load cuda/12.4.1
bash rendering/LIC/build_lic_gpu.sh

# CPU driver, if you need it instead
bash rendering/LIC/build_lic.sh
```

Produces `./ParaFlow_lic_gpu` (or `./ParaFlow_lic`) at the repo root. If you
pull new changes or edit anything under `ParaFlow/OSUFlow/`, rerun the first
`build.sh` command before rebuilding the LIC driver.

## Run one dataset

The launch flags differ per cluster — the binary this repo builds by default is
for **OSC Cardinal**.

### OSC Cardinal (MVAPICH + spack NetCDF)

```bash
export LD_LIBRARY_PATH=/apps/spack/0.21/cardinal/linux-rhel9-sapphirerapids/netcdf-cxx4/gcc/12.3.0/mvapich/3.0/4.3.1-tgr36hp/lib64:${LD_LIBRARY_PATH}
srun -n <nproc> --ntasks-per-node=4 --gpus-per-task=1 --mpi=pmi2 \
  ./ParaFlow_lic_gpu conf/lic_streamline_gpu.yaml
```

`LD_LIBRARY_PATH` is required (`build_lic_gpu.sh` links `libnetcdf-cxx4.so` with
no rpath; `jobs/osc_lic_batch_gpu.sh` already sets it). `<nproc>` must match
`nproc`/`nblocks` in the config. `--mpi=pmi2` is required — this binary links
Slurm's PMI client, not `libpmix` (confirm with `ldd ./ParaFlow_lic_gpu | grep
pmi`); with `--mpi=pmix` every rank becomes a singleton `rank=0 size=1` world
and aborts at the first cross-rank op with `Invalid rank ...`.

### NERSC Perlmutter (Cray MPICH + Cray NetCDF)

Rebuild for the Cray toolchain first (the default binary is Cardinal-linked and
won't run here):

```bash
./switch.sh nersc
OSUFLOW_ENABLE_CUDA=1 CUDA_ARCH=80 bash build.sh   # A100 = arch 80
bash rendering/LIC/build_lic_gpu.sh
```

Then launch — no `LD_LIBRARY_PATH` export (the NetCDF module sets it) and no
`--mpi` flag (Cray MPICH speaks Perlmutter's default PMI natively):

```bash
export MPICH_GPU_SUPPORT_ENABLED=0      # binary not linked with libmpi_gtl_cuda
srun -n <nproc> --ntasks-per-node=4 --gpus-per-task=1 --gpu-bind=closest \
  ./ParaFlow_lic_gpu conf/lic_streamline_gpu.yaml
```

## Run several datasets in one job

`conf/lic_batch_gpu.yaml` is a manifest — a list of config file paths.
`jobs/osc_lic_batch_gpu.sh` runs `ParaFlow_lic_gpu` on each one in turn:

```bash
bash jobs/osc_lic_batch_gpu.sh conf/lic_batch_gpu.yaml
```

To add a dataset: copy `conf/lic_streamline_gpu.yaml`, change `datafile_src`
(and `meshfile_src` if the mesh differs), and add the new file's path to the
manifest.

## Config options you'll actually touch (`conf/*.yaml`)

- **`datafile_src`** / **`meshfile_src`**: the velocity data and mesh files to render.
- **`lic_direction`**: `forward` | `backward` | `both`.
- **`lic_width`** / **`lic_height`**: output image resolution.
- **`nproc`** / **`nblocks`** / **`graph_partition_indices`**: how many
  ranks/GPUs to use — must match each other and your `srun -n`.
- **`streamline_storage`**: segments are gathered to rank 0 over MPI (no rank
  writes a `.bin` during tracing). `memory` (default) drops them once folded —
  only the final PNG; `disk` writes the reassembled trajectories as `<gid>.bin`
  in the run folder for later reuse (e.g. `scripts/plot/read_traces.py`).
- **`outputDir`**: where results land, in a subfolder named after the mode,
  direction, and dataset — so different datasets/directions don't overwrite
  each other, but two runs of the *same* dataset/direction will.
