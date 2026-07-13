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

```bash
srun -n <nproc> --ntasks-per-node=4 --gpus-per-task=1 --mpi=pmi2 \
  ./ParaFlow_lic_gpu conf/lic_streamline_gpu.yaml
```

`<nproc>` must match `nproc`/`nblocks` in the config. `--mpi=pmi2` is
required — without it Slurm won't form a single MPI job across ranks.

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
- **`streamline_storage`**: `memory` (default) keeps everything in RAM and
  writes nothing but the final PNG; `disk` also keeps the raw per-block
  trajectory files for later reuse.
- **`outputDir`**: where results land, in a subfolder named after the mode,
  direction, and dataset — so different datasets/directions don't overwrite
  each other, but two runs of the *same* dataset/direction will.
