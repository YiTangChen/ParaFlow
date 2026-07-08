# ParaFlow

ParaFlow is a distributed particle-tracing engine for **MPAS-Ocean** model output.
It integrates **streamlines** and **pathlines** through the time-varying velocity
field of an unstructured global-ocean mesh, scaling across nodes with MPI.

**Architecture**
- **MPI + [DIY](https://github.com/diatomic/diy)** — block-parallel decomposition, one mesh partition per rank.
- **OSUFlow** — RK4 field-line integration, with an optional **CUDA** backend.

The CUDA tracing kernels are factored into a standalone submodule,
**[mpaso-gpu-kernels](https://github.com/Jiaxin-yyjx/mpaso-gpu-kernels)** — a
framework-agnostic (plain-array) core shared with the MOPS project. ParaFlow adds a thin host uploader (`MPASOGPUTracer`) that flattens its OSUFlow mesh into that core.

## Execution modes

ParaFlow support two build modes. The CUDA build additionally selects CPU vs GPU per run via `useGPU`.

| Mode | Build | `useGPU` | Tracing |
|------|-------|----------|---------|
| **CPU** | `./build.sh` | `false` | RK4 on CPU cores |
| **CPU + GPU** | `OSUFLOW_ENABLE_CUDA=1 ./build.sh` | `true` | RK4 on GPU (auto-fallback to CPU if no device) |

Runs use **one block per rank**: keep `nproc == nblocks ==` the partition's `part.N`,
and launch with `srun -n N`.

## Requirements

MPI (Cray MPICH / MVAPICH) · CMake ≥ 3.18 · C++17 · NetCDF-C / NetCDF-C++ ·
CUDA ≥ 11 (GPU mode only).

## Quick start

ParaFlow has **two submodules** (`diy`, `mpaso-gpu-kernels`), so clone recursively:

    git clone --recurse-submodules git@github.com:YiTangChen/ParaFlow.git
    cd ParaFlow
    # cloned flat already?  
    ->  git submodule update --init --recursive

The `mpaso-gpu-kernels` submodule is compiled only in the CUDA build.

<!-- ## Running

    # 1 · seeds + output folders  (default 10000)
    python3 scripts/setup/gen_random_seeds.py [N]

    # 2 · subdomain / ocean mask  (once; used for masking + plots)
    srun -n N ./DrawSubdomain conf/<machine>_<res>_drawsubdomain.yaml

    # 3 · trace  (N must equal the config's nproc)
    srun -n N ./ParaFlow_streamline conf/<cfg>.yaml
    srun -n N ./ParaFlow_pathline   conf/<cfg>.yaml

    # 4 · reassemble ranks + render on the ocean mask
    jobs/plot_traces.sh -->


## Main setup flow

Run the following steps in order.

### 1. Generate seeds and output folders

Run this first.

```bash
python3 scripts/setup/gen_random_seeds.py
```

The script accepts an optional seed count:

```bash
python3 scripts/setup/gen_random_seeds.py 100
python3 scripts/setup/gen_random_seeds.py 10000
```

If no argument is given, it generates `10000` seeds by default.

This creates `seeds/random_seeds.bin` and the required output folders such as `seeds/`, `drawSubdomain/`, `png/`, `streamlines/`, and `pathlines/`.

### 2. Select the machine configuration

Run this once at the beginning, or anytime you switch machines.

```bash
./switch.sh {nersc|ascend|cardinal}
```

This updates:

- `build.sh`
- `CMakeLists.txt`
- `ParaFlow/CMakeLists.txt`
- `ParaFlow/OSUFlow/CMakeLists.txt`

### 3. Clean old CMake files

```bash
./clean.sh
```

### 4. Build everything

```bash
./build.sh
```

For GPU compile:
```bash
OSUFLOW_ENABLE_CUDA=1 ./build.sh
```

### 5. Generate the subdomain output

```bash
srun -n N ./DrawSubdomain conf/<machine>_<res>_drawsubdomain.yaml
```

### 6. Plot the subdomain

```bash
python3 scripts/plot/drawsubdomain.py
```

After steps 1-6 are done, the project is ready for the trace runs below.

## Run a trace

Once setup (steps 1–6) is done, trace **streamlines** or **pathlines** — they are
independent, so run either or both. Match `srun -n N` to the config's `nproc`
(= `nblocks` = the partition's `part.N`):

​```bash
srun -n N ./ParaFlow_streamline conf/<cfg>.yaml   # streamlines
srun -n N ./ParaFlow_pathline   conf/<cfg>.yaml   # pathlines
​```

On NERSC, submit the batch scripts instead — e.g. `sbatch jobs/nersc_lowres_gpu.sh`
(each pairs with its like-named `conf/` file).

## Plot

​```bash
jobs/plot_traces.sh   # reassemble per-rank output + render it on the ocean mask
​```

`plot_traces.sh` plots **pathlines** by default; to plot **streamlines**, swap the two
commented lines at the top of the script (point it at `streamlines/` instead of `pathlines/`).
