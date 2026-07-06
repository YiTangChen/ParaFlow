# ParaFlow

ParaFlow is a parallel particle-tracing / flow-visualization code for **MPAS-Ocean** data.
It generates **streamlines** and **pathlines** over an unstructured ocean mesh using:

- **MPI + DIY** for block-parallel domain decomposition,
- **OSUFlow** for the field-line integrators (RK4), with an **optional CUDA GPU backend**,
- **yaml-cpp** for run configuration.

It runs on HPC clusters (NERSC Perlmutter, OSC ascend/cardinal) under Slurm via `srun`.

## Repository layout

```
ParaFlow/                     ← repo root
├─ ParaFlow_streamline.cpp    ← thin driver → GenStreamLines()
├─ ParaFlow_pathline.cpp      ← thin driver → GenPathLines()
├─ DrawSubdomain.cpp          ← driver: generates the subdomain / ocean-mask data
├─ CMakeLists.txt             ← top-level build
├─ build.sh / clean.sh / switch.sh   ← build + machine-select (run from root)
│
├─ ParaFlow/                  ← core library (libParaFlow.a)
│  ├─ ParaFlow.cpp / .hpp     ←   main tracing logic + config parsing
│  ├─ block.hpp               ←   DIY block definition
│  ├─ timing.hpp / utils.hpp  ←   instrumentation + helpers
│  ├─ OSUFlow/                ←   vendored flow library (CPU + GPU/CUDA integrators)
│  ├─ diy/                    ←   git submodule (block-parallel framework)
│  └─ yaml-cpp/               ←   vendored config parser
│
├─ scripts/                   ← Python tooling (run from repo root)
│  ├─ setup/                  ←   gen_random_seeds, gen_partition, check_atlantic_seeds
│  ├─ plot/                   ←   read_traces, plot_*, drawsubdomain, convert_traces_to_vtk
│  └─ profiling/              ←   parse_timing, plot_timing, parse_mem
├─ jobs/                      ← Slurm launch scripts (nersc_*, osc_run.sh, plot_traces.sh)
├─ conf/                      ← YAML run configs (see naming below)
└─ rendering/ , LIC/          ← visualization / line-integral-convolution
```

> Python and job scripts use **CWD-relative** paths — always run them from the repo
> root (e.g. `python3 scripts/plot/plot_traces.py ...`), not from inside their folder.

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
srun -n 256 -N 12 --ntasks-per-node=24 --distribution=block:cyclic ./DrawSubdomain conf/osc_drawsubdomain.yaml
```

### 6. Plot the subdomain

```bash
python3 scripts/plot/drawsubdomain.py
```

After steps 1-6 are done, the project is ready for the trace runs below.

## Optional runs

The following parts do not need to follow a strict global order.

- `7` and `8` belong together for the streamline workflow.
- `9` and `10` belong together for the pathline workflow.
- You can run either workflow depending on what you want to generate.

## Streamline workflow

### 7. Run streamline tracing

```bash
srun -n 256 -N 12 --ntasks-per-node=24 --distribution=block:cyclic ./ParaFlow_streamline conf/osc_ParaFlow_streamline.yaml
```

### 8. Plot streamline results

```bash
jobs/plot_traces.sh
```

## Pathline workflow

### 9. Run pathline tracing

```bash
srun -n 256 -N 12 --ntasks-per-node=24 --distribution=block:cyclic ./ParaFlow_pathline conf/osc_ParaFlow_pathline.yaml
```

### 10. Plot pathline results

```bash
jobs/plot_traces.sh
```

## Note

`jobs/plot_traces.sh` currently uses the pathline plotting command by default.

If you want to plot streamlines instead, edit `jobs/plot_traces.sh`, uncomment the first two lines below, and comment out the last line:

```bash
# python3 scripts/plot/read_traces.py 256 streamlines
# python3 scripts/plot/plot_traces.py streamlines/streamlines.bin drawSubdomain streamlines_map.png

# python3 scripts/plot/read_traces.py 256 pathlines
python3 scripts/plot/plot_traces.py pathlines/pathlines.bin drawSubdomain pathlines_map.png
```

## Config naming (`conf/`)

- **NERSC** (canonical, 256 ranks): `nersc_{highres,lowres}_{cpu,gpu}.yaml` (pathline) +
  `nersc_{highres,lowres}_drawsubdomain.yaml`. highres = `18to6v3` 3.7M-cell mesh;
  lowres = `LowRes` climatology mesh; `gpu` sets `useGPU: true`.
- **OSC:** `osc_ParaFlow_{pathline,streamline}[_gpu].yaml`, `osc_drawsubdomain.yaml`.
- Each NERSC job script pairs with its like-named config
  (`jobs/nersc_highres_cpu.sh` → `conf/nersc_highres_cpu.yaml`).

## Clone with submodules

This repository uses `ParaFlow/diy` as a git submodule.

When cloning the repository for the first time, use:

```bash
git clone --recurse-submodules git@github.com:YiTangChen/ParaFlow.git
```

If you already cloned the repository without submodules, run:

```bash
git submodule update --init --recursive
```
