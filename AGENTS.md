# AGENTS

Root guide for the **ParaFlow** repository. Nested `AGENTS.md` files (e.g. `debug/AGENTS.md`)
take precedence within their own subtree.

## What this is

ParaFlow is a parallel particle-tracing / flow-visualization code for **MPAS-Ocean** data.
It generates **streamlines** and **pathlines** over an unstructured ocean mesh using:

- **MPI + DIY** for block-parallel domain decomposition,
- **OSUFlow** for the field-line integrators (RK4), with an **optional CUDA GPU backend**,
- **yaml-cpp** for run configuration.

It runs on HPC clusters (NERSC Perlmutter, OSC ascend/cardinal) under Slurm via `srun`.

## Repository layout

```
ParaFlow/                     ← repo root (this file)
├─ ParaFlow_streamline.cpp    ← thin driver: build ParaFlow → GenStreamLines()
├─ ParaFlow_pathline.cpp      ← thin driver: build ParaFlow → GenPathLines()
├─ DrawSubdomain.cpp          ← driver: generates the subdomain / ocean-mask data
├─ CMakeLists.txt             ← top-level build
├─ build.sh / clean.sh / switch.sh   ← build + machine-select (run from root)
│
├─ ParaFlow/                  ← CORE LIBRARY (libParaFlow.a)
│  ├─ ParaFlow.cpp / .hpp     ←   main tracing logic + config parsing
│  ├─ gpu_parity_tests.hpp    ←   CPU-vs-GPU parity/smoke tests (debug-only; textual
│  │                          ←   include of ParaFlow.cpp, dormant unless OSUFLOW_GPU_SMOKE set)
│  ├─ block.hpp               ←   DIY block definition
│  ├─ timing.hpp / utils.hpp  ←   instrumentation + helpers
│  ├─ OSUFlow/                ←   vendored flow library (CPU integrators + GPU/CUDA/)
│  ├─ diy/                    ←   git submodule (block-parallel framework)
│  └─ yaml-cpp/               ←   vendored config parser
│
├─ scripts/                   ← Python tooling (run from repo root)
│  ├─ setup/                  ←   gen_random_seeds, gen_partition, check_atlantic_seeds
│  ├─ plot/                   ←   read_traces, plot_*, drawsubdomain, convert_traces_to_vtk
│  └─ profiling/              ←   parse_timing, plot_timing, parse_mem
├─ jobs/                      ← Slurm launch scripts
│  ├─ nersc_{highres,lowres}_{cpu,gpu}.sh   ← NERSC pathline (256 ranks, sbatch)
│  ├─ nersc_drawsubdomain.sh  ←   NERSC subdomain job
│  ├─ osc_run.sh              ←   OSC (ascend/cardinal) launcher
│  └─ plot_traces.sh          ←   reassemble + plot helper
├─ conf/                      ← YAML run configs (see naming below)
│  └─ debug/                  ←   debug_* + gpu_parity_* configs
├─ NERSC_job/                 ← archive of the full NERSC rank sweep (16–256); gitignored
├─ debug/                     ← CPU/GPU parity debugging workspace (own AGENTS.md; gitignored)
└─ rendering/ , LIC/          ← visualization / line-integral-convolution
```

> Python scripts and job scripts use **CWD-relative** data paths — always invoke them
> from the repo root (e.g. `python3 scripts/plot/plot_traces.py ...`, `jobs/plot_traces.sh`),
> not from inside their folder.

**Generated, gitignored (do not commit):**
`streamlines/`, `pathlines/`, `png/`, `output/`, `block_outputs/`, `seeds/`,
`logs/`, `results/`, `partition/`, `drawSubdomain/`, all build artifacts
(`CMakeFiles/`, `CMakeCache.txt`, `*.a`, top-level executables), and `debug/`.

## Build

```bash
./switch.sh {nersc|ascend|cardinal}   # select machine (edits build.sh + CMakeLists)
./clean.sh                            # clear stale CMake files
./build.sh                            # CPU build
OSUFLOW_ENABLE_CUDA=1 ./build.sh      # GPU build (CUDA_ARCH=90 to override, default 80)
```

Build order: `OSUFlow` → `yaml-cpp` → `ParaFlow/` (libParaFlow.a) → top-level executables.

## Run (Slurm)

```bash
python3 scripts/setup/gen_random_seeds.py [N]   # default 10000; also creates output dirs

# NERSC: submit a batch job (each script cd's to repo root, reads its conf/ file)
sbatch jobs/nersc_highres_cpu.sh      # or _gpu / nersc_lowres_{cpu,gpu} / nersc_drawsubdomain
# OSC (ascend/cardinal): interactive / batch launcher
bash jobs/osc_run.sh

# Or run a binary directly (executables live at repo root):
srun -n 256 -N 12 --ntasks-per-node=24 --distribution=block:cyclic \
     ./ParaFlow_pathline conf/<cfg>.yaml     # or ./ParaFlow_streamline / ./DrawSubdomain
```

**Config naming (`conf/`):**
- **NERSC** (canonical, 256 ranks): `nersc_{highres,lowres}_{cpu,gpu}.yaml`
  (pathline) + `nersc_{highres,lowres}_drawsubdomain.yaml`.
  highres = `18to6v3` 3.7M-cell mesh; lowres = `LowRes` climatology mesh; `gpu` sets `useGPU: true`.
- **OSC default:** `osc_ParaFlow_{pathline,streamline}[_gpu].yaml`, `osc_drawsubdomain.yaml`.
- **`conf/debug/`:** `debug_*` (1-rank / small-block) and `gpu_parity_*` configs.
- Each NERSC job script pairs with its like-named config (`jobs/nersc_highres_cpu.sh` →
  `conf/nersc_highres_cpu.yaml`). The full 16–256 rank sweep is archived in `NERSC_job/`.

## Analysis / plotting (Python — run from repo root)

- **`scripts/setup/`:** `gen_random_seeds.py`, `gen_partition.py`, `check_atlantic_seeds.py`
- **`scripts/plot/`:** `read_traces.py`, `convert_traces_to_vtk.py`, `plot_traces.py`,
  `plot_streamlines.py`, `plot_pathlines.py`, `plot_traces_latlon.py`, `drawsubdomain.py`
- **`scripts/profiling/`:** `parse_timing.py`, `plot_timing.py`, `parse_mem.py`

These scripts have no cross-imports (references between them are docstring notes only),
so they can be moved/renamed freely as long as callers are updated.

## Work rules

1. **Read related code and files first**, then begin revisions.
2. **State the expected effect** of a change before making it.
3. Keep changes **small and tightly scoped** to the current task; do not touch unrelated files.
4. Use **English** for all code, commands, and file names.
5. **Do not modify vendored code** (`ParaFlow/diy`, `ParaFlow/OSUFlow`, `ParaFlow/yaml-cpp`)
   unless the task is explicitly about them.
6. **Never commit** build artifacts or generated outputs (respect `.gitignore`).
7. `debug/` is the parity/debug sandbox and has its own `AGENTS.md` — defer to it there.

## Known context

- Three RK4 integrators exist (CPU full-scan reference, CPU walk, GPU). **CPU/GPU parity**
  is a live concern; the GPU 1-ring cell walk was a past divergence source, fixed with an
  iterative march. See `debug/GPU_CPU_PARITY_DEBUG_STORY.md` and related debug docs.
