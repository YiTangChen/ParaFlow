# Parallel Workflow

## Main setup flow

Run the following steps in order.

### 1. Generate seeds and output folders

Run this first.

```bash
python3 gen_random_seeds.py
```

The script accepts an optional seed count:

```bash
python3 gen_random_seeds.py 100
python3 gen_random_seeds.py 10000
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

### 5. Generate the subdomain output

```bash
srun -n 256 -N 12 --ntasks-per-node=24 --distribution=block:cyclic ./DrawSubdomain conf/osc_drawsubdomain.yaml
```

### 6. Plot the subdomain

```bash
python3 drawsubdomain.py
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
./plot_traces.sh
```

## Pathline workflow

### 9. Run pathline tracing

```bash
srun -n 256 -N 12 --ntasks-per-node=24 --distribution=block:cyclic ./ParaFlow_pathline conf/osc_ParaFlow_pathline.yaml
```

### 10. Plot pathline results

```bash
./plot_traces.sh
```

## Note

`plot_traces.sh` currently uses the pathline plotting command by default.

If you want to plot streamlines instead, edit `plot_traces.sh`, uncomment the first two lines below, and comment out the last line:

```bash
# python3 read_traces.py 256 streamlines
# python3 plot_traces.py streamlines/streamlines.bin drawSubdomain streamlines_map.png

# python3 read_traces.py 256 pathlines
python3 plot_traces.py pathlines/pathlines.bin drawSubdomain pathlines_map.png
```

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
