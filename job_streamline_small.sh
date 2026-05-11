#!/bin/bash

# OSC short streamline validation jobs for the asynchronous timing design.
# Runs 1-seed and 100-seed cases with maxsteps*dt = 100*20s = 2000s.

#SBATCH -N 1
#SBATCH --ntasks=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --account=PAS0027
#SBATCH --mem=64G
#SBATCH -J StreamlineSmall

set -euo pipefail

export LD_LIBRARY_PATH=/apps/spack/0.21/ascend/linux-rhel9-zen2/netcdf-cxx4/oneapi/2024.1.0/mvapich/3.0/4.3.1-rxfah3c/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

unset SLURM_CPUS_PER_TASK
unset SLURM_TRES_PER_TASK

mkdir -p logs/small_async results/small_async

srun --export=ALL -n 2 --cpus-per-task=1 --distribution=block:cyclic ./ParaFlow_streamline conf/small_streamline_s1.yaml \
  > logs/small_async/b002_s0000001_run1.stdout \
  2> logs/small_async/b002_s0000001_run1.stderr

srun --export=ALL -n 2 --cpus-per-task=1 --distribution=block:cyclic ./ParaFlow_streamline conf/small_streamline_s100.yaml \
  > logs/small_async/b002_s0000100_run1.stdout \
  2> logs/small_async/b002_s0000100_run1.stderr

python3 parse_timing.py --results-dir=results/small_async --reset logs/small_async/*.stderr
