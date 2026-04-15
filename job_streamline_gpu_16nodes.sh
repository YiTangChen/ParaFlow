#!/bin/bash
#SBATCH --job-name=StreamlineGPU16
#SBATCH --account=PAS0027
#SBATCH --nodes=16
#SBATCH --ntasks=16
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=1
#SBATCH --mem=256G
#SBATCH --time=04:00:00
#SBATCH --output=osc_streamline_gpu.sbatch.%j.out
#SBATCH --error=osc_streamline_gpu.sbatch.%j.err

set -euo pipefail

cd /users/PAS2171/chen10522/OSUFlow/Parallel

export LD_LIBRARY_PATH=/apps/spack/0.21/ascend/linux-rhel9-zen2/netcdf-cxx4/oneapi/2024.1.0/mvapich/3.0/4.3.1-rxfah3c/lib:${LD_LIBRARY_PATH:-}

unset OSUFLOW_GPU_SMOKE
export OSUFLOW_GPU_STREAMLINE_CHUNK_STEPS=${OSUFLOW_GPU_STREAMLINE_CHUNK_STEPS:-43200}

export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

unset SLURM_TRES_PER_TASK

srun --export=ALL \
  -N 16 \
  -n 16 \
  --ntasks-per-node=1 \
  --cpus-per-task=8 \
  --gpus-per-node=1 \
  --distribution=block:cyclic \
  ./ParaFlow_streamline conf/osc_ParaFlow_streamline_gpu.yaml \
  2>&1 | tee osc_streamline_gpu.log
