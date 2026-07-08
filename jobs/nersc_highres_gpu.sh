#!/bin/bash

#SBATCH --nodes=64
#SBATCH --ntasks=256
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=16
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -J paraflow_nersc_highres_gpu
#SBATCH --mail-user=yang.5039@osu.edu
#SBATCH --mail-type=ALL
#SBATCH -A m4259
#SBATCH -t 24:0:0
#SBATCH -o /global/homes/y/yang5039/ParaFlow/logs/nersc_highres_gpu_%j.out
#SBATCH -e /global/homes/y/yang5039/ParaFlow/logs/nersc_highres_gpu_%j.err

# Layout: 64 GPU nodes * 4 A100 GPUs/node = 256 GPUs.
# Run 256 MPI ranks total: 4 ranks/node, 1 A100/rank, 16 CPUs/rank.
# HighRes mesh: 3.7M cells, 80 levels.
# Pair with nersc_highres_cpu.sh for CPU vs GPU comparison.

export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
# Binary not linked with libmpi_gtl_cuda; disable GPU-aware MPI to avoid GTL abort
export MPICH_GPU_SUPPORT_ENABLED=0

cd /global/homes/y/yang5039/ParaFlow

time srun --ntasks=256 \
          --ntasks-per-node=4 \
          --cpus-per-task=16 \
          --cpu-bind=cores \
          --gpus-per-task=1 \
          --gpu-bind=closest \
          ./ParaFlow_pathline conf/nersc_highres_gpu.yaml
