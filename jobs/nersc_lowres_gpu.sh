#!/bin/bash

#SBATCH -N 64
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -J paraflow_nersc_lowres_gpu
#SBATCH --mail-user=yang.5039@osu.edu
#SBATCH --mail-type=ALL
#SBATCH -A m4259
#SBATCH -t 24:0:0
#SBATCH -o /global/homes/y/yang5039/ParaFlow/logs/nersc_lowres_gpu_%j.out
#SBATCH -e /global/homes/y/yang5039/ParaFlow/logs/nersc_lowres_gpu_%j.err

# Layout: 64 GPU nodes, 4 ranks/node, 1 A100/rank, 16 CPUs/rank (256 ranks total)
# Cells/block: ~925   Seeds/block: ~39 (10K seeds)
# Pair with nersc_lowres_cpu.sh for CPU vs GPU comparison at same rank count
# NOTE: 64 GPU nodes is large — run lowres_16_gpu.sh first to validate

export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
# Binary not linked with libmpi_gtl_cuda; disable GPU-aware MPI to avoid GTL abort
export MPICH_GPU_SUPPORT_ENABLED=0

cd /global/homes/y/yang5039/ParaFlow

time srun -n 256 -c 16 --gpus-per-task=1 --gpu-bind=map_gpu:0,1,2,3 ./ParaFlow_pathline conf/nersc_lowres_gpu.yaml
