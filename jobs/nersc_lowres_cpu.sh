#!/bin/bash

#SBATCH -N 16
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -J paraflow_nersc_lowres_cpu
#SBATCH --mail-user=yang.5039@osu.edu
#SBATCH --mail-type=ALL
#SBATCH -A m4259
#SBATCH -t 24:0:0
#SBATCH -o /global/homes/y/yang5039/ParaFlow/logs/nersc_lowres_cpu_%j.out
#SBATCH -e /global/homes/y/yang5039/ParaFlow/logs/nersc_lowres_cpu_%j.err

# Layout: 16 CPU nodes, 16 ranks/node, 16 CPUs/rank (fills node with hyperthreading)
# Cells/block: ~925   Seeds/block: ~39 (10K seeds)
# Pair with nersc_lowres_gpu.sh for CPU vs GPU comparison at same rank count
# At ~39 seeds/block, GPU overhead likely dominates — expect CPU to win here

export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

cd /global/homes/y/yang5039/ParaFlow

time srun -n 256 -c 16 --cpu_bind=cores ./ParaFlow_pathline conf/nersc_lowres_cpu.yaml
