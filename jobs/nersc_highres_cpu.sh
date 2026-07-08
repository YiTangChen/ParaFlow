#!/bin/bash

#SBATCH --nodes=16
#SBATCH --ntasks=256
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=16
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -J paraflow_nersc_highres_cpu
#SBATCH --mail-user=yang.5039@osu.edu
#SBATCH --mail-type=ALL
#SBATCH -A m4259
#SBATCH -t 24:0:0
#SBATCH -o /global/homes/y/yang5039/ParaFlow/logs/nersc_highres_cpu_%j.out
#SBATCH -e /global/homes/y/yang5039/ParaFlow/logs/nersc_highres_cpu_%j.err

# Layout: 16 CPU nodes, 256 MPI ranks total, 16 ranks/node, 16 CPUs/rank → 32 GB/rank
# Memory sweep: nblocks=256 point (pair with 16/32/64/128 for scaling analysis)

export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

cd /global/homes/y/yang5039/ParaFlow

time srun --ntasks=256 \
          --ntasks-per-node=16 \
          --cpus-per-task=16 \
          --cpu-bind=cores \
          ./ParaFlow_pathline conf/nersc_highres_cpu.yaml
