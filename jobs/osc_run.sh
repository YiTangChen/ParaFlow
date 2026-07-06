#!/bin/bash
#SBATCH --job-name=TestTime
#SBATCH --ntasks=16
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=02:00:00
#SBATCH --account=PAS0027
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G

export LD_LIBRARY_PATH=/apps/spack/0.21/ascend/linux-rhel9-zen2/netcdf-cxx4/oneapi/2024.1.0/mvapich/3.0/4.3.1-rxfah3c/lib:$LD_LIBRARY_PATH

# Clear stale SLURM CPU env vars inherited from shell to prevent cpus-per-task conflicts
unset SLURM_CPUS_PER_TASK
unset SLURM_TRES_PER_TASK

# srun --export=ALL --distribution=block:cyclic ./DrawSubdomain conf/osc_drawsubdomain.yaml
srun --export=ALL --cpus-per-task=1 --distribution=block:cyclic ./ParaFlow_streamline conf/osc_ParaFlow_streamline.yaml

srun -n 32 --distribution=block:cyclic ./ParaFlow_streamline conf/osc_ParaFlow_streamline.yaml 1>logs/b${032}_s${100}_run${1}.stdout 2>logs/b${032}_s${100}_run${1}.stderr
# srun --export=ALL --ntasks=16 --ntasks-per-node=4 --gpus-per-node=4 --gpu-bind=closest ./DrawSubdomain config_drawsubdomain_normal.yaml

# srun --export=ALL --ntasks=16 --ntasks-per-node=4 ./ParaFlow_streamline conf/ParaFlow_streamline.yaml

# srun -n 256 -N 12 --ntasks-per-node=24 --distribution=block:cyclic ./DrawSubdomain conf/osc_drawsubdomain.yaml

# srun -n 256 -N 12 --ntasks-per-node=24 --distribution=block:cyclic ./ParaFlow_streamline conf/osc_ParaFlow_streamline.yaml

# srun -n 256 -N 12 --ntasks-per-node=24 --distribution=block:cyclic ./ParaFlow_pathline conf/osc_ParaFlow_pathline.yaml