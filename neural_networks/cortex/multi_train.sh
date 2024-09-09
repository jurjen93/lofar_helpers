#!/bin/bash
#SBATCH --job-name=cortex_multi_node
#SBATCH --partition=gpu              
#SBATCH --time=00:30:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH --cpus-per-task=18
#SBATCH --output=out/multi_cortex%A_%a.out
#SBATCH --gpu-bind=None

set -e

cd ~/projects/lofar_helpers/neural_networks

source venv/bin/activate
module load 2023

export MASTER_PORT=$(expr 10000 + $(echo -n $SLURM_JOBID | tail -c 4))
export MASTER_ADDR=$(scontrol show hostnames "$SLURM_JOB_NODELIST" | head -n 1)
echo "MASTER_ADDR:MASTER_PORT="${MASTER_ADDR}:${MASTER_PORT}

srun python train_nn_multi.py