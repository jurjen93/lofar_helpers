#!/bin/bash
#SBATCH --job-name=cortex_multi_node
#SBATCH --partition=gpu
#SBATCH --time=00:30:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH --output=out/multi_cortex%A_%a.out
#SBATCH --gpu-bind=None

set -e

cd ~/projects/lofar_helpers/neural_networks


module load 2023
# module load NCCL/2.18.3-GCCcore-12.3.0-CUDA-12.1.1
# module load PyTorch/2.1.2-foss-2023a-CUDA-12.1.1
# module load libjpeg-turbo/2.1.5.1-GCCcore-12.3.0
# module load torchvision/0.16.0-foss-2023a-CUDA-12.1.1
source venv/bin/activate

export MASTER_PORT=$(expr 10000 + $(echo -n $SLURM_JOBID | tail -c 4))
export MASTER_ADDR=$(scontrol show hostnames "$SLURM_JOB_NODELIST" | head -n 1)
# export NCCL_SOCKET_IFNAME='eno2np0' # Change for a100
echo "MASTER_ADDR:MASTER_PORT="${MASTER_ADDR}:${MASTER_PORT}

NCCL_DEBUG=INFO srun python train_nn_multi.py