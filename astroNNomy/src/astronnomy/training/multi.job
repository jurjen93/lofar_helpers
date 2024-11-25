#!/bin/bash
#SBATCH --job-name=pyg-multinode-tutorial # identifier for the job listings
#SBATCH --output=pyg-multinode.log        # outputfile
#SBATCH --partition=gpu              # ADJUST this to your system
# SBATCH -t 00:10:00
#SBATCH --time=10:00
#SBATCH --nodes=2 
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4 
#SBATCH --cpus-per-task=18



cd ~/lofar_helpers/neural_networks

source venv/bin/activate

# module load 2023
# module load PyTorch/2.1.2-foss-2023a-CUDA-12.1.1
# module load torchvision/0.16.0-foss-2023a-CUDA-12.1.1

export MASTER_PORT=$(expr 10000 + $(echo -n $SLURM_JOBID | tail -c 4))
export MASTER_ADDR=$(scontrol show hostnames "$SLURM_JOB_NODELIST" | head -n 1)
echo "MASTER_ADDR:MASTER_PORT="${MASTER_ADDR}:${MASTER_PORT}


srun python test2.py