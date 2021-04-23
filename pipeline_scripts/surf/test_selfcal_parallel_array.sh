#!/bin/bash
#SBATCH -c 2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jurjendejong@strw.leidenuniv.nl
#SBATCH --job-name=array
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err
#SBATCH --array=1-5%2

SOURCE=$1 #L626678
TO=/project/lofarvwf/Share/jdejong/output/${SOURCE}

SCRIPT_FOLDER=/home/lofarvwf-jdejong/scripts/lofar_helpers/pipeline_scripts/surf

echo "Task number ${SLURM_ARRAY_TASK_ID}"
echo $SLURM_ARRAY_TASK_ID
sleep 60
srun /home/lofarvwf-jdejong/scripts/monitor.sh