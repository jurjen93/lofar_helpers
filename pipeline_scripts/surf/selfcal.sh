#!/bin/bash
#SBATCH -c 20
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jurjendejong@strw.leidenuniv.nl
#SBATCH --array=1-6

SOURCE=$1 #L626678
TO=/project/lofarvwf/Share/jdejong/output/${SOURCE}/selfcal
N=$2 #box number
echo "${TO}/box_${N}.${SLURM_ARRAY_TASK_ID}"
until [[ -d "${TO}/box_${N}.${SLURM_ARRAY_TASK_ID}" ]]
do
  sleep 5
  echo "NOT EXISTING"
done
sleep 5
srun ${TO}/box_${N}.${SLURM_ARRAY_TASK_ID}/command.sh