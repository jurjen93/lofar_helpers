#!/bin/bash
#SBATCH -c 20
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jurjendejong@strw.leidenuniv.nl
#SBATCH --array=1-6

SOURCE=$1 #L626678
TO=/project/lofarvwf/Share/jdejong/output/${SOURCE}/selfcal/box_${N}
N=$2 #box number

until [[ -f "${TO}/selfcal/box_${N}.${SLURM_ARRAY_TASK_ID}/command.txt" ]]
do
  sleep 60
done
COMMAND="cat ${TO}/selfcal/box_${N}.${SLURM_ARRAY_TASK_ID}/command.txt"
srun ${COMMAND}