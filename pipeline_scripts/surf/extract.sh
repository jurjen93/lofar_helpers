#!/bin/bash
#SBATCH -c 32
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jurjendejong@strw.leidenuniv.nl
#SBATCH --array=1-200%10

FIELD=$1 #L626678
TO=/project/lofarvwf/Share/jdejong/output/${FIELD}
SCRIPT_FOLDER=/home/lofarvwf-jdejong/scripts

SING_IMAGE=/project/lofarvwf/Software/lofar_sksp_fedora27_ddf_slurmfix.sif
SING_BIND=/project/lofarvwf/Share/jdejong

#GET LAST BOX NUMBER
TOTAL_BOXES=$(ls -dq ${TO}/boxes/box*.reg | wc -l)

#START EXTRACT
echo "-----STARTED EXTRACT-----"
cd ${TO}/extract || { echo "Missing path"; exit 1; }
if [[ ! ${SLURM_ARRAY_TASK_ID} > ${TOTAL_BOXES} ]]
then
  singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/sub-sources-outside-region.py -b ${TO}/boxes/box_${SLURM_ARRAY_TASK_ID}.reg --overwriteoutput -p box_${SLURM_ARRAY_TASK_ID}
  echo "Extracted box_${BOX}"
else
  :
fi
echo "-----FINISHED EXTRACT-----"