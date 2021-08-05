#!/bin/bash
#SBATCH -c 32
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jurjendejong@strw.leidenuniv.nl
#SBATCH --array=2-11


FIELD=$1 #L626678
TO=/project/lofarvwf/Share/jdejong/output/${FIELD}
SCRIPT_FOLDER=/home/lofarvwf-jdejong/scripts

SING_IMAGE=/project/lofarvwf/Software/lofar_sksp_fedora27_ddf_slurmfix.sif
SING_BIND=/project/lofarvwf/Share/jdejong

#GET LAST BOX NUMBER
TOTAL_BOXES=$(ls -dq ${TO}/boxes/box*.reg | wc -l)

#START EXTRACT
echo "-----STARTED EXTRACT-----"
if [[ ! ${SLURM_ARRAY_TASK_ID} > ${TOTAL_BOXES} ]]
then
  mkdir ${TO}/extract/box_${SLURM_ARRAY_TASK_ID}
  cp ${TO}/extract/data_archive.tar.gz ${TO}/extract/box_${SLURM_ARRAY_TASK_ID}/
  cd ${TO}/extract/box_${SLURM_ARRAY_TASK_ID} || { echo "Missing path"; exit 1; }
  tar -xvf data_archive.tar.gz
  START="$(date -u +%s)"
  singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/sub-sources-outside-region.py -b ${TO}/boxes/box_${SLURM_ARRAY_TASK_ID}.reg --overwriteoutput -p box_${SLURM_ARRAY_TASK_ID}
  END="$(date -u +%s)"
  echo "Extracted box_${BOX}"
  echo "Extracted in $((${END}-${START})) seconds" > ${TO}/extract/finished/box_${SLURM_ARRAY_TASK_ID}.txt
else
  :
fi
echo "-----FINISHED EXTRACT-----"