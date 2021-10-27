#!/bin/bash
#SBATCH -c 32
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jurjendejong@strw.leidenuniv.nl
#SBATCH --array=1-82
#SBATCH --job-name=extract

SOURCE=$1
TO=/project/lofarvwf/Share/jdejong/output/${SOURCE}
SCRIPT_FOLDER=/home/lofarvwf-jdejong/scripts

SING_IMAGE=/project/lofarvwf/Software/lofar_sksp_fedora27_ddf_slurmfix.sif
SING_BIND=/project/lofarvwf/Share/jdejong

#GET LAST BOX NUMBER
TOTAL_BOXES=$(ls -dq ${TO}/boxes/box*.reg | wc -l)

#START EXTRACT
echo "-----STARTED EXTRACT-----"
if [[ ! ${SLURM_ARRAY_TASK_ID} -gt ${TOTAL_BOXES} ]]
then
  mkdir ${TO}/extract/box_${SLURM_ARRAY_TASK_ID}
  cp ${TO}/extract/data_archive.tar.gz ${TO}/extract/box_${SLURM_ARRAY_TASK_ID}/
  cd ${TO}/extract/box_${SLURM_ARRAY_TASK_ID} || { echo "Missing path"; exit 1; }
  tar -xvf data_archive.tar.gz
  rm -rf data_archive.tar.gz
  singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/sub-sources-outside-region.py -b ${TO}/boxes/box_${SLURM_ARRAY_TASK_ID}.reg --overwriteoutput -p box_${SLURM_ARRAY_TASK_ID}
  echo "Extracted box_${SLURM_ARRAY_TASK_ID}"
  echo "Selfcal box_${SLURM_ARRAY_TASK_ID} finished" > ${TO}/finished/box_${SLURM_ARRAY_TASK_ID}.txt
else
  :
fi

echo "#####CLEANING UP#####"
source ${SCRIPT_FOLDER}/lofar_helpers/cleaning/cleaning_extract_per_box.sh ${SLURM_ARRAY_TASK_ID}

echo "-----FINISHED EXTRACT-----"