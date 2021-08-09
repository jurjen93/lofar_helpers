#!/bin/bash
#SBATCH -c 32
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jurjendejong@strw.leidenuniv.nl

FIELD=$1
BOX=$2
TO=/project/lofarvwf/Share/jdejong/output/${FIELD}
SCRIPT_FOLDER=/home/lofarvwf-jdejong/scripts

SING_IMAGE=/project/lofarvwf/Software/lofar_sksp_fedora27_ddf_slurmfix.sif
SING_BIND=/project/lofarvwf/Share/jdejong

#GET LAST BOX NUMBER
TOTAL_BOXES=$(ls -dq ${TO}/boxes/box*.reg | wc -l)

#START EXTRACT
echo "-----STARTED EXTRACT-----"
if [[ ! ${BOX} > ${TOTAL_BOXES} ]]
then
  mkdir ${TO}/extract/box_${BOX}
  cp ${TO}/extract/data_archive.tar.gz ${TO}/extract/box_${BOX}/
  cd ${TO}/extract/box_${BOX} || { echo "Missing path"; exit 1; }
  tar -xvf data_archive.tar.gz
  START="$(date -u +%s)"
  singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/sub-sources-outside-region.py -b ${TO}/boxes/box_${BOX}.reg --overwriteoutput -p box_${BOX}
  END="$(date -u +%s)"
  echo "Extracted box_${BOX}"
  echo "Extracted in $((${END}-${START})) seconds" > ${TO}/extract/finished/box_${BOX}.txt
else
  :
fi
echo "-----FINISHED EXTRACT-----"