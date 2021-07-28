#!/bin/bash
#SBATCH -c 1
#SBATCH --mail-type=END,FAIL
#SBATCH --wait

#THIS SCRIPT IS WRITTEN FOR SLURM ON SURFSARA
echo "----------START----------"

SOURCE=$1 #L626678
TO=/project/lofarvwf/Share/jdejong/output/${SOURCE}
SCRIPT_FOLDER=/home/lofarvwf-jdejong/scripts/lofar_helpers

SING_IMAGE=/home/lofarvwf-jdejong/singularities/pill-latest.simg
SING_BIND=/project/lofarvwf/Share/jdejong

#CREATE BOXES
echo "Create boxes..."
#singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/make_boxes.py -f ${TO}/extract/image_full_ampphase_di_m.NS.app.restored.fits -l ${TO}
TOTAL_BOXES=20
echo "Succesfully created boxes..."

mkdir ${TO}/test

#EXTRACT
echo "EXTRACT STARTED"
sbatch ${SCRIPT_FOLDER}/pipeline_scripts/surf/extract_test.sh L626678 &
wait &

#SELFCAL
echo "SELFCAL STARTED"
for ((N=1;N<=${TOTAL_BOXES};N++))
do
  echo "SELFCAL ${N}"
  until [[ -f ${TO}/test/test_${N}.txt ]]
  do
    sleep 5
    echo "Still waiting for ${N}"
  done
  echo "selfcal_${N}" > ${TO}/test/test_selfcal_${N}.txt
done
wait

echo "----------END----------"