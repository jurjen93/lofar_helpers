#!/bin/bash
#SBATCH -c 1
#SBATCH --mail-type=END,FAIL

#THIS SCRIPT IS WRITTEN FOR SLURM ON SURFSARA
echo "----------START----------"

SOURCE=$1 #L626678
TO=/project/lofarvwf/Share/jdejong/output/${SOURCE}
SCRIPT_FOLDER=/home/lofarvwf-jdejong/scripts/lofar_helpers

SING_IMAGE=/home/lofarvwf-jdejong/singularities/pill-latest.simg
SING_BIND=/project/lofarvwf/Share/jdejong

#MOVE NEEDED FILES
#echo "Moving files to ${TO}/extract and untar..."
#cp -r /project/lofarvwf/Share/jdejong/data/${SOURCE}/data_archive.tar.gz ${TO}/extract
#rm -r /project/lofarvwf/Share/jdejong/data/${SOURCE}/data_archive.tar.gz
#echo "Succesfully finished moving files..."

#cd ${TO}/extract
#tar -zxvf data_archive.tar.gz
#echo"Untarred succesfuly..."

#CREATE BOXES
echo "Create boxes..."
#singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/make_boxes.py -f ${TO}/extract/image_full_ampphase_di_m.NS.app.restored.fits -l ${TO}
TOTAL_BOXES=20
echo "Succesfully created boxes..."

mkdir ${TO}/test

#EXTRACT
echo "EXTRACT STARTED"
sbatch ${SCRIPT_FOLDER}/pipeline_scripts/surf/extract_test.sh L626678 &

#SELFCAL
echo "SELFCAL STARTED"
for ((N=1;N<=${TOTAL_BOXES};N++))
do
  echo "SELFCAL ${N}"
  while [[ ! -f ${TO}/test/text_${N}.txt ]]
  do
    sleep 5
    echo "Still waiting for ${N}"
  done &
  echo "selfcal_${N}" > ${TO}/test/test_selfcal_${N}.txt
  exit
done
wait

echo "----------END----------"