#!/bin/bash
#SBATCH -c 1
#SBATCH --mail-type=END,FAIL

echo "----------START RECALIBRATING----------"

FIELD=Abell399-401
TO=/project/lofarvwf/Share/jdejong/output/${FIELD}
SCRIPT_FOLDER=/home/lofarvwf-jdejong/scripts/lofar_helpers

SING_IMAGE=/home/lofarvwf-jdejong/singularities/pill-latest.simg
SING_BIND=/project/lofarvwf/Share/jdejong

#MOVE NEEDED FILES
#echo "Moving files to ${TO}/extract and untar..."
#cp -r /project/lofarvwf/Share/jdejong/data/${FIELD}/data_archive.tar.gz ${TO}/extract
#rm -r /project/lofarvwf/Share/jdejong/data/${FIELD}/data_archive.tar.gz
#echo "Succesfully finished moving files..."

#cd ${TO}/extract
#tar -zxvf data_archive.tar.gz
#echo"Untarred succesfuly..."

#CREATE FILES
mkdir ${TO}

#CREATE BOXES
echo "Create boxes..."
singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/make_boxes.py -f ${TO}/extract/image_full_ampphase_di_m.NS.app.restored.fits -l ${TO}
TOTAL_BOXES=$(ls -dq ${TO}/boxes/box*.reg | wc -l)
echo "Succesfully created boxes..."

#EXTRACT WITH PARALLEL ARRAY
echo "There are ${TOTAL_BOXES} boxes to extract"
mkdir ${TO}/extract
sbatch ${SCRIPT_FOLDER}/pipeline_scripts/surf/extract.sh ${FIELD} &
wait &

#SELFCAL
mkdir ${TO}/selfcal
for ((N=1;N<=${TOTAL_BOXES};N++))
do
  until [[ -f ${TO}/extract/${FIELD}_box_${N}.dysco.sub.shift.avg.weights.ms.archive0 ]]
  do
    sleep 60
  done
  sbatch ${SCRIPT_FOLDER}/pipeline_scripts/surf/selfcal_per_box.sh ${FIELD} ${N} &
done
wait

#MERGE ALL H5 FILES
singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/lofar_helpers/pipeline_scripts/surf/merge_selfcals.py -d ${TO}/selfcal

echo "----------END RECALIBRATING----------"