#!/bin/bash
#SBATCH -c 1
#SBATCH --mail-type=END,FAIL

echo "----------START RECALIBRATING----------"

FIELD=$1
TO=/project/lofarvwf/Share/jdejong/output/${FIELD}
SCRIPT_FOLDER=/home/lofarvwf-jdejong/scripts/lofar_helpers

SING_IMAGE=/home/lofarvwf-jdejong/singularities/pill-latest.simg
SING_BIND=/project/lofarvwf/Share/jdejong,

#CREATE FILES [SHOULD HAVE ALREADY BEEN DONE]
#mkdir ${TO}

#CREATE BOXES
echo "Create boxes..."
singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/make_boxes.py -f ${TO}/extract/image_full_ampphase_di_m.NS.app.restored.fits -l ${TO} -ac 2.5
TOTAL_BOXES=$(ls -dq ${TO}/boxes/box*.reg | wc -l)
if [[ ${TOTAL_BOXES} = 0 ]]; then
  echo "Boxes selection failed, see slurm output."
  exit
fi
echo "Succesfully created boxes..."

#EXTRACT WITH PARALLEL ARRAY
echo "There are ${TOTAL_BOXES} boxes to extract"
mkdir ${TO}/extract
sbatch ${SCRIPT_FOLDER}/pipeline_scripts/test_scripts/extract_test_1_box.sh ${FIELD} &
wait &

#SELFCAL
mkdir ${TO}/selfcal
for ((N=1;N<=2;N++))
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

#MOVE H5 SOLUTION DONE ON STRW
#srun ${SCRIPT_FOLDER}/move_files/move_result/move_result_selfcal_surf-strw.sh /project/lofarvwf/Share/jdejong/output/${FIELD}/selfcal/all_directions.h5 /net/tussenrijn/data2/jurjendejong/${FIELD}

echo "----------END RECALIBRATING----------"