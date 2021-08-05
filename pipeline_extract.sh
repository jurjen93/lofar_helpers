#!/bin/bash
#SBATCH -c 1
#SBATCH --mail-type=END,FAIL

echo "----------START RECALIBRATING----------"

FIELD=$1
TO=/project/lofarvwf/Share/jdejong/output/${FIELD}
SCRIPT_FOLDER=/home/lofarvwf-jdejong/scripts/lofar_helpers

SING_IMAGE=/home/lofarvwf-jdejong/singularities/pill-latest.simg
SING_BIND=/project/lofarvwf/Share/jdejong

#CREATE FILES [SHOULD HAVE ALREADY BEEN DONE]
#mkdir ${TO}

#CREATE BOXES
echo "Create boxes..."
singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/make_boxes.py -f ${TO}/extract/image_full_ampphase_di_m.NS.app.restored.fits -l ${TO} -ac 2.5
rm ${TO}/source_file.csv && rm ${TO}/excluded_sources.csv
TOTAL_BOXES=$(ls -dq ${TO}/boxes/box*.reg | wc -l)
if [[ ${TOTAL_BOXES} = 0 ]]; then
  echo "Boxes selection failed, see slurm output."
  exit
fi
echo "Succesfully created boxes..."

#EXTRACT WITH PARALLEL ARRAY
echo "There are ${TOTAL_BOXES} boxes to extract"
mkdir ${TO}/extract
mkdir ${TO}/extract/finished
sbatch ${SCRIPT_FOLDER}/pipeline_scripts/surf/extract.sh ${FIELD} &
wait &

echo "----------END RECALIBRATING----------"