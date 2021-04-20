#!/bin/bash
#SBATCH -c 10

SOURCE=$1 #L626678
TO=/project/lofarvwf/Share/jdejong/output/${SOURCE}

SCRIPT_FOLDER=/home/lofarvwf-jdejong/scripts/lofar_helpers/pipeline_scripts

START_N=1
END_N=$(ls -dq ${TO}/extract/*.dysco.sub.shift.avg.weights.ms.archive* | wc -l)

#PARAMETER SWEEP
for i in 5 6 7; do
  sh SCRIPT_FOLDER/selfcal_per_box.sh ${SOURCE} ${i} &
done
wait