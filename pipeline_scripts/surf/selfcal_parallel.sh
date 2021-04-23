#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jurjendejong@strw.leidenuniv.nl

SOURCE=$1 #L626678
TO=/project/lofarvwf/Share/jdejong/output/${SOURCE}

SCRIPT_FOLDER=/home/lofarvwf-jdejong/scripts/lofar_helpers/pipeline_scripts/surf

#START_N=1
#END_N=$(ls -dq ${TO}/extract/L626678/*.dysco.sub.shift.avg.weights.ms.archive* | wc -l)

START_N=2
END_N=4

#PARAMETER SWEEP
for ((i=${START_N};i<=${END_N};i++)); do
    sbatch ${SCRIPT_FOLDER}/selfcal_per_box.sh ${SOURCE} ${i}
    sleep 1
done
wait
