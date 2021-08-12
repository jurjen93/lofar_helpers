#!/bin/bash

BOX=$1

FOLDER=/project/lofarvwf/Share/jdejong/output/A399/selfcal_test2
RESULT=${FOLDER}/box_${BOX}_result
DATAFOLDER=${FOLDER}/box_${BOX}

mkdir ${RESULT}
mv ${RESULT}/image_*.png ${DATAFOLDER}
mv ${RESULT}/box_*.dysco.sub.shift.avg.weights.ms.archive*.avg ${DATAFOLDER}
mv ${RESULT}/plotlosoto* ${DATAFOLDER}