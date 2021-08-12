#!/bin/bash

BOX=$1

FOLDER=/project/lofarvwf/Share/jdejong/output/A399/selfcal_test2
RESULT=${FOLDER}/box_${BOX}_result
DATAFOLDER=${FOLDER}/box_${BOX}

mkdir ${RESULT}
mv ${DATAFOLDER}/image_*.png ${RESULT}
mv ${DATAFOLDER}/box_*.dysco.sub.shift.avg.weights.ms.archive*.avg ${RESULT}
mv ${DATAFOLDER}/plotlosoto* ${RESULT}