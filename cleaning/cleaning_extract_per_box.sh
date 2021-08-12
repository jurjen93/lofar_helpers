#!/bin/bash

BOX=$1

FOLDER=/project/lofarvwf/Share/jdejong/output/A399/extract
RESULT=${FOLDER}/box_${BOX}_result
DATAFOLDER=${FOLDER}/box_${BOX}

mkdir ${RESULT}
mv ${DATAFOLDER}/*box_${BOX}.dysco.sub.shift.avg.weights.ms.archive* ${RESULT}
rm -rf ${DATAFOLDER} # BE VERY CAREFUL HERE
mv ${RESULT} ${DATAFOLDER}