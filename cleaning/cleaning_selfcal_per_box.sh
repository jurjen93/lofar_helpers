#!/bin/bash

BOX=$1

FOLDER=/project/lofarvwf/Share/jdejong/output/A399/selfcal_test2
#RESULT=${FOLDER}/box_${BOX}_result
DATAFOLDER=${FOLDER}/box_${BOX}

rm ${DATAFOLDER}/*.ddfcache
rm ${DATAFOLDER}/*.py
#rm ${DATAFOLDER}/merged_selfcalcyle000*.h5
#rm ${DATAFOLDER}/merged_selfcalcyle001*.h5
#rm ${DATAFOLDER}/merged_selfcalcyle002*.h5
#rm ${DATAFOLDER}/merged_selfcalcyle003*.h5
#rm ${DATAFOLDER}/merged_selfcalcyle004*.h5
#rm ${DATAFOLDER}/merged_selfcalcyle005*.h5
#rm ${DATAFOLDER}/merged_selfcalcyle006*.h5

#mkdir ${RESULT}
#mv ${DATAFOLDER}/image_*.png ${RESULT}
#mv ${DATAFOLDER}/*box_*.dysco.sub.shift.avg.weights.ms.archive*.avg ${RESULT}
#mv ${DATAFOLDER}/plotlosoto* ${RESULT}
#mv ${DATAFOLDER}/merged_selfcalcyle007*.h5 ${RESULT} # note that this is for max 8 cycles
#mv ${RESULT} ${DATAFOLDER}