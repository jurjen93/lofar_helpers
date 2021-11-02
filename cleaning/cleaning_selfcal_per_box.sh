#!/bin/bash

BOX=$1

FOLDER=/project/lofarvwf/Share/jdejong/output/A399/selfcal
#RESULT=${FOLDER}/box_${BOX}_result
DATAFOLDER=${FOLDER}/box_${BOX}

rm -rf ${DATAFOLDER}/*.ddfcache
rm ${DATAFOLDER}/*.py
rm ${DATAFOLDER}/merged_selfcalcyle000*.h5
rm ${DATAFOLDER}/merged_selfcalcyle001*.h5
rm ${DATAFOLDER}/merged_selfcalcyle002*.h5
rm ${DATAFOLDER}/merged_selfcalcyle003*.h5
rm ${DATAFOLDER}/merged_selfcalcyle004*.h5
rm ${DATAFOLDER}/merged_selfcalcyle005*.h5
rm ${DATAFOLDER}/merged_selfcalcyle006*.h5
rm -rf ${DATAFOLDER}/*.ddfcache
rm -rf ${DATAFOLDER}/Abell399*box_${BOX}.dysco.sub.shift.avg.weights.ms*
rm -rf ${DATAFOLDER}/*templatejones.h5
rm -rf ${DATAFOLDER}/final_merge*.h5
rm -rf ${DATAFOLDER}/box_${BOX}.tar.gz
rm ${DATAFOLDER}/antennaconstraint.p
rm ${DATAFOLDER}/image_001*
rm ${DATAFOLDER}/image_002*
rm ${DATAFOLDER}/image_003*
rm ${DATAFOLDER}/image_004*
rm ${DATAFOLDER}/image_005*
rm ${DATAFOLDER}/image_006*
rm -rf ${DATAFOLDER}/scalarcomplexgain*
rm -rf ${DATAFOLDER}/tecandphase*
rm ${DATAFOLDER}/smoothness*
rm ${DATAFOLDER}/solint.p
rm ${DATAFOLDER}/soltypecycles.p
rm ${DATAFOLDER}/nchan.p
rm ${DATAFOLDER}/mslist.txt
rm ${DATAFOLDER}/last_DDFacet.obj
rm ${DATAFOLDER}/last_MyCasapy2BBS.obj
rm ${DATAFOLDER}/losotobeam.parset
rm ${DATAFOLDER}/losoto_flag_apgrid.parset
rm ${DATAFOLDER}/*brutalRestored*

#mkdir ${RESULT}
#mv ${DATAFOLDER}/image_*.png ${RESULT}
#mv ${DATAFOLDER}/*box_*.dysco.sub.shift.avg.weights.ms.archive*.avg ${RESULT}
#mv ${DATAFOLDER}/plotlosoto* ${RESULT}
#mv ${DATAFOLDER}/merged_selfcalcyle007*.h5 ${RESULT} # note that this is for max 8 cycles
#mv ${RESULT} ${DATAFOLDER}