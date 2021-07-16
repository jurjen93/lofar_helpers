#!/bin/bash

SOURCE=$1 #L626678
TO=/net/tussenrijn/data2/jurjendejong/${SOURCE}
SCRIPT_FOLDER=/net/rijn/data2/rvweeren/LoTSS_ClusterCAL
N=$2 #box number

#SING_IMAGE=/home/lofarvwf-jdejong/singularities/pill-latest.simg
SING_IMAGE=/net/lofar1/data1/sweijen/software/LOFAR/singularity/lofar_sksp_fedora31_ddf.sif
SING_BIND=/tmp,/dev/shm,/disks/paradata,/data1,/net/lofar1,/net/rijn,/net/nederrijn/,/net/bovenrijn,/net/botlek,/net/para10,/net/lofar2,/net/lofar3,/net/lofar4,/net/lofar5,/net/lofar6,/net/lofar7,/disks/ftphome,/net/krommerijn,/net/voorrijn,/net/achterrijn,/net/tussenrijn,/net/ouderijn,/net/nieuwerijn,/net/lofar8,/net/lofar9,/net/rijn8,/net/rijn7,/net/rijn5,/net/rijn4,/net/rijn3,/net/rijn2

singularity exec -B ${SING_BIND} ${SING_IMAGE} CleanSHM.py

echo "-----STARTED SELFCAL-----"
echo "Started selfcal for box_${N}"
cd ${TO}/selfcal/
mkdir ${TO}/selfcal/box_${N}
cp -r ${TO}/extract/*box_${N}.dysco.sub.shift.avg.weights.ms.archive0 ${TO}/selfcal/
singularity exec -B ${SING_BIND} ${SING_IMAGE} DPPP msin=${TO}/selfcal/Abell399-401_box_${N}.dysco.sub.shift.avg.weights.ms.archive0 msout.storagemanager=dysco msout=${TO}/selfcal/box_${N}/box_${N}.dysco.sub.shift.avg.weights.ms.archive0.goodtimes msin.ntimes=1500 steps=[]
cd ${TO}/selfcal/box_${N}
singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/runwscleanLBautoR.py -b ${TO}/boxes/box_${N}.reg --auto --imager=DDFACET box_${N}.dysco.sub.shift.avg.weights.ms.archive0.goodtimes --helperscriptspath ${SCRIPT_FOLDER}/
echo "Finished selfcal for box_${N}"
echo "-----FINISHED SELFCAL-----"