#!/bin/bash

#THIS SCRIPT IS WRITTEN FOR SLURM ON SURFSARA
echo "----------START----------"

SOURCE=$1 #L626678
TO=/project/lofarvwf/Share/jdejong/output/${SOURCE}
SCRIPT_FOLDER=/home/lofarvwf-jdejong/scripts
N=$2 #box number

SING_IMAGE=/home/lofarvwf-jdejong/singularities/pill-latest.simg
SING_BIND=/project/lofarvwf/Share/jdejong

#SELFCAL

echo "-----STARTED SELFCAL-----"
echo "Started selfcal for box_${N}"
cd ${TO}/selfcal/
mkdir ${TO}/selfcal/box_${N}
cp -r ${TO}/extract/*box_${N}.dysco.sub.shift.avg.weights.ms.archive0 ${TO}/selfcal/
singularity exec -B ${SING_BIND} ${SING_IMAGE} DPPP msin=${TO}/selfcal/Abell399-401_box_${N}.dysco.sub.shift.avg.weights.ms.archive0 msout.storagemanager=dysco msout=${TO}/selfcal/box_${N}/box_${N}.dysco.sub.shift.avg.weights.ms.archive0.goodtimes msin.ntimes=1500 steps=[]
cd ${TO}/selfcal/box_${N}
singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/runwscleanLBautoR.py -b ${TO}/boxes/box_${N}.reg --auto box_${N}.dysco.sub.shift.avg.weights.ms.archive0.goodtimes
# singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/runwscleanLBauto.py -b ${TO}/boxes/box_${N}.reg --forwidefield --smoothnessconstraint-slow=5.0 --usemodeldataforsolints --slow-soltype=complexgain --usewgridder --avgfreqstep=2 --avgtimestep=2 --imager=DDFACET --useaoflagger --noarchive box_${N}.dysco.sub.shift.avg.weights.ms.archive0.goodtimes
echo "Finished selfcal for box_${N}"
echo "-----FINISHED SELFCAL-----"

#merge selfcals
singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/lofar_helpers/merge_selfcals.py -d ${TO}/selfcal

echo "----------END----------"