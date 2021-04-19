#!/bin/bash
#SBATCH -c 10

#THIS SCRIPT IS WRITTEN FOR SLURM ON SURFSARA
echo "----------START----------"

SOURCE=$1 #L626678
TO=/project/lofarvwf/Share/jdejong/output/${SOURCE}
SCRIPT_FOLDER=/home/lofarvwf-jdejong/scripts
i=$2#box number

SING_IMAGE=/home/lofarvwf-jdejong/singularities/pill-latest.simg
SING_BIND=/project/lofarvwf/Share/jdejong

#SELFCAL

echo "There are ${END_N} boxes ready for selfcal"

echo "-----STARTED SELFCAL-----"
cd ${TO}/selfcal/
mkdir ${TO}/selfcal/box_${i}
cp -r ${TO}/extract/*box_${i}.dysco.sub.shift.avg.weights.ms.archive0 ${TO}/selfcal/
singularity exec -B ${SING_BIND} ${SING_IMAGE} DPPP msin=${TO}/selfcal/Abell399-401_box_${i}.dysco.sub.shift.avg.weights.ms.archive0 msout.storagemanager=dysco msout=${TO}/selfcal/box_${i}/box_${i}.dysco.sub.shift.avg.weights.ms.archive0.goodtimes msin.ntimes=1500 steps=[]
cd ${TO}/selfcal/box_${i}
singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/runwscleanLBautoR.py -b ${TO}/boxes/box_${i}.reg --auto box_${i}.dysco.sub.shift.avg.weights.ms.archive0.goodtimes
# singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/runwscleanLBauto.py -b ${TO}/boxes/box_${i}.reg --forwidefield --smoothnessconstraint-slow=5.0 --usemodeldataforsolints --slow-soltype=complexgain --usewgridder --avgfreqstep=2 --avgtimestep=2 --imager=DDFACET --useaoflagger --noarchive box_${i}.dysco.sub.shift.avg.weights.ms.archive0.goodtimes
echo "Finished selfcal for box_${i}"
echo "-----FINISHED SELFCAL-----"

#merge selfcals
singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/lofar_helpers/merge_selfcals.py -d ${TO}/selfcal

echo "----------END----------"