#!/bin/bash
#SBATCH -c 10

#THIS SCRIPT IS WRITTEN FOR SLURM ON SURFSARA
echo "----------START----------"

SOURCE=$1 #L626678
TO=/project/lofarvwf/Share/jdejong/output/${SOURCE}
SCRIPT_FOLDER=/home/lofarvwf-jdejong/scripts

SING_IMAGE=/home/lofarvwf-jdejong/singularities/pill-latest.simg
SING_BIND=/project/lofarvwf/Share/jdejong

#start box number
START_N=1

#SELFCAL
mkdir ${TO}/selfcal

END_N=$(ls -dq ${TO}/extract/*.dysco.sub.shift.avg.weights.ms.archive* | wc -l)
echo "There are ${END_N} boxes ready for selfcal"
END_N=1

echo "-----STARTED SELFCAL-----"
cd ${TO}/selfcal/
for ((i=${START_N};i<=${END_N};i++)); do
mkdir ${TO}/selfcal/box_${i}
cp -r ${TO}/extract/*box_${i}.dysco.sub.shift.avg.weights.ms.archive0 ${TO}/selfcal/
singularity exec -B ${SING_BIND} ${SING_IMAGE} DPPP msin=${TO}/selfcal/Abell399-401_box_${i}.dysco.sub.shift.avg.weights.ms.archive0 msout.storagemanager=dysco msout=${TO}/selfcal/box_${i}/box_${i}.dysco.sub.shift.avg.weights.ms.archive0.goodtimes msin.ntimes=1500 steps=[]
cd ${TO}/selfcal/box_${i}
singularity exec -B ${SING_BIND} ${SING_IMAGE}  python ${SCRIPT_FOLDER}/runwscleanLBautoR.py -b ${TO}/boxes/box_${i}.reg --auto box_${i}.dysco.sub.shift.avg.weights.ms.archive0.goodtimes
# singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/runwscleanLBauto.py -b ${TO}/boxes/box_${i}.reg --forwidefield --smoothnessconstraint-slow=5.0 --usemodeldataforsolints --slow-soltype=complexgain --usewgridder --avgfreqstep=2 --avgtimestep=2 --imager=DDFACET --useaoflagger --noarchive box_${i}.dysco.sub.shift.avg.weights.ms.archive0.goodtimes
echo "Finished selfcal for box_${i}"
done
echo "-----FINISHED SELFCAL-----"

#merge selfcals
singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/lofar_helpers/pipeline_scripts/surf/merge_selfcals.py -d ${TO}/selfcal

echo "----------END----------"