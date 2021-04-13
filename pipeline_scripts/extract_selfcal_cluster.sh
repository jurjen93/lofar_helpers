#!/bin/bash

#SBATCH -c 1 -t 60:00

echo "---START----" 

SOURCE=$1 #L626678
TO=/project/lofarvwf/Share/jdejong/output/${SOURCE}
SCRIPT_FOLDER=/home/lofarvwf-jdejong/scripts

#singularity exec -B /tmp,/dev/shm,/disks/paradata,/data1,/net/lofar1,/net/rijn,/net/nederrijn/,/net/bovenrijn,/net/botlek,/net/para10,/net/lofar2,/net/lofar3,/net/lofar4,/net/lofar5,/net/lofar6,/net/lofar7,/disks/ftphome,/net/krommerijn,/net/voorrijn,/net/achterrijn,/net/tussenrijn,/net/ouderijn,/net/nieuwerijn,/net/lofar8,/net/lofar9 /net/lofar1/data1/sweijen/software/LOFAR/singularity/lofar_sksp_fedora27_ddf.sif
SING_IMAGE=/home/lofarvwf-jdejong/singularities/pill-latest.simg
SING_BIND=/project/lofarvwf/Share/jdejong

#start box number
START_N=0

#Clean cash
# singularity exec -B ${SING_BIND} ${SING_IMAGE} python CleanSHM.py

#REORGANIZE FILES
cp ~/scripts/lofar_helpers/h5_merger.py ~/scripts

#MAKE DIRECTORIES
if test -f "${TO}"; then
  echo "${TO} exists"
else
  mkdir ${TO}
  echo "Created ${TO}"
mkdir ${TO}/extract
echo "Created ${TO}/extract"

#MOVE NEEDED FILES
echo "Moving files to ${FOLDER}/extract..."
cp -r /project/lofarvwf/Share/jdejong/data/${SOURCE}/DDS3_full_*.01_merged.npz ${TO}/extract
cp -r /project/lofarvwf/Share/jdejong/data/${SOURCE}/DDS3_full_*.01_smoothed.npz ${TO}/extract
cp -r /project/lofarvwf/Share/jdejong/data/${SOURCE}/*_uv.pre-cal_*.pre-cal.ms.archive ${TO}/extract
cp -r /project/lofarvwf/Share/jdejong/data/${SOURCE}/SOLSDIR ${TO}/extract
cp -r /project/lofarvwf/Share/jdejong/data/${SOURCE}/image_dirin_SSD_m.npy.ClusterCat.npy ${TO}/extract
cp -r /project/lofarvwf/Share/jdejong/data/${SOURCE}/image_full_ampphase_di_m.NS.DicoModel ${TO}/extract
cp -r /project/lofarvwf/Share/jdejong/data/${SOURCE}/image_full_ampphase_di_m.NS.mask01.fits ${TO}/extract
cp -r /project/lofarvwf/Share/jdejong/data/${SOURCE}/image_full_ampphase_di_m.NS.tessel.reg ${TO}/extract
cp -r /project/lofarvwf/Share/jdejong/data/${SOURCE}/image_full_ampphase_di_m.NS.app.restored.fits ${TO}/extract
echo "Succesfully finished moving files..."

cd ${TO}/extract

#CREATE BOXES
echo "Create boxes..."
python3 ${SCRIPT_FOLDER}/lofar_helpers/make_boxes.py -f ${TO}/extract/image_full_ampphase_di_m.NS.app.restored.fits -l ${TO} -i false
echo "Succesfully created boxes..."

END_N=$(ls -dq ${TO}/boxes/box*.reg | wc -l)

#EXTRACT
echo "Started extract..."
for ((i=${START_N};i<=${END_N};i++)); do 
sinularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/sub-sources-outside-region.py -b ${TO}/boxes/box_${i}.reg --overwriteoutput -p box_${i} 
echo "Extracted box_${i}"
done
echo "Finished extract..."

#SELFCAL
mkdir ${TO}/selfcal

echo "There are ${START_N} boxes"
END_N=$(ls -dq ${TO}/extract/*.dysco.sub.shift.avg.weights.ms.archive* | wc -l)
echo "There are ${START_N} boxes ready for selfcal"

cd ${TO}/selfcal/
for ((i=${START_N};i<=${END_N};i++)); do
mkdir ${TO}/selfcal/box_${i}
echo "Created ${TO}/selfcal/box_${i}"
cp -r ${TO}/extract/*box_${i}.dysco.sub.shift.avg.weights.ms.archive0 ${TO}/selfcal/
echo "Moved box_${i}*.ms.archive0 to ${TO}/selfcal/"
cp -r ${TO}/extract/boxes/box_${i}.reg ${TO}/selfcal/box_${i}/
echo "Extracting goodtimes from box_${i}"
singularity exec -B ${SING_BIND} ${SING_IMAGE} DPPP msin=${TO}/selfcal/Abell399-401_box_${i}.dysco.sub.shift.avg.weights.ms.archive0 msout.storagemanager=dysco msout=${TO}/selfcal/box_${i}/box_${i}.dysco.sub.shift.avg.weights.ms.archive0.goodtimes msin.ntimes=1500 steps=[]
echo "Started selfcal box_${i}"
cd ${TO}/selfcal/box_${i}
singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/runwscleanLBautoR.py -b ${TO}/selfcal/box_${i}/box_${i}.reg --auto box_${i}.dysco.sub.shift.avg.weights.ms.archive0.goodtimes
# singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/runwscleanLBauto.py -b ${TO}/selfcal/box_${i}/box_${i}.reg --forwidefield --smoothnessconstraint-slow=5.0 --usemodeldataforsolints --slow-soltype=complexgain --usewgridder --avgfreqstep=2 --avgtimestep=2 --imager=DDFACET --useaoflagger --noarchive box_${i}.dysco.sub.shift.avg.weights.ms.archive0.goodtimes
echo "Finished selfcal for box_${i}"
done

#merge selfcals
singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/lofar_helpers/merge_selfcals.py -d ${TO}/selfcal

echo "---END---"