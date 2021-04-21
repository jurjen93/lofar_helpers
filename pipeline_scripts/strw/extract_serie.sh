#!/bin/bash

#THIS SCRIPT IS WRITTEN FOR STRW SERVERS

SOURCE=$1
FOLDER=/net/tussenrijn/data2/jurjendejong/${SOURCE}
FROM=$2 #/disks/paradata/shimwell/LoTSS-DR2/archive_other/L626678
SING_IMAGE=/net/rijn/data2/rvweeren/data/pill-latest.simg
SING_BIND=/tmp,/dev/shm,/disks/paradata,/data1,/net/lofar1,/net/rijn,/net/nederrijn/,/net/bovenrijn,/net/botlek,/net/para10,/net/lofar2,/net/lofar3,/net/lofar4,/net/lofar5,/net/lofar6,/net/lofar7,/disks/ftphome,/net/krommerijn,/net/voorrijn,/net/achterrijn,/net/tussenrijn,/net/ouderijn,/net/nieuwerijn,/net/lofar8,/net/lofar9,/net/rijn8,/net/rijn7,/net/rijn5,/net/rijn4,/net/rijn3,/net/rijn2

#Clean cash
singularity exec -B ${SING_BIND} ${SING_IMAGE} CleanSHM.py

#MAKE DIRECTORIES
mkdir ${FOLDER}/extract
echo "Created ${FOLDER}/extract"

#MOVE NEEDED FILES
echo "Moving files to ${FOLDER}/extract..."
cp -r ${FROM}/DDS3_full_*.01_merged.npz ${FOLDER}/extract
cp -r ${FROM}/DDS3_full_*.01_smoothed.npz ${FOLDER}/extract
cp -r ${FROM}/*_uv.pre-cal_*.pre-cal.ms.archive ${FOLDER}/extract
cp -r ${FROM}/SOLSDIR ${FOLDER}/extract
cp -r ${FROM}/image_dirin_SSD_m.npy.ClusterCat.npy ${FOLDER}/extract
cp -r ${FROM}/image_full_ampphase_di_m.NS.DicoModel ${FOLDER}/extract
cp -r ${FROM}/image_full_ampphase_di_m.NS.mask01.fits ${FOLDER}/extract
cp -r ${FROM}/image_full_ampphase_di_m.NS.tessel.reg ${FOLDER}/extract
cp -r ${FROM}/image_full_ampphase_di_m.NS.app.restored.fits ${FOLDER}/extract
echo "Succesfully finished moving files..."

#CREATE BOXES
echo "Create boxes..."
python3 ~/Documents/python/lofar_boxes/make_boxes.py -f ${FOLDER}/extract/image_full_ampphase_di_m.NS.app.restored.fits -l ${FOLDER}
echo "Succesfully created boxes..."

START_N=1
END_N=$(ls -dq ${FOLDER}/boxes/box*.reg | wc -l)

#EXTRACT
echo "Started extract..."
for ((i=${START_N};i<=${END_N};i++)); do
singularity exec -B ${SING_BIND} ${SING_IMAGE} python /net/rijn/data2/rvweeren/LoTSS_ClusterCAL/sub-sources-outside-region.py -b ${FOLDER}/boxes/box_${i}.reg --overwriteoutput -p box_${i}
echo "Extracted box_${i}"
done
echo "Finished extract..."