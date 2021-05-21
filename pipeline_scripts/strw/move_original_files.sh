#!/bin/bash

#THIS SCRIPT IS WRITTEN FOR STRW SERVERS

FROM=$1#/disks/paradata/shimwell/LoTSS-DR2/archive_other/L626678
TO=$2#/net/tussenrijn/data2/jurjendejong/L626678/extract
SING_IMAGE=/net/rijn/data2/rvweeren/data/pill-latest.simg
SING_BIND=/tmp,/dev/shm,/disks/paradata,/data1,/net/lofar1,/net/rijn,/net/nederrijn/,/net/bovenrijn,/net/botlek,/net/para10,/net/lofar2,/net/lofar3,/net/lofar4,/net/lofar5,/net/lofar6,/net/lofar7,/disks/ftphome,/net/krommerijn,/net/voorrijn,/net/achterrijn,/net/tussenrijn,/net/ouderijn,/net/nieuwerijn,/net/lofar8,/net/lofar9,/net/rijn8,/net/rijn7,/net/rijn5,/net/rijn4,/net/rijn3,/net/rijn2

#Clean cash
singularity exec -B ${SING_BIND} ${SING_IMAGE} CleanSHM.py

#MAKE DIRECTORIES
mkdir ${TO}
echo "Created ${TO}"

#MOVE NEEDED FILES
echo "Moving files to ${TO}..."
cp -r ${FROM}/DDS3_full_*.01_merged.npz ${TO}
cp -r ${FROM}/DDS3_full_*.01_smoothed.npz ${TO}
cp -r ${FROM}/*_uv.pre-cal_*.pre-cal.ms.archive ${TO}
cp -r ${FROM}/SOLSDIR ${TO}
cp -r ${FROM}/image_dirin_SSD_m.npy.ClusterCat.npy ${TO}
cp -r ${FROM}/image_full_ampphase_di_m.NS.DicoModel ${TO}
cp -r ${FROM}/image_full_ampphase_di_m.NS.mask01.fits ${TO}
cp -r ${FROM}/image_full_ampphase_di_m.NS.tessel.reg ${TO}
cp -r ${FROM}/image_full_ampphase_di_m.NS.app.restored.fits ${TO}
echo "Succesfully finished moving files..."