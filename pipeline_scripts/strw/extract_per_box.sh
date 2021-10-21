#!/bin/bash

FIELD=$1
BOX=$2
TO=/net/nieuwerijn/data2/jurjendejong/${FIELD}/${BOX}
SCRIPT_FOLDER=/home/lofarvwf-jdejong/scripts

SING_IMAGE=/net/lofar1/data1/sweijen/software/LOFAR/singularity/lofar_sksp_fedora31_ddf.sif
SING_BIND=/tmp,/dev/shm,/disks/paradata,/data1,/net/lofar1,/net/rijn,/net/nederrijn/,/net/bovenrijn,/net/botlek,/net/para10,/net/lofar2,/net/lofar3,/net/lofar4,/net/lofar5,/net/lofar6,/net/lofar7,/disks/ftphome,/net/krommerijn,/net/voorrijn,/net/achterrijn,/net/tussenrijn,/net/ouderijn,/net/nieuwerijn,/net/lofar8,/net/lofar9,/net/rijn8,/net/rijn7,/net/rijn5,/net/rijn4,/net/rijn3,/net/rijn2

#START EXTRACT
echo "-----STARTED EXTRACT-----"
if [[ ! ${BOX} -gt ${TOTAL_BOXES} ]]
then
  mkdir ${TO}/extract/box_${BOX}
  cp ${TO}/extract/data_archive.tar.gz ${TO}/extract/box_${BOX}/
  cd ${TO}/extract/box_${BOX} || { echo "Missing path"; exit 1; }
  tar -xvf data_archive.tar.gz
  rm data_archive.tar.gz
  singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/sub-sources-outside-region.py -b ${TO}/boxes/box_${BOX}.reg --overwriteoutput -p box_${BOX}
  echo "Extracted box_${BOX}"
  echo "Selfcal box_${BOX} finished" > ${TO}/finished/box_${BOX}.txt
else
  :
fi
echo "-----FINISHED EXTRACT-----"