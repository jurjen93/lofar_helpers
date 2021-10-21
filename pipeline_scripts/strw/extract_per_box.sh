#!/bin/bash

FIELD=$1
BOX=$2
TO=/net/nieuwerijn/data2/jurjendejong/${BOX}
SCRIPT_FOLDER=/net/rijn/data2/rvweeren/LoTSS_ClusterCAL

SING_IMAGE=/net/lofar1/data1/sweijen/software/LOFAR/singularity/lofar_sksp_fedora31_ddf.sif
SING_BIND=/tmp,/dev/shm,/disks/paradata,/data1,/net/lofar1,/net/rijn,/net/nederrijn/,/net/bovenrijn,/net/botlek,/net/para10,/net/lofar2,/net/lofar3,/net/lofar4,/net/lofar5,/net/lofar6,/net/lofar7,/disks/ftphome,/net/krommerijn,/net/voorrijn,/net/achterrijn,/net/tussenrijn,/net/ouderijn,/net/nieuwerijn,/net/lofar8,/net/lofar9,/net/rijn8,/net/rijn7,/net/rijn5,/net/rijn4,/net/rijn3,/net/rijn2

#START EXTRACT
singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/sub-sources-outside-region.py -b box_${BOX}.reg --overwriteoutput -p box_${BOX}

echo "-----FINISHED EXTRACT-----"