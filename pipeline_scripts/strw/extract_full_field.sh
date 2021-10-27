#!/bin/bash

#EXTRACT THE FULL FIELD ON 2.5 DEGREES

SING_IMAGE=/net/lofar1/data1/sweijen/software/LOFAR/singularity/lofar_sksp_fedora31_ddf.sif
SING_BIND=/tmp,/dev/shm,/disks/paradata,/data1,/net/lofar1,/net/rijn,/net/nederrijn/,/net/bovenrijn,/net/botlek,/net/para10,/net/lofar2,/net/lofar3,/net/lofar4,/net/lofar5,/net/lofar6,/net/lofar7,/disks/ftphome,/net/krommerijn,/net/voorrijn,/net/achterrijn,/net/tussenrijn,/net/ouderijn,/net/nieuwerijn,/net/lofar8,/net/lofar9,/net/rijn8,/net/rijn7,/net/rijn5,/net/rijn4,/net/rijn3,/net/rijn2
SCRIPT_FOLDER=/home/jurjendejong/scripts/lofar_helpers

#singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/move_files/move_extract_files.py --frm /net/rijn/data2/jdejong/A399_DEEP --to /net/nieuwerijn/data2/jurjendejong/A399/fullextracted
cd /net/nieuwerijn/data2/jurjendejong/A399/fullextracted
singularity exec -B ${SING_BIND} ${SING_IMAGE} python /net/rijn/data2/rvweeren/LoTSS_ClusterCAL/sub-sources-outside-region.py -b extracted.reg --overwriteoutput --noconcat --nophaseshift --adjustboxrotation=False -p extr
