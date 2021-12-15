#!/bin/bash

SING_IMAGE=/net/lofar1/data1/sweijen/software/LOFAR/singularity/lofar_sksp_fedora31_ddf.sif
SING_BIND=/tmp,/dev/shm,/disks/paradata,/data1,/net/lofar1,/net/rijn,/net/nederrijn/,/net/bovenrijn,/net/botlek,/net/para10,/net/lofar2,/net/lofar3,/net/lofar4,/net/lofar5,/net/lofar6,/net/lofar7,/disks/ftphome,/net/krommerijn,/net/voorrijn,/net/achterrijn,/net/tussenrijn,/net/ouderijn,/net/nieuwerijn,/net/lofar8,/net/lofar9,/net/rijn8,/net/rijn7,/net/rijn5,/net/rijn4,/net/rijn3,/net/rijn2
SCRIPT_FOLDER=/home/jurjendejong/scripts/lofar_helpers
EXTRACT_TO=/net/tussenrijn/data2/jurjendejong/A399_extracted

mkdir ${EXTRACT_TO}
mv /net/tussenrijn/data2/jurjendejong/extractedregion.reg ${EXTRACT_TO}
singularity exec -B ${SING_BIND} /net/rijn/data2/rvweeren/data/pill-latest.simg python ${SCRIPT_FOLDER}/move_files/move_extract_files.py --frm /net/rijn/data2/jdejong/A399_DEEP --to ${EXTRACT_TO}
cd ${EXTRACT_TO}
singularity exec -B ${SING_BIND} ${SING_IMAGE} python /net/rijn/data2/rvweeren/LoTSS_ClusterCAL/sub-sources-outside-region.py --boxfile extractedregion.reg --freqavg=2 --timeavg=2 --overwriteoutput --nophaseshift --adjustboxrotation=False --prefixname extr
rm -rf ${EXTRACT_TO}/L6*.ms
rm -rf ${EXTRACT_TO}/L6*.ms.archive
#singularity exec -B ${SING_BIND} ${SING_IMAGE} python /home/jurjendejong/scripts/lofar_helpers/supporting_scripts/make_new_dicomodel.py