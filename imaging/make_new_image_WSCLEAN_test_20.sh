#!/bin/bash

N=$1
NMITER=$2

SING_BIND=/tmp,/dev/shm,/disks/paradata,/data1,/net/lofar1,/net/rijn,/net/nederrijn/,/net/bovenrijn,/net/botlek,/net/para10,/net/lofar2,/net/lofar3,/net/lofar4,/net/lofar5,/net/lofar6,/net/lofar7,/disks/ftphome,/net/krommerijn,/net/voorrijn,/net/achterrijn,/net/tussenrijn,/net/ouderijn,/net/nieuwerijn,/net/lofar8,/net/lofar9,/net/rijn8,/net/rijn7,/net/rijn5,/net/rijn4,/net/rijn3,/net/rijn2
SING_IMAGE=/net/rijn/data2/rvweeren/data/pill-latestJune2021.simg
SING_IMAGE_WSCLEAN=/net/lofar1/data1/sweijen/software/LOFAR/singularity/test/test_wsclean_facet_fix_sep30.sif
TO=/net/nieuwerijn/data2/jurjendejong/Abell399-401_${N}_20
FROM=/net/tussenrijn/data2/jurjendejong/A399_extracted_avg
H5=all_directions${N}.h5
MS=Abell399-401_extr.dysco.sub.shift.avg.weights.ms.archive${N}.avg.goodtimes

singularity exec -B ${SING_BIND} ${SING_IMAGE} CleanSHM.py
mkdir -p ${TO}
cp ${FROM}/${H5} ${TO}
cp -r ${FROM}/${MS} ${TO}
singularity exec -B ${SING_BIND} ${SING_IMAGE} python /net/rijn/data2/rvweeren/LoTSS_ClusterCAL/ds9facetgenerator.py --h5 ${TO}/${H5} --DS9regionout ${TO}/tess.reg --imsize 3000 --ms ${TO}/${MS}
singularity exec -B ${SING_BIND} ${SING_IMAGE_WSCLEAN} python ~/scripts/lofar_helpers/imaging/make_new_image_WSCLEAN_test_20.py --N ${N} --nmiter ${NMITER}