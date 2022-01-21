#!/bin/bash

#input
H5=$1
MS=$2

NAME=image_test_A399_cleanbridge

SING_BIND=/tmp,/dev/shm,/disks/paradata,/data1,/net/lofar1,/net/rijn,/net/nederrijn/,/net/bovenrijn,/net/botlek,/net/para10,/net/lofar2,/net/lofar3,/net/lofar4,/net/lofar5,/net/lofar6,/net/lofar7,/disks/ftphome,/net/krommerijn,/net/voorrijn,/net/achterrijn,/net/tussenrijn,/net/ouderijn,/net/nieuwerijn,/net/lofar8,/net/lofar9,/net/rijn8,/net/rijn7,/net/rijn5,/net/rijn4,/net/rijn3,/net/rijn2
SING_IMAGE=/net/rijn/data2/rvweeren/data/pill-latestJune2021.simg
SING_IMAGE_WSCLEAN=/net/lofar1/data1/sweijen/software/LOFAR/singularity/test/idgtest.sif
SING_IMAGE_P2=/net/lofar1/data1/sweijen/software/LOFAR/singularity/lofar_sksp_fedora31_ddf.sif

TO=/net/${HOSTNAME%%.*}/data2/jurjendejong/Abell399-401_cleanbridge
FROM=/net/tussenrijn/data2/jurjendejong/A399_extracted_avg

singularity exec -B ${SING_BIND} ${SING_IMAGE_WSCLEAN} \
wsclean \
-size 6000 6000 \
-use-wgridder \
-no-update-model-required \
-reorder \
-channels-out 6 \
-weight briggs -0.5 \
-weighting-rank-filter 3 \
-clean-border 1 \
-parallel-reordering 6 \
-padding 1.2 \
-fits-mask ${NAME}_compact-MFS-image.fits.mask.fits \
-pol i \
-name ${NAME}_compactmask \
-scale 1.5arcsec \
-niter 50000 \
-mgain 0.8 \
-fit-beam \
-multiscale \
-join-channels \
-multiscale-max-scales 4 \
-nmiter 9 \
-log-time \
-multiscale-scale-bias 0.7 \
-facet-regions tess.reg \
-minuv-l 2100.0 \
-parallel-gridding 6 \
-fit-spectral-pol 3 \
-apply-facet-solutions ${H5// /,} amplitude000,phase000 \
${MS}

#predict
singularity exec -B ${SING_BIND} ${SING_IMAGE_WSCLEAN} \
wsclean \
-size 1500 1500 \
-channels-out 6 \
-padding 1.2 \
-predict \
-name ${NAME}_compactmask \
${MS}

#subtract
singularity exec -B ${SING_BIND} ${SING_IMAGE} \
python ~/scripts/lofar_helpers/supporting_scripts/substract_mscols.py --ms ${MS} --colname DIFFUSE_SUB