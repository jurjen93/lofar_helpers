#!/bin/bash

SING_BIND=/tmp,/dev/shm,/disks/paradata,/data1,/net/lofar1,/net/rijn,/net/nederrijn/,/net/bovenrijn,/net/botlek,/net/para10,/net/lofar2,/net/lofar3,/net/lofar4,/net/lofar5,/net/lofar6,/net/lofar7,/disks/ftphome,/net/krommerijn,/net/voorrijn,/net/achterrijn,/net/tussenrijn,/net/ouderijn,/net/nieuwerijn,/net/lofar8,/net/lofar9,/net/rijn8,/net/rijn7,/net/rijn5,/net/rijn3,/net/rijn2
SING_IMAGE=/net/rijn/data2/rvweeren/data/pill-latestJune2021.simg
SING_IMAGE_WSCLEAN=/net/lofar1/data1/sweijen/software/LOFAR/singularity/test/idgtest_23_02_2022.sif
SING_IMAGE_P2=/net/lofar1/data1/sweijen/software/LOFAR/singularity/lofar_sksp_fedora31_ddf.sif

TO=/net/${HOSTNAME%%.*}/data2/jurjendejong/Abell399-401_cleanbridge_300kpc
FROM=/net/rijn5/data2/jurjendejong/A399_extracted_avg
H5=all_directions0_freqav2.h5
MS=Abell399-401_extr.dysco.sub.shift.avg.weights.ms.archive0.avg.goodtimes
TESS=tessupdate.reg
NAME=image_test_A399_cleanbridge

singularity exec -B ${SING_BIND} ${SING_IMAGE} CleanSHM.py

cd ${TO}

#predict
singularity exec -B ${SING_BIND} ${SING_IMAGE_WSCLEAN} \
wsclean \
-size 6000 6000 \
-use-wgridder \
-channels-out 6 \
-padding 1.2 \
-predict \
-name ${NAME}_compactmask \
-apply-facet-solutions ${H5} amplitude000,phase000 \
-facet-regions ${TESS} \
${MS}

#subtract
singularity exec -B ${SING_BIND} ${SING_IMAGE} \
python ~/scripts/lofar_helpers/supporting_scripts/substract_mscols.py --ms ${MS} --colname DIFFUSE_SUB

#make final image
singularity exec -B ${SING_BIND} ${SING_IMAGE_WSCLEAN} \
wsclean \
-size 1500 1500 \
-use-wgridder \
-no-update-model-required \
-reorder \
-channels-out 6 \
-weight briggs -0.5 \
-weighting-rank-filter 3 \
-clean-border 1 \
-parallel-reordering 6 \
-padding 1.2 \
-auto-mask 2.5 \
-auto-threshold 0.5 \
-pol i \
-name ${NAME} \
-scale 6arcsec \
-niter 50000 \
-mgain 0.8 \
-fit-beam \
-multiscale \
-join-channels \
-multiscale-max-scales 10 \
-nmiter 8 \
-log-time \
-multiscale-scale-bias 0.7 \
-facet-regions ${TESS} \
-parallel-gridding 6 \
-fit-spectral-pol 3 \
-taper-gaussian 60arcsec \
-data-column DIFFUSE_SUB \
-apply-facet-solutions ${H5} amplitude000,phase000 \
${MS}