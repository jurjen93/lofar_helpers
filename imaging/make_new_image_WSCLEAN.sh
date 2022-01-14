#!/bin/bash

#parameters
SING_BIND=/tmp,/dev/shm,/disks/paradata,/data1,/net/lofar1,/net/rijn,/net/nederrijn/,/net/bovenrijn,/net/botlek,/net/para10,/net/lofar2,/net/lofar3,/net/lofar4,/net/lofar5,/net/lofar6,/net/lofar7,/disks/ftphome,/net/krommerijn,/net/voorrijn,/net/achterrijn,/net/tussenrijn,/net/ouderijn,/net/nieuwerijn,/net/lofar8,/net/lofar9,/net/rijn8,/net/rijn7,/net/rijn5,/net/rijn4,/net/rijn3,/net/rijn2
SING_IMAGE=/net/rijn/data2/rvweeren/data/pill-latestJune2021.simg
SING_IMAGE_WSCLEAN=/net/lofar1/data1/sweijen/software/LOFAR/singularity/test/test_wsclean_facet_fix_sep30.sif
TO=/net/nieuwerijn/data2/jurjendejong/Abell399-401
FROM=/net/tussenrijn/data2/jurjendejong/A399_extracted_avg
H5_1=all_directions0.h5
H5_2=all_directions1.h5
MS_1=Abell399-401_extr.dysco.sub.shift.avg.weights.ms.archive0.avg.goodtimes
MS_2=Abell399-401_extr.dysco.sub.shift.avg.weights.ms.archive1.avg.goodtimes

#cache
singularity exec -B ${SING_BIND} ${SING_IMAGE} CleanSHM.py

#make directory
mkdir -p ${TO}

#copy files
cp ${FROM}/${H5_1} ${TO} && \
cp ${FROM}/${H5_2} ${TO} && \
cp -r ${FROM}/${MS_1} ${TO} && \
cp -r ${FROM}/${MS_2} ${TO} && wait

#aoflagger
aoflagger ${TO}/${MS_1} && \
aoflagger ${TO}/${MS_2}

cd ${TO}

#make facet
cp ${FROM}/tessupdate.reg ${TO} && wait
#singularity exec -B ${SING_BIND} ${SING_IMAGE} python /net/rijn/data2/rvweeren/LoTSS_ClusterCAL/ds9facetgenerator.py \
#--h5 ${TO}/${H5} \
#--DS9regionout ${TO}/tess.reg \
#--imsize 6000 \
#--ms ${TO}/${MS_1}

#run wsclean
singularity exec -B ${SING_BIND} ${SING_IMAGE_WSCLEAN} wsclean \
-data-column DATA \
-use-wgridder \
-update-model-required \
-reorder \
-weight briggs \
-0.5 \
-weighting-rank-filter 3 \
-clean-border 1 \
-parallel-reordering 5 \
-padding 1.2 \
-auto-mask 2.5 \
-auto-threshold 0.5 \
-pol i \
-niter 150000 \
-mgain 0.7 \
-fit-beam \
-multiscale \
-channels-out 6 \
-fit-spectral-pol 3 \
-join-channels \
-multiscale-max-scales 10 \
-log-time \
-multiscale-scale-bias 0.7 \
-parallel-deconvolution 1600 \
-parallel-gridding 5 \
-facet-regions ${TO}/tessupdate.reg \
-apply-facet-solutions ${H5_1} ${H5_2} amplitude000,phase000 \
-name image_test_A399 \
-size 6000 6000 \
-scale 1.5arcsec \
-nmiter 11 \
-minuv-l 80.0 \
${MS_1} ${MS_2}