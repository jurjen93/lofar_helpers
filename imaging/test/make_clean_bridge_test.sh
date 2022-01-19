#!/bin/bash

N=$1
NMITER=$2

NAME=image_test_A399_cleanbridge

SING_BIND=/tmp,/dev/shm,/disks/paradata,/data1,/net/lofar1,/net/rijn,/net/nederrijn/,/net/bovenrijn,/net/botlek,/net/para10,/net/lofar2,/net/lofar3,/net/lofar4,/net/lofar5,/net/lofar6,/net/lofar7,/disks/ftphome,/net/krommerijn,/net/voorrijn,/net/achterrijn,/net/tussenrijn,/net/ouderijn,/net/nieuwerijn,/net/lofar8,/net/lofar9,/net/rijn8,/net/rijn7,/net/rijn5,/net/rijn4,/net/rijn3,/net/rijn2
SING_IMAGE=/net/rijn/data2/rvweeren/data/pill-latestJune2021.simg
SING_IMAGE_WSCLEAN=/net/lofar1/data1/sweijen/software/LOFAR/singularity/test/idgtest.sif
TO=/net/${HOSTNAME%%.*}/data2/jurjendejong/Abell399-401_${N}_cleanbridge
FROM=/net/tussenrijn/data2/jurjendejong/A399_extracted_avg
H5=all_directions${N}.h5
MS=Abell399-401_extr.dysco.sub.shift.avg.weights.ms.archive${N}.avg.goodtimes

singularity exec -B ${SING_BIND} ${SING_IMAGE} CleanSHM.py

#make dir
mkdir -p ${TO}

#copy files
cp ${FROM}/${H5} ${TO} && wait
cp -r ${FROM}/${MS} ${TO} && wait

#make shorter time axis
singularity exec -B ${SING_BIND} ${SING_IMAGE} \
python /home/jurjendejong/scripts/lofar_helpers/supporting_scripts/flag_time.py \
--time_flag 0 300 \
-msin ${TO}/${MS} \
-msout ${TO}/${MS}.test && wait
rm -rf ${TO}/${MS}

singularity exec -B ${SING_BIND} ${SING_IMAGE} \
python ~/scripts/lofar_helpers/h5_merger.py \
--h5_tables ${TO}/${H5} \
--h5_out ${TO}/short_${H5} \
--ms ${TO}/${MS}.test && wait
rm -rf ${TO}/${H5}

#aoflagger
singularity exec -B ${SING_BIND} ${SING_IMAGE} \
aoflagger ${TO}/${MS}.test && wait

cd ${TO}

#make tesselation
singularity exec -B ${SING_BIND} ${SING_IMAGE} \
python /net/rijn/data2/rvweeren/LoTSS_ClusterCAL/ds9facetgenerator.py \
--h5 ${TO}/short_${H5} \
--DS9regionout ${TO}/tess.reg \
--imsize 6000 \
--ms ${TO}/${MS}.test

# make first image
singularity exec -B ${SING_BIND} ${SING_IMAGE_WSCLEAN} \
wsclean \
-size 1500 1500 \
-use-wgridder \
-no-update-model-required \
-reorder \
-channels-out 2 \
-weight briggs -0.5 \
-weighting-rank-filter 3 \
-clean-border 1 \
-parallel-reordering 6 \
-padding 1.2 \
-auto-mask 2.5 \
-auto-threshold 0.5 \
-pol i \
-name ${NAME}_compact \
-scale 6arcsec \
-niter 50000 \
-mgain 0.8 \
-fit-beam \
-join-channels \
-nmiter ${NMITER} \
-log-time \
-facet-regions tess.reg \
-minuv-l 2000.0 \
-parallel-gridding 6 \
-fit-spectral-pol 3 \
-apply-facet-solutions short_${H5} amplitude000,phase000 \
${MS}.test

##mask compact objects
#singularity exec -B ${SING_BIND} ${SING_IMAGE} \
#python /net/para10/data1/shimwell/software/killmsddf/new-install/DDFacet/SkyModel/MakeMask.py \
#--Th=3.0 \
#--RestoredIm=${NAME}_compact-MFS-image.fits
#
#singularity exec -B ${SING_BIND} ${SING_IMAGE_WSCLEAN} \
#wsclean \
#-size 1500 1500 \
#-use-wgridder \
#-no-update-model-required \
#-reorder \
#-channels-out 2 \
#-weight briggs -0.5 \
#-weighting-rank-filter 3 \
#-clean-border 1 \
#-parallel-reordering 6 \
#-padding 1.2 \
#-fits-mask ${NAME}_compact-MFS-image.fits.mask.fits \
#-pol i \
#-name ${NAME}_compactmask \
#-scale 6arcsec \
#-niter 50000 \
#-mgain 0.8 \
#-fit-beam \
#-multiscale \
#-join-channels \
#-multiscale-max-scales 4 \
#-nmiter ${NMITER} \
#-log-time \
#-multiscale-scale-bias 0.7 \
#-facet-regions tess.reg \
#-minuv-l 2000.0 \
#-parallel-gridding 6 \
#-fit-spectral-pol 3 \
#-apply-facet-solutions short_${H5} amplitude000,phase000 \
#${MS}.test
#
##predict
#singularity exec -B ${SING_BIND} ${SING_IMAGE_WSCLEAN} \
#wsclean \
#-size 1500 1500 \
#-channels-out 6 \
#-padding 1.2 \
#-predict \
#-name ${NAME}_compactmask
#
##subtract
#singularity exec -B ${SING_BIND} ${SING_IMAGE} \
#python ~/scripts/lofar_helpers/supporting_scripts/substract_mscols.py --ms ${MS}.test --colname DIFFUSE_SUB
#
##make final image
#singularity exec -B ${SING_BIND} ${SING_IMAGE_WSCLEAN} \
#wsclean \
#-size 1500 1500 \
#-use-wgridder \
#-no-update-model-required \
#-reorder \
#-channels-out 2 \
#-weight briggs -0.5 \
#-weighting-rank-filter 3 \
#-clean-border 1 \
#-parallel-reordering 6 \
#-padding 1.2 \
#-auto-mask 2.5 \
#-auto-threshold 0.5 \
#-pol i \
#-name ${NAME} \
#-scale 6arcsec \
#-niter 50000 \
#-mgain 0.8 \
#-fit-beam \
#-multiscale \
#-join-channels \
#-multiscale-max-scales 10 \
#-nmiter ${NMITER} \
#-log-time \
#-multiscale-scale-bias 0.7 \
#-facet-regions tess.reg \
#-parallel-gridding 6 \
#-fit-spectral-pol 3 \
#-taper-gaussian 60arcsec \
#-data-column DIFFUSE_SUB \
#-apply-facet-solutions short_${H5} amplitude000,phase000 \
#${MS}.test