#!/bin/bash

#input
#H5='all_directions0.h5 all_directions1.h5 all_directions2.h5 all_directions3.h5 all_directions4.h5 all_directions5.h5'
MS='Abell399-401_extr.dysco.sub.shift.avg.weights.ms.archive0.avg.goodtimes Abell399-401_extr.dysco.sub.shift.avg.weights.ms.archive1.avg.goodtimes Abell399-401_extr.dysco.sub.shift.avg.weights.ms.archive2.avg.goodtimes Abell399-401_extr.dysco.sub.shift.avg.weights.ms.archive3.avg.goodtimes Abell399-401_extr.dysco.sub.shift.avg.weights.ms.archive4.avg.goodtimes Abell399-401_extr.dysco.sub.shift.avg.weights.ms.archive5.avg.goodtimes'
H5='all_directions0.h5 all_directions1.h5 all_directions2.h5 all_directions3.h5 all_directions4.h5 all_directions5.h5'

#parameters
SING_BIND=/tmp,/dev/shm,/disks/paradata,/data1,/net/lofar1,/net/rijn,/net/nederrijn/,/net/bovenrijn,/net/botlek,/net/para10,/net/lofar2,/net/lofar3,/net/lofar4,/net/lofar5,/net/lofar6,/net/lofar7,/disks/ftphome,/net/krommerijn,/net/voorrijn,/net/achterrijn,/net/tussenrijn,/net/ouderijn,/net/nieuwerijn,/net/lofar8,/net/lofar9,/net/rijn8,/net/rijn7,/net/rijn5,/net/rijn4,/net/rijn3,/net/rijn2
SING_IMAGE=/net/rijn/data2/rvweeren/data/pill-latestJune2021.simg
SING_IMAGE_WSCLEAN=/net/lofar1/data1/sweijen/software/LOFAR/singularity/test/idgtest_23_02_2022.sif
TO=/net/${HOSTNAME%%.*}/data2/jurjendejong/Abell399-401_cleanbridge_all_300kpc
FROM=/net/rijn5/data2/jurjendejong/A399_extracted_avg
TESS=tessupdate.reg
NAME=final_image_A399

#check if directory exists
if [[ -f ${TO} ]]
then
  echo "${TO} exists. Exit script"
  exit 0
fi


#cache
singularity exec -B ${SING_BIND} ${SING_IMAGE} CleanSHM.py


#make directory
mkdir -p ${TO}


#copy files
for H in ${H5}
do
  cp ${FROM}/${H} ${TO}
done


#aoflagger
for M in ${MS}
do
  cp -r ${FROM}/${M} ${TO} && wait
  singularity exec -B ${SING_BIND} ${SING_IMAGE} aoflagger ${TO}/${M}
done


cd ${TO}


#make facet
cp ${FROM}/${TESS} ${TO} && wait


# make first image
singularity exec -B ${SING_BIND} ${SING_IMAGE_WSCLEAN} wsclean \
-use-wgridder \
-no-update-model-required \
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
-log-time \
-parallel-deconvolution 1600 \
-parallel-gridding 5 \
-facet-regions ${TESS} \
-apply-facet-solutions ${H5// /,} amplitude000,phase000 \
-name ${NAME}_compact \
-size 6000 6000 \
-scale 1.5arcsec \
-nmiter 7 \
-minuv-l 943.0 \
${MS} > logcompact.txt


#mask compact objects
singularity exec -B ${SING_BIND} ${SING_IMAGE_P2} \
python /net/para10/data1/shimwell/software/killmsddf/new-install/DDFacet/SkyModel/MakeMask.py \
--Th=3.0 \
--RestoredIm=${NAME}_compact-MFS-image.fits


singularity exec -B ${SING_BIND} ${SING_IMAGE_WSCLEAN} wsclean \
-use-wgridder \
-no-update-model-required \
-reorder \
-weight briggs \
-0.5 \
-weighting-rank-filter 3 \
-clean-border 1 \
-parallel-reordering 5 \
-padding 1.2 \
-fits-mask ${NAME}_compact-MFS-image.fits.mask.fits \
-pol i \
-niter 150000 \
-mgain 0.7 \
-fit-beam \
-multiscale \
-channels-out 6 \
-fit-spectral-pol 3 \
-join-channels \
-log-time \
-parallel-deconvolution 1600 \
-parallel-gridding 5 \
-facet-regions ${TESS} \
-apply-facet-solutions ${H5// /,} amplitude000,phase000 \
-name ${NAME}_compactmask \
-size 6000 6000 \
-scale 1.5arcsec \
-nmiter 9 \
-minuv-l 943.0 \
${MS} > logcompactmask.txt


rm -rf ${NAME}_compactmask-0000-model.fits
rm -rf ${NAME}_compactmask-0001-model.fits
rm -rf ${NAME}_compactmask-0002-model.fits
rm -rf ${NAME}_compactmask-0003-model.fits
rm -rf ${NAME}_compactmask-0004-model.fits
rm -rf ${NAME}_compactmask-0005-model.fits
mv ${NAME}_compactmask-0000-model-pb.fits ${NAME}_compactmask-0000-model.fits
mv ${NAME}_compactmask-0001-model-pb.fits ${NAME}_compactmask-0001-model.fits
mv ${NAME}_compactmask-0002-model-pb.fits ${NAME}_compactmask-0002-model.fits
mv ${NAME}_compactmask-0003-model-pb.fits ${NAME}_compactmask-0003-model.fits
mv ${NAME}_compactmask-0004-model-pb.fits ${NAME}_compactmask-0004-model.fits
mv ${NAME}_compactmask-0005-model-pb.fits ${NAME}_compactmask-0005-model.fits


cp -r ${TO} ${TO}_backup && wait


#predict
singularity exec -B ${SING_BIND} ${SING_IMAGE_WSCLEAN} \
wsclean \
-use-wgridder \
-size 6000 6000 \
-channels-out 6 \
-padding 1.2 \
-predict \
-parallel-gridding 5 \
-name ${NAME}_compactmask \
-apply-facet-solutions ${H5// /,} amplitude000,phase000 \
-facet-regions ${TESS} \
${MS} > logpredict.txt


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
-nmiter 5 \
-log-time \
-multiscale-scale-bias 0.7 \
-facet-regions ${TESS} \
-parallel-gridding 6 \
-fit-spectral-pol 3 \
-taper-gaussian 60arcsec \
-data-column DIFFUSE_SUB \
-apply-facet-solutions ${H5// /,} amplitude000,phase000 \
${MS} > logfinal.txt
