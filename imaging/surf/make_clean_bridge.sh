#!/bin/bash
#SBATCH --partition=infinite
#SBATCH -c 24
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jurjendejong@strw.leidenuniv.nl


#input
H5=$1
MS=$2

NAME=image_test_A399_cleanbridge

SING_BIND=/project/lofarvwf/Share/jdejong,/home/lofarvwf-jdejong/scripts
SING_IMAGE=/home/lofarvwf-jdejong/singularities/pill-latest.simg
SING_IMAGE_WSCLEAN=/home/lofarvwf-jdejong/singularities/idgtest.sif
SING_IMAGE_P2=/home/lofarvwf-jdejong/singularities/lofar_sksp_fedora31_ddf.sif

TO=/project/lofarvwf/Share/jdejong/output/A399/imaging/Abell399-401_cleanbridge
FROM=/project/lofarvwf/Share/jdejong/output/A399/imaging/A399_extracted_avg

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
cp ${FROM}/tess60.reg ${TO} && wait

# make first image
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
-auto-mask 2.5 \
-auto-threshold 0.5 \
-pol i \
-name ${NAME}_compact \
-scale 1.5arcsec \
-niter 50000 \
-mgain 0.8 \
-fit-beam \
-join-channels \
-nmiter 7 \
-log-time \
-facet-regions tess.reg \
-minuv-l 2100.0 \
-parallel-gridding 6 \
-fit-spectral-pol 3 \
-apply-facet-solutions ${H5// /,} amplitude000,phase000 \
${MS}

#mask compact objects
singularity exec -B ${SING_BIND} ${SING_IMAGE_P2} \
python /home/lofarvwf-jdejong/scripts/MakeMask.py \
--Th=3.0 \
--RestoredIm=${NAME}_compact-MFS-image.fits

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
-nmiter 10 \
-log-time \
-multiscale-scale-bias 0.7 \
-facet-regions tess.reg \
-parallel-gridding 6 \
-fit-spectral-pol 3 \
-taper-gaussian 60arcsec \
-data-column DIFFUSE_SUB \
-apply-facet-solutions ${H5// /,} amplitude000,phase000 \
${MS}