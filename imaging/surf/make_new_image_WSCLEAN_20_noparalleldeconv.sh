#!/bin/bash
#SBATCH --partition=infinite
#SBATCH -c 24
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jurjendejong@strw.leidenuniv.nl

#input
H5=$1
MS=$2

SING_BIND=/project/lofarvwf/Share/jdejong,/home/lofarvwf-jdejong/scripts
SING_IMAGE=/home/lofarvwf-jdejong/singularities/pill-latest.simg
SING_IMAGE_WSCLEAN=/home/lofarvwf-jdejong/singularities/idgtest.sif

TO=/project/lofarvwf/Share/jdejong/output/A399/imaging/Abell399-401_20_$(echo "$H5" | tr -cd ' ' | wc -c)
FROM=/project/lofarvwf/Share/jdejong/output/A399/imaging/A399_extracted_avg
TESS=tess60.reg

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

#run wsclean
singularity exec -B ${SING_BIND} ${SING_IMAGE_WSCLEAN} \
wsclean \
-size 3000 3000 \
-use-wgridder \
-no-update-model-required \
-reorder \
-weight briggs -0.5 \
-weighting-rank-filter 3 \
-clean-border 1 \
-parallel-reordering 6 \
-padding 1.2 \
-auto-mask 2.5 \
-auto-threshold 0.5 \
-pol i \
-name image_test_A399_20 \
-scale 3arcsec \
-niter 50000 \
-mgain 0.8 \
-fit-beam \
-multiscale \
-channels-out 6 \
-join-channels \
-multiscale-max-scales 10 \
-nmiter 10 \
-log-time \
-multiscale-scale-bias 0.7 \
-facet-regions ${TESS} \
-parallel-gridding 6 \
-fit-spectral-pol 3 \
-taper-gaussian 20arcsec \
-apply-facet-solutions ${H5// /,} amplitude000,phase000 \
${MS}