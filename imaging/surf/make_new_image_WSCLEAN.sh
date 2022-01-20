#!/bin/bash
#SBATCH --partition=infinite
#SBATCH -c 30
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jurjendejong@strw.leidenuniv.nl


#EXAMPLE ~/scripts/lofar_helpers/imaging/make_new_image_WSCLEAN.sh \
# 'all_directions0.h5 all_directions1.h5' \
# 'Abell399-401_extr.dysco.sub.shift.avg.weights.ms.archive0.avg.goodtimes Abell399-401_extr.dysco.sub.shift.avg.weights.ms.archive1.avg.goodtimes'

#input
H5=$1
MS=$2

SING_BIND=/project/lofarvwf/Share/jdejong,/home/lofarvwf-jdejong/scripts
SING_IMAGE=/home/lofarvwf-jdejong/singularities/pill-latest.simg
SING_IMAGE_WSCLEAN=/home/lofarvwf-jdejong/singularities/idgtest.sif

TO=/project/lofarvwf/Share/jdejong/output/A399/imaging/Abell399-401
FROM=/project/lofarvwf/Share/jdejong/output/A399/imaging/jurjendejong/A399_extracted_avg

echo "Running on ${HOSTNAME%%.*}"

#check if directory exists
if [[ -f ${TO} ]]
then
  echo "${TO} exists. Make Abell399-401_new"
  TO=/net/${HOSTNAME%%.*}/data2/jurjendejong/Abell399-401_new
fi

#cache
singularity exec -B ${SING_BIND} ${SING_IMAGE} CleanSHM.py

#make directory
mkdir -p ${TO}
cd ${FROM}

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
cp ${FROM}/tessupdate.reg ${TO} && wait
#singularity exec -B ${SING_BIND} ${SING_IMAGE} python /net/rijn/data2/rvweeren/LoTSS_ClusterCAL/ds9facetgenerator.py \
#--h5 ${TO}/${H5} \
#--DS9regionout ${TO}/tessupdate.reg \
#--imsize 6000 \
#--ms ${TO}/${MS}

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
-facet-regions tessupdate.reg \
-apply-facet-solutions ${H5// /,} amplitude000,phase000 \
-name image_test_A399 \
-size 6000 6000 \
-scale 1.5arcsec \
-nmiter 10 \
${MS}