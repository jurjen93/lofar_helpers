#!/bin/bash

#input
H5=$1
MS=$2

#parameters
SING_BIND=/tmp,/dev/shm,/disks/paradata,/data1,/net/lofar1,/net/rijn,/net/nederrijn/,/net/bovenrijn,/net/botlek,/net/para10,/net/lofar2,/net/lofar3,/net/lofar4,/net/lofar5,/net/lofar6,/net/lofar7,/disks/ftphome,/net/krommerijn,/net/voorrijn,/net/achterrijn,/net/tussenrijn,/net/ouderijn,/net/nieuwerijn,/net/lofar8,/net/lofar9,/net/rijn8,/net/rijn7,/net/rijn5,/net/rijn4,/net/rijn3,/net/rijn2
SING_IMAGE=/net/rijn/data2/rvweeren/data/pill-latestJune2021.simg
SING_IMAGE_WSCLEAN=/net/lofar1/data1/sweijen/software/LOFAR/singularity/test/test_wsclean_facet_fix_sep30.sif
TO=/net/nieuwerijn/data2/jurjendejong/Abell399-401_60
FROM=/net/tussenrijn/data2/jurjendejong/A399_extracted_avg

#cache
singularity exec -B ${SING_BIND} ${SING_IMAGE} CleanSHM.py

#check if directory exists
if [[ -f ${TO} ]]
then
  echo "${TO} exists. Exit script"
  exit 0
fi

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
  cp ${FROM}/${M} ${TO}
  aoflagger ${TO}/${M}
done

cd ${TO}

#make facet
cp ${FROM}/tess.reg ${TO} && wait
#singularity exec -B ${SING_BIND} ${SING_IMAGE} python /net/rijn/data2/rvweeren/LoTSS_ClusterCAL/ds9facetgenerator.py \
#--h5 ${TO}/${H5} \
#--DS9regionout ${TO}/tess.reg \
#--imsize 6000 \
#--ms ${TO}/${MS_1}

#run wsclean
singularity exec -B ${SING_BIND} ${SING_IMAGE_WSCLEAN} \
wsclean \
-size 1500 1500 \
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
-name image_test_A399_60 \
-scale 6arcsec \
-niter 50000 \
-mgain 0.8 \
-fit-beam \
-multiscale \
-channels-out 6 \
-join-channels \
-multiscale-max-scales 10 \
-nmiter 11 \
-log-time \
-multiscale-scale-bias 0.7 \
-facet-regions ${TO}/tessupdate.reg \
-parallel-gridding 6 \
-fit-spectral-pol 3 \
-taper-gaussian 60arcsec \
-apply-facet-solutions ${H5} amplitude000,phase000 \
-minuv-l 80.0 \
${TO}/${MS}

lines = ['wsclean -size 1500 1500 -use-wgridder -no-update-model-required -reorder -weight briggs -0.5 -weighting-rank-filter 3 ' \
        '-clean-border 1 -parallel-reordering 6 -padding 1.2 -auto-mask 2.5 -auto-threshold 0.5 -pol i -name image_test_A399_60 ' \
        '-scale 6arcsec -niter 50000 -mgain 0.8 -fit-beam -multiscale -channels-out 6 -join-channels -multiscale-max-scales 10 -nmiter ' + nmiter + \
        ' -log-time -multiscale-scale-bias 0.7 -facet-regions '+TO+'/'+FACET+ \
        ' -minuv-l 80.0' \
        '-parallel-gridding 6 -fit-spectral-pol 3 -taper-gaussian 60arcsec ' \
        '-apply-facet-solutions '+TO+'/'+H5+' amplitude000,phase000 '+ TO + \
        '/'+MS]