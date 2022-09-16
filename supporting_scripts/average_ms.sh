#!/bin/bash
#SBATCH -c 1

SING_IMAGE=/home/lofarvwf-jdejong/singularities/pill-latest.simg
SING_BIND=/project/lofarvwf/Share/jdejong

MS_IN=$1

singularity exec -B ${SING_BIND} ${SING_IMAGE} DPPP msin=${MS_IN} msout=avg_${MS_IN} steps=[av] msout.storagemanager=dysco steps=[av] av.type=averager av.freqstep=16 av.timestep=16