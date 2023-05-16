#!/bin/bash
#SBATCH -c 10 --constraint=amd

#MSLIST WITH PATH TO MS
MSLIST=$1

#SINGULARITY
if [[ "$HOSTNAME" == *"surfsara.nl" ]]; then
    BIND=$PWD,/home/lofarvwf-jdejong/scripts
    SIMG=/project/lofarvwf/Software/singularity/lofar_sksp_v4.1.0_znver2_znver2_noavx512_aocl3_cuda_ddf.sif
elif [[ "$HOSTNAME" == *"leidenuniv.nl" ]]; then
    BIND=$PWD,/net/rijn,/net/tussenrijn,/net/rijn1,/net/rijn2,/net/rijn3,/net/rijn4,/net/rijn5,/net/rijn6,/net/rijn7,/net/rijn8,/net/rijn9,/net/rijn10
    SIMG=/net/achterrijn/data1/sweijen/software/containers/lofar_sksp_v4.0.2_x86-64_cascadelake_cascadelake_avx512_mkl_cuda_ddf.sif
else
    echo "Host not recognized, please set paths"
    exit 0
fi

mkdir -p h5output

while read -r MS; do
#  singularity exec -B $BIND $SIMG ./phasediff_inttest.sh ${MS} 1min &
#  singularity exec -B $BIND $SIMG ./phasediff_inttest.sh ${MS} 2min &
#  singularity exec -B $BIND $SIMG ./phasediff_inttest.sh ${MS} 4min
#  singularity exec -B $BIND $SIMG ./phasediff_inttest.sh ${MS} 8min &
  singularity exec -B $BIND $SIMG ./phasediff_inttest.sh ${MS} 10min
#  singularity exec -B $BIND $SIMG ./phasediff_inttest.sh ${MS} 12min
#  singularity exec -B $BIND $SIMG ./phasediff_inttest.sh ${MS} 15min &
#  singularity exec -B $BIND $SIMG ./phasediff_inttest.sh ${MS} 20min
done <$MSLIST

singularity exec -B $BIND $SIMG python phasediff_output.py --h5 h5output/*.h5
