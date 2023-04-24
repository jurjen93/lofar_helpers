#!/bin/bash
#SBATCH -c 10 --constraint=amd

#MSLIST WITH PATH TO MS
MSLIST=$1

#SINGULARITY
BIND=$PWD
SIMG=/net/achterrijn/data1/sweijen/software/containers/lofar_sksp_v4.0.2_x86-64_cascadelake_cascadelake_avx512_mkl_cuda_ddf.sif
#Spider --> /project/lofarvwf/Software/singularity/lofar_sksp_v4.0.2_znver2_znver2_noavx512_ddf_10_02_2023.sif
#Leiden --> /net/achterrijn/data1/sweijen/software/containers/lofar_sksp_v4.0.2_x86-64_cascadelake_cascadelake_avx512_mkl_cuda_ddf.sif

mkdir -p h5output

while read -r MS; do
  singularity exec -B $BIND $SIMG ./scalarphasediff.sh ${MS}
done <$MSLIST

singularity exec -B $BIND $SIMG python scalarphasediff_output.py --h5 h5output/*.h5
