#!/bin/bash
#SBATCH -c 10 --constraint=amd

MSLIST=$1

export SIMG=/project/lofarvwf/Software/singularity/lofar_sksp_v4.0.2_znver2_znver2_noavx512_ddf_10_02_2023.sif

mkdir -p h5output

while read -r MS; do
  singularity exec -B $PWD,/project,/home/lofarvwf-jdejong/scripts $SIMG ./scalarphasediff.sh ${MS}
done <$MSLIST

singularity exec -B $PWD,/project,/home/lofarvwf-jdejong/scripts $SIMG python scalarphasediff_output.py --h5 h5output/*.h5
