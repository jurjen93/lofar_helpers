#!/bin/bash

BIND=$PWD,/home/lofarvwf-jdejong/scripts,/project/lofarvwf/Share/jdejong,/project/lotss/Public/jdejong,/project/lofarvwf
SIMG=/project/lofarvwf/Software/singularity/flocs_v5.0.0_znver2_znver2_aocl4_cuda.sif

FACET=$1

MAJOR=1.51
MINOR=1.00
PA=108.85

singularity exec -B $BIND $SIMG \
wsclean -beam-shape $MAJOR $MINOR $PA -restore facet_${FACET}-MFS-residual.fits facet_${FACET}-MFS-model.fits facet_${FACET}-MFS-image.restored.fits
