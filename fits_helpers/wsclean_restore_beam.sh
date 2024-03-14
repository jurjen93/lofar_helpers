#!/bin/bash

BIND=$PWD,/home/lofarvwf-jdejong/scripts,/project/lofarvwf/Share/jdejong,/project/lotss/Public/jdejong,/project/lofarvwf
SIMG=/project/lofarvwf/Software/singularity/flo4.6.0_znver2_znver2.sif

FACET=$1

MAJOR=0.6
MINOR=0.65

singularity exec -B $BIND $SIMG \
wsclean -beam-shape 0.42 0.35 0 -restore facet_${FACET}-MFS-residual.fits facet_${FACET}-MFS-model.fits facet_${FACET}-MFS-image.restored.fits