#!/bin/bash
#SBATCH -c 20
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jurjendejong@strw.leidenuniv.nl

BOX=$1
SING_IMAGE=/home/lofarvwf-jdejong/singularities/lofar_sksp_fedora31_ddf.sif
SING_BIND=/project/lofarvwf/Share/jdejong,/home/lofarvwf-jdejong/scripts
SCRIPT_FOLDER=/home/lofarvwf-jdejong/scripts/lofar_helpers
TO=/project/lofarvwf/Share/jdejong/output/A399/selfcal

START="$(date -u +%s)"
singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/pipeline_scripts/surf/selfcal_A399.py --box ${BOX}
END="$(date -u +%s)"
echo "Selfcal in $((${END}-${START})) seconds" > ${TO}/finished/box_${BOX}.txt