#!/bin/bash
#SBATCH -c 20
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jurjendejong@strw.leidenuniv.nl
#SBATCH --job-name=array
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err
#SBATCH --array=1-10

SING_IMAGE=/home/lofarvwf-jdejong/singularities/lofar_sksp_fedora31_ddf.sif
SING_BIND=/project/lofarvwf/Share/jdejong,/home/lofarvwf-jdejong/scripts
SCRIPT_FOLDER=/home/lofarvwf-jdejong/scripts/lofar_helpers
TO=/project/lofarvwf/Share/jdejong/output/A399/selfcal

mkdir -p ${TO} && mkdir -p ${TO}/finished

START="$(date -u +%s)"
singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/pipeline_scripts/surf/selfcal_A399.py --box ${SLURM_ARRAY_TASK_ID}
END="$(date -u +%s)"
echo "Selfcal in $((${END}-${START})) seconds" > ${TO}/finished/box_${SLURM_ARRAY_TASK_ID}.txt