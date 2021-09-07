#!/bin/bash
#SBATCH -c 4
#SBATCH -p normal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jurjendejong@strw.leidenuniv.nl
#SBATCH --job-name=merge
#SBATCH --array=0-5

SOURCE=$1
SING_IMAGE=/home/lofarvwf-jdejong/singularities/lofar_sksp_fedora31_ddf_fixed.sif
SING_BIND=/project/lofarvwf/Share/jdejong,/home/lofarvwf-jdejong/scripts
SCRIPT_FOLDER=/home/lofarvwf-jdejong/scripts/lofar_helpers
SELFCAL=/project/lofarvwf/Share/jdejong/output/${SOURCE}/selfcal

cd ${SELFCAL} && singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/merge_selfcals.py -d ${SELFCAL} -a ${SLURM_ARRAY_TASK_ID}