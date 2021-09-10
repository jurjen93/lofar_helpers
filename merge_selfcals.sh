#!/bin/bash
#SBATCH -c 4
#SBATCH -p normal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jurjendejong@strw.leidenuniv.nl
#SBATCH --job-name=merge
#SBATCH --array=0-5%1

SOURCE=$1
SING_IMAGE=/home/lofarvwf-jdejong/singularities/pill-latest.simg
SING_BIND=/project/lofarvwf/Share/jdejong,/home/lofarvwf-jdejong/scripts
SCRIPT_FOLDER=/home/lofarvwf-jdejong/scripts/lofar_helpers
SELFCAL=/project/lofarvwf/Share/jdejong/output/${SOURCE}/selfcal

cd ${SELFCAL} && singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/merge_correct.py -d ${SELFCAL} -a ${SLURM_ARRAY_TASK_ID} -del 31