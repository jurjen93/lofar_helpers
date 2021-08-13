#!/bin/bash
#SBATCH -c 20
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jurjendejong@strw.leidenuniv.nl
#SBATCH --job-name=array
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err
#SBATCH --array=11-82

SOURCE=$1
SING_IMAGE=/home/lofarvwf-jdejong/singularities/lofar_sksp_fedora31_ddf.sif
SING_BIND=/project/lofarvwf/Share/jdejong,/home/lofarvwf-jdejong/scripts
SCRIPT_FOLDER=/home/lofarvwf-jdejong/scripts/lofar_helpers
TO=/project/lofarvwf/Share/jdejong/output/${SOURCE}/selfcal

# get CPU ID currently used by self
echo "CPU ID:"
CPU_ID=$(cat /proc/self/stat)
echo $CPU_ID | gawk '{print $39}'
# find all CPU ID's that are allocated (reserved) by self
echo "CPU ALLOC:"
CPU_ALLOC=$(cat /proc/self/status | grep 'Cpus_allowed_list')
echo $CPU_ALLOC

mkdir -p ${TO}/box_${SLURM_ARRAY_TASK_ID} && mkdir -p ${TO}/finished

START="$(date -u +%s)"
singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/pipeline_scripts/surf/prepare_data_${SOURCE}.py --box ${SLURM_ARRAY_TASK_ID}
wait
cd /project/lofarvwf/Share/jdejong/output/${SOURCE}/selfcal/box_${SLURM_ARRAY_TASK_ID}
singularity exec -B ${SING_BIND} ${SING_IMAGE} python /home/lofarvwf-jdejong/scripts/runwscleanLBautoR.py -b /project/lofarvwf/Share/jdejong/output/${SOURCE}/boxes/box_${SLURM_ARRAY_TASK_ID}.reg --auto --imager=DDFACET --helperscriptspath=/home/lofarvwf-jdejong/scripts --autofrequencyaverage-calspeedup='True' *box_${SLURM_ARRAY_TASK_ID}.dysco.sub.shift.avg.weights.ms.archive*
END="$(date -u +%s)"
echo "Selfcal in $((${END}-${START})) seconds" > ${TO}/finished/box_${SLURM_ARRAY_TASK_ID}.txt