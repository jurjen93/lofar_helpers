#!/bin/bash
#SBATCH -c 24
#SBATCH -p normal
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jurjendejong@strw.leidenuniv.nl
#SBATCH --job-name=selfcal
#SBATCH --array=1-82

SOURCE=$1
SING_IMAGE=/home/lofarvwf-jdejong/singularities/lofar_sksp_fedora31_ddf_fixed.sif
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

singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/pipeline_scripts/surf/prepare_data_${SOURCE}.py --box ${SLURM_ARRAY_TASK_ID}
wait
cd /project/lofarvwf/Share/jdejong/output/${SOURCE}/selfcal/box_${SLURM_ARRAY_TASK_ID}
singularity exec -B ${SING_BIND} ${SING_IMAGE} python /home/lofarvwf-jdejong/scripts/runwscleanLBautoR.py -b /project/lofarvwf/Share/jdejong/output/${SOURCE}/boxes/box_${SLURM_ARRAY_TASK_ID}.reg --auto --imager=DDFACET --helperscriptspath=/home/lofarvwf-jdejong/scripts --autofrequencyaverage-calspeedup --useaoflagger --uvmin=750 --tecfactorsolint=1.5 --gainfactorsolint=2.0 *box_${SLURM_ARRAY_TASK_ID}.dysco.sub.shift.avg.weights.ms.archive*
echo "Selfcal box_${SLURM_ARRAY_TASK_ID} finished" > ${TO}/finished/box_${SLURM_ARRAY_TASK_ID}.txt