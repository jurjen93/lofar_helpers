#!/bin/bash
#SBATCH -c 24
#SBATCH -p normal
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jurjendejong@strw.leidenuniv.nl
#SBATCH --job-name=selfcal
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --constraint=intel

BOX=$1
SING_IMAGE=/home/lofarvwf-jdejong/singularities/lofar_sksp_fedora31_ddf_fixed.sif
SING_BIND=/project/lofarvwf/Share/jdejong,/home/lofarvwf-jdejong/scripts
SCRIPT_FOLDER=/home/lofarvwf-jdejong/scripts/lofar_helpers
TO=/project/lofarvwf/Share/jdejong/output/A399/selfcal

# get CPU ID currently used by self
echo "CPU ID:"
CPU_ID=$(cat /proc/self/stat)
echo $CPU_ID | gawk '{print $39}'
# find all CPU ID's that are allocated (reserved) by self
echo "CPU ALLOC:"
CPU_ALLOC=$(cat /proc/self/status | grep 'Cpus_allowed_list')
echo $CPU_ALLOC

singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/pipeline_scripts/surf/prepare_data_A399.py --box ${BOX}
cd /project/lofarvwf/Share/jdejong/output/A399/selfcal/box_${BOX}
singularity exec -B ${SING_BIND} ${SING_IMAGE} python /home/lofarvwf-jdejong/scripts/runwscleanLBautoR.py -b /project/lofarvwf/Share/jdejong/output/A399/boxes/box_${BOX}.reg --auto --imager=DDFACET --helperscriptspath=/home/lofarvwf-jdejong/scripts --autofrequencyaverage-calspeedup --useaoflagger --uvmin=750 --tecfactorsolint=1.5 --gainfactorsolint=2.0 *box_${BOX}.dysco.sub.shift.avg.weights.ms.archive*
echo "Selfcal box_${BOX} finished" > ${TO}/finished/box_${BOX}.txt