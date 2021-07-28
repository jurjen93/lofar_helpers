#!/bin/bash
#SBATCH -c 20
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jurjendejong@strw.leidenuniv.nl 

FIELD=$1 #L626678
TO=/project/lofarvwf/Share/jdejong/output/${FIELD}
SCRIPT_FOLDER=/home/lofarvwf-jdejong/scripts
N=$2 #box number

SING_IMAGE=/home/lofarvwf-jdejong/singularities/lofar_sksp_fedora31_ddf.sif
SING_BIND=/project/lofarvwf/Share/jdejong,/home/lofarvwf-jdejong/scripts

echo "-----START SELFCAL-----"
echo "Started selfcal for box_${N}"
cd ${TO}/selfcal/
mkdir ${TO}/selfcal/box_${N}
cp -r ${TO}/extract/*box_${N}.dysco.sub.shift.avg.weights.ms.archive0 ${TO}/selfcal/
singularity exec -B ${SING_BIND} ${SING_IMAGE} DPPP msin=${TO}/selfcal/Abell399-401_box_${N}.dysco.sub.shift.avg.weights.ms.archive0 msout.storagemanager=dysco msout=${TO}/selfcal/box_${N}/box_${N}.dysco.sub.shift.avg.weights.ms.archive0.goodtimes msin.ntimes=1500 steps=[]
rm -rf ${TO}/selfcal/Abell399-401_box_${N}.dysco.sub.shift.avg.weights.ms.archive0
cd ${TO}/selfcal/box_${N}
singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/runwscleanLBautoR.py -b ${TO}/boxes/box_${N}.reg --auto --imager=DDFACET --helperscriptspath=${SCRIPT_FOLDER}/ --autofrequencyaverage-calspeedup='True' box_${N}.dysco.sub.shift.avg.weights.ms.archive0.goodtimes
echo "Finished selfcal for box_${N}"
echo "-----END SELFCAL-----"