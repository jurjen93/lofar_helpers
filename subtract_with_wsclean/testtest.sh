#!/bin/bash
#SBATCH -c 3
#SBATCH --job-name=subtract
#SBATCH --array=0-24

#SINGULARITY SETTINGS
SING_BIND=$( python3 $HOME/parse_settings.py --BIND )
SIMG=$( python3 $HOME/parse_settings.py --SIMG )

pattern="*.ms"
MS_FILES=( $pattern )
MS=${MS_FILES[${SLURM_ARRAY_TASK_ID}]}

singularity exec -B $BIND $SIMG DP3 \
msin.missingdata=True \
msin.datacolumn=SUBTRACT_DATA \
msin.orderms=False \
msout.storagemanager=dysco \
msout.writefullreslag=False \
ps.type=phaseshifter \
ps.phasecenter=[243.58963deg,55.61411deg] \
beam.type=applybeam \
beam.direction=[] \
beam.updateweights=True \
ac.type=applycal \
ac.parmdb=/project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/imaging/DD_1.2/test_facet/merged_L686962.h5 \
ac.correction=fulljones \
ac.soltab=[amplitude000,phase000] \
ac.direction=[243.58963deg,55.61411deg] \
avg.type=averager \
avg.freqstep=4 \
avg.timestep=4 \
steps=[ps,beam,ac,avg] \
msin=${MS} \
msout=sub1.2asec_avg_applycal_sub6asec_L686962_SB001_uv_12CFFDDA9t_164MHz.pre-cal.ms.sub.shift.avg.ms