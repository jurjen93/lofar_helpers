#!/bin/bash
#SBATCH -c 10 --constraint=amd

#SINGULARITY
BIND=$( python3 $HOME/parse_settings.py --BIND )
SIMG=$( python3 $HOME/parse_settings.py --SIMG )
FACETSELFCAL=$( python3 $HOME/parse_settings.py --facet_selfcal )

#input MS
DIR=$1

CURDIR=$PWD

mkdir -p ${DIR}
cd ${DIR}

#cp -r /project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/ddcal/alldirs/*${DIR}*.ms .
#
##PRE-AVERAGE
#singularity exec -B $BIND $SIMG DP3 \
#msin=L686962_${DIR}.ms \
#msin.orderms=False \
#msin.missingdata=True \
#msin.datacolumn=DATA \
#msout=L686962_${DIR}_spd.ms  \
#msout.storagemanager=dysco \
#msout.writefullresflag=False \
#steps=[avg] \
#avg.type=averager \
#avg.freqresolution=390.56kHz \
#avg.timeresolution=60
#
##PRE-AVERAGE
#singularity exec -B $BIND $SIMG DP3 \
#msin=L769393_${DIR}.ms \
#msin.orderms=False \
#msin.missingdata=True \
#msin.datacolumn=DATA \
#msout=L769393_${DIR}_spd.ms  \
#msout.storagemanager=dysco \
#msout.writefullresflag=False \
#steps=[avg] \
#avg.type=averager \
#avg.freqresolution=390.56kHz \
#avg.timeresolution=60
#
##PRE-AVERAGE
#singularity exec -B $BIND $SIMG DP3 \
#msin=L798074_${DIR}.ms \
#msin.orderms=False \
#msin.missingdata=True \
#msin.datacolumn=DATA \
#msout=L798074_${DIR}_spd.ms  \
#msout.storagemanager=dysco \
#msout.writefullresflag=False \
#steps=[avg] \
#avg.type=averager \
#avg.freqresolution=390.56kHz \
#avg.timeresolution=60
#
##PRE-AVERAGE
#singularity exec -B $BIND $SIMG DP3 \
#msin=L816272_${DIR}.ms \
#msin.orderms=False \
#msin.missingdata=True \
#msin.datacolumn=DATA \
#msout=L816272_${DIR}_spd.ms  \
#msout.storagemanager=dysco \
#msout.writefullresflag=False \
#steps=[avg] \
#avg.type=averager \
#avg.freqresolution=390.56kHz \
#avg.timeresolution=60
#
#rm -rf L??????_${DIR}.ms

#GET PHASEDIFF SCORES
singularity exec -B $BIND $SIMG python ${FACETSELFCAL} \
-i phasediff_${DIR} \
--forwidefield \
--phaseupstations='core' \
--skipbackup \
--uvmin=20000 \
--soltype-list="['scalarphasediff']" \
--solint-list="['10min']" \
--nchan-list="[6]" \
--docircular \
--uvminscalarphasediff=0 \
--stop=1 \
--soltypecycles-list="[0]" \
--imsize=1600 \
--skymodelpointsource=1.0 \
--stopafterskysolve \
--helperscriptspath=/home/lofarvwf-jdejong/scripts/lofar_facet_selfcal \
--helperscriptspathh5merge=/home/lofarvwf-jdejong/scripts/lofar_helpers \
*${DIR}_spd.ms

#BIG CLEAN UP
rm -rf *.fits
rm -rf *.p
rm -rf tmpfile
#rm *.h5
rm -rf *.avg
rm -rf *.avg.phaseup
rm -rf *.parset
rm BLsmooth.py
rm lin2circ.py
rm lib_multiproc.py
rm h5_merger.py
rm polconv.py
rm plot_tecandphase.py
rm VLASS_dyn_summary.php
rm vlass_search.py
rm -rf __pycache__
rm *.log
rm merged*.h5
rm *.scbackup
rm *templatejones.h5
rm *.png

cd ../