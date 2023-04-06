#!/bin/bash

#input MS
MS=$1

mkdir -p h5output

#PRE-AVERAGE
DP3 \
msin=${MS} \
msin.orderms=False \
msin.missingdata=True \
msin.datacolumn=DATA \
msout=${MS}.out \
msout.storagemanager=dysco \
msout.writefullresflag=False \
steps=[avg] \
avg.type=averager \
avg.freqresolution=390.56kHz \
avg.timeresolution=60


#GET SCALARPHASEDIFF SCORES
python /home/lofarvwf-jdejong/scripts/lofar_facet_selfcal/facetselfcal.py \
-i scalarphasediff_${MS} \
--forwidefield \
--phaseupstations='core' \
--msinnchan=120 \
--avgfreqstep=2 \
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
--helperscriptspath=/home/lofarvwf-jdejong/scripts/lofar_facet_selfcal \
--helperscriptspathh5merge=/home/lofarvwf-jdejong/scripts/lofar_helpers \
--stopafterskysolve \
${MS}.out


#BIG CLEAN UP
mv scalarphasediff*.h5 h5output/scalarphasediff_${MS}.h5
rm -rf *.fits
rm -rf *.p
rm -rf tmpfile
rm *.h5
rm -rf *.avg
rm -rf *.avg.phaseup
