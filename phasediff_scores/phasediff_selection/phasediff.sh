#!/bin/bash

#input MS
MS=$1

#SCRIPT PATHS
FACETSELFCAL=$( python3 $HOME/parse_settings.py --facet_selfcal )
LOFARFACETSELFCAL=$( python3 $HOME/parse_settings.py --lofar_facet_selfcal )
LOFARHELPERS=$( python3 $HOME/parse_settings.py --lofar_helpers )

IFS='/' read -ra MSS <<< "$MS"
MSOUT=spd_${MSS[-1]}

CURDIR=$PWD

if [[ -f ${MS} ]]
then
  echo "ERROR: ${MS} NEEDS TO HAVE ABSOLUTE PATHS"
  exit 0
fi

#PRE-AVERAGE
DP3 \
msin=${MS} \
msin.orderms=False \
msin.missingdata=True \
msin.datacolumn=DATA \
msout=${MSOUT} \
msout.storagemanager=dysco \
msout.writefullresflag=False \
steps=[avg] \
avg.type=averager \
avg.freqresolution=390.56kHz \
avg.timeresolution=60

#GET PHASEDIFF H5
python $FACETSELFCAL \
-i phasediff \
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
--helperscriptspath=$LOFARFACETSELFCAL \
--helperscriptspathh5merge=$LOFARHELPERS \
--stopafterskysolve \
${MSOUT}

#BIG CLEAN UP
#mv *phasediff*.h5 ${CURDIR}/h5output
rm -rf *.fits
rm -rf *.p
rm -rf tmpfile
#rm *.h5
rm -rf *.avg
rm -rf *.phaseup
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