#!/bin/bash

#input MS
MS=$1

IFS='/' read -ra MSS <<< "$MS"
MSOUT=spd_${MSS[-1]}

CURDIR=$PWD

#re="L[0-9][0-9][0-9][0-9][0-9][0-9]"
#if [[ $MS =~ $re ]]; then OBSERVATION=${BASH_REMATCH}; fi
#
#re="P[0-9][0-9][0-9][0-9][0-9]"
#if [[ $MS =~ $re ]]; then DIR=${BASH_REMATCH}; fi

#mkdir -p h5output
#mkdir -p ${OBSERVATION}_${DIR}
#cd ${OBSERVATION}_${DIR}


if [[ -f ${MS} ]]
then
  echo "ERROR: ${MS} NEEDS TO HAVE ABSOLUTE PATHS"
  exit 0
fi

FACETSELFCAL=/home/lofarvwf-jdejong/scripts/lofar_facet_selfcal/facetselfcal.py
LOFARFACETSELFCAL=/home/lofarvwf-jdejong/scripts/lofar_facet_selfcal
LOFARHELPERS=/home/lofarvwf-jdejong/scripts/lofar_helpers

#if [[ "$HOSTNAME" == *"surfsara.nl" ]]; then
#  FACETSELFCAL=/home/lofarvwf-jdejong/scripts/lofar_facet_selfcal/facetselfcal.py
#  LOFARFACETSELFCAL=/home/lofarvwf-jdejong/scripts/lofar_facet_selfcal
#  LOFARHELPERS=/home/lofarvwf-jdejong/scripts/lofar_helpers
#elif [[ "$HOSTNAME" == *"leidenuniv.nl" ]]; then
#  FACETSELFCAL=/net/rijn/data2/rvweeren/LoTSS_ClusterCAL/facetselfcal.py
#  LOFARFACETSELFCAL=/net/rijn/data2/rvweeren/LoTSS_ClusterCAL
#  LOFARHELPERS=/net/rijn/data2/rvweeren/LoTSS_ClusterCAL
#elif [ -f "facetselfcal.py" ]; then
#  FACETSELFCAL=facetselfcal.py
#  LOFARFACETSELFCAL=.
#  LOFARHELPERS=.
#elif : >/dev/tcp/8.8.8.8/53; then
#  wget https://raw.githubusercontent.com/rvweeren/lofar_facet_selfcal/main/facetselfcal.py
#  wget https://raw.githubusercontent.com/rvweeren/lofar_facet_selfcal/main/lib_multiproc.py
#  wget https://raw.githubusercontent.com/rvweeren/lofar_facet_selfcal/main/plot_tecandphase.py
#  wget https://raw.githubusercontent.com/rvweeren/lofar_facet_selfcal/main/lin2circ.py
#  wget https://raw.githubusercontent.com/rvweeren/lofar_facet_selfcal/main/BLsmooth.py
#  wget https://raw.githubusercontent.com/rvweeren/lofar_facet_selfcal/main/polconv.py
#  wget https://raw.githubusercontent.com/rvweeren/lofar_facet_selfcal/main/vlass_search.py
#  wget https://raw.githubusercontent.com/rvweeren/lofar_facet_selfcal/main/VLASS_dyn_summary.php
#  wget https://raw.githubusercontent.com/rvweeren/lofar_facet_selfcal/main/find_solint.py
#  wget https://raw.githubusercontent.com/rvweeren/lofar_facet_selfcal/main/ds9facetgenerator.py
#  wget https://raw.githubusercontent.com/rvweeren/lofar_facet_selfcal/main/default_StokesV.lua
#  wget https://raw.githubusercontent.com/jurjen93/lofar_helpers/master/h5_merger.py
#  FACETSELFCAL=facetselfcal.py
#  LOFARFACETSELFCAL=.
#  LOFARHELPERS=.
#else
#  # FACETSELFCAL=<SET_PATH_TO_FACETSELFCAL>
#  # LOFARFACETSELFCAL=<SET_PATH_TO_FACETSELFCAL_FOLDER>
#  # LOFARHELPERS=<SET_PATH_TO_LOFAR_HELPERS_FOLDER>
#  echo "ERROR: NO INTERNET CONNECTION AND NO SCRIPTS IN FOLDER"
#  exit 0
#fi

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