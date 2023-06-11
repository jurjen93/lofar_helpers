#!/bin/bash

#input MS
MS=$1

CURDIR=$PWD

re="L[0-9][0-9][0-9][0-9][0-9][0-9]"
if [[ $MS =~ $re ]]; then OBSERVATION=${BASH_REMATCH}; fi

re="P[0-9][0-9][0-9][0-9][0-9]"
if [[ $MS =~ $re ]]; then DIR=${BASH_REMATCH}; fi

mkdir -p h5output
mkdir -p ${OBSERVATION}_${DIR}
cd ${OBSERVATION}_${DIR}

if [[ -f ${MS} ]]
then
    echo "${MS} needs to have an absolute path"
    exit 0
fi

#PRE-AVERAGE
DP3 \
msin=${MS} \
msin.orderms=False \
msin.missingdata=True \
msin.datacolumn=DATA \
msout=${OBSERVATION}_${DIR}_spd.ms  \
msout.storagemanager=dysco \
msout.writefullresflag=False \
steps=[avg] \
avg.type=averager \
avg.freqresolution=390.56kHz \
avg.timeresolution=60

#GET PHASEDIFF SCORES
if [[ "$HOSTNAME" == *"surfsara.nl" ]]; then
    python /home/lofarvwf-jdejong/scripts/lofar_facet_selfcal/facetselfcal.py \
    -i phasediff_${OBSERVATION}_${DIR} \
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
    --helperscriptspath=/home/lofarvwf-jdejong/scripts/lofar_facet_selfcal \
    --helperscriptspathh5merge=/home/lofarvwf-jdejong/scripts/lofar_helpers \
    --stopafterskysolve \
    ${OBSERVATION}_${DIR}_spd.ms
elif [[ "$HOSTNAME" == *"leidenuniv.nl" ]]; then
    python /net/rijn/data2/rvweeren/LoTSS_ClusterCAL/facetselfcal.py \
    -i phasediff_${OBSERVATION}_${DIR} \
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
    ${OBSERVATION}_${DIR}_spd.ms
else
    echo "Host not recognized, please set paths"
    exit 0
fi



#BIG CLEAN UP
mv phasediff*.h5 ${CURDIR}/h5output
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