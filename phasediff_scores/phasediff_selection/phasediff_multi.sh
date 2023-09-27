#!/bin/bash

#SCRIPT PATHS
FACETSELFCAL=$( python3 $HOME/parse_settings.py --facet_selfcal )
LOFARFACETSELFCAL=$( python3 $HOME/parse_settings.py --lofar_facet_selfcal )
LOFARHELPERS=$( python3 $HOME/parse_settings.py --lofar_helpers )

#TODO: NO PRE-AVERAGE

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
*.ms

cd ../