#!/usr/bin/env bash

# ---- Input ---- #

MS=$1
CHANNELS=$2
TIME_INTERVALS=$3

# ---- Run WSClean ---- #

echo "Running WSClean with:"
echo "  MS          = $MS"
echo "  Channels    = $CHANNELS"
echo "  Intervals   = $INTERVALS"
echo "------------------------------------------------------------"

wsclean \
-no-update-model-required \
-minuv-l 1500.0 \
-size 1024 1024 \
-reorder \
-weight briggs -1.5 \
-parallel-reordering 4 \
-mgain 0.6 \
-data-column DATA \
-channels-out ${CHANNELS} \
-fit-spectral-pol 3 \
-join-channels \
-auto-mask 2.5 \
-auto-threshold 0.5 \
-pol i \
-gridder wgridder \
-wgridder-accuracy 0.0001 \
-use-differential-lofar-beam \
-facet-beam-update 120 \
-name transient_imaging \
-scale 0.075arcsec \
-nmiter 9 \
-niter 15000 \
-intervals-out ${TIME_INTERVALS} \
$MS

echo "------------------------------------------------------------"
echo "Imaging complete."