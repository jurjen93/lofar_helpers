#!/bin/bash

cd /home/jurjen/Documents/Python/lofar_helpers

OUTFILE=ptp_results_bridgebubble

mkdir -p ${OUTFILE}

for N in `seq 0 4 54`
do
    rm -rf ptp_dir && wait

    F=fits/60rudnick.fits
    CELLSIZE=$((27+${N}))

    python analysis/ppt.py -radio1 ${F} -limits $(($((1+$RANDOM % 10))+1190)) $(($((1+$RANDOM % 10))+1263)) $(($((1+$RANDOM % 10))+1672)) $(($((1+$RANDOM % 10))+1800)) -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -sys1 0.1 -excluderegion regions/excluderegionsbridge.reg
    mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE}/grid_${CELLSIZE}x${CELLSIZE}_ds9_image.reg ${OUTFILE}/gridbridge_rudnick_${CELLSIZE}.reg
    mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE}/grid_${CELLSIZE}x${CELLSIZE}_results.fits ${OUTFILE}/bridge_results_rudnick_${CELLSIZE}.fits
    python analysis/make_corr_images.py -filein ${OUTFILE}/bridge_results_rudnick_${CELLSIZE}.fits -fileout ${OUTFILE}/bridgecorr_rudnick_${CELLSIZE}.png -no_y -noisefits ${F} > ${OUTFILE}/bridgeresults_rudnick_${CELLSIZE}.txt
    echo $N

    rm -rf ptp_dir && wait

    F=fits/60cleanbridge_300kpc.fits
    CELLSIZESUB=$((13+${N}/2))

    python analysis/ppt.py -radio1 ${F} -limits $(($((1+$RANDOM % 10))+595)) $(($((1+$RANDOM % 10))+631)) $(($((1+$RANDOM % 10))+836)) $(($((1+$RANDOM % 10))+900)) -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZESUB} -sys1 0.1 -excluderegion regions/excluderegionsbridge.reg
    mv ptp_dir/grid_${CELLSIZESUB}x${CELLSIZESUB}/grid_${CELLSIZESUB}x${CELLSIZESUB}_ds9_image.reg ${OUTFILE}/gridbridge_cb_${CELLSIZE}.reg
    mv ptp_dir/grid_${CELLSIZESUB}x${CELLSIZESUB}/grid_${CELLSIZESUB}x${CELLSIZESUB}_results.fits ${OUTFILE}/bridge_results_cb_${CELLSIZE}.fits
    python analysis/make_corr_images.py -filein ${OUTFILE}/bridge_results_cb_${CELLSIZE}.fits -fileout ${OUTFILE}/bridgecorr_cb_${CELLSIZE}.png -no_y -noisefits ${F} > ${OUTFILE}/bridgeresults_cb_${CELLSIZE}.txt
    echo $CELLSIZE
done