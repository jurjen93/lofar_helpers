#!/bin/bash

cd /home/jurjen/Documents/Python/lofar_helpers

F=fits/60rudnick.fits #fits/60rudnick.fits
CELLSIZE=25



F=fits/60cleanbridge_200kpc.fits
CELLSIZE=13

#A399

OUTFILE=ptp_results_a399_trail
mkdir ${OUTFILE}

rm -rf ptp_dir && wait
python analysis/ppt.py -radio1 ${F} -limits $(($((1+$RANDOM % 10))+725)) $(($((1+$RANDOM % 10))+520)) $(($((1+$RANDOM % 10))+945)) $(($((1+$RANDOM % 10))+740)) -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -sys1 0.2 -excluderegion regions/excluderegionshalo_trail.reg
mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/A399_trail && wait
mv ptp_dir/A399_trail/grid_${CELLSIZE}x${CELLSIZE}_ds9_image.reg ${OUTFILE}/gridA399_cb_${N}.reg
python analysis/make_corr_images.py -filein ptp_dir/A399_trail/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout ${OUTFILE}/A399corr_cb_${N}.png -no_y -noisefits ${F} > ${OUTFILE}/a399results_cb_${N}.txt
mv ptp_dir/A399_trail/grid_${CELLSIZE}x${CELLSIZE}_results.fits ${OUTFILE}/A399_results_cb_${N}.fits

OUTFILE=ptp_results_a399_extended
mkdir ${OUTFILE}

rm -rf ptp_dir && wait
python analysis/ppt.py -radio1 ${F} -limits $(($((1+$RANDOM % 10))+725)) $(($((1+$RANDOM % 10))+520)) $(($((1+$RANDOM % 10))+945)) $(($((1+$RANDOM % 10))+740)) -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -sys1 0.2 -excluderegion regions/excluderegionshalo_extended.reg
mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/A399_extended && wait
mv ptp_dir/A399_extended/grid_${CELLSIZE}x${CELLSIZE}_ds9_image.reg ${OUTFILE}/gridA399_cb_${N}.reg
python analysis/make_corr_images.py -filein ptp_dir/A399_extended/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout ${OUTFILE}/A399corr_cb_${N}.png -no_y -noisefits ${F} > ${OUTFILE}/a399results_cb_${N}.txt
mv ptp_dir/A399_extended/grid_${CELLSIZE}x${CELLSIZE}_results.fits ${OUTFILE}/A399_results_cb_${N}.fits

OUTFILE=ptp_results_a399
mkdir ${OUTFILE}

rm -rf ptp_dir && wait
python analysis/ppt.py -radio1 ${F} -limits $(($((1+$RANDOM % 10))+725)) $(($((1+$RANDOM % 10))+520)) $(($((1+$RANDOM % 10))+945)) $(($((1+$RANDOM % 10))+740)) -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -sys1 0.2 -excluderegion regions/excluderegionshalo.reg
mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/A399 && wait
mv ptp_dir/A399/grid_${CELLSIZE}x${CELLSIZE}_ds9_image.reg ${OUTFILE}/gridA399_cb_${N}.reg
python analysis/make_corr_images.py -filein ptp_dir/A399/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout ${OUTFILE}/A399corr_cb_${N}.png -no_y -noisefits ${F} > ${OUTFILE}/a399results_cb_${N}.txt
mv ptp_dir/A399/grid_${CELLSIZE}x${CELLSIZE}_results.fits ${OUTFILE}/A399_results_cb_${N}.fits