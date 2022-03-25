#!/bin/bash

cd /home/jurjen/Documents/Python/lofar_helpers

F=fits/60cleanbridge_200kpc.fits #fits/60rudnick.fits
CELLSIZE=13
rm -rf ptp_dir && wait

#A399
#python analysis/ppt.py -radio1 fits/A399.fits -limits 0 0 180 180 -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -sys1 0.2 -excluderegion regions/excluderegionshalo.reg
#mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/A399 && wait
#mv ptp_dir/A399/grid_${CELLSIZE}x${CELLSIZE}_ds9_image.reg regions/gridA399_cut.reg
#python analysis/make_corr_images.py -filein ptp_dir/A399/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout analysis/A399corr.png -no_y -noisefits ${F} > a399results_cut.txt
#mv analysis/A399corr.png analysis/A399corr_cut.png
#mv ptp_dir/A399/grid_${CELLSIZE}x${CELLSIZE}_results.fits analysis/A399_results_cut.fits
#
##A401
#python analysis/ppt.py -radio1 fits/A401.fits -limits 0 0 180 180 -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -sys1 0.2 -excluderegion regions/excluderegionshalo.reg
#mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/A401 && wait
#mv ptp_dir/A401/grid_${CELLSIZE}x${CELLSIZE}_ds9_image.reg regions/gridA401_cut.reg
#python analysis/make_corr_images.py -filein ptp_dir/A401/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout analysis/A401corr.png -no_y -noisefits ${F} > a401results_cut.txt
#mv analysis/A401corr.png analysis/A401corr_cut.png
#mv ptp_dir/A401/grid_${CELLSIZE}x${CELLSIZE}_results.fits analysis/A401_results_cut.fits

#Bridge
python analysis/ppt.py -radio1 fits/Bridge.fits -limits 0 0 400 400 -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -sys1 0.2 -excluderegion regions/excluderegionsbridge.reg
mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/bridge && wait
mv ptp_dir/bridge/grid_${CELLSIZE}x${CELLSIZE}_ds9_image.reg regions/gridbridge_cut.reg
python analysis/make_corr_images.py -filein ptp_dir/bridge/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout analysis/bridgecorr.png -no_y -noisefits ${F} > bridgeresults_cut.txt
mv analysis/bridgecorr.png analysis/bridgecorr_cut.png
mv ptp_dir/bridge/grid_${CELLSIZE}x${CELLSIZE}_results.fits analysis/bridge_results_cut.fits