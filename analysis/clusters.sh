#!/bin/bash

cd /home/jurjen/Documents/Python/lofar_helpers

F=fits/60rudnick.fits #fits/60rudnick.fits
CELLSIZE=35
rm -rf ptp_dir && wait

#A399
python analysis/ppt.py -radio1 ${F} -limits 1450 1040 1890 1480 -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -sys1 0.2 -excluderegion regions/excluderegionshalo.reg
mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/A399 && wait
mv ptp_dir/A399/grid_${CELLSIZE}x${CELLSIZE}_ds9_image.reg regions/gridA399_rudnick.reg
python analysis/make_corr_images.py -filein ptp_dir/A399/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout analysis/A399corr.png -no_y > a399results_rudnick.txt
mv analysis/A399corr.png analysis/A399corr_rudnick.png
mv ptp_dir/A399/grid_${CELLSIZE}x${CELLSIZE}_results.fits analysis/A399_results_rudnick.fits

#A401
python analysis/ppt.py -radio1 ${F} -limits 1200 1700 1430 2050 -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -sys1 0.2 -excluderegion regions/excluderegionshalo.reg
mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/A401 && wait
mv ptp_dir/A401/grid_${CELLSIZE}x${CELLSIZE}_ds9_image.reg regions/gridA401_rudnick.reg
python analysis/make_corr_images.py -filein ptp_dir/A401/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout analysis/A401corr.png -no_y > a401results_rudnick.txt
mv analysis/A401corr.png analysis/A401corr_rudnick.png
mv ptp_dir/A401/grid_${CELLSIZE}x${CELLSIZE}_results.fits analysis/A401_results_rudnick.fits

#Bridge
python analysis/ppt.py -radio1 ${F} -limits 1190 1263 1672 1800 -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -sys1 0.2 -excluderegion regions/excluderegionsbridge.reg
mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/bridge && wait
mv ptp_dir/bridge/grid_${CELLSIZE}x${CELLSIZE}_ds9_image.reg regions/gridbridge_rudnick.reg
python analysis/make_corr_images.py -filein ptp_dir/bridge/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout analysis/bridgecorr.png -no_y > bridgeresults_rudnick.txt
mv analysis/bridgecorr.png analysis/bridgecorr_rudnick.png
mv ptp_dir/bridge/grid_${CELLSIZE}x${CELLSIZE}_results.fits analysis/bridge_results_rudnick.fits

F=fits/60cleanbridge_200kpc.fits
CELLSIZE=17
rm -rf ptp_dir && wait

#A399
python analysis/ppt.py  -radio1 ${F} -limits 725 520 945 740 -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -sys1 0.2 -excluderegion regions/excluderegionshalo.reg
mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/A399 && wait
mv ptp_dir/A399/grid_${CELLSIZE}x${CELLSIZE}_ds9_image.reg regions/gridA399_cb.reg
python analysis/make_corr_images.py -filein ptp_dir/A399/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout analysis/A399corr.png -no_y > a399results_cb.txt
mv analysis/A399.png analysis/A399corr_cb.png
mv ptp_dir/A399/grid_${CELLSIZE}x${CELLSIZE}_results.fits analysis/A399_results_cb.fits

#A401
python analysis/ppt.py  -radio1 ${F} -limits 600 850 715 1025 -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -sys1 0.2 -excluderegion regions/excluderegionshalo.reg
mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/A401 && wait
mv ptp_dir/A401/grid_${CELLSIZE}x${CELLSIZE}_ds9_image.reg regions/gridA401_cb.reg
python analysis/make_corr_images.py -filein ptp_dir/A401/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout analysis/A401corr.png -no_y > a401results_cb.txt
mv analysis/A401.png analysis/A401corr_cb.png
mv ptp_dir/A401/grid_${CELLSIZE}x${CELLSIZE}_results.fits analysis/A401_results_cb.fits

#Bridge
python analysis/ppt.py  -radio1 ${F} -limits 595 631 836 900 -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -sys1 0.2 -excluderegion regions/excluderegionsbridge.reg
mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/bridge && wait
mv ptp_dir/bridge/grid_${CELLSIZE}x${CELLSIZE}_ds9_image.reg regions/gridbridge_cb.reg
python analysis/make_corr_images.py -filein ptp_dir/bridge/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout analysis/bridgecorr.png -no_y > bridgeresults_cb.txt
mv analysis/bridgecorr.png analysis/bridgecorr_cb.png
mv ptp_dir/bridge/grid_${CELLSIZE}x${CELLSIZE}_results.fits analysis/bridge_results_cb.fits
