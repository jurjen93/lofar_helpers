#!/bin/bash

cd /home/jurjen/Documents/Python/lofar_helpers

F=fits/60rudnick.fits #fits/60rudnick.fits
CELLSIZE=25
rm -rf ptp_dir && wait

#A399
python analysis/ppt.py -radio1 ${F} -limits 1450 1040 1890 1480 -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -sys1 0.2 -excluderegion regions/excluderegionshalo.reg
mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/A399 && wait
mv ptp_dir/A399/grid_${CELLSIZE}x${CELLSIZE}_ds9_image.reg ptp_results/gridA399_rudnick.reg
python analysis/make_corr_images.py -filein ptp_dir/A399/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout ptp_results/A399corr_rudnick.png -no_y -noisefits ${F} > ptp_results/a399results_rudnick.txt
mv ptp_dir/A399/grid_${CELLSIZE}x${CELLSIZE}_results.fits ptp_results/A399_results_rudnick.fits

#A401
python analysis/ppt.py -radio1 ${F} -limits 1200 1700 1430 2050 -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -sys1 0.2 -excluderegion regions/excluderegionshalo.reg
mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/A401 && wait
mv ptp_dir/A401/grid_${CELLSIZE}x${CELLSIZE}_ds9_image.reg ptp_results/gridA401_rudnick.reg
python analysis/make_corr_images.py -filein ptp_dir/A401/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout ptp_results/A401corr_rudnick.png -no_y -noisefits ${F} > ptp_results/a401results_rudnick.txt
mv ptp_dir/A401/grid_${CELLSIZE}x${CELLSIZE}_results.fits ptp_results/A401_results_rudnick.fits

#Bridge
python analysis/ppt.py -radio1 ${F} -limits 1190 1263 1672 1800 -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -sys1 0.2 -excluderegion regions/excluderegionsbridge.reg
mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/bridge && wait
mv ptp_dir/bridge/grid_${CELLSIZE}x${CELLSIZE}_ds9_image.reg ptp_results/gridbridge_rudnick.reg
python analysis/make_corr_images.py -filein ptp_dir/bridge/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout ptp_results/bridgecorr_rudnick.png -no_y -noisefits ${F} > ptp_results/bridgeresults_rudnick.txt
mv ptp_dir/bridge/grid_${CELLSIZE}x${CELLSIZE}_results.fits ptp_results/bridge_results_rudnick.fits

F=fits/60cleanbridge_200kpc.fits
CELLSIZE=13
rm -rf ptp_dir && wait

A399
python analysis/ppt.py -radio1 ${F} -limits 725 520 945 740 -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -sys1 0.2 -excluderegion regions/excluderegionshalo.reg
mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/A399 && wait
mv ptp_dir/A399/grid_${CELLSIZE}x${CELLSIZE}_ds9_image.reg ptp_results/gridA399_cb.reg
python analysis/make_corr_images.py -filein ptp_dir/A399/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout ptp_results/A399corr_cb.png -no_y -noisefits ${F} > ptp_results/a399results_cb.txt
mv ptp_dir/A399/grid_${CELLSIZE}x${CELLSIZE}_results.fits ptp_results/A399_results_cb.fits

#A401
python analysis/ppt.py -radio1 ${F} -limits 600 850 715 1025 -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -sys1 0.2 -excluderegion regions/excluderegionshalo.reg
mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/A401 && wait
mv ptp_dir/A401/grid_${CELLSIZE}x${CELLSIZE}_ds9_image.reg ptp_results/gridA401_cb.reg
python analysis/make_corr_images.py -filein ptp_dir/A401/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout ptp_results/A401corr_cb.png -no_y -noisefits ${F} > ptp_results/a401results_cb.txt
mv ptp_dir/A401/grid_${CELLSIZE}x${CELLSIZE}_results.fits ptp_results/A401_results_cb.fits

#Bridge
python analysis/ppt.py -radio1 ${F} -limits 595 631 836 900 -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -sys1 0.2 -excluderegion regions/excluderegionsbridge.reg
mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/bridge && wait
mv ptp_dir/bridge/grid_${CELLSIZE}x${CELLSIZE}_ds9_image.reg ptp_results/gridbridge_cb.reg
python analysis/make_corr_images.py -filein ptp_dir/bridge/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout ptp_results/bridgecorr_cb.png -no_y -noisefits ${F} > ptp_results/bridgeresults_cb.txt
mv ptp_dir/bridge/grid_${CELLSIZE}x${CELLSIZE}_results.fits ptp_results/bridge_results_cb.fits

#Bridge
python analysis/ppt.py -radio1 ${F} -limits 595 631 836 900 -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -sys1 0.2 -excluderegion regions/excluderegionsbridge_ext.reg
mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/bridge && wait
mv ptp_dir/bridge/grid_${CELLSIZE}x${CELLSIZE}_ds9_image.reg ptp_results/gridbridge_cb_ext.reg
python analysis/make_corr_images.py -filein ptp_dir/bridge/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout ptp_results/bridgecorr_cb_ext.png -no_y -noisefits ${F} > ptp_results/bridgeresults_cb_ext.txt
mv ptp_dir/bridge/grid_${CELLSIZE}x${CELLSIZE}_results.fits ptp_results/bridge_results_cb_ext.fits

#python analysis/make_corr_mix.py -obj A399
#python analysis/make_corr_mix.py -obj A401
#python analysis/make_corr_mix.py -obj bridge