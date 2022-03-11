#!/bin/bash

cd /home/jurjen/Documents/Python/lofar_helpers

CELLSIZE=35

F=$1 #fits/60rudnick.fits

rm -rf ptp_dir && wait

#A399
python analysis/ppt.py -radio1 ${F} -limits 1450 1040 1890 1480 -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -sys1 0.2 -excluderegion regions/excluderegionshalo.reg
mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/A399 && wait
python analysis/make_corr_images.py -filein ptp_dir/A399/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout analysis/A399corr.png > a399results.txt


#A401
python analysis/ppt.py -radio1 ${F} -limits 1200 1700 1430 2050 -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -sys1 0.2 -excluderegion regions/excluderegionshalo.reg
mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/A401 && wait
python analysis/make_corr_images.py -filein ptp_dir/A401/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout analysis/A401corr.png > a401results.txt


#Bridge
python analysis/ppt.py -radio1 ${F} -limits 1190 1263 1672 1800 -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -sys1 0.2 -excluderegion regions/excluderegionsbridge.reg
mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/bridge && wait
python analysis/make_corr_images.py -filein ptp_dir/bridge/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout analysis/bridgecorr.png > bridgeresults.txt

#python analysis/ppt.py  -radio1 fits/20median.fits -limits 719 300 1147 677 -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -sys1 0.2 -excluderegion regions/excluderegionshalo.reg
#mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/A399 && wait
#python analysis/make_corr_images.py -filein ptp_dir/A399/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout analysis/A399corr.png > a399results.txt


#A401
#python analysis/ppt.py  -radio1 fits/20median.fits -limits 515 1000 650 1230 -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -sys1 0.2 -excluderegion regions/excluderegionshalo.reg
#mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/A401 && wait
#python analysis/make_corr_images.py -filein ptp_dir/A401/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout analysis/A401corr.png > a401results.txt


#Bridge
#python analysis/ppt.py  -radio1 fits/20median.fits -limits 500 490 950 1080 -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -sys1 0.2 -excluderegion regions/excluderegionsbridge.reg
#mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/bridge && wait
#python analysis/make_corr_images.py -filein ptp_dir/bridge/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout analysis/bridgecorr.png > bridgeresults.txt