#!/bin/bash

cd /home/jurjen/Documents/Python/lofar_helpers

CELLSIZE=35

#A399
python analysis/ppt.py  -radio1 fits/20median.fits -limits 1490 1079 1699 1308 -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits
mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/A399
python analysis/make_corr_images.py -filein ptp_dir/A399/grid_${CELLSIZE}x${CELLSIZE}/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout analysis/A399corr.png

#A401
python analysis/ppt.py  -radio1 fits/20median.fits -limits 1216 1779 1388 1991 -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits
mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/A401
python analysis/make_corr_images.py -filein ptp_dir/A401/grid_${CELLSIZE}x${CELLSIZE}/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout analysis/A401corr.png


#Bridge
#python analysis/ppt.py  -radio1 fits/20median.fits -limits 1190 1263 1672 1776 -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize 51 -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -excluderegion excluderegion.reg
#mv ptp_dir/grid_51x51 ptp_dir/bridge
#python analysis/make_corr_images.py -filein ptp_dir/bridge/grid_51x51_results.fits -fileout analysis/bridgecorr.png