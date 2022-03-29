#!/bin/bash

cd /home/jurjen/Documents/Python/lofar_helpers

F=fits/60rudnick.fits #fits/60rudnick.fits
CELLSIZE=25
OUTFILE=ptp_results_${CELLSIZE}

mkdir ${OUTFILE}

for N in {84..100}
do
    rm -rf ptp_dir && wait
    #A399
     python analysis/ppt.py -radio1 ${F} -limits $(($((1+$RANDOM % 10))+1450)) $(($((1+$RANDOM % 10))+1040)) $(($((1+$RANDOM % 10))+1890)) $(($((1+$RANDOM % 10))+1480)) -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -sys1 0.2 -excluderegion regions/excluderegionshalo.reg
    mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/A399 && wait
    mv ptp_dir/A399/grid_${CELLSIZE}x${CELLSIZE}_ds9_image.reg ${OUTFILE}/gridA399_rudnick_${N}.reg
    python analysis/make_corr_images.py -filein ptp_dir/A399/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout ${OUTFILE}/A399corr_rudnick_${N}.png -no_y -noisefits ${F} > ${OUTFILE}/a399results_rudnick_${N}.txt
    mv ptp_dir/A399/grid_${CELLSIZE}x${CELLSIZE}_results.fits ${OUTFILE}/A399_results_rudnick_${N}.fits

    #A401
    python analysis/ppt.py -radio1 ${F} -limits $(($((1+$RANDOM % 10))+1200)) $(($((1+$RANDOM % 10))+1700)) $(($((1+$RANDOM % 10))+1430)) $(($((1+$RANDOM % 10))+2050)) -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -sys1 0.2 -excluderegion regions/excluderegionshalo.reg
    mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/A401 && wait
    mv ptp_dir/A401/grid_${CELLSIZE}x${CELLSIZE}_ds9_image.reg ${OUTFILE}/gridA401_rudnick_${N}.reg
    python analysis/make_corr_images.py -filein ptp_dir/A401/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout ${OUTFILE}/A401corr_rudnick_${N}.png -no_y -noisefits ${F} > ${OUTFILE}/a401results_rudnick_${N}.txt
    mv ptp_dir/A401/grid_${CELLSIZE}x${CELLSIZE}_results.fits ${OUTFILE}/A401_results_rudnick_${N}.fits

    #Bridge
    python analysis/ppt.py -radio1 ${F} -limits $(($((1+$RANDOM % 10))+1190)) $(($((1+$RANDOM % 10))+1263)) $(($((1+$RANDOM % 10))+1672)) $(($((1+$RANDOM % 10))+1800)) -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -sys1 0.2 -excluderegion regions/excluderegionsbridge.reg
    mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/bridge && wait
    mv ptp_dir/bridge/grid_${CELLSIZE}x${CELLSIZE}_ds9_image.reg ${OUTFILE}/gridbridge_rudnick_${N}.reg
    python analysis/make_corr_images.py -filein ptp_dir/bridge/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout ${OUTFILE}/bridgecorr_rudnick_${N}.png -no_y -noisefits ${F} > ${OUTFILE}/bridgeresults_rudnick_${N}.txt
    mv ptp_dir/bridge/grid_${CELLSIZE}x${CELLSIZE}_results.fits ${OUTFILE}/bridge_results_rudnick_${N}.fits
done

F=fits/60cleanbridge_200kpc.fits
CELLSIZE=13

for N in {1..100}
do
    #A399
    rm -rf ptp_dir && wait
    python analysis/ppt.py -radio1 ${F} -limits $(($((1+$RANDOM % 10))+725)) $(($((1+$RANDOM % 10))+520)) $(($((1+$RANDOM % 10))+945)) $(($((1+$RANDOM % 10))+740)) -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -sys1 0.2 -excluderegion regions/excluderegionshalo.reg
    mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/A399 && wait
    mv ptp_dir/A399/grid_${CELLSIZE}x${CELLSIZE}_ds9_image.reg ${OUTFILE}/gridA399_cb_${N}.reg
    python analysis/make_corr_images.py -filein ptp_dir/A399/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout ${OUTFILE}/A399corr_cb_${N}.png -no_y -noisefits ${F} > ${OUTFILE}/a399results_cb_${N}.txt
    mv ptp_dir/A399/grid_${CELLSIZE}x${CELLSIZE}_results.fits ${OUTFILE}/A399_results_cb_${N}.fits

    #A401
    python analysis/ppt.py -radio1 ${F} -limits $(($((1+$RANDOM % 10))+600)) $(($((1+$RANDOM % 10))+850)) $(($((1+$RANDOM % 10))+715)) $(($((1+$RANDOM % 10))+1025)) -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -sys1 0.2 -excluderegion regions/excluderegionshalo.reg
    mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/A401 && wait
    mv ptp_dir/A401/grid_${CELLSIZE}x${CELLSIZE}_ds9_image.reg ${OUTFILE}/gridA401_cb_${N}.reg
    python analysis/make_corr_images.py -filein ptp_dir/A401/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout ${OUTFILE}/A401corr_cb_${N}.png -no_y -noisefits ${F} > ${OUTFILE}/a401results_cb_${N}.txt
    mv ptp_dir/A401/grid_${CELLSIZE}x${CELLSIZE}_results.fits ${OUTFILE}/A401_results_cb_${N}.fits

    #Bridge
    python analysis/ppt.py -radio1 ${F} -limits $(($((1+$RANDOM % 10))+595)) $(($((1+$RANDOM % 10))+631)) $(($((1+$RANDOM % 10))+836)) $(($((1+$RANDOM % 10))+900)) -xsou fits/mosaic_a399_a401.fits -xbkg fits/mosaic_a399_a401_bkg.fits -xexp fits/mosaic_a399_a401_exp.fits -cellsize ${CELLSIZE} -y fits/a401_curdecmaps_0.2_1.5s_sz.fits -sys1 0.2 -excluderegion regions/excluderegionsbridge.reg
    mv ptp_dir/grid_${CELLSIZE}x${CELLSIZE} ptp_dir/bridge && wait
    mv ptp_dir/bridge/grid_${CELLSIZE}x${CELLSIZE}_ds9_image.reg ${OUTFILE}/gridbridge_cb_${N}.reg
    python analysis/make_corr_images.py -filein ptp_dir/bridge/grid_${CELLSIZE}x${CELLSIZE}_results.fits -fileout ${OUTFILE}/bridgecorr_cb_${N}.png -no_y -noisefits ${F} > ${OUTFILE}/bridgeresults_cb_${N}.txt
    mv ptp_dir/bridge/grid_${CELLSIZE}x${CELLSIZE}_results.fits ${OUTFILE}/bridge_results_cb_${N}.fits
done
