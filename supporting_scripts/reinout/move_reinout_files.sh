#!/bin/bash

FOLDER=/net/rijn/data2/rvweeren/LoTSS_ClusterCAL

cp -r ${FOLDER}/lin2circ.py ~/scripts/reinout/ ~/scripts/lofar_helpers/supporting_scripts/reinout
cp -r ${FOLDER}/plot_tecandphase.py ~/scripts/reinout/ ~/scripts/lofar_helpers/supporting_scripts/reinout
cp -r ${FOLDER}/BLsmooth.py ~/scripts/reinout/ ~/scripts/lofar_helpers/supporting_scripts/reinout
cp -r ${FOLDER}/lib_multiproc.py ~/scripts/reinout/ ~/scripts/lofar_helpers/supporting_scripts/reinout
cp -r ${FOLDER}/runwscleanLBautoR.py ~/scripts/reinout/ ~/scripts/lofar_helpers/supporting_scripts/reinout
cp -r ${FOLDER}/sub-sources-outside-cluster.py ~/scripts/reinout/ ~/scripts/lofar_helpers/supporting_scripts/reinout