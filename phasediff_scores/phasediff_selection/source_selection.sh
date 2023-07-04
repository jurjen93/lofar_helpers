#!/bin/bash
#SBATCH -c 10 --constraint=amd

#MSLIST TEXT FILE WITH PATH TO MS
MSLIST=$1

#SCRIPT FOLDER
LOFAR_HELPERS=/home/lofarvwf-jdejong/scripts/lofar_helpers
SCRIPT_DIR=${LOFAR_HELPERS}/phasediff_scores/phasediff_selection

#SINGULARITY
BIND=$( python3 $HOME/parse_settings.py --BIND ) # SEE --> https://github.com/jurjen93/lofar_vlbi_helpers/blob/main/parse_settings.py
SIMG=$( python3 $HOME/parse_settings.py --SIMG )

#phasediff output folder
mkdir -p phasediff_h5s

#RUN MS FROM MS LIST
while read -r MS; do
  mkdir ${MS}_folder
  mv ${MS} ${MS}_folder
  cd ${MS}_folder
  chmod 755 ${SCRIPT_DIR}/*
  singularity exec -B $BIND $SIMG source ${SCRIPT_DIR}/phasediff.sh ${MS}
  mv scalarphasediff0*phaseup.h5 ../phasediff_h5s
  mv ${MS} ../
  cd ../
  rm -rf ${MS}_folder
done <$MSLIST

#RETURN SCORES
singularity exec -B $BIND $SIMG python ${SCRIPT_DIR}/phasediff_output.py --h5 phasediff_h5s/*.h5
