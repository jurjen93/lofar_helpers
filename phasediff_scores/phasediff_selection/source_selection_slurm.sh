#!/bin/bash
#SBATCH -c 10 --constraint=amd --array=0-30 -t 1:00:00

#MSLIST TEXT FILE WITH PATH TO MS
PLIST=$1

#SINGULARITY
BIND=$( python3 $HOME/parse_settings.py --BIND ) # SEE --> https://github.com/jurjen93/lofar_vlbi_helpers/blob/main/parse_settings.py
SIMG=$( python3 $HOME/parse_settings.py --SIMG )
#SCRIPT FOLDER
LOFAR_HELPERS=$( python3 $HOME/parse_settings.py --lofar_helpers )
SCRIPT_DIR=$PWD

#phasediff output folder
mkdir -p phasediff_h5s

DIR=$( awk NR==${SLURM_ARRAY_TASK_ID} $PLIST )

echo "$DIR"
mkdir $DIR
cp -r *${DIR}*.ms $DIR
cd $DIR
chmod 755 ${SCRIPT_DIR}/*
singularity exec -B $BIND $SIMG ${SCRIPT_DIR}/phasediff_multi.sh
mv scalarphasediff0*phaseup.h5 ../phasediff_h5s
cd ../
#rm -rf $DIR

#RETURN SCORES
singularity exec -B $BIND $SIMG python ${SCRIPT_DIR}/phasediff_output.py --h5 phasediff_h5s/*.h5
