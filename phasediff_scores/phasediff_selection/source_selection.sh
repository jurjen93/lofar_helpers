#!/bin/bash
#SBATCH -c 10 --constraint=amd

#MSLIST WITH PATH TO MS
MSLIST=$1

#SINGULARITY
BIND=$( python3 $HOME/parse_settings.py --BIND )
SIMG=$( python3 $HOME/parse_settings.py --SIMG )

mkdir -p h5output

while read -r MS; do
  singularity exec -B $BIND $SIMG ./phasediff.sh ${MS}
done <$MSLIST

singularity exec -B $BIND $SIMG python phasediff_output.py --h5 h5output/*.h5
