#!/bin/bash
#SBATCH -c 10 --constraint=amd

#MSLIST TEXT FILE WITH PATH TO MS
MSLIST=$1

#SINGULARITY
BIND=$( python3 $HOME/parse_settings.py --BIND ) # SEE --> https://github.com/jurjen93/lofar_vlbi_helpers/blob/main/parse_settings.py
SIMG=$( python3 $HOME/parse_settings.py --SIMG )
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

mkdir -p h5output

while read -r MS; do
  singularity exec -B $BIND $SIMG SCRIPT_DIR/phasediff.sh ${MS}
done <$MSLIST

singularity exec -B $BIND $SIMG python SCRIPT_DIR/phasediff_output.py --h5 h5output/*.h5
