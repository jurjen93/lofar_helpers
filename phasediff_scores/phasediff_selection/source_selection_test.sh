#!/bin/bash
#SBATCH -c 10 --constraint=amd

#MSLIST WITH PATH TO MS
DIRLIST=$1

while read -r DIR; do
  sbatch ./phasediff.sh ${DIR}
done <$DIRLIST

#python phasediff_output.py --h5 P?????/*.h5