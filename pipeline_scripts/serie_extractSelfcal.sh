#!/bin/bash
#SBATCH -c 1 -t 60:00

SOURCE=$1

#extract
sh ./extract.sh ${SOURCE}
#selfcal
sh ./selfcal.sh ${SOURCE}