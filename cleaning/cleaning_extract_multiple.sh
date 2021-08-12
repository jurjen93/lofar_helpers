#!/bin/bash

BOXSTART=$1
BOXEND=$2

for ((N=${BOXSTART};N<=${BOXEND};N++))
do
  source /home/lofarvwf-jdejong/scripts/lofar_helpers/cleaning/cleaning_extract_per_box.sh ${N}
done