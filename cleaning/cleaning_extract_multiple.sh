#!/bin/bash

BOXSTART=$1
BOXEND=$2

for ((N=${BOXSTART};N<=${BOXEND};N++))
do
  ./cleaning_extract_per_box.sh ${N}
done