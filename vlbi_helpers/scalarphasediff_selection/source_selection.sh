#!/bin/bash

MSLIST=$1

mkdir -p h5output

while read -r MS; do
  source ./scalarphasediff.sh ${MS}
done <$MSLIST

python scalarphasediff_output.py --h5 h5output/*.h5
