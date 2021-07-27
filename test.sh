#!/bin/bash


for ((N=1;N<=20;N++))
do
  until [[ -f /home/jurjen/Documents/Python/lofar_helpers/test_${N}.txt ]]
  do
    sleep 2
    echo "Waiting for ${N}"
  done &
  echo "FINISHED ${N}"
done
wait