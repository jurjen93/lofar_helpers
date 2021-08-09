#!/bin/bash

number=$1

TOTAL_BOXES=$(ls -dq boxes/box*.reg | wc -l)

if [[ ! ${number} > ${TOTAL_BOXES} ]]
then
  echo "HOI"
fi

if [[ ! ${number} -gt ${TOTAL_BOXES} ]]
then
  echo "HOI2"
fi
