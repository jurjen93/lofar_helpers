#!/bin/bash

FOLDER=$1
NAME=$2

REGEX='([^\/]+$)'

for BOXFOLDER in ${FOLDER}/box*
do
  [[ ${BOXFOLDER} =~ ${REGEX} ]] # $pat must be unquoted
  BOX=${BASH_REMATCH[1]}
  rclone --config=~/macaroon/maca_lofarvwf.conf mkdir maca_lofarvwf:/disk/${NAME}/extract/${BOX}
  rclone --config=~/macaroon/maca_lofarvwf.conf copy ${FOLDER}/${BOX}/ maca_lofarvwf:/disk/${NAME}/extract/${BOX} -P
done