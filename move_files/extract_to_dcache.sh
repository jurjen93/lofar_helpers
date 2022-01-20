#!/bin/bash

BOXES=$1/box*
NAME=$2

REGEX='([^\/]+$)'
echo ${BOXES}
for BOXFOLDER in ${BOXES}
do
  [[ ${BOXFOLDER} =~ ${REGEX} ]] # $pat must be unquoted
  BOX=${BASH_REMATCH[1]}
  rclone --config=/home/lofarvwf-jdejong/macaroon/maca_lofarvwf.conf mkdir maca_lofarvwf:/disk/${NAME}/extract/${BOX}
  rclone --config=/home/lofarvwf-jdejong/macaroon/maca_lofarvwf.conf copy ${FOLDER}/${BOX}/ maca_lofarvwf:/disk/${NAME}/extract/${BOX} -P
done