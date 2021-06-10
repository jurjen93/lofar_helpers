#!/bin/bash

#Clean up wisdom files if you change python version in singularity

DIR="~/.fftw_wisdom/*"

if [[ -f ${DIR} ]]; then
  rm -rf ~/.fftw_wisdom/*
else
  echo "${DIR} not found. Nothing deleted."