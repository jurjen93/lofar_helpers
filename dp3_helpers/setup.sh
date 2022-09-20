#!/bin/bash

#RUN THIS SCRIPT IN THE SAME FOLDER AS polconv.py

BASEDIR=$(realpath $0)

export PYTHONPATH="${BASEDIR%/*}/":$PYTHONPATH
