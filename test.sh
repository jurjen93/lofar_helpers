#!/bin/bash

test=test

while [[ ! -f /home/jurjen/Documents/Python/lofar_helpers/${test}.txt ]]
  do
    sleep 5
    echo "test"
  done