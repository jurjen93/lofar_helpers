#!/bin/bash

taskset 1 cat /proc/self/stat | gawk '{print 39}'

echo "CPU ID:"
CPU_ID=$(cat /proc/self/stat)
echo $CPU_ID | gawk '{print $39}'
# find all CPU ID's that are allocated (reserved) by self
echo "CPU ALLOC:"
CPU_ALLOC=$(cat /proc/self/status | grep 'Cpus_allowed_list')
echo $CPU_ALLOC