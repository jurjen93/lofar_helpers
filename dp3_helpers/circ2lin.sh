#!/bin/bash

MSIN=$1
MSOUT=$2

DPPP \
msin=$MSIN \
msout=$MSOUT \
msout.overwrite=1 \
steps=[pystep] \
pystep.python.module=polconv \
pystep.python.class=PolConv \
pystep.type=PythonDPPP \
pystep.circ2lin=1