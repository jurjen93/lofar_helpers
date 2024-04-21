#!/bin/bash

L=$1

DP3 \
msin=concat_${L}.ms \
msout=apply_concat_${L}.ms \
msin.datacolumn=DATA \
msout.storagemanager=dysco \
beam_dir.type=applybeam \
beam_dir.direction=[243.04908752441406deg,55.42097854614258deg] \
beam_dir.updateweights=True \
ac0.type=applycal \
ac0.parmdb=../merged_${L}.h5 \
ac0.correction=amplitude000 \
ac1.type=applycal \
ac1.parmdb=../merged_${L}.h5 \
ac1.correction=phase000 \
beam_center.type=applybeam \
beam_center.direction=[] \
beam_center.updateweights=True \
steps=[beam_dir,ac0,ac1,bl,beam_center,avg] \
avg.timestep=3 \
avg.freqstep=5 \
avg.type=averager \
bl.type=filter \
bl.baseline='[CR]S*&&'
