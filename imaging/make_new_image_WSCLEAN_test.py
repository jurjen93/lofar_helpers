"""
Example:
    python ~/scripts/lofar_helpers/make_new_image.py
        -from /disks/paradata/shimwell/LoTSS-DR2/archive_other/L626678
        -to /net/tussenrijn/data2/jurjendejong/L626678/result
"""

__author__ = "Jurjen de Jong (jurjendejong@strw.leidenuniv.nl)"

import os
from glob import glob
from os import path
import casacore.tables as ct
import tables
import numpy as np

TO='/net/nieuwerijn/data2/jurjendejong/L626678_WSCLEAN'
FROM='/net/nieuwerijn/data2/jurjendejong/L626678_WSCLEAN'

#CREATE DESTINATION DIRECTORY IF NOT EXISTS
if not path.exists(TO):
    os.system('mkdir {LOCATION}'.format(LOCATION=TO))
    print('Created {LOCATION}'.format(LOCATION=TO))

#CLEAN CACHE
os.system('CleanSHM.py')

#----------------------------------------------------------------------------------------------------------------------
CUTFREQS = [5021107868.011121, 5021107864.005561]

for MS in glob(FROM+'/L626678_concat.ms'):
    t = ct.table(MS)
    time = t.getcol('TIME')[0]
    t.close()
    if not (time in CUTFREQS and '127' in MS):
        print('Making goodtimes for'+MS)
        os.system("python /home/jurjendejong/scripts/lofar_helpers/supporting_scripts/flag_time.py -tf 0 3000 -msin " + MS + " -msout " + TO + '/' + MS.split('/')[-1] + '.goodtimes')

#----------------------------------------------------------------------------------------------------------------------

#MAKE LIST WITH MEASUREMENT SETS
os.system('ls -1d '+TO+'/L626678_concat.ms.goodtimes > '+TO+'/big-mslist.txt'.format(LOCATION=TO))

#----------------------------------------------------------------------------------------------------------------------

#MAKE WSCLEAN COMMAND
with open('/'.join(__file__.split('/')[0:-1])+'/WSCLEAN_scripts/wsclean.txt') as f:
    lines = [l.replace('\n','') for l in f.readlines()]

#RUN DDF COMMAND
print('Running WSCLEAN COMMAND')
os.system(' '.join(['cd', TO, '&&'] + lines))
print('Finished making new image')