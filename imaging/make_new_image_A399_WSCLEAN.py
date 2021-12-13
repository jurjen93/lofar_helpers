"""
NOT FULLY FINISHED. PROBLEM WITH FEEDING WSCLEAN MULTIPLE SOLUTIONS.
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
#starting times for measurement sets that have to be cutted for freq
CUTFREQS = [5021107868.011121, 5021107864.005561]

for n, L6 in ['L626678', 'L626692', 'L626706', 'L632229', 'L632511', 'L632525']:
    os.system('cd '+TO+' && python /net/rijn/data2/rvweeren/LoTSS_ClusterCAL/ds9facetgenerator.py --h5 all_directions0.h5 --DS9regionout '+
              TO+'/tess.reg --imsize 14000 --plottesselation '+'--ms '+glob(FROM+'/extr*.ms')[0])

for MS in glob(FROM+'/extr*.ms'):
    t = ct.table(MS)
    time = t.getcol('TIME')[0]
    t.close()
    if not (time in CUTFREQS and '127' in MS):
        print('Making goodtimes for'+MS)
        os.system("python /home/jurjendejong/scripts/lofar_helpers/supporting_scripts/flag_time.py -tf 0 3000 -msin " + MS + " -msout " + TO + '/' + MS.split('/')[-1] + '.goodtimes')

#----------------------------------------------------------------------------------------------------------------------

#MAKE LIST WITH MEASUREMENT SETS
os.system('ls -1d '+TO+'/extr*.ms.goodtimes > '+TO+'/big-mslist.txt'.format(LOCATION=TO))

#----------------------------------------------------------------------------------------------------------------------

#MAKE WSCLEAN COMMAND
with open('/'.join(__file__.split('/')[0:-1])+'/WSCLEAN_scripts/wsclean.txt') as f:
    lines = [l.replace('\n','') for l in f.readlines()]
    lines += ['-facet-regions '+TO+'/tess.reg',
    '-apply-facet-solutions '+TO+'/all_directions*.h5 amplitude000,phase000']
    lines += [TO+'/A399_concat.ms']

#RUN DDF COMMAND
print('Running WSCLEAN COMMAND')
print(' '.join(['cd', TO, '&&'] + lines))
os.system(' '.join(['cd', TO, '&&'] + lines))
print('Finished making new image')