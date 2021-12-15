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

os.system('cd '+TO+' && python /net/rijn/data2/rvweeren/LoTSS_ClusterCAL/ds9facetgenerator.py --h5 all_directions0.h5 --DS9regionout '+
          TO+'/tess.reg --imsize 14000 --plottesselation '+'--ms '+FROM+'/Abell399-401_extr*.ms.archive0')


for MS in glob(FROM+'/Abell399-401_extr*.ms.archive*'):
    if ct.table(MS).getcol('TIME')[0] in CUTFREQS:
        print('Cutting freq for ' + MS)
        os.system("python /home/lofarvwf-jdejong/scripts/lofar_helpers/supporting_scripts/flag_freq.py -ff='[15..19]' -msin " + MS+" -msout " + TO + '/' + MS.split('/')[-1] + '.goodfreq')
        os.system("python /home/lofarvwf-jdejong/scripts/lofar_helpers/supporting_scripts/flag_time.py -tf 0 3000 -msin " + MS + " -msout " + TO + '/' + MS.split('/')[-1] + '.goodfreq.goodtimes')
        os.system("rm -rf " + TO + '/' + MS.split('/')[-1] + '.goodfreq')
    else:
        print('Cutting time for '+MS)
        os.system("python /home/lofarvwf-jdejong/scripts/lofar_helpers/supporting_scripts/flag_time.py -tf 0 3000 -msin " + MS + " -msout " + TO + '/' + MS.split('/')[-1] + '.goodtimes')

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