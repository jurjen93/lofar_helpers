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

TO='/net/tussenrijn/data2/jurjendejong/L626678_WSCLEAN_avg'
FROM='/net/tussenrijn/data2/jurjendejong/A399_extracted_avg'

#CREATE DESTINATION DIRECTORY IF NOT EXISTS
if not path.exists(TO):
    os.system('mkdir {LOCATION}'.format(LOCATION=TO))
    print('Created {LOCATION}'.format(LOCATION=TO))

#CLEAN CACHE
os.system('CleanSHM.py')

#----------------------------------------------------------------------------------------------------------------------


os.system('cd '+TO+' && python /net/rijn/data2/rvweeren/LoTSS_ClusterCAL/ds9facetgenerator.py --h5 all_directions0.h5 --DS9regionout '+
          TO+'/tess.reg --imsize 6000 '+'--ms '+FROM+'/Abell399-401_extr.dysco.sub.shift.avg.weights.ms.archive0.avg')

os.system('cp -r '+FROM+'/Abell399-401_extr.dysco.sub.shift.avg.weights.ms.archive0.avg.goodtimes'+' '+TO+' && wait')

CUTFREQS = [5021107868.011121, 5021107864.005561]

# for MS in glob(FROM+'/L626678_concat.ms'):
#     t = ct.table(MS)
#     time = t.getcol('TIME')[0]
#     t.close()
#     if not (time in CUTFREQS and '127' in MS):
#         print('Making goodtimes for'+MS)
#         os.system("python /home/jurjendejong/scripts/lofar_helpers/supporting_scripts/flag_time.py -tf 0 3000 -msin " + MS + " -msout " + TO + '/' + MS.split('/')[-1] + '.goodtimes')


# for MS in glob(FROM+'/Abell399-401_extr*.ms.archive0.avg'):
#     if ct.table(MS).getcol('TIME')[0] in CUTFREQS:
#         print('Cutting freq for ' + MS)
#         os.system("python /home/jurjendejong/scripts/lofar_helpers/supporting_scripts/flag_freq.py -ff='[15..19]' -msin " + MS+" -msout " + TO + '/' + MS.split('/')[-1] + '.goodfreq')
#         os.system("python /home/jurjendejong/scripts/lofar_helpers/supporting_scripts/flag_time.py -tf 0 1500 -msin " + MS + " -msout " + TO + '/' + MS.split('/')[-1] + '.goodfreq.goodtimes')
#         os.system("rm -rf " + TO + '/' + MS.split('/')[-1] + '.goodfreq')
#     else:
#         print('Cutting time for '+MS)
#         os.system("python /home/jurjendejong/scripts/lofar_helpers/supporting_scripts/flag_time.py -tf 0 1500 -msin " + MS + " -msout " + TO + '/' + MS.split('/')[-1] + '.goodtimes')


#----------------------------------------------------------------------------------------------------------------------

#MAKE LIST WITH MEASUREMENT SETS
# os.system('ls -1d '+TO+'/L626678_concat.ms.goodtimes > '+TO+'/big-mslist.txt'.format(LOCATION=TO))

#----------------------------------------------------------------------------------------------------------------------

#MAKE WSCLEAN COMMAND
#TODO: SINGULARITY: /net/lofar1/data1/sweijen/software/LOFAR/singularity/test/test_wsclean_facet_fix_sep30.sif
with open('/'.join(__file__.split('/')[0:-1])+'/WSCLEAN_scripts/wsclean.txt') as f:
    lines = [l.replace('\n', '') for l in f.readlines()]
    lines += ['-facet-regions '+TO+'/tess.reg']
    lines += ['-apply-facet-solutions '+TO+'/all_directions0.h5 amplitude000,phase000']
    lines += ['-name image_test_L626678']
    lines += [TO+'/Abell399-401_extr.dysco.sub.shift.avg.weights.ms.archive0.avg.goodtimes']

cmd = ' '.join(['cd', TO, '&&'] + lines)
#RUN DDF COMMAND
print('Running WSCLEAN COMMAND')
print(cmd)
os.system(cmd)
print('Finished making new image')