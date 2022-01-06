"""
Example:
    python ~/scripts/lofar_helpers/make_new_image.py
        -from /disks/paradata/shimwell/LoTSS-DR2/archive_other/L626678
        -to /net/tussenrijn/data2/jurjendejong/L626678/result
"""

__author__ = "Jurjen de Jong (jurjendejong@strw.leidenuniv.nl)"

import os
import sys
from glob import glob
from os import path
import casacore.tables as ct
import tables
import numpy as np
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--N', type=str, help='archive number', required=True)
args = parser.parse_args()

N = args.N
if N!='4':
    MS = 'Abell399-401_extr.dysco.sub.shift.avg.weights.ms.archive' + N + '.avg.goodtimes'
else:
    MS = 'Abell399-401_extr.dysco.sub.shift.avg.weights.ms.archive' + N + '.avg.goodfreq.goodtimes'
H5 = 'all_directions'+N+'.h5'

TO='/net/nieuwerijn/data2/jurjendejong/Abell399-401_' + N+'_20arcsec'
FROM='/net/tussenrijn/data2/jurjendejong/A399_extracted_avg'

#CREATE DESTINATION DIRECTORY IF NOT EXISTS
if not path.exists(TO):
    os.system('mkdir {LOCATION}'.format(LOCATION=TO))
    print('Created {LOCATION}'.format(LOCATION=TO))

#CLEAN CACHE
os.system('CleanSHM.py')

#----------------------------------------------------------------------------------------------------------------------


# os.system('cd '+TO+' && python /net/rijn/data2/rvweeren/LoTSS_ClusterCAL/ds9facetgenerator.py --h5 all_directions0.h5 --DS9regionout '+
#           TO+'/tess.reg --imsize 6000 '+'--ms '+FROM+'/Abell399-401_extr.dysco.sub.shift.avg.weights.ms.archive0.avg')

os.system('cp ' + FROM + '/tess.reg ' + TO)
os.system('cp ' + FROM + '/'+H5+' ' + TO + ' && wait')
os.system('cp -r '+FROM+'/'+MS+' ' + TO + ' && wait')

f = open(FROM+'/tess.reg')
tess = f.read()
f.close()
H = tables.open_file(FROM+'/'+H5)
if len(H.root.sol000.phase000.dir[:])!=len(tess.split('polygon'))-1:
    sys.exit('ERROR: H5 and tess.reg do not match')

CUTFREQS = [5021107868.011121, 5021107864.005561]

#----------------------------------------------------------------------------------------------------------------------

#MAKE WSCLEAN COMMAND
#TODO: SINGULARITY: /net/lofar1/data1/sweijen/software/LOFAR/singularity/test/test_wsclean_facet_fix_sep30.sif
with open('/'.join(__file__.split('/')[0:-1])+'/WSCLEAN_scripts/wsclean.txt') as f:
    lines = [l.replace('\n', '') for l in f.readlines()]
    lines += ['-facet-regions '+TO+'/tess.reg']
    lines += ['-apply-facet-solutions '+TO+'/'+H5+' amplitude000,phase000']
    lines += ['-name image_test_L626678']
    lines += ['-size 3000 3000']
    lines += ['-scale 3arcsec']
    lines += [TO+'/'+MS]

os.system('aoflagger '+TO+'/'+MS+' && wait')

cmd = ' '.join(['cd', TO, '&&'] + lines)
#RUN DDF COMMAND
print('Running WSCLEAN COMMAND')
print(cmd)
os.system(cmd)
print('Finished making new image')