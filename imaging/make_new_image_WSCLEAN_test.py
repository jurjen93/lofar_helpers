__author__ = "Jurjen de Jong (jurjendejong@strw.leidenuniv.nl)"

import os
import sys
import tables
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--N', type=str, help='archive number', required=True)
parser.add_argument('--nmiter', type=str, default=None, help='max major iterations')
parser.add_argument('--h5', type=str, default=None, help='max major iterations')
parser.add_argument('--ms', type=str, default=None, help='max major iterations')

args = parser.parse_args()

N = args.N
if not args.nmiter:
    nmiter = '8'
else:
    nmiter = args.nmiter

MS = args.ms
H5 = args.h5
FACET = 'tessupdate.reg'

TO='/net/nieuwerijn/data2/jurjendejong/Abell399-401_' + N
FROM='/net/tussenrijn/data2/jurjendejong/A399_extracted_avg'


f = open(TO+'/'+FACET)
tess = f.read()
f.close()
H = tables.open_file(TO+'/'+H5)
if len(H.root.sol000.phase000.dir[:])!=len(tess.split('polygon'))-1:
    sys.exit('ERROR: H5 and tess.reg do not match')

#----------------------------------------------------------------------------------------------------------------------

#MAKE WSCLEAN COMMAND
with open('/'.join(__file__.split('/')[0:-1])+'/WSCLEAN_scripts/wsclean.txt') as f:
    lines = [l.replace('\n', '') for l in f.readlines()]
    lines += ['-facet-regions '+TO+'/'+FACET]
    lines += ['-apply-facet-solutions '+TO+'/'+H5+' amplitude000,phase000']
    lines += ['-name image_test_L626678']
    lines += ['-size 6000 6000']
    lines += ['-scale 1.5arcsec']
    lines += ['-nmiter '+nmiter]
    lines += [TO+'/'+MS]

os.system('aoflagger '+TO+'/'+MS+' && wait')

cmd = ' '.join(['cd', TO, '&&'] + lines)
#RUN DDF COMMAND
print('Running WSCLEAN COMMAND')
print(cmd)
os.system(cmd + ' > '+TO+'/log.txt')
print('Finished making new image')