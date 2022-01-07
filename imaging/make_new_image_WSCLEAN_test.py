__author__ = "Jurjen de Jong (jurjendejong@strw.leidenuniv.nl)"

import os
import sys
import tables
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--N', type=str, help='archive number', required=True)
args = parser.parse_args()

N = args.N
MS = 'Abell399-401_extr.dysco.sub.shift.avg.weights.ms.archive' + N + '.avg.goodtimes'
H5 = 'all_directions'+N+'.h5'

TO='/net/nieuwerijn/data2/jurjendejong/Abell399-401_' + N
FROM='/net/tussenrijn/data2/jurjendejong/A399_extracted_avg'


f = open(TO+'/tess.reg')
tess = f.read()
f.close()
H = tables.open_file(TO+'/'+H5)
if len(H.root.sol000.phase000.dir[:])!=len(tess.split('polygon'))-1:
    sys.exit('ERROR: H5 and tess.reg do not match')

#----------------------------------------------------------------------------------------------------------------------

#MAKE WSCLEAN COMMAND
with open('/'.join(__file__.split('/')[0:-1])+'/WSCLEAN_scripts/wsclean.txt') as f:
    lines = [l.replace('\n', '') for l in f.readlines()]
    lines += ['-facet-regions '+TO+'/tess.reg']
    lines += ['-apply-facet-solutions '+TO+'/'+H5+' amplitude000,phase000']
    lines += ['-name image_test_L626678']
    lines += ['-size 6000 6000']
    lines += ['-scale 1.5arcsec']
    lines += [TO+'/'+MS]

os.system('aoflagger '+TO+'/'+MS+' && wait')

cmd = ' '.join(['cd', TO, '&&'] + lines)
#RUN DDF COMMAND
print('Running WSCLEAN COMMAND')
print(cmd)
os.system(cmd + ' > '+TO+'/log.txt')
print('Finished making new image')