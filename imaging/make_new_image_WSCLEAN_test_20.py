__author__ = "Jurjen de Jong (jurjendejong@strw.leidenuniv.nl)"

import os
import sys
import tables
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--N', type=str, help='archive number', required=True)
parser.add_argument('--nmiter', type=str, default=None, help='max major iterations')
args = parser.parse_args()

N = args.N
if not args.nmiter:
    nmiter = '8'
else:
    nmiter = args.nmiter

MS = 'Abell399-401_extr.dysco.sub.shift.avg.weights.ms.archive' + N + '.avg.goodtimes'
H5 = 'all_directions'+N+'.h5'

TO='/net/nieuwerijn/data2/jurjendejong/Abell399-401_' + N + '_20'
FROM='/net/tussenrijn/data2/jurjendejong/A399_extracted_avg'


f = open(TO+'/tess.reg')
tess = f.read()
f.close()
H = tables.open_file(TO+'/'+H5)
if len(H.root.sol000.phase000.dir[:])!=len(tess.split('polygon'))-1:
    sys.exit('ERROR: H5 and tess.reg do not match')

#----------------------------------------------------------------------------------------------------------------------

#MAKE WSCLEAN COMMAND
lines = ['wsclean -size 3000 3000 -use-wgridder -no-update-model-required -reorder -weight briggs -0.5 -weighting-rank-filter 3 ' \
        '-clean-border 1 -parallel-reordering  6 -padding 1.2 -auto-mask 2.5 -auto-threshold 0.5 -pol i -name image_test_L626678_rvw20 ' \
        '-scale 3arcsec -niter 50000 -mgain 0.8 -fit-beam -multiscale -channels-out 6 -join-channels -multiscale-max-scales 7 -nmiter ' + nmiter +\
        ' -log-time -multiscale-scale-bias 0.7 -facet-regions '+TO+'/' \
        'tess.reg ' \
        '-parallel-gridding 6 -fit-spectral-pol 3 -taper-gaussian 20arcsec -parallel-deconvolution 1600 ' \
        '-apply-facet-solutions '+TO+'/'+H5+' amplitude000,phase000 '+ TO + \
        '/'+MS]

os.system('aoflagger '+TO+'/'+MS+' && wait')

cmd = ' '.join(['cd', TO, '&&'] + lines)
#RUN DDF COMMAND
print('Running WSCLEAN COMMAND')
print(cmd)
os.system(cmd + ' > '+TO+'/log.txt')
print('Finished making new image')