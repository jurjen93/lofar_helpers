import os
import sys
sys.path.insert(1, '~/scripts/lofar_helpers/supporting_scripts')
from supporting_scripts.get_DDS3 import get_DDS3
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--frm', type=str, help='from where')
parser.add_argument('--to', type=str, help='to where')
args = parser.parse_args()

# FROM --> /disks/paradata/shimwell/LoTSS-DR2/archive_other/L626678
# TO --> /project/lofarvwf/Share/jdejong/data/L626678

files = [D.split('/')[-1] for D in get_DDS3(args.frm)] + \
        ['*.ms.archive',
         'image_full_ampphase_di_m.NS.tessel.reg',
         'image_full_ampphase_di_m.NS.mask01.fits',
         'image_full_ampphase_di_m.NS.DicoModel',
         'image_dirin_SSD_m.npy.ClusterCat.npy',
         'SOLSDIR']

command_1 = 'tar -C {FROM} -cvzf /net/tussenrijn/data2/jurjendejong/data_archive.tar.gz --absolute-names '.format(FROM=args.frm) + ' '.join(files)
os.system(command_1)
command_2 = 'scp -r /net/tussenrijn/data2/jurjendejong/data_archive.tar.gz lofarvwf-jdejong@spider.surfsara.nl:{TO}'.format(TO=args.to)
os.system(command_2)
