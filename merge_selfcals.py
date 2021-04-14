from h5_merger import merge_h5
from glob import glob
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-d', '--directory', type=str, help='directory path')
args = parser.parse_args()

h5_files = []
for box in glob('{directory}/box_*'.format(directory=args.directory)):
  try:
    h5_files.append(sorted(glob('{box}/merged_selfcalcyle*_*.ms.*h5'.format(box=box)))[-1])
  except:
    print("No merged_selfcal* in {box}".format(box=box.split('/')[-1]))

merge_h5(h5_out='all_directions.h5',
         h5_files=h5_files)