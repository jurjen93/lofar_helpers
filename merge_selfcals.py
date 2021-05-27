from h5_merger import merge_h5
from glob import glob
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-d', '--directory', type=str, help='directory path')
parser.add_argument('-del', '--exclude_boxes', help='Exclude the following boxes (numbers only)')
args = parser.parse_args()


excluded_boxes=args.exclude_boxes
if type(excluded_boxes)==str:
    excluded_boxes = excluded_boxes.split(',')

excluded_boxes = ['box_'+n for n in excluded_boxes]

h5_files = []
for box in sorted(glob('{directory}/box_*'.format(directory=args.directory)))[0:5]:
    try:
        last_merged = sorted(glob('{box}/merged_selfcalcyle*_*.ms.*h5'.format(box=box)))[-1]
        print(last_merged)
        if any(box in last_merged for box in excluded_boxes):
            print('Exclude '+last_merged.split('/')[-1])
        else:
            h5_files.append(last_merged)
    except:
        print("No merged_selfcal* in {box}".format(box=box.split('/')[-1]))

merge_h5(h5_out='all_directions.h5',
         h5_files=h5_files)