__author__ = "Jurjen de Jong (jurjendejong@strw.leidenuniv.nl)"

from h5_merger import merge_h5
from glob import glob
from argparse import ArgumentParser, ArgumentTypeError

def str2bool(v):
    v = str(v)
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise ArgumentTypeError('Boolean value expected.')

parser = ArgumentParser()
parser.add_argument('-d', '--directory', type=str, help='directory path')
parser.add_argument('-del', '--exclude_boxes', help='Exclude the following boxes (numbers only)')
parser.add_argument('-nd', '--make_new_direction', type=str2bool, nargs='?', const=True, default=True, help='make new directions')
args = parser.parse_args()

def get_digits(x):
    return int(''.join([d for d in x if d.isdigit()]))

if args.exclude_boxes:
    excluded_boxes=args.exclude_boxes
    if type(excluded_boxes)==str:
        excluded_boxes = excluded_boxes.split(',')

    excluded_boxes = ['box_'+n for n in excluded_boxes]
else:
    excluded_boxes = []

h5_tables = []
boxes_h5_list = glob('{directory}/box_*'.format(directory=args.directory))
boxes_h5_list.sort(key=lambda x: get_digits(x))
for box in boxes_h5_list:
    try:
        last_merged = sorted(glob('{box}/merged_selfcalcyle*_*.ms.*h5'.format(box=box)))[-1]
        print(last_merged)
        if any(box in last_merged for box in excluded_boxes):
            print('Exclude '+last_merged.split('/')[-1])
        else:
            h5_tables.append(last_merged)
    except:
        print("No merged_selfcal* in {box}".format(box=box.split('/')[-1]))

merge_h5(h5_out='all_directions.h5',
         h5_tables=h5_tables,
         make_new_direction=args.make_new_direction)