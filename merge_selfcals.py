__author__ = "Jurjen de Jong (jurjendejong@strw.leidenuniv.nl)"

from h5_merger import merge_h5
from glob import glob
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-d', '--directory', type=str, help='directory path')
parser.add_argument('-a', '--archive', type=str, help='archive number', default='0')
parser.add_argument('-del', '--exclude_boxes', help='Exclude the following boxes (numbers only)')
# parser.add_argument('-nd', '--make_new_direction', type=str2bool, nargs='?', const=True, default=True, help='make new directions')
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

boxes_h5_list = glob('{directory}/box_*'.format(directory=args.directory))
boxes_h5_list.sort(key=lambda x: get_digits(x))
# number_of_measurements = len(glob(boxes_h5_list[0]+'/'+'tecandphase0_selfcalcyle000*'))

for box in boxes_h5_list[59:61]:
    merge_h5(h5_out='{box}/final_merge_{n}.h5'.format(box=box, n=str(args.archive)),
             h5_tables=[sorted(glob('{box}/tecandphase*_*.ms.archive{n}*h5'.format(box=box, n=args.archive)))[-1],
                        sorted(glob('{box}/scalarcomplexgain*_*.ms.archive{n}*h5'.format(box=box, n=args.archive)))[-1]],
             h5_time_freq=sorted(glob('{box}/merged_selfcalcyle*_*.ms.archive{n}*h5'.format(box=boxes_h5_list[0], n=args.archive)))[-1],
             merge_all_in_one=True)

merge_h5(h5_out='all_directions{n}.h5'.format(n=str(args.archive)),
         h5_tables=glob('box_*/final_merge_{n}.h5'.format(n=str(args.archive))))

merge_h5(h5_out='all_directions{n}_wrong.h5'.format(n=str(args.archive)),
         h5_tables=[sorted(glob('{box}/merged_selfcal*_*.ms.archive{n}*h5'.format(box=box, n=args.archive)))[-1] for box in boxes_h5_list[59:61]])