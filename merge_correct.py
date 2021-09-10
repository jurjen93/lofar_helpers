__author__ = "Jurjen de Jong (jurjendejong@strw.leidenuniv.nl)"

from h5_merger import merge_h5
from glob import glob
from argparse import ArgumentParser
import tables
import numpy as np
import os

parser = ArgumentParser()
parser.add_argument('-d', '--directory', type=str, help='directory path')
parser.add_argument('-a', '--archive', type=str, help='archive number', default='0')
parser.add_argument('-del', '--exclude_boxes', help='Exclude the following boxes (numbers only)')
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

for box in boxes_h5_list:
    h5out = '{box}/final_merge_{n}.h5'.format(box=box, n=str(args.archive))
    os.system('rm '+h5out)
    tecandphase1 = sorted(glob('{box}/tecandphase1*_*.ms.archive{n}*h5'.format(box=box, n=args.archive)))
    tecandphase0 = sorted(glob('{box}/tecandphase0*_*.ms.archive{n}*h5'.format(box=box, n=args.archive)))
    scalarcomplexgain2 = sorted(glob('{box}/scalarcomplexgain2*_*.ms.archive{n}*h5'.format(box=box, n=args.archive)))
    h5_files = []
    if len(tecandphase1)>1:
        h5_files.append(tecandphase1[-1])
    if len(tecandphase0)>1:
        h5_files.append(tecandphase0[-1])
    if len(scalarcomplexgain2)>1:
        h5_files.append(scalarcomplexgain2[-1])
    merge_h5(h5_out=h5out,
             h5_tables=h5_files,
             h5_time_freq=sorted(glob('{box}/merged_selfcalcyle*_*.ms.archive{n}*h5'.format(box=boxes_h5_list[0], n=args.archive)))[-1],
             merge_all_in_one=True)
    T = tables.open_file(h5out, 'r+')
    T.root.sol000.source._f_remove()
    alt_files = []
    if len(tecandphase1)>1:
        openfile = tecandphase1[0]
    elif len(tecandphase0)>1:
        openfile = tecandphase0[0]
    elif len(scalarcomplexgain2)>1:
        openfile = scalarcomplexgain2[0]
    else:
        openfile=''
    if openfile!='':
        H = tables.open_file(openfile, 'r+')
        new_source = np.array(T.root.sol000.source[:], dtype=[('name', 'S128'), ('dir', '<f4', (2,))])
        H.close()
        T.create_table(T.root.sol000, 'source', new_source, "Source names and directions")
    T.close()

merge_h5(h5_out='all_directions{n}.h5'.format(n=str(args.archive)),
         h5_tables=sorted(glob('box_*/final_merge_{n}.h5'.format(n=str(args.archive)))))

# final_boxes = []
# for box in boxes_h5_list:
#     try:
#         final_boxes.append(sorted(glob('{box}/merged_selfcal*.ms.archive{n}*h5'.format(box=box, n=str(args.archive))))[-1])
#     except:
#         print('Issues with finding:')
#         print('{box}/merged_selfcal*.ms.archive{n}*h5'.format(box=box, n=str(args.archive)))
#
# merge_h5(h5_out='all_directions{n}_wrong.h5'.format(n=str(args.archive)),
#          h5_tables=final_boxes)