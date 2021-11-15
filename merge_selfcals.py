"""
Merge selfcals together.
We do an extra check for coordinates, as we have seen wrong coordinates.
We therefore make a final_merge_[N].h5 that we will merge.
"""

__author__ = "Jurjen de Jong (jurjendejong@strw.leidenuniv.nl)"

from h5_merger import merge_h5
from glob import glob
from argparse import ArgumentParser
import tables
import numpy as np
import os
import sys

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

merged_boxes = []

for box in boxes_h5_list:
    print(box)
    h5out = '{box}/final_merge_{n}.h5'.format(box=box, n=str(args.archive))
    h5merge = sorted(glob('{box}/merged_selfcalcyle*_*.ms.archive{n}*h5'.format(box=box, n=args.archive)))[-1]
    merged_boxes.append(h5merge)
    os.system('rm '+h5out)
    os.system('cp '+h5merge+' '+h5out)

    tecandphase1 = sorted(glob('{box}/tecandphase1*_*.ms.archive{n}*h5'.format(box=box, n=args.archive)))
    tecandphase0 = sorted(glob('{box}/tecandphase0*_*.ms.archive{n}*h5'.format(box=box, n=args.archive)))
    scalarcomplexgain2 = sorted(glob('{box}/scalarcomplexgain2*_*.ms.archive{n}*h5'.format(box=box, n=args.archive)))

    T = tables.open_file(h5out, 'r+')
    T.root.sol000.source._f_remove()
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
        new_source = np.array(H.root.sol000.source[:], dtype=[('name', 'S128'), ('dir', '<f4', (2,))])
        H.close()
        T.create_table(T.root.sol000, 'source', new_source, "Source names and directions")
    T.close()

final_merge = glob('box_*/final_merge_{n}.h5'.format(n=str(args.archive)))
final_merge.sort(key=lambda x: get_digits(x))
merge_h5(h5_out='all_directions{n}.h5'.format(n=str(args.archive)),
         h5_tables=final_merge)

merged_boxes.sort(key=lambda x: get_digits(x))
H = tables.open_file('all_directions{n}.h5'.format(n=str(args.archive)))
for n, h5 in enumerate(merged_boxes):
    T = tables.open_file(h5)
    if H.root.sol000.phase000.val[0,0,0,n,0]!=T.root.sol000.phase000.val[0,0,0,0,0] or \
        H.root.sol000.amplitude000.val[0,0,0,n,0]!=T.root.sol000.amplitude000.val[0,0,0,0,0]:
        sys.exit('check direction '+str(n))
    T.close()
H.close()