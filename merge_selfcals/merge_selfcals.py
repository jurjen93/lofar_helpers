"""
Merge selfcals together.
We do an extra check for coordinates, as we have seen wrong coordinates.
We therefore make a final_merge_[N].h5 that we will merge.
"""

__author__ = "Jurjen de Jong (jurjendejong@strw.leidenuniv.nl)"


from glob import glob
from argparse import ArgumentParser
import tables
import numpy as np
import os
import sys
import re
sys.path.append('/home/lofarvwf-jdejong/scripts/lofar_helpers')
from h5_merger import merge_h5

parser = ArgumentParser()
parser.add_argument('-d', '--directory', type=str, help='directory path')
parser.add_argument('-a', '--archive', type=str, help='archive number', default='0')
parser.add_argument('-del', '--exclude_boxes', help='Exclude the following boxes (numbers only)')
parser.add_argument('-in', '--include_boxes', help='Include only the following boxes (numbers only)')
args = parser.parse_args()

def get_digits(x):
    return int(''.join([d for d in x if d.isdigit()]))

def filter_box_N(boxnumber):
    return bool(re.match('^box_[0-9]+$', boxnumber))

excluded_boxes, included_boxes = [], []

if args.exclude_boxes and args.include_boxes:
    sys.exit('You need to choose between --exclude_boxes or --include_boxes')
elif args.exclude_boxes:
    excluded_boxes=args.exclude_boxes
    if type(excluded_boxes)==str:
        excluded_boxes = excluded_boxes.split(',')
    excluded_boxes = ['box_'+n for n in excluded_boxes]
elif args.include_boxes:
    included_boxes=args.include_boxes
    if type(included_boxes)==str:
        included_boxes = included_boxes.split(',')
    included_boxes = ['box_'+n for n in included_boxes]

boxes_h5_list = [b for b in glob('{directory}/box_*'.format(directory=args.directory)) if filter_box_N(b.split('/')[-1])]

if included_boxes:
    boxes_h5_list = [b for b in boxes_h5_list if b.split('/')[-1] in included_boxes]
elif excluded_boxes:
    boxes_h5_list = [b for b in boxes_h5_list if not b.split('/')[-1] in excluded_boxes]

boxes_h5_list.sort(key=lambda x: get_digits(x))

merged_boxes = []
final_merge = []
for box in boxes_h5_list:
    print(box)
    h5out = '{box}/final_merge_{n}.h5'.format(box=box, n=str(args.archive))
    h5merge = sorted(glob('{box}/merged_selfcalcyle*_*.ms.archive{n}*h5'.format(box=box, n=args.archive)))[-1]
    merged_boxes.append(h5merge)
    final_merge.append(h5out)
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

final_merge.sort(key=lambda x: get_digits(x))
merge_h5(h5_out='all_directions{n}.h5'.format(n=str(args.archive)),
         h5_tables=final_merge)

merged_boxes.sort(key=lambda x: get_digits(x))
H = tables.open_file('all_directions{n}.h5'.format(n=str(args.archive)))
for n, h5 in enumerate(merged_boxes):
    T = tables.open_file(h5)
    num_elements = 1
    for i in T.root.sol000.phase000.val.shape:
        num_elements*=i
    phase_elements = H.root.sol000.phase000.val[:,:,:,n,:]==T.root.sol000.phase000.val[:,:,:,0,:]
    amplitude_elements = H.root.sol000.amplitude000.val[:,:,:,n,:]==T.root.sol000.amplitude000.val[:,:,:,0,:]
    if np.sum(phase_elements)!=num_elements:
        print('Problems with following phase elements in direction '+str(n)+':\n'+str(np.argwhere(np.invert(phase_elements))))
        sys.exit('ERROR: CHECK DIRECTION '+str(n))
    if np.sum(amplitude_elements) != num_elements:
        print('Problems with following phase elements '+str(n)+':\n' + str(np.argwhere(np.invert(amplitude_elements))))
        sys.exit('ERROR: CHECK DIRECTION '+str(n))

    T.close()
H.close()