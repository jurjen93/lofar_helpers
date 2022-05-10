import numpy as np
import os
from astropy.wcs import WCS
from astropy.io import fits
from argparse import ArgumentParser
from glob import glob
import re
import sys


parser = ArgumentParser()
parser.add_argument('-f', '--fits', type=str, help='fitsfile name')
parser.add_argument('--path', type=str, help='data path', default='.')
parser.add_argument('-in', '--include_boxes', help='Include only the following boxes (numbers only)')
parser.add_argument('--imager', default='DDFACET')
args = parser.parse_args()

def filter_box_N(boxnumber):
    return bool(re.match('^box_[0-9]+$', boxnumber))

fits_file = args.fits

boxes = [b for b in glob(args.path+'/box_*') if filter_box_N(b.split('/')[-1])]

if args.include_boxes:
    boxes = [b for b in boxes if b.split('/')[-1] in ['box_'+n for n in args.include_boxes.split(',')]]

os.system('cp ' + fits_file + ' ' + fits_file.replace('.fits', '_new.fits') + ' && wait')
fits_file = fits_file.replace('.fits', '_new.fits')

hdu = fits.open(fits_file, mode='update')[0]
hdu.data[0, 0, :, :] = np.zeros(hdu.data.shape)
wcs = WCS(hdu.header, naxis=2)
for b in boxes:
    if args.imager.lower()=='ddfacet':
        im = sorted(glob(b+'/image_00*.app.restored.fits'))[-1]
    elif args.imager.lower()=='wsclean':
        im = sorted(glob(b + '/image_00*-MFS-image.fits'))[-1]
    else:
        sys.exit('ERROR: Choose --imager=WSCLEAN or --imager=DDFACET')
    print('Merge:\n'+im)
    hdu_box = fits.open(im)[0]
    box_data = hdu_box.data[0][0]
    wcs_box = WCS(hdu_box.header, naxis=2)
    x, y = wcs.world_to_pixel(wcs_box.pixel_to_world(0, 0))
    shape = box_data.shape
    imdata = hdu.data[0, 0, :, :]
    imdata[int(y):int(y+shape[1]), int(x):int(x+shape[0])] = box_data
    hdu.data[0, 0, :, :] = imdata


hdu.writeto(fits_file, overwrite=True)