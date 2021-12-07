import numpy as np
import os
from astropy.wcs import WCS
from astropy.io import fits
from argparse import ArgumentParser
from glob import glob
import re


parser = ArgumentParser()
parser.add_argument('-f', '--fits', type=str, help='fitsfile name')
parser.add_argument('--path', type=str, help='data path', default='.')
args = parser.parse_args()

def filter_box_N(boxnumber):
    return bool(re.match('^box_[0-9]+$', boxnumber))

fits_file = args.fits
boxes = [b for b in glob(args.path+'/box_*') if filter_box_N(b.split('/')[-1])]

os.system('cp ' + fits_file + ' ' + fits_file.replace('.fits', '_new.fits') + ' && wait')
fits_file = fits_file.replace('.fits', '_new.fits')

hdu = fits.open(fits_file, mode='update')[0]
hdu.data[0, 0, :, :] = np.zeros(hdu.data.shape)
wcs = WCS(hdu.header, naxis=2)
for b in boxes:
    im = sorted(glob(b+'/image_00*.app.restored.fits'))[-1]
    hdu_box = fits.open(im)[0]
    box_data = hdu_box.data[0][0]
    wcs_box = WCS(hdu_box.header, naxis=2)
    x, y = wcs.world_to_pixel(wcs_box.pixel_to_world(0, 0))
    shape = box_data.shape
    imdata = hdu.data[0, 0, :, :]
    imdata[int(y+shape[1]/10):int(y+0.9*shape[1]), int(x+shape[0]/10):int(x+0.9*shape[0])] = \
        box_data[int(box_data.shape[1]/10):int(box_data.shape[1]*0.9), int(box_data.shape[0]/10):int(box_data.shape[0]*0.9)]
    hdu.data[0, 0, :, :] = imdata


hdu.writeto(fits_file, overwrite=True)