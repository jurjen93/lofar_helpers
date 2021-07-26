#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import astropy
import sys
from astropy.io import fits
from astropy.wcs import WCS
import pyregion


def flatten(f):
    """ Flatten a fits file so that it becomes a 2D image. Return new header and data """

    naxis=f[0].header['NAXIS']
    if naxis<2:
        raise RadioError('Can\'t make map from this')
    if naxis==2:
        return fits.PrimaryHDU(header=f[0].header,data=f[0].data)

    w = WCS(f[0].header)
    wn=WCS(naxis=2)
    
    wn.wcs.crpix[0]=w.wcs.crpix[0]
    wn.wcs.crpix[1]=w.wcs.crpix[1]
    wn.wcs.cdelt=w.wcs.cdelt[0:2]
    wn.wcs.crval=w.wcs.crval[0:2]
    wn.wcs.ctype[0]=w.wcs.ctype[0]
    wn.wcs.ctype[1]=w.wcs.ctype[1]
    
    header = wn.to_header()
    header["NAXIS"]=2
    copy=('EQUINOX','EPOCH','BMAJ', 'BMIN', 'BPA', 'RESTFRQ', 'TELESCOP', 'OBSERVER')
    for k in copy:
        r=f[0].header.get(k)
        if r is not None:
            header[k]=r

    slice=[]
    for i in range(naxis,0,-1):
        if i<=2:
            slice.append(np.s_[:],)
        else:
            slice.append(0)
        
    hdu = fits.PrimaryHDU(header=header,data=f[0].data[tuple(slice)])
    return hdu

#Check input
try:
    fits_filename = sys.argv[1]
except IndexError:
    print("Please pass all filenames (e.g.: 'insert_highres.py input.fits model.fits region.reg')")
    sys.exit()
    
try:
    model_filename = sys.argv[2]
except IndexError:
    print("Please pass all filenames (e.g.: 'insert_highres.py input.fits model.fits region.reg')")
    sys.exit()
    
try:
    region_filename = sys.argv[3]
except IndexError:
    print("Please pass all filenames (e.g.: 'insert_highres.py input.fits model.fits region.reg')")
    sys.exit()
    
try:
    out_filename = sys.argv[4]
except IndexError:
    out_filename = sys.argv[1]
    
    
#Open files
hdu_i=fits.open(fits_filename)
hduflat_i = flatten(hdu_i)

hdu_m=fits.open(model_filename)

#Get array dimensions
input_size = hdu_i[0].data[0][0].shape
model_size = hdu_m[0].data[0][0].shape

#Match model to input size
model_array = np.zeros(input_size)
idx_x = int(input_size[0]/2-model_size[0]/2)
idx_y = int(input_size[1]/2-model_size[1]/2)
model_array[idx_x:idx_x+model_size[0], idx_y:idx_y+model_size[1]] = hdu_m[0].data[0][0]

#Create mask
r = pyregion.open(region_filename)
mask = r.get_mask(hdu=hduflat_i)

#Set everything within the region to the model values
hdu_i[0].data[0][0][np.where(mask == True)] = model_array[np.where(mask == True)]

#Store the new model
hdu_i.writeto(out_filename,overwrite=True)