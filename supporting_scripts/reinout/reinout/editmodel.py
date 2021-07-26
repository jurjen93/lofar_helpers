#!/usr/bin/env python
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
    print("Please pass all filenames (e.g.: 'editmodel.py model.fits inner.reg outer.reg')")
    sys.exit()
    
try:
    inner_region = sys.argv[2]
except IndexError:
    print("Please pass all filenames (e.g.: 'editmodel.py model.fits inner.reg outer.reg')")
    sys.exit()
    
try:
    outer_region = sys.argv[3]
except IndexError:
    print("Please pass all filenames (e.g.: 'editmodel.py model.fits inner.reg outer.reg')")
    sys.exit()
    
try:
    out_filename = sys.argv[4]
except IndexError:
    out_filename = sys.argv[1]
    
    
#Open file
hdu=fits.open(fits_filename)
hduflat = flatten(hdu)

hdu_2=fits.open(fits_filename)
hduflat_2 = flatten(hdu_2)

#Create masks
r = pyregion.open(outer_region)
manualmask_outer = r.get_mask(hdu=hduflat)

r = pyregion.open(inner_region)
manualmask_inner = r.get_mask(hdu=hduflat)

#Set everything within the outer region to 0
hdu[0].data[0][0][np.where(manualmask_outer == True)] = 0.0

#Set everything within the inner region back to original value
hdu[0].data[0][0][np.where(manualmask_inner == True)] = hdu_2[0].data[0][0][np.where(manualmask_inner == True)]

#Store the new model
hdu.writeto(out_filename,overwrite=True)