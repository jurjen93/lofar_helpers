import os
from astropy.io import fits
import pyregion
import numpy as np
from astropy.wcs import WCS

def flatten(f):
    """ Flatten a fits file so that it becomes a 2D image. Return new header and data """

    naxis = f[0].header['NAXIS']
    # if naxis < 2:
    #     raise RadioError('Can\'t make map from this')
    if naxis == 2:
        return fits.PrimaryHDU(header=f[0].header, data=f[0].data)

    w = WCS(f[0].header)
    wn = WCS(naxis=2)

    wn.wcs.crpix[0] = w.wcs.crpix[0]
    wn.wcs.crpix[1] = w.wcs.crpix[1]
    wn.wcs.cdelt = w.wcs.cdelt[0:2]
    wn.wcs.crval = w.wcs.crval[0:2]
    wn.wcs.ctype[0] = w.wcs.ctype[0]
    wn.wcs.ctype[1] = w.wcs.ctype[1]

    header = wn.to_header()
    header["NAXIS"] = 2
    copy = ('EQUINOX', 'EPOCH', 'BMAJ', 'BMIN', 'BPA', 'RESTFRQ', 'TELESCOP', 'OBSERVER')
    for k in copy:
        r = f[0].header.get(k)
        if r is not None:
            header[k] = r

    slice = []
    for i in range(naxis, 0, -1):
        if i <= 2:
            slice.append(np.s_[:], )
        else:
            slice.append(0)

    hdu = fits.PrimaryHDU(header=header, data=f[0].data[tuple(slice)])
    return hdu

def mask_except_region(infilename,ds9region,outfilename):

    hdu=fits.open(infilename)
    hduflat = flatten(hdu)
    map=hdu[0].data

    r = pyregion.open(ds9region)
    manualmask = r.get_mask(hdu=hduflat)
    hdu[0].data[0][0][np.where(manualmask == True)] = 0.0
    hdu.writeto(outfilename,overwrite=True)

    return

mask='dicoMask.fits'

mask_except_region('image_full_ampphase_di_m.NS.mask01.fits', 'extracted.reg', mask)

indico='image_full_ampphase_di_m.NS.DicoModel'
outdico='extr.DicoModel'

os.system("MaskDicoModel.py --MaskName=%s --InDicoModel=%s --OutDicoModel=%s"%(mask, indico, outdico))