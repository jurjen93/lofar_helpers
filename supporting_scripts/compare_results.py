import pyregion
from astropy.io import fits
import numpy as np
from astropy.wcs import WCS
import sys
from argparse import ArgumentParser
from glob import glob
from astropy.nddata import Cutout2D

def flatten(f):
    """ Flatten a fits file so that it becomes a 2D image. Return new header and data """

    naxis = f[0].header['NAXIS']
    if naxis < 2:
        sys.exit('Can\'t make map from this')
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


def mask_region(infilename,ds9region):

    hdu=fits.open(infilename)
    hduflat = flatten(hdu)

    r = pyregion.open(ds9region)
    manualmask = r.get_mask(hdu=hduflat)
    # hdu[0].data[0][0][np.where(manualmask == True)] = 0.0
    # hdu.writeto(outfilename,overwrite=True)

    print(manualmask)

def make_cutout(image_data, pos: tuple = None, size: tuple = (1000, 1000), hdu: str = None):
    """
    Make cutout from your image with this method.
    pos (tuple) -> position in pixels
    size (tuple) -> size of your image in pixel size, default=(1000,1000)
    """
    wcs = WCS(hdu.header, naxis=2)
    out = Cutout2D(
        data=image_data,
        position=pos,
        size=size,
        wcs=wcs,
        mode='partial'
    )

    return out.data, out.wcs


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--surf_im', type=str, help='fits images produced on surfsara', required=True)
    parser.add_argument('--int_restored', type=str, nargs='+', help='*int.restored.fits file', required=True)
    parser.add_argument('--image_final', type=str, nargs='+', help='final image', required=True)
    parser.add_argument('--app_restored', type=str, nargs='+', help='*app.restored.fits image', required=True)
    parser.add_argument('--region_files', type=str, nargs='+', help='region files location', required=True)
    args = parser.parse_args()

    print(args.surf_im)
    print(args.region_files)
    surf_images = glob(args.surf_im)
    region_files = glob(args.region_files)
    int_restored = args.int_restored
    image_final = args.image_final
    app_restored = args.app_restored

    for reg in region_files:
        mask_region(image_final, region_files)
        break