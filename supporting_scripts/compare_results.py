import pyregion
from astropy.io import fits
import numpy as np
from astropy.wcs import WCS
import sys
from argparse import ArgumentParser
from glob import glob
from astropy.nddata import Cutout2D


def make_cutout(fitsfile: str = None, region: str = None ):

    hdu = fits.open(fitsfile)
    image_data = hdu[0].data
    header = hdu.header
    hdu.close()

    r = pyregion.open(region)
    ra  = r[0].coord_list[0]
    dec = r[0].coord_list[1]
    sizex = r[0].coord_list[2]
    sizey = r[0].coord_list[3]

    wcs = WCS(header, naxis=2)
    out = Cutout2D(
        data=image_data,
        position=(ra, dec),
        size=(sizex, sizey),
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
    parser.add_argument('--region_files', type=str, help='region files location', required=True)
    args = parser.parse_args()

    surf_images = glob(args.surf_im)
    region_files = glob(args.region_files)
    int_restored = args.int_restored
    image_final = args.image_final
    app_restored = args.app_restored

    for reg in region_files:
        print(image_final)
        make_cutout(image_final, reg)
        break