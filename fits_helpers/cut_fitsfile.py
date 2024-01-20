"""With this script you can cut a region from a fits file"""

import pyregion
from astropy.io import fits
from numpy import nan
from argparse import ArgumentParser


def parse_args():
    """
    Command line argument parser
    :return: parsed arguments
    """
    parser = ArgumentParser(description='Cut fits file with region file')
    parser.add_argument('--fits_input', help='fits input file', required=True, type=str)
    parser.add_argument('--fits_output', help='fits output file', required=True, type=str)
    parser.add_argument('--region', help='region file', required=True, type=str)
    return parser.parse_args()


def main():
    """ Main function"""
    args = parse_args()

    fitsfile = args.fits_input
    regionfile = args.region
    outputfits = args.fits_output

    hdu = fits.open(fitsfile)

    header = hdu[0].header

    r = pyregion.open(regionfile).as_imagecoord(header=header)
    mask = r.get_mask(hdu=hdu[0], shape=(header["NAXIS1"], header["NAXIS2"])).astype(int)
    imagedata = hdu[0].data.reshape(header["NAXIS1"], header["NAXIS2"])
    imagedata *= mask
    imagedata[imagedata == 0] = nan

    hdu = fits.PrimaryHDU(header=header, data=imagedata)
    hdu.writeto(outputfits, overwrite=True)


if __name__ == '__main__':
    main()
