from argparse import ArgumentParser
from astropy.io import fits
import pyregion

__author__ = "Jurjen de Jong"

def parse_args():
    """
    Command line argument parser
    :return: parsed arguments
    """
    parser = ArgumentParser(description='Mask fits file with region file')
    parser.add_argument('fits_input', help='fits input file', type=str)
    parser.add_argument('--output_name', help='fits output file', type=str)
    parser.add_argument('--region', help='region file', required=True, type=str)
    return parser.parse_args()


def main():
    """ Main function"""
    args = parse_args()

    hdu = fits.open(args.fits_input)
    header = hdu[0].header

    r = pyregion.open(args.region).as_imagecoord(header=header)
    mask = r.get_mask(hdu=hdu[0], shape=(header["NAXIS2"], header["NAXIS1"])).astype(int)

    imagedata = hdu[0].data.reshape(header["NAXIS2"], header["NAXIS1"])
    imagedata[mask] = 0

    hdu = fits.PrimaryHDU(header=header, data=imagedata)
    hdu.writeto(args.output_name, overwrite=True)


if __name__ == '__main__':
    main()
