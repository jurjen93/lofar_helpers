"""With this script you can make a cutout with a ds9 region file"""

import pyregion
from astropy.io import fits
from argparse import ArgumentParser
import numpy as np


def crop_nan_boundaries(fits_in, fits_out):
    """
    Crop nan boundaries

    input:
        - fits_in: input fits file
        - fits_out: output fits file
    """

    with fits.open(fits_in) as hdul:
        image_data = hdul[0].data
        header = hdul[0].header

    mask = ~np.isnan(image_data)
    non_nan_indices = np.where(mask)

    ymin, ymax = non_nan_indices[0].min(), non_nan_indices[0].max()
    xmin, xmax = non_nan_indices[1].min(), non_nan_indices[1].max()

    print(f"Original shape {image_data.shape}")
    print(ymin, ymax)
    print(xmin, xmax)

    cropped_image = image_data[ymin:ymax + 1, xmin:xmax + 1]

    header['NAXIS1'] = cropped_image.shape[1]
    header['NAXIS2'] = cropped_image.shape[0]
    header['CRPIX1'] -= xmin
    header['CRPIX2'] -= ymin

    print(f"New shape {cropped_image.shape}")

    hdu = fits.PrimaryHDU(cropped_image, header=header)
    hdu.writeto(fits_out, overwrite=True)


def parse_args():
    """
    Command line argument parser
    :return: parsed arguments
    """
    parser = ArgumentParser(description='Cut fits file with region file')
    parser.add_argument('fits_input', help='fits input file', type=str)
    parser.add_argument('--output_name', help='fits output file', type=str)
    parser.add_argument('--region', help='region file', required=True, type=str)
    return parser.parse_args()


def main():
    """ Main function"""
    args = parse_args()

    fitsfile = args.fits_input
    regionfile = args.region
    outputfits = args.output_name
    if outputfits is None:
        outputfits = fitsfile

    hdu = fits.open(fitsfile)

    header = hdu[0].header

    r = pyregion.open(regionfile).as_imagecoord(header=header)
    mask = r.get_mask(hdu=hdu[0], shape=(header["NAXIS2"], header["NAXIS1"])).astype(int)
    imagedata = hdu[0].data.reshape(header["NAXIS2"], header["NAXIS1"])
    imagedata *= mask
    imagedata[imagedata == 0] = np.nan

    hdu = fits.PrimaryHDU(header=header, data=imagedata)
    hdu.writeto(outputfits, overwrite=True)

    crop_nan_boundaries(outputfits, outputfits)


if __name__ == '__main__':
    main()
