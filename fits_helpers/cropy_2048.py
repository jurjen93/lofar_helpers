"""
Crop image to 2048x2048 (for neural network)
"""

from astropy.io import fits
import numpy as np
from argparse import ArgumentParser


def crop_fits_image(input_filename, output_filename, center=None):
    """
    Crops a FITS image to 2048x2048 pixels.

    Parameters:
    - input_filename: str, path to the input FITS file.
    - output_filename: str, path to the output cropped FITS file.
    - center: tuple (x_center, y_center), the pixel coordinates of the center of the crop.
              If None, the image will be cropped from the center of the input image.
    """
    # Open the FITS file
    with fits.open(input_filename) as hdul:
        data = hdul[0].data.squeeze()
        header = hdul[0].header

        # Get the dimensions of the image
        y_size, x_size = data.shape

        # Determine the center of the crop
        if center is None:
            x_center = x_size // 2
            y_center = y_size // 2
        else:
            x_center, y_center = center

        # Calculate the starting and ending indices for the crop
        x_start = int(x_center - 1024)
        x_end = int(x_center + 1024)
        y_start = int(y_center - 1024)
        y_end = int(y_center + 1024)

        # Ensure the indices are within the bounds of the image
        x_start = max(0, x_start)
        y_start = max(0, y_start)
        x_end = min(x_size, x_end)
        y_end = min(y_size, y_end)

        # Crop the image data
        cropped_data = data[y_start:y_end, x_start:x_end]

        # Update header to reflect new image size
        header['NAXIS1'] = cropped_data.shape[1]
        header['NAXIS2'] = cropped_data.shape[0]

        # Adjust the reference pixel if present
        if 'CRPIX1' in header:
            header['CRPIX1'] -= x_start
        if 'CRPIX2' in header:
            header['CRPIX2'] -= y_start

        # Write the cropped image to the output file
        hdu = fits.PrimaryHDU(data=np.array([[cropped_data]]), header=header)
        hdu.writeto(output_filename, overwrite=True)


def parse_args():
    """
    Command line argument parser
    :return: parsed arguments
    """
    parser = ArgumentParser(description='Crop image to 2048x2048')
    parser.add_argument('--fits_input', help='fits input file', required=True, type=str)
    parser.add_argument('--fits_output', help='fits output file', required=True, type=str)
    return parser.parse_args()

def main():
    """ Main function"""
    args = parse_args()
    crop_fits_image(args.fits_input, args.fits_output)

if __name__ == '__main__':
    main()