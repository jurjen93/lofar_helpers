from astropy.io import fits
from radio_beam import Beam
from astropy.convolution import convolve_fft, Gaussian2DKernel
from argparse import ArgumentParser
import astropy.units as u

def parse_args():
    """
    Command line argument parser

    :return: parsed arguments
    """

    parser = ArgumentParser(description='Convolve radio map with other beam')
    parser.add_argument('--maj', type=str, help='beam size major axis (arcsec)')
    parser.add_argument('--min', type=str, help='beam size minor axis (arcsec)')
    parser.add_argument('--fits_in', type=str, help='fits image input name')
    parser.add_argument('--fits_out', type=str, help='fits image output name')
    return parser.parse_args()

def main():

    args = parse_args()

    # Load your FITS file
    fits_file = args.fits
    hdul = fits.open(fits_file)

    # Get data and header
    map = hdul[0].data
    header = hdul[0].header

    maj = args.maj
    min = args.min

    # Calculate the kernel size needed for the convolution
    target_beam = Beam(major=maj, minor=min, pa=0, default_unit=u.arcsec)
    kernel = Gaussian2DKernel(target_beam.as_kernel(header['CDELT2']*u.degree))

    # Convolve the map
    convolved_map = convolve_fft(map, kernel)

    # Update the header or modify as needed
    header['BMIN'] = min/3600 #degrees
    header['BMAJ'] = maj/3600 #degrees
    header['HISTORY'] = 'Convolved to {}x{} arcsec resolution'.format(maj, min)

    # Create a new HDU
    hdu = fits.PrimaryHDU(convolved_map, header)

    # Save to a new file
    hdu.writeto('convolved_radio_map.fits', overwrite=True)
