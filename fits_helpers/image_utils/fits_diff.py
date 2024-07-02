from astropy.io import fits

def get_fits_diff(fits_file_1, fits_file_2, out_name='diff.fits'):
    """
    Get subtracted difference between two equally sized fits files
    """

    f1 = fits.open(fits_file_1)[0]
    f2 = fits.open(fits_file_2)[0]

    hdu = fits.PrimaryHDU(header=f1.header, data=f1.data - f2.data)
    hdu.writeto(out_name, overwrite=True)
