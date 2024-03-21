"""
IN CASA:

PATH=path_to_fits_files
targetres = True/False

execfile("casa_deconvolve.py")
"""

from glob import glob
from numpy import mean, sqrt
import os

# path to fits files
PATH='.'

def calc_conv_beam(fwhm_orig_major, fwhm_orig_minor, fwhm_new_major, fwhm_new_minor):
    """
    Calculate the FWHM of the convolution beam needed to convolve an image from
    an original beam to a new beam.

    Parameters:
    - fwhm_orig_major: FWHM of the original beam's major axis (arcseconds)
    - fwhm_orig_minor: FWHM of the original beam's minor axis (arcseconds)
    - fwhm_new_major: FWHM of the desired new beam's major axis (arcseconds)
    - fwhm_new_minor: FWHM of the desired new beam's minor axis (arcseconds)

    Returns:
    - fwhm_conv_major: FWHM of the convolution beam's major axis (arcseconds)
    - fwhm_conv_minor: FWHM of the convolution beam's minor axis (arcseconds)
    """
    # Calculate squared FWHM values for the original and new beams
    fwhm_orig_major_squared = fwhm_orig_major ** 2
    fwhm_orig_minor_squared = fwhm_orig_minor ** 2
    fwhm_new_major_squared = fwhm_new_major ** 2
    fwhm_new_minor_squared = fwhm_new_minor ** 2

    # Calculate squared FWHM for the convolution beam (ensure it's non-negative)
    fwhm_conv_major_squared = max(fwhm_new_major_squared - fwhm_orig_major_squared, 0)
    fwhm_conv_minor_squared = max(fwhm_new_minor_squared - fwhm_orig_minor_squared, 0)

    # Calculate the FWHM for the convolution beam
    fwhm_conv_major = sqrt(fwhm_conv_major_squared)
    fwhm_conv_minor = sqrt(fwhm_conv_minor_squared)

    return fwhm_conv_major, fwhm_conv_minor

beam_major = []
beam_minor = []
beam_pa = []

for fts in glob(PATH+'/*.fits'):
    beam_major.append(imhead(fts, mode='get', hdkey='beammajor')['value'])
    beam_minor.append(imhead(fts, mode='get', hdkey='beamminor')['value'])
    beam_pa.append(imhead(fts, mode='get', hdkey='beampa')['value'])

final_beam_major = max(beam_major)+0.001
final_beam_minor = max(beam_minor)+0.001
final_beam_pa = mean(beam_pa)

print('Final beam major axis: '+str(final_beam_major))
print('Final beam minor axis: '+str(final_beam_minor))
print('Final beam pa: '+str(final_beam_pa))

for fts in glob(PATH+'/*.fits'):
    print('Smooth '+fts)
    # b_maj, b_min = calc_conv_beam(imhead(fts, mode='get', hdkey='beammajor')['value'],
    #                               imhead(fts, mode='get', hdkey='beamminor')['value'],
    #                               final_beam_major,
    #                               final_beam_minor)
    # print('Smooth with '+str(round(b_maj, 3))+' X '+str(round(b_min, 3)))
    imsmooth(imagename=fts,
             kernel='gauss',
             major=str(final_beam_major)+'arcsec',
             minor=str(final_beam_minor)+'arcsec',
             targetres=True,
             pa=str(final_beam_pa)+'deg',
             outfile=fts.replace('.fits', '')+'_newbeam.tmp')
    exportfits(imagename=fts.replace('.fits', '')+'_newbeam.tmp',
               fitsimage=fts.replace('.fits', '')+'_newbeam.fits')
    os.system('rm -rf '+fts.replace('.fits', '')+'_newbeam.tmp')
    print('Created '+fts.replace('.fits', '')+'_newbeam.fits')
