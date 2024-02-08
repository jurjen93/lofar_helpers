"""
Use this script if you want to do a primary beam correction with wsclean in case you forgot to turn on the primary beam
settings and don't want to run the entire imaging.

Run the script in your wsclean output folder.
"""

from astropy.io import fits
from glob import glob
import numpy as np
import sys
import os
import re


def has_num(string):
    """
    Verify if string has number
    """
    return bool(re.search(r'\d', string))


def has_char(string):
    """
    Verify if string has letter
    """
    return bool(re.search('[a-zA-Z]', string))


def split_num_char(string):
    """
    split characters and numbers from string
    """
    return ((' '.join(re.split('(\d+)', string)).strip().
            replace(' .','.').replace('. ','.')).
            replace(' _','_').replace('_ ','_'))


def make_beam_images(cmd):
    """
    Make beam images with wsclean, by running 1 major and minor cycle

    :param cmd: wsclean command line from history
    """

    hist_split = cmd.split()

    for idx, element in enumerate(hist_split):
        if element=='-nmiter':
            cmd = cmd.replace('-nmiter ' + hist_split[idx + 1], '-nmiter 0 -no-reorder -no-dirty -apply-primary-beam -use-differential-lofar-beam')
        if element=='-niter':
            cmd = cmd.replace('-niter ' + hist_split[idx + 1], '-niter 0')
        if element=='-name':
            cmd = cmd.replace('-name ' + hist_split[idx + 1], '-name beam')
        if element=='-baseline-averaging':
            cmd = cmd.replace('-baseline-averaging ' + hist_split[idx + 1], '')
        # if '/imaging/' in element and 'split_facets2' in element:
        #     cmd = cmd.replace(element, element.split('/imaging/')[-1])
        if element=='-':
            cmd = cmd.replace(' - -', ' -')
        if has_char(element) and has_num(element) and '.ms' not in element and '-' in element \
                and hist_split[idx-1]!='-taper-gaussian' and hist_split[idx-1]!='-scale':
            cmd = cmd.replace(element, split_num_char(element))
        if element[-3:]=='.ms':
            cmd = cmd.replace(element, '')
    cmd = cmd.replace('- ', ' ').replace('  ', ' ').strip()
    cmd += ' *.ms'

    if 'facet_14-scale' in hist_split:
        cmd = 'wsclean -gridder wgridder -no-update-model-required -minuv-l 80.0 -size 45000 45000 -weighting-rank-filter 3 '\
              '-reorder -weight briggs -1.5 -parallel-reordering 4 -mgain 0.75 -data-column DATA -auto-mask 2.5 '\
              '-auto-threshold 1.0 -pol i -name beam -scale 0.1arcsec -niter 0 -multiscale-scale-bias 0.6'\
              ' -parallel-deconvolution 2600 -multiscale -multiscale-max-scales 9 -nmiter 0 -no-reorder -no-dirty '\
              '-apply-primary-beam -use-differential-lofar-beam -mem 75 -channels-out 6 -join-channels -fit-spectral-pol 3 '\
              '-local-rms -local-rms-window 50 -taper-gaussian 0.3asec *.ms'
    if 'facet_20-scale' in hist_split:
        cmd = ('wsclean -gridder wgridder -no-update-model-required -minuv-l 80.0 -size 45000 45000 -weighting-rank-filter 3 '
               '-reorder -weight briggs -1.5 -parallel-reordering 2 -mgain 0.75 -data-column DATA -auto-mask 2.5 '
               '-auto-threshold 1.0 -pol i -name beam -scale 0.1arcsec -niter 0 -multiscale-scale-bias 0.6 '
               '-parallel-deconvolution 2600 -multiscale -multiscale-max-scales 9 -nmiter 0 -no-reorder -no-dirty '
               '-apply-primary-beam -use-differential-lofar-beam -mem 30 -channels-out 6 -join-channels -fit-spectral-pol 3 '
               '-local-rms -local-rms-window 50 -taper-gaussian 0.3asec *.ms')
    if 'facet_29' in hist_split:
        cmd = ('wsclean -gridder wgridder -no-update-model-required -minuv-l 80.0 -size 43500 43500 -weighting-rank-filter 3 '
               '-reorder -parallel-reordering 2 -weight briggs -1.5 -mgain 0.75 -data-column DATA -auto-mask 2.5 '
               '-auto-threshold 1.0 -pol i -name beam -scale 0.1arcsec -niter 0 -multiscale-scale-bias 0.6 '
               '-parallel-deconvolution 2600 -multiscale -multiscale-max-scales 9 -nmiter 0 -no-reorder -no-dirty '
               '-apply-primary-beam -use-differential-lofar-beam -parallel-gridding 2 -channels-out 6 -mem 30 -join-channels -fit-spectral-pol 3 '
               '-local-rms -local-rms-window 50 -taper-gaussian 0.3asec -mem 50 *.ms')

    print(cmd)
    if len(glob('apply*.ms'))>0:
        os.system('mkdir -p beamrun && mv apply*.ms beamrun && cd beamrun && ' + cmd + ' > wsclean.txt && cd ../ && mv beamrun/*beam-*.fits .')
    else:
        os.system('mkdir -p beamrun && mv *.ms beamrun && cd beamrun && ' + cmd + ' > wsclean.txt && cd ../ && mv beamrun/*beam-*.fits .')

def get_history(fitsfile):
    """
    Get history information from fits file

    :param fitsfile: fits file input

    :return: history
    """
    history = (str(fitsfile[0].header["HISTORY"]).replace('\n','').
               replace('mssub', 'ms sub')).replace('log-time','')
    for i in range(10):
        history = history.replace('size'+str(i), 'size '+str(i))
    return history


def main():
    """
    Main function
    """

    image_file = glob('*-MFS-image.fits')[0]

    if len(glob(image_file.replace('-image.fits', '-image-pb.fits'))) > 0:
        sys.exit("ERROR: primary beam corrected image already exists, see: "
                 + glob(image_file.replace('-image.fits', '-image-pb.fits'))[0])
    f = fits.open(image_file)

    history = get_history(f)
    make_beam_images(history)

    b0 = fits.open(glob('*beam-0.fits')[0])[0].data
    b15 = fits.open(glob('*beam-15.fits')[0])[0].data

    hdu = fits.PrimaryHDU(header=f[0].header, data=np.divide(f[0].data, np.sqrt(np.multiply(b0, b15))))
    hdu.writeto(image_file.replace('-image.fits', '-image-pb.fits'), overwrite=True)


if __name__ == '__main__':
    main()
