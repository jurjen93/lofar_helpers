#!/usr/bin/env python

# Perform point-to-point analysis (radio/X or radio/radio/X) and save the results in a fits table. If (i) the X-ray counts < 0, (ii) the cell is not totally inside the X-ray FoV, and (ii) the radio flux density is < 0, it saves NaNs.
#
# Written in python3 by Andrea Botteon (botteon@strw.leidenuniv.nl).
#
# version 1.0       2021/05/24: complete
# version 1.1       2021/08/22: alpha_err was not computed correctly, fixed ln_freq21_ratio (np.log instead of np.log10)
#
# To do:
# -Reproject using montage, faster?
# -Calc cellsize using beam info?
# -Analysis in annuli
# -Convolve images at the same resolution
# -Summarize stuff at the end of the script
# -lLogging?
# -Compared to the ciao/bash version the code is slower and there is a difference on evaluation of the pixel factions
# -Make the code parallel with pool?
# -Split sys error and stat error?
# -Remove the reading/deleting print while reading regions
# -Accept XMM images
# -Change name of ptp_dir

import os
import sys
import numpy as np
import glob
import argparse
import re
from past.utils import old_div
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table, Column
from reproject import reproject_interp
import pyregion
import warnings

warnings.filterwarnings('ignore')


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

    hdu = fits.PrimaryHDU(header=header, data=f[0].data[slice])
    return hdu


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def natural_sort(l):
    # from https://stackoverflow.com/questions/4836710/is-there-a-built-in-function-for-string-natural-sort
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def calc_beamarea_fromfits(f):
    # Given a fitsfile this calculates the beamarea in pixels

    hdu = fits.open(f)

    bmaj = hdu[0].header['BMAJ']
    bmin = hdu[0].header['BMIN']

    beammaj = bmaj / (2.0 * (2 * np.log(2)) ** 0.5)  # Convert to sigma
    beammin = bmin / (2.0 * (2 * np.log(2)) ** 0.5)  # Convert to sigma
    pixarea = abs(hdu[0].header['CDELT1'] * hdu[0].header['CDELT2'])

    beamarea = 2 * np.pi * 1.0 * beammaj * beammin  # Note that the volume of a two dimensional gaus$
    beamarea_pix = beamarea / pixarea

    return beamarea_pix

def calc_beamarea(hdu):
    # Given a fitsfile this calculates the beamarea in pixels

    bmaj = hdu[0].header['BMAJ']
    bmin = hdu[0].header['BMIN']

    beammaj = bmaj / (2.0 * (2 * np.log(2)) ** 0.5)  # Convert to sigma
    beammin = bmin / (2.0 * (2 * np.log(2)) ** 0.5)  # Convert to sigma
    pixarea = abs(hdu[0].header['CDELT1'] * hdu[0].header['CDELT2'])

    beamarea = 2 * np.pi * 1.0 * beammaj * beammin  # Note that the volume of a two dimensional gaus$
    beamarea_pix = beamarea / pixarea

    return beamarea_pix


def print_beam(hdu):
    # print the beam from the header

    bmaj = hdu[0].header['BMAJ']
    bmin = hdu[0].header['BMIN']
    bpa = hdu[0].header['BPA']

    return bmaj, bmin, bpa


def findrms(mIn, maskSup=1e-7):
    """
    find the rms of an array, from Cycil Tasse/kMS
    """
    m = mIn[np.abs(mIn) > maskSup]
    rmsold = np.std(m)
    diff = 1e-1
    cut = 3.
    bins = np.arange(np.min(m), np.max(m), (np.max(m) - np.min(m)) / 30.)
    med = np.median(m)
    for i in range(10):
        ind = np.where(np.abs(m - med) < rmsold * cut)[0]
        rms = np.std(m[ind])
        if np.abs(old_div((rms - rmsold), rmsold)) < diff: break
        rmsold = rms
    return rms


def add_mask(hdu, ds9region, flat_it=True):
    if flat_it == True:
        hduflat = flatten(hdu)
        # map=hdu[0].data
        w = WCS(hduflat.header)
    else:
        hduflat = hdu

    r = pyregion.open(ds9region)
    manualmask = r.get_mask(hdu=hduflat)

    return manualmask


def get_image_info(infile, xray=False, y=False):
    hdu = fits.open(infile)
    hdr = hdu[0].header

    imsize1 = hdu[0].header['NAXIS1']
    imsize2 = hdu[0].header['NAXIS2']
    print('The image size is: {} pixel x {} pixel'.format(imsize1, imsize2))

    pixsize = abs(hdu[0].header['CDELT1'] * 3600.)  # in arcsec
    print('The pixel scale is: {:.2f} arcsec'.format(pixsize))

    if xray == True:
        return pixsize

    if y == True:
        rmsy = np.nanstd(np.ndarray.flatten(hdu[0].data))
        return pixsize, rmsy

    beamarea_pix = calc_beamarea(hdu)
    print('The beam area is: {:.2f} pixel^2'.format(beamarea_pix))  # in pixel^2
    beamarea_arcsec = beamarea_pix * pixsize * pixsize
    print('The beam area is: {:.2f} arcsec^2'.format(beamarea_arcsec))  # in arcsec^2

    rms_jyb = findrms(np.ndarray.flatten(hdu[0].data))
    print('The rms is: {:.3f} mJy/beam'.format(rms_jyb * 1e3))  # in mJy/beam
    rms_jyarcsec = rms_jyb / beamarea_arcsec
    print('The rms is: {:.3f} muJy/arcsec^2'.format(rms_jyarcsec * 1e6))  # in muJy/arcsec^2

    freq = hdu[0].header['CRVAL3']  # in Hz
    print('The frequency is {:.2f} MHz'.format(freq / 1e6))

    hdu.close()

    return pixsize, beamarea_pix, rms_jyb, rms_jyarcsec, freq


def reproject_hdu(hdu1, hdu2, outputname):
    array, _ = reproject_interp(hdu2, hdu1.header)
    fits.writeto(outputname, array, hdu1.header, overwrite=True)

    return


# def reproject_montage(footprint, infile):
#     '''
#     not used nor tested. it may be faster
#     footprint = 'a2255c1.75im2340r-0.5t5u60-MFS-image_15x14.fits'
#     infile = 'a2255_0.5-2.0_thresh_time.expmap.gz'
#     '''
#     import montage_wrapper as montage
#     montage.commands.mGetHdr(footprint, 'in_header')
#     montage.wrappers.reproject(infile, 'infile_reprojected.fits', header='in_header', exact_size=True)
#     return


def make_ds9_image_region_header(regionfile):
    with open(regionfile, 'w') as f:
        f.write('# Region file format: DS9 version 4.1\n')
        f.write(
            'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        f.write('image\n')

    return


def make_grid(regionfile, xmin, ymin, xmax, ymax, size):
    make_ds9_image_region_header(regionfile)

    counter = 0

    j = xmin
    while j <= xmax:
        i = ymin
        while i <= ymax:
            with open(regionfile, 'a') as f:
                f.write('box({},{},{},{},0)\n'.format(j, i, size, size))

            xaf_region = 'xaf_' + str(counter) + '.reg'
            make_ds9_image_region_header(xaf_region)
            with open(xaf_region, 'a') as f:
                f.write('box({},{},{},{},0)\n'.format(j, i, size, size))

            counter += 1
            i += size
        j += size

    print('Created: ' + regionfile)

    return


def calc_sum_pix(hduflat, xaf_region, cellsize,
                 reproj_corr_factor=1.0):  # merge with function below, need to pass hduflatten
    # reporject dont conserve sum! https://reproject.readthedocs.io/en/stable/celestial.html
    # as reprojection does not conserve flux, we need to use a correction (new pixel area/old pixel area) when we compute the sum

    print('Reading ' + xaf_region)

    xaf_mask = add_mask(hduflat, xaf_region, flat_it=False)
    xaf_sum = np.sum(hduflat.data[np.where(xaf_mask == True)])
    xaf_median = np.median(hduflat.data[np.where(xaf_mask == True)])
    xaf_Npix = np.count_nonzero(hduflat.data[np.where(xaf_mask == True)])

    if reproj_corr_factor != 1.0:
        xaf_sum = xaf_sum * reproj_corr_factor

    return xaf_sum, xaf_Npix, xaf_median


def make_alpha_map(hduflat, xaf_region, value):
    xaf_mask = add_mask(hduflat, xaf_region, flat_it=False)
    hduflat.data[xaf_mask] = value

    return


def add_info_to_table(intable):
    '''
    Add info computed of the radio images to the header of the table with the results.
    This info will be used by the plotting script.
    '''
    hdul = fits.open(intable)
    hdr = hdul[0].header

    hdr['freq1'] = (radio1_freq, 'Frequency [Hz]')
    hdr['rms1asec'] = (radio1_rms_jyarcsec, 'rms [Jy/arcsec^2]')
    hdr['rms1beam'] = (radio1_rms_jyb, 'rms [Jy/beam]')
    hdr['sys1'] = (sys1, 'Systematic error')

    hdul.writeto(intable, overwrite=True)
    hdul.close()

    return


# PARSER

parser = argparse.ArgumentParser(
    description='Perform point-to-point analysis (radio/X or radio/radio/X) and save the results in a fits table. If (i) the X-ray counts < 0, (ii) the cell is not totally inside the X-ray FoV, and (ii) the radio flux density is < 0, it saves NaNs.')
parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
# REQUIRED
required.add_argument('-radio1', help='First radio image (low freq).', type=str, required=True)
required.add_argument('-limits',
                      help='Space-separted list of image coordinates defining the region used to draw the grid in the -radio1 image, in format xmin, ymin, xmax, ymax.',
                      nargs=4, type=int, required=True)
required.add_argument('-xsou', help='X-ray source image (in count units).', type=str, required=True)
required.add_argument('-xbkg', help='X-ray background image (in count units).', type=str, required=True)
required.add_argument('-xexp', help='X-ray exposure image (in second units).', type=str, required=True)
required.add_argument('-cellsize',
                      help='Size of the cell, it should be AN ODD NUMER smaller than the beam! (in pixel units).',
                      type=int, required=True)
# OPTIONAL
optional.add_argument('-median', help='take median value')
optional.add_argument('-y', help='y-map (planck)')
optional.add_argument('-rms1', help='Noise of first radio image (in Jy/b). If not provided, it is internally computed.',
                      type=float, required=False, default=0.0)
optional.add_argument('-sys1', help='Systematic error on the first radio image (default=0).', type=float,
                      required=False, default=0.0)
optional.add_argument('-excluderegion', help='ds9 region file with the regions to exclude from the grid.', type=str,
                      required=False, default=None)
optional.add_argument('-skip_reproj',
                      help='Skip the reprojection of X-ray and secondary radio images (True/False, default=False).',
                      type=str2bool, required=False, default=False)

# python ptp_analysis.py  -radio1 a2255c1.75im2340r-0.5t5u60-MFS-image_15x14.fits -limits 1052 895 1532 1374 -xsou a2255_0.5-2.0_thresh_nopt.img.gz -xbkg a2255_0.5-2.0_thresh.bkgmap.gz -xexp a2255_0.5-2.0_thresh_time.expmap.gz -cellsize 35  -skip_reproj T -sys1 0.2 -sys2 0.05 -radio2 A2255_25CM_BEAM.image.regridded.aligned.fits -excluderegion GB_sources_ds9.reg

# python ptp_analysis_nan.py  -radio1 b3-sub-15arcsec.fits -limits 195 181 285 272 -xsou gordo_0.5-2.0_thresh_nopt.img.gz -xbkg gordo_0.5-2.0_thresh.bkgmap.gz -xexp gordo_0.5-2.0_thresh_time.expmap.gz -cellsize 9 -radio2 b4-sub-15arcsec.fits -excluderegion gordo_mask_ds9image.reg -skip_reproj F

args = vars(parser.parse_args())

radio1 = args['radio1']
median = args['median']
ymap = args['y']
rms1 = args['rms1']
sys1 = args['sys1']
excluderegion = args['excluderegion']
limitslist = args['limits']
xsou = args['xsou']
xbkg = args['xbkg']
xexp = args['xexp']
cellsize = args['cellsize']

xmin, ymin, xmax, ymax = limitslist[0], limitslist[1], limitslist[2], limitslist[3]

ptp_dir = './ptp_dir/'
grid_dir = ptp_dir + 'grid_{}x{}'.format(cellsize, cellsize)
grid_name = 'grid_{}x{}_ds9_image.reg'.format(cellsize, cellsize)
result_name = 'grid_{}x{}_results.fits'.format(cellsize, cellsize)

# check that the cellsize is an odd number
if (cellsize % 2) != 1:
    print('The cellsize must be an odd number. Exiting...')
    sys.exit()

if not os.path.isdir(ptp_dir):
    os.system('mkdir ' + ptp_dir)
if not os.path.isdir(grid_dir):
    os.system('mkdir ' + grid_dir)

# Read info from files
print('Information of ' + radio1)
radio1_pixscale, radio1_beam_area_pix, radio1_rms_jyb, radio1_rms_jyarcsec, radio1_freq = get_image_info(radio1)

radio1_hdu = fits.open(radio1)
radio1_hduflat = flatten(radio1_hdu)

if rms1 != 0.0:
    print('You set -rms1, so using the input value for the noise, ie:')
    radio1_rms_jyb = rms1
    radio1_rms_jyarcsec = radio1_rms_jyb / (radio1_beam_area_pix * radio1_pixscale * radio1_pixscale)
    print('The rms is: {:.3f} mJy/beam'.format(radio1_rms_jyb * 1e3))  # in mJy/beam
    print('The rms is: {:.3f} muJy/arcsec^2'.format(radio1_rms_jyarcsec * 1e6))  # in muJy/arcsec^2

print('Information of ' + xsou)
xray_pixscale = get_image_info(xsou, xray=True)

if ymap:
    y_pixscale, rmsy = get_image_info(ymap, y=True)


if args['skip_reproj'] == False:
    xsou_hdu = fits.open(xsou)
    xbkg_hdu = fits.open(xbkg)
    xexp_hdu = fits.open(xexp)
    xsou_hduflat = flatten(xsou_hdu)
    xbkg_hduflat = flatten(xbkg_hdu)
    xexp_hduflat = flatten(xexp_hdu)

    reproject_hdu(radio1_hduflat, xsou_hduflat, 'xsou_reproj.fits.gz')
    reproject_hdu(radio1_hduflat, xbkg_hduflat, 'xbkg_reproj.fits.gz')
    reproject_hdu(radio1_hduflat, xexp_hduflat, 'xexp_reproj.fits.gz')

    xsou_hdu.close()
    xbkg_hdu.close()
    xexp_hdu.close()

    if ymap:
        y_hdu = fits.open(ymap)
        y_hduflat = flatten(y_hdu)
        reproject_hdu(radio1_hduflat, y_hduflat, 'y_reproj.fits.gz')
        y_hdu.close()


print('IMAGE REPROJECTED, BE CAREFUL ABOUT THE DIFFERENCE IN PIXELSCALE')

print('Information of xsou_reproj.fits.gz')
xray_reproj_pixscale = get_image_info('xsou_reproj.fits.gz', xray=True)

xray_reproj_corr_factor = (xray_reproj_pixscale / xray_pixscale) ** 2

if ymap:
    y_reproj_pixscale, rmsy = get_image_info('y_reproj.fits.gz', y=True)

    y_reproj_corr_factor = (y_reproj_pixscale / y_pixscale) **2

# create grid & single xaf regions into grid_dir directory
os.chdir(grid_dir)
make_grid(grid_name, xmin, ymin, xmax, ymax, cellsize)
os.chdir('../../')

xaf_region_list = natural_sort(glob.glob(grid_dir + '/xaf_*.reg'))

xsou_reproj_hdu = fits.open('xsou_reproj.fits.gz')
xbkg_reproj_hdu = fits.open('xbkg_reproj.fits.gz')
xexp_reproj_hdu = fits.open('xexp_reproj.fits.gz')
xsou_reproj_hduflat = flatten(xsou_reproj_hdu)
xbkg_reproj_hduflat = flatten(xbkg_reproj_hdu)
xexp_reproj_hduflat = flatten(xexp_reproj_hdu)
xexp_reproj_hduflat.data = np.where(xexp_reproj_hduflat.data==0., 1., xexp_reproj_hduflat.data)

if ymap:
    y_reproj_hdu = fits.open('y_reproj.fits.gz')
    y_reproj_hduflat = flatten(y_reproj_hdu)
    y_reproj_hduflat.data = np.where(y_reproj_hduflat.data==0., 1., y_reproj_hduflat.data)

# apply exclude_mask
if excluderegion != None:
    exclude_mask = add_mask(radio1_hduflat, excluderegion, flat_it=False)
    radio1_hduflat.data[exclude_mask] = 0.0
    exclude_mask = add_mask(xexp_reproj_hduflat, excluderegion, flat_it=False)
    xexp_reproj_hduflat.data[exclude_mask] = 0.0
    if ymap:
        exclude_mask = add_mask(y_reproj_hduflat, excluderegion, flat_it=False)
        y_reproj_hduflat.data[exclude_mask] = 0.0


# prepare table
if os.path.isfile(grid_dir + '/' + result_name):
    os.system('rm ' + grid_dir + '/' + result_name)

if ymap:
    tab_init = {'region_name': Column(dtype='S8'),
                'xray_sb': Column(dtype=float, unit='count/s/arcsec^2'),
                'xray_sb_err': Column(dtype=float, unit='count/s/arcsec^2'),
                'radio1_sb': Column(dtype=float, unit='Jy/arcsec^2'),
                'radio1_sb_err': Column(dtype=float, unit='Jy/arcsec^2'),
                'radio1_fluxdensity': Column(dtype=float, unit='Jy'),
                'radio1_fluxdensity_err': Column(dtype=float, unit='Jy'),
                'y_sb': Column(dtype=float, unit='Compton-y/arcsec^2'),
                'y_sb_err': Column(dtype=float, unit='Compton-y/arcsec^2'),
                'y_fluxdensity': Column(dtype=float, unit='Compton-y'),
                'y_fluxdensity_err': Column(dtype=float, unit='Compton-y'),
                }
    t = Table(tab_init, names=(
        'region_name', 'xray_sb', 'xray_sb_err', 'radio1_sb', 'radio1_sb_err', 'radio1_fluxdensity',
        'radio1_fluxdensity_err', 'y_sb', 'y_sb_err', 'y_fluxdensity', 'y_fluxdensity_err'))

else:
    tab_init = {'region_name': Column(dtype='S8'),
                'xray_sb': Column(dtype=float, unit='count/s/arcsec^2'),
                'xray_sb_err': Column(dtype=float, unit='count/s/arcsec^2'),
                'radio1_sb': Column(dtype=float, unit='Jy/arcsec^2'),
                'radio1_sb_err': Column(dtype=float, unit='Jy/arcsec^2'),
                'radio1_fluxdensity': Column(dtype=float, unit='Jy'),
                'radio1_fluxdensity_err': Column(dtype=float, unit='Jy'),
                }
    t = Table(tab_init, names=(
    'region_name', 'xray_sb', 'xray_sb_err', 'radio1_sb', 'radio1_sb_err', 'radio1_fluxdensity', 'radio1_fluxdensity_err'))

# start looping on regions
for xaf_region in xaf_region_list:

    # RADIO1
    print('Doing radio1 image:')
    radio1_sum, radio1_Npix, _ = calc_sum_pix(radio1_hduflat, xaf_region, cellsize)

    if radio1_Npix / (cellsize * cellsize) != 1:
        print('Removing region {}: it contains pixels with zero value'.format(xaf_region))
        #        os.system('rm -rf ' + xaf_region)
        #        continue
        radio1_fluxdensity = np.nan
        radio1_fluxdensity_err = np.nan
        radio1_sb = np.nan
        radio1_sb_err = np.nan
    elif radio1_sum <= 0:
        print('Removing region {}: the sum of pixel values is negative'.format(xaf_region))
        #        os.system('rm -rf ' + xaf_region)
        #        continue
        radio1_fluxdensity = np.nan
        radio1_fluxdensity_err = np.nan
        radio1_sb = np.nan
        radio1_sb_err = np.nan
    else:
        radio1_fluxdensity = radio1_sum / radio1_beam_area_pix
        radio1_fluxdensity_err = np.sqrt(
            (radio1_rms_jyb * np.sqrt(radio1_Npix / radio1_beam_area_pix)) ** 2 + (sys1 * radio1_fluxdensity) ** 2)
        radio1_sb = radio1_fluxdensity / (radio1_Npix * radio1_pixscale * radio1_pixscale)
        radio1_sb_err = radio1_fluxdensity_err / (radio1_Npix * radio1_pixscale * radio1_pixscale)

    # XRAY
    print('Doing xray images:')
    xsou_sum, xsou_Npix, _ = calc_sum_pix(xsou_reproj_hduflat, xaf_region, cellsize,
                                          reproj_corr_factor=xray_reproj_corr_factor)
    xbkg_sum, xbkg_Npix, _ = calc_sum_pix(xbkg_reproj_hduflat, xaf_region, cellsize,
                                          reproj_corr_factor=xray_reproj_corr_factor)
    xexp_sum, xexp_Npix, xexp_median = calc_sum_pix(xexp_reproj_hduflat, xaf_region, cellsize,
                                                  reproj_corr_factor=xray_reproj_corr_factor)
    xray_net_count = (xsou_sum - xbkg_sum)

    if xexp_Npix / (cellsize * cellsize) != 1:
        # THIS IS NOT OK FOR XMM, A SOLUTION CUOLD BE PLACE 1s WHERE THE EXPOSURE IS 0s DUE TO CCD GAPS (THEY ARE SMALL PORTIONS)
        print('Removing region {}: it contains pixels with zero value'.format(xaf_region))
        #        os.system('rm -rf ' + xaf_region)
        #        continue
        xray_count_rate = np.nan
        xray_count_rate_err = np.nan
        xray_sb = np.nan
        xray_sb_err = np.nan
    elif xray_net_count <= 0:
        print('Removing region {}: it has <= 0 net counts'.format(xaf_region))
        #        os.system('rm -rf ' + xaf_region)
        #        continue
        xray_count_rate = np.nan
        xray_count_rate_err = np.nan
        xray_sb = np.nan
        xray_sb_err = np.nan
    else:
        xray_count_rate = (xsou_sum - xbkg_sum) / xexp_median
        xray_count_rate_err = np.sqrt(xsou_sum + xbkg_sum) / xexp_median
        xray_sb = xray_count_rate / ((cellsize * cellsize) * xray_reproj_pixscale * xray_reproj_pixscale)
        xray_sb_err = xray_count_rate_err / ((cellsize * cellsize) * xray_reproj_pixscale * xray_pixscale)


    print('Doing y-map (SZ-effect)')
    #Ymap
    y_sum, y_Npix, y_median = calc_sum_pix(y_reproj_hduflat, xaf_region, cellsize,
                                          reproj_corr_factor=y_reproj_corr_factor)
    if y_Npix / (cellsize * cellsize) != 1:
        print('Removing region {}: it contains pixels with zero value'.format(xaf_region))
        #            os.system('rm -rf ' + xaf_region)
        #            continue
        y_fluxdensity = np.nan
        y_fluxdensity_err = np.nan
        y_sb = np.nan
        y_sb_err = np.nan
    elif y_sum <= 0:
        print('Removing region {}: the sum of pixel values is negative'.format(xaf_region))
        #            os.system('rm -rf ' + xaf_region)
        #            continue
        y_fluxdensity = np.nan
        y_fluxdensity_err = np.nan
        y_sb = np.nan
        y_sb_err = np.nan
    else:
        y_fluxdensity = y_sum
        y_fluxdensity_err = rmsy * np.sqrt(y_Npix)
        y_sb = y_fluxdensity / (y_Npix * y_reproj_pixscale * y_reproj_pixscale)
        y_sb_err = y_fluxdensity_err / (y_Npix * y_reproj_pixscale * y_reproj_pixscale)

    if ymap:
        t.add_row(
            [xaf_region.replace(grid_dir + '/', '').replace('.reg', ''), xray_sb, xray_sb_err, radio1_sb, radio1_sb_err,
             radio1_fluxdensity, radio1_fluxdensity_err, y_sb, y_sb_err, y_fluxdensity,
             y_fluxdensity_err])
    else:
        t.add_row(
            [xaf_region.replace(grid_dir + '/', '').replace('.reg', ''), xray_sb, xray_sb_err, radio1_sb, radio1_sb_err,
             radio1_fluxdensity, radio1_fluxdensity_err])


t.write(grid_dir + '/' + result_name, format='fits')

# add info in the result table
hdul = fits.open(grid_dir + '/' + result_name)
hdr = hdul[0].header

hdr['freq1'] = (radio1_freq, 'Frequency [Hz]')
hdr['rms1asec'] = (radio1_rms_jyarcsec, 'rms [Jy/arcsec^2]')
hdr['rms1beam'] = (radio1_rms_jyb, 'rms [Jy/beam]')
hdr['sys1'] = (sys1, 'Systematic error')

hdul.writeto(grid_dir + '/' + result_name, overwrite=True)
hdul.close()

# add_info_to_table(ptp_dir + '/' + result_name) #not used because I would need to pass many variables


radio1_hdu.close()
xsou_reproj_hdu.close()
xbkg_reproj_hdu.close()
xexp_reproj_hdu.close()

#python ptp_analysis.py  -radio1 image_test_L626678_rvw20-MFS-image.fits -limits 600 500 950 1025 -xsou a399-401/mosaic_a399_a401.fits.gz -xbkg a399-401/mosaic_a399_a401_bkg.fits.gz -xexp a399-401/mosaic_a399_a401_exp.fits.gz -cellsize 35