"""
This script is meant to mosaic together facets.
The algorithm uses feathering for overlapping facets.

For this you will need:
- facet fits files
- region files

Make sure that if you sort the facet fits files and the region files on string name, the fits and region correspond with
each other on index.
For example, it is best to name the fits and region files:
facet_1.fits, facet_2.fits, ... facet_10.fits, ...
region_1.reg, region_2.reg, ... region_10.reg, ...

Also update the coordinates and object name below the script if you image another field than ELAIS-N1
"""

from __future__ import print_function
import sys
from argparse import ArgumentParser
from reproj_test import reproject_interp_chunk_2d
from auxcodes import flatten
from shapely import geometry
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
import pyregion
import astropy.units as u
from astropy.coordinates import SkyCoord
from past.utils import old_div
from astropy.wcs.utils import skycoord_to_pixel
import gc

plt.style.use('ggplot')

def make_header(fitsfile, fullpixsize):
    """
    Make header

    :param fitsfile: fits file input
    :param fullpixsize: full wide-field pixel size

    :return: header
    """
    hdu = fits.open(fitsfile)
    himsize = fullpixsize // 2
    header = fits.Header()
    header['BITPIX'] = -32
    header['NAXIS'] = 2
    header['WCSAXES'] = 2
    header['NAXIS1'] = 2 * himsize
    header['NAXIS2'] = 2 * himsize
    header['CTYPE1'] = 'RA---SIN'
    header['CTYPE2'] = 'DEC--SIN'
    header['CUNIT1'] = 'deg'
    header['CUNIT2'] = 'deg'
    header['CRPIX1'] = himsize
    header['CRPIX2'] = himsize
    header['CRVAL1'] = CRVAL1
    header['CRVAL2'] = CRVAL2
    header['CDELT1'] = -hdu[0].header['CDELT2']
    header['CDELT2'] = hdu[0].header['CDELT2']
    header['LATPOLE'] = header['CRVAL2']
    header['BMAJ'] = hdu[0].header['BMAJ']
    header['BMIN'] = hdu[0].header['BMIN']
    header['BPA'] = hdu[0].header['BPA']
    header['TELESCOPE'] = 'LOFAR'
    header['OBSERVER'] = 'LoTSS'
    header['BUNIT'] = 'JY/BEAM'
    header['BSCALE'] = 1.0
    header['BZERO'] = 0
    header['BTYPE'] = 'Intensity'
    header['OBJECT'] = OBJECT_NAME
    return header


def get_polygon_center(regionfile):
    """
    get polygon center
    :param regionfile: region file
    :return: polygon center
    """
    regname = regionfile
    regionfile = open(regionfile, 'r')
    lines = regionfile.readlines()
    regionfile.close()
    try:
        polygon = lines[4]
    except IndexError:
        polygon = lines[3]
    polygon = polygon.replace(' # text={'+regname.replace('.reg','')+'}\n','')
    polyp = [float(p) for p in polygon.replace('polygon(', '').replace(')', '').replace('\n', '').split(',')]
    poly_geo = geometry.Polygon(tuple(zip(polyp[0::2], polyp[1::2])))
    return SkyCoord(f'{poly_geo.centroid.x}deg', f'{poly_geo.centroid.y}deg', frame='icrs')


def get_distance_weights(center, arr, wcsheader):
    """
    Get weights based on center polygon to coordinates in array

    :param center: center polygon
    :param arr: numpy array (for shape)
    :param wcsheader: header
    :return: weights from center polygon
    """

    pix_center = [int(c) for c in skycoord_to_pixel(center, wcsheader)]
    pix_pos = np.array(np.where(np.ones(arr.shape))).T.astype(np.float32)
    pix_dist = np.linalg.norm(np.subtract(pix_center, pix_pos), axis=1).astype(np.float32)
    weights = np.divide(wcsheader.to_header()['CRPIX1'], pix_dist).astype(np.float32)
    # remove infinite
    weights[np.where(weights == np.inf)[0].squeeze()] = weights[np.where(weights != np.inf)].max()
    return weights.reshape(arr.shape).T.astype(np.float32)


def rms(image_data):
    """
    from Cyril Tasse/kMS

    :param image_data: image data array
    :return: rms (noise measure)
    """

    maskSup = 1e-7
    m = image_data[np.abs(image_data)>maskSup]
    rmsold = np.std(m)
    diff = 1e-1
    cut = 3.
    med = np.median(m)
    for _ in range(10):
        ind = np.where(np.abs(m - med) < rmsold*cut)[0]
        rms = np.std(m[ind])
        if np.abs(old_div((rms-rmsold), rmsold)) < diff: break
        rmsold = rms
    print(f'Noise : {str(round(rms * 1000, 4))} {u.mJy/u.beam}')
    return rms


def make_image(image_data=None, hdu=None, save=None, cmap: str = 'CMRmap', header=None):
    """
    Image your data with this method.
    image_data -> insert your image_data or plot full image
    cmap -> choose your preferred cmap
    """

    RMS = 1.76e-05 #TODO: TEMPORARY
    vmin = RMS
    vmax = RMS*20

    # if hdu is None:
    #     wcs = WCS(header, naxis=2)
    # else:
    #     wcs = WCS(hdu[0].header, naxis=2)


    plt.figure(figsize=(7, 10), dpi=200)
    # plt.subplot(projection=wcs)
    plt.subplot()
    # WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=wcs)

    # image_data[image_data == np.inf] = np.nan
    # image_data[image_data == 0] = np.nan

    im = plt.imshow(image_data, origin='lower', cmap=cmap)

    im.set_norm(PowerNorm(vmin=0, vmax=vmax, gamma=1 / 2))

    # plt.xlabel('Right Ascension (J2000)', size=14)
    # plt.ylabel('Declination (J2000)', size=14)
    # plt.tick_params(axis='both', which='major', labelsize=12)
    plt.xticks([])
    plt.yticks([])

    plt.grid(False)
    plt.grid('off')
    plt.grid(None)
    plt.tight_layout()

    if save is not None:
        plt.savefig(save, dpi=250, bbox_inches='tight')
    else:
        plt.show()


def parse_arg():
    """
    Command line argument parser
    """

    parser = ArgumentParser(description='make wide-field by combining facets')
    parser.add_argument('--resolution', help='resolution in arcsecond', required=True, type=float)
    parser.add_argument('--fits', type=str, nargs='+', help='facets to merge')
    parser.add_argument('--regions', type=str, nargs='+', help='regions corresponding to facets')
    return parser.parse_args()


def main():
    """
    Main function
    """

    args = parse_arg()

    resolution = args.resolution
    facets = sorted(args.fits)
    regions = sorted(args.regions)

    if resolution == 0.3:
        pixelscale = 0.1  # arcsec
    elif resolution == 0.6:
        pixelscale = 0.2
    elif resolution == 1.2:
        pixelscale = 0.4
    else:
        sys.exit('ERROR: only use resolution 0.3 or 1.2')

    fullpixsize = int(2.5 * 3600 / pixelscale)

    header_new = make_header(facets[0], fullpixsize)
    print(header_new)

    xsize = header_new['NAXIS1']
    ysize = header_new['NAXIS2']

    isum = np.zeros([ysize, xsize], dtype="float32")
    weights = np.zeros_like(isum, dtype="float32")
    fullmask = np.zeros_like(isum, dtype=bool)

    for n, facet in enumerate(facets):
        print(facet, regions[n])

        hdu = fits.open(facet)
        hduflatten = flatten(hdu)

        imagedata, _ = reproject_interp_chunk_2d(hduflatten, header_new, hdu_in=0, parallel=True)
        imagedata = imagedata.astype(np.float32)

        # trying to clear cache? Not sure if it works..
        hduflatten = None
        del hduflatten

        reg = regions[n]
        polycenter = get_polygon_center(reg)
        r = pyregion.open(reg).as_imagecoord(header=header_new)
        mask = r.get_mask(hdu=hdu[0], shape=(header_new["NAXIS1"], header_new["NAXIS2"])).astype(np.int16)
        hdu.close()

        make_image(mask*imagedata, None, facet+'.png', 'CMRmap', header_new)

        fullmask |= ~np.isnan(imagedata)
        facetweight = get_distance_weights(polycenter, mask, WCS(header_new, naxis=2)) * mask
        # facetweight = mask
        facetweight[(~np.isfinite(imagedata)) | (~np.isfinite(facetweight)) | (imagedata == 0)] = 0  # so we can add
        imagedata *= facetweight
        imagedata[~np.isfinite(imagedata)] = 0  # so we can add
        isum += imagedata

        # trying to clear cache? Not sure if it works..
        imagedata = None
        del imagedata
        weights += facetweight

        # trying to clear cache? Not sure if it works..
        facetweight = None
        del facetweight
        gc.collect()

        make_image(isum, None, facet.split('/')[-1]+'full.png', 'CMRmap', header_new)

    print('FACET COMPLETED ... ADDING CORRECT WEIGHTS ...')

    isum /= weights
    isum[isum == np.inf] = np.nan
    isum[~fullmask] = np.nan

    make_image(isum, None, 'finalfaceted.png', 'CMRmap', header_new)

    hdu = fits.PrimaryHDU(header=header_new, data=isum)

    hdu.writeto('full-mosaic.fits', overwrite=True)


if __name__ == '__main__':
    #ELAIS-N1 --> UPDATE THIS FOR OWN USE
    CRVAL1=-117.25
    CRVAL2=54.95
    OBJECT_NAME='ELAIS-N1'

    main()
