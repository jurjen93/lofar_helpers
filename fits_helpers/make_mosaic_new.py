"""
(Alternative to make_mosaic.py, meant for doing the mosaicing in parallel)

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
# from reproject import reproject_interp
from auxcodes import flatten
from shapely import geometry
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm, PowerNorm
import pyregion
from astropy.visualization.wcsaxes import WCSAxes
import astropy.units as u
from astropy.coordinates import SkyCoord
from past.utils import old_div
from astropy.wcs.utils import skycoord_to_pixel

plt.style.use('ggplot')

def make_header(fitsfile, fullpixsize):
    """
    Make header

    :param fitsfile: fits file input
    :param fullpixsize: full wide-field pixel size

    :return: header
    """
    hdu = update_header(fits.open(fitsfile))
    himsize = fullpixsize // 2
    header = fits.Header()
    for k in list(hdu[0].header.keys()):
        if k=='HISTORY' or k=='COMMENT':
            continue
        header[k] = hdu[0].header[k]
    header['NAXIS'] = 2
    header['NAXIS1'] = 2 * himsize
    header['NAXIS2'] = 2 * himsize
    header['CRPIX1'] = himsize
    header['CRPIX2'] = himsize
    header['CRVAL1'] = CRVAL1
    header['CRVAL2'] = CRVAL2
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

def get_array_coordinates(pix_array, wcsheader):
    """
    Get coordinates from pixel

    :param pix_array: array with pixel coordinates
    :param wcsheader: wcs header
    :return: array with coordinates from pixel array
    """
    pixarray = np.argwhere(pix_array)
    return wcsheader.pixel_to_world(pixarray[:, 0], pixarray[:, 1], 0, 0)[0]


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

    if hdu is None:
        wcs = WCS(header, naxis=2)
    else:
        wcs = WCS(hdu[0].header, naxis=2)


    fig = plt.figure(figsize=(7, 10), dpi=200)
    plt.subplot(projection=wcs)
    WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=wcs)

    # image_data[image_data == np.inf] = np.nan
    # image_data[image_data == 0] = np.nan

    im = plt.imshow(image_data, origin='lower', cmap=cmap)

    im.set_norm(PowerNorm(vmin=0, vmax=vmax, gamma=1 / 2))

    plt.xlabel('Right Ascension (J2000)', size=14)
    plt.ylabel('Declination (J2000)', size=14)
    plt.tick_params(axis='both', which='major', labelsize=12)

    plt.grid(False)
    plt.grid('off')

    if save is not None:
        plt.savefig(save, dpi=250, bbox_inches='tight')
    else:
        plt.show()


def update_header(hdu):
    """Important to prevent mismatches in wcs headers"""

    for _ in range(100):
        for n, h in enumerate(hdu):
            max_d = len(h.data.shape)
            for key in h.header.keys():
                if key[-1].isdigit():
                    if int(float(key[-1]))>2:
                        h.header.remove(key)
                        print(key)

            hdu[n] = h

    return hdu


def parse_arg():
    """
    Command line argument parser
    """

    parser = ArgumentParser(description='make wide-field by combining facets')
    parser.add_argument('--resolution', help='resolution in arcsecond', required=True, type=float)
    parser.add_argument('--fits', type=str, nargs='+', help='facets to merge')
    parser.add_argument('--regions', type=str, nargs='+', help='regions corresponding to facets')
    parser.add_argument('--skip_reproject', action='store_true', help='Skip reproject')
    return parser.parse_args()


def reproject(fitsfile, header, region, skip):
    """
    Reproject fits file with new header and region file
    """

    hdu = update_header(fits.open(fitsfile))
    if not skip:
        hduflatten = flatten(hdu)
        imagedata, _ = reproject_interp_chunk_2d(hduflatten, header, parallel=True)
        del hduflatten
    else:
        imagedata = hdu[0].data
    imagedata = imagedata.astype(np.float32)
    polycenter = get_polygon_center(region)
    r = pyregion.open(region).as_imagecoord(header=header)
    mask = r.get_mask(hdu=hdu[0], shape=(header["NAXIS1"], header["NAXIS2"])).astype(np.int16)
    hdu.close()
    hdu = fits.PrimaryHDU(header=header, data=mask*imagedata)
    hdu.writeto(fitsfile.replace('.fits', '.reproject.fits'), overwrite=True)
    make_image(mask*imagedata, None, fitsfile.replace('.fits', 'full.png'), 'CMRmap', header)
    del imagedata

    facetweight = get_distance_weights(polycenter, mask, WCS(header, naxis=2)) * mask
    del mask
    facetweight[~np.isfinite(facetweight)] = 0  # so we can add
    hdu = fits.PrimaryHDU(header=header, data=facetweight)
    hdu.writeto(fitsfile.replace('.fits', '.weights.fits'), overwrite=True)

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

    for n, facet in enumerate(facets):
        print(facet)
        reproject(facet,  header_new, regions[n], args.skip_reproject)



if __name__ == '__main__':
    #ELAIS-N1
    CRVAL1=-117.25
    CRVAL2=54.95
    OBJECT_NAME='ELAIS-N1'

    main()