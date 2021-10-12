"""
python ~/scripts/lofar_helpers/supporting_scripts/compare_results.py
--surf_im='/net/tussenrijn/data2/jurjendejong/A399/box_images/*.fits'
--int_restored='/net/rijn/data2/jdejong/A399_DEEP/image_full_ampphase_di_m.NS.int.restored.fits'
--image_final='/net/tussenrijn/data2/jurjendejong/A399/result/image_final.app.restored.fits'
--app_restored='/net/rijn/data2/jdejong/A399_DEEP/image_full_ampphase_di_m.NS.app.restored.fits'
--region_files='A399/boxes/*.reg'
"""

import pyregion
from astropy.io import fits
import numpy as np
from astropy.wcs import WCS
from argparse import ArgumentParser
from glob import glob
from astropy.nddata import Cutout2D
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import os

def make_cutout(fitsfile: str = None, region: str = None ):

    hdu = fits.open(fitsfile)
    image_data = hdu[0].data[0][0]
    header = hdu[0].header

    r = pyregion.open(region)
    ra  = r[0].coord_list[0]
    dec = r[0].coord_list[1]
    sizex = int(r[0].coord_list[2]/hdu[0].header['CDELT2'])
    sizey = int(r[0].coord_list[3]/hdu[0].header['CDELT2'])

    hdu.close()

    wcs = WCS(header, naxis=2)

    posx, posy = wcs.all_world2pix([[ra, dec]], 0)[0]

    out = Cutout2D(
        data=image_data,
        position=(int(posx), int(posy)),
        size=(sizex, sizey),
        wcs=wcs,
        mode='partial'
    )

    return out.data, out.wcs

def findrms(mIn,maskSup=1e-7):
    """
    find the rms of an array, from Cycil Tasse/kMS
    """
    m=mIn[np.abs(mIn)>maskSup]
    rmsold=np.std(m)
    diff=1e-1
    cut=3.
    bins=np.arange(np.min(m),np.max(m),(np.max(m)-np.min(m))/30.)
    med=np.median(m)
    for i in range(10):
        ind=np.where(np.abs(m-med)<rmsold*cut)[0]
        rms=np.std(m[ind])
        if np.abs((rms-rmsold)/rmsold)<diff: break
        rmsold=rms
    return rms

def make_image(image_final, app_restored, int_restored, surf_im, region, name):
    image_final, wcs = make_cutout(image_final, region)
    imagenoise_final = findrms(image_final)

    app_restored, wcs = make_cutout(app_restored, region)
    imagenoise_app = findrms(app_restored)

    int_restored, wcs = make_cutout(int_restored, region)
    imagenoise_int = findrms(int_restored)

    surf_im, wcs = make_cutout(surf_im, region)
    imagenoise_surf = findrms(surf_im)

    imagenoise = max(imagenoise_final, imagenoise_app, imagenoise_int, imagenoise_surf)

    print(imagenoise)

    fig, axs = plt.subplots(2, 2, figsize=(10, 10))

    axs[0, 0].imshow(image_final, norm=SymLogNorm(linthresh=imagenoise, vmin=imagenoise/10, vmax=16*imagenoise),
               origin='lower', cmap='bone')
    axs[0, 0].set_title('image_final')

    axs[0, 1].imshow(app_restored, norm=SymLogNorm(linthresh=imagenoise, vmin=imagenoise/10, vmax=16*imagenoise),
               origin='lower', cmap='bone')
    axs[0, 1].set_title('app_restored')

    axs[1, 0].imshow(int_restored, norm=SymLogNorm(linthresh=imagenoise, vmin=imagenoise/10, vmax=16*imagenoise),
               origin='lower', cmap='bone')
    axs[1, 0].set_title('int_restored')

    axs[1, 1].imshow(surf_im, norm=SymLogNorm(linthresh=imagenoise, vmin=imagenoise/10, vmax=16*imagenoise),
               origin='lower', cmap='bone')
    axs[1, 1].set_title('surf_im')

    fig.suptitle(region)

    plt.savefig(name+'/'+region.split('/')[-1].split('.')[0]+'.png')

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--surf_im', type=str, help='fits images produced on surfsara', required=True)
    parser.add_argument('--int_restored', type=str, help='*int.restored.fits file', required=True)
    parser.add_argument('--image_final', type=str, help='final image', required=True)
    parser.add_argument('--app_restored', type=str, help='*app.restored.fits image', required=True)
    parser.add_argument('--region_files', type=str, help='region files location', required=True)
    args = parser.parse_args()

    surf_images = glob(args.surf_im)
    region_files = glob(args.region_files)
    int_restored = args.int_restored
    image_final = args.image_final
    app_restored = args.app_restored

    outputfile = '/'.join(image_final.split('/')[0:-2])+'/comparisons'

    os.system('mkdir -p '+outputfile)


    for reg in region_files:
        try:
            surf_image = '/'.join(args.surf_im.split('/')[0:-1])+'/'+reg.split('/')[-1].split('.')[0]+'.fits'
            print(surf_image)
            make_image(image_final, app_restored, int_restored, surf_image, reg, outputfile)
        except:
            print(reg+' does not exist')
