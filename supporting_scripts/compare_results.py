import pyregion
from astropy.io import fits
import numpy as np
from astropy.wcs import WCS
import sys
from argparse import ArgumentParser
from glob import glob
from astropy.nddata import Cutout2D
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm


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

def make_image(image_final, app_restored, int_restored, surf_im, region):
    image_final, wcs = make_cutout(image_final, region)
    imagenoise_final = findrms(image_final)

    app_restored, wcs = make_cutout(app_restored, region)
    imagenoise_app = findrms(app_restored)

    int_restored, wcs = make_cutout(int_restored, region)
    imagenoise_int = findrms(int_restored)

    surf_im, wcs = make_cutout(surf_im, region)
    imagenoise_surf = findrms(surf_im)

    imagenoise = max(imagenoise_final, imagenoise_app, imagenoise_int, imagenoise_surf)


    # plt.figure(figsize=(10, 10))
    # plt.subplot(projection=wcs)
    # plt.imshow(data, norm=SymLogNorm(linthresh=imagenoise*6, vmin=imagenoise, vmax=16*imagenoise),
    #            origin='lower', cmap='CMRmap')
    # plt.xlabel('Galactic Longitude')
    # plt.ylabel('Galactic Latitude')
    # plt.savefig('test.png')

    fig, axs = plt.subplots(2, 2)
    axs[0, 0].imshow(imagenoise_final, norm=SymLogNorm(linthresh=imagenoise*6, vmin=imagenoise, vmax=16*imagenoise),
               origin='lower', cmap='CMRmap')
    axs[0, 0].set_title('imagenoise_final')

    axs[0, 1].imshow(imagenoise_app, norm=SymLogNorm(linthresh=imagenoise*6, vmin=imagenoise, vmax=16*imagenoise),
               origin='lower', cmap='CMRmap')
    axs[0, 1].set_title('imagenoise_app')

    axs[1, 0].imshow(int_restored, norm=SymLogNorm(linthresh=imagenoise*6, vmin=imagenoise, vmax=16*imagenoise),
               origin='lower', cmap='CMRmap')
    axs[1, 0].set_title('int_restored')

    axs[1, 1].imshow(surf_im, norm=SymLogNorm(linthresh=imagenoise*6, vmin=imagenoise, vmax=16*imagenoise),
               origin='lower', cmap='CMRmap')
    axs[1, 1].set_title('surf_im')

    fig.suptitle(region)

    plt.savefig('test.png')

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

    for reg in region_files:
        print(image_final)
        surf_image = args.surf_im.replace('*','')+'/'+reg.split('.')[0].split('/')[-1]+'.fits'
        make_image(image_final, app_restored, int_restored, surf_image, reg)
        break