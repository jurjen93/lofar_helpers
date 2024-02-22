from glob import glob
import numpy as np
from astropy.io import fits
import astropy.units as u
from past.utils import old_div

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

print("0.3''")
for fts in sorted(glob('/project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/imaging/split_facets2/final_img_0.3_rms_old/facet_*.fits')):
    print(fts)
    f = fits.open(fts)[0]
    print(f"{f.header['BMIN']*3600}'' x {f.header['BMAJ']*3600}''")
    data = f.data
    pixelsize = f.header['CDELT1'] ** 2
    print(f'Size: {len(data[data == data]) * pixelsize} degree^2')
    rms(data)
    print()

print("0.6''")
for fts in sorted(glob('/project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/imaging/split_facets2/final_img_0.6_rms_old/facet_*.fits')):
    print(fts)
    f = fits.open(fts)[0]
    print(f"{f.header['BMIN']*3600}'' x {f.header['BMAJ']*3600}''")
    data = f.data
    pixelsize = f.header['CDELT1'] ** 2
    print(f'Size: {len(data[data == data]) * pixelsize} degree^2')
    rms(data)
    print()

print("1.2''")
for fts in sorted(glob('/project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/imaging/split_facets2/final_img_1.2_rms_old/facet_?.fits') +
                  glob('/project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/imaging/split_facets2/final_img_1.2_rms_old/facet_??.fits')):
    print(fts)
    f = fits.open(fts)[0]
    print(f"{round(f.header['BMIN']*3600, 2)}'' x {round(f.header['BMAJ']*3600,2)}''")
    pixelsize = f.header['CDELT1'] ** 2
    data = f.data
    print(f'Size: {len(data[data == data]) * pixelsize} degree^2')
    rms(data)
    print()
