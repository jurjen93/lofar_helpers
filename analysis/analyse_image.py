from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from past.utils import old_div



def findrms(mIn, maskSup=1e-7):
    """
    from Cyril Tasse/kMS
    """

    m = mIn[np.abs(mIn)>maskSup]
    rmsold = np.std(m)
    diff = 1e-1
    cut = 3.
    med = np.median(m)
    for _ in range(10):
        ind = np.where(np.abs(m - med) < rmsold*cut)[0]
        rms = np.std(m[ind])
        if np.abs(old_div((rms-rmsold), rmsold)) < diff: break
        rmsold = rms
    return rms

# final_noise=0
# for i in range(6):
#     f = fits.open(f'/net/nieuwerijn/data2/jurjendejong/Abell399-401_indimages/archive{i}.fits')
#     noise = findrms(np.ndarray.flatten(f[0].data))
#     print(f'archive{i}')
#     print('Noise : '+str(round(noise*1000,2))+' mJy/beam')
#     final_noise+=noise**2
#     f.close()
# print('Final predicted noise : '+str(round(np.sqrt(final_noise)*1000/6, 2))+' mJy/beam')

f = fits.open(f'../image_test_A399-MFS-image.fits')
noise = findrms(np.ndarray.flatten(f[0].data))
print('Noise : ' + str(round(noise * 1000, 2)) + ' mJy/beam')


