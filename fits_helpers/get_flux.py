from astropy.io import fits
from astropy import units as u
import numpy as np

def findrms(mIn, maskSup=1e-7):
    """
    find the rms of an array, from Cycil Tasse/kMS
    """
    try:
        m = mIn[np.abs(mIn) > maskSup]
        rmsold = np.std(m)
        diff = 1e-1
        cut = 3.
        med = np.median(m)
        for i in range(10):
            ind = np.where(np.abs(m - med) < rmsold * cut)[0]
            rms = np.std(m[ind])
            if np.abs((rms - rmsold) / rmsold) < diff: break
            rmsold = rms
        return rms
    except ValueError:
        return np.nan


def getflux(filename: str):
    """
    Gef flux from fits file
    Args:
        filename: Input file name

    Returns: Flux density
    """

    with fits.open(filename) as hdu:

        header = hdu[0].header
        data = hdu[0].data
        rms = findrms(data)

        data = data[data > 5 * rms]  # filter noise/non-detections
        bmaj = header['BMAJ'] * u.deg
        bmin = header['BMIN'] * u.deg
        beam_area = (bmaj * bmin * np.pi / (4 * np.log(2))).to(u.steradian)

        try:
            cdelt1 = abs(header['CDELT1']) * u.deg
            cdelt2 = abs(header['CDELT2']) * u.deg
        except KeyError:
            cdelt1 = abs(header['CD1_1']) * u.deg
            cdelt2 = abs(header['CD2_2']) * u.deg

        pixel_area = (cdelt1 * cdelt2).to(u.steradian)
        conversion_factor = beam_area / pixel_area

        return float(np.sum(data) / conversion_factor.value)