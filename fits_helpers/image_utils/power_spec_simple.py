import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft2, fftshift
from astropy.io import fits


def spectrum(fts):
    """Make power spectrum"""
    f = fits.open(fts)
    img = f[0].data
    fft_result = fft2(img)
    fft_result_shifted = fftshift(fft_result)
    magnitude_spectrum = np.abs(fft_result_shifted)
    plt.figure(figsize=(12, 6))
    plt.subplot(121)
    plt.imshow(np.log1p(img[0][0]), cmap='RdBu_r')
    plt.title('Original Image')
    plt.subplot(122)
    plt.imshow(np.log1p(magnitude_spectrum[0][0]), cmap='RdBu_r')
    plt.title('2D Fourier Transform')
    plt.show()
