from analyse_image import Imaging
from past.utils import old_div
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
from astropy.nddata import Cutout2D
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from reproject import reproject_interp
import string
import sys
from astropy.modeling.models import Gaussian2D
from astropy.convolution import convolve, Gaussian2DKernel
from matplotlib.ticker import LogLocator, LogFormatterSciNotation as LogFormatter
import os


Image = Imaging(f'../fits/60arcsec_full.fits')
Image.taper(60/80)
Image.make_image()

Image = Imaging(f'../fits/80arcsec_full.fits')
Image.make_image()
