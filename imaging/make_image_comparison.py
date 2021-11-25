import numpy as np
import os
from astropy.wcs import WCS
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy import units as u
from argparse import ArgumentParser
import warnings
from pandas import DataFrame, concat, read_csv
import csv
import sys
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable

import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm

class Imaging:
    def __init__(self, fits_file: str = None, vmin: float = None, vmax: float = None):
        self.hdu = fits.open(fits_file)[0]
        self.image_data = self.hdu.data
        while len(self.image_data.shape) != 2:
            self.image_data = self.image_data[0]
        self.wcs = WCS(self.hdu.header, naxis=2)
        self.header = self.wcs.to_header()

        if vmin is None:
            self.vmin = np.nanstd(self.image_data)
        else:
            self.vmin = vmin
        if vmax is None:
            self.vmax = np.nanstd(self.image_data)*25
        else:
            self.vmax = vmax


    def imaging(self, image_data=None, cmap: str = 'CMRmap'):
        """
        Image your data with this method.
        image_data -> insert your image_data or plot full image
        cmap -> choose your preferred cmap
        """

        if image_data is None:
            image_data = self.image_data
        try:
            plt.figure(figsize=(10, 10))
            plt.subplot(projection=self.wcs)
            plt.imshow(image_data, norm=SymLogNorm(linthresh=self.vmin/20, vmin=self.vmin/50, vmax=self.vmax), origin='lower', cmap=cmap)
            plt.xlabel('Galactic Longitude')
            plt.ylabel('Galactic Latitude')
            plt.show()
        except:
            print('Error making images with matplotlib. Images will not be made.')

        return self

    def make_cutout(self, pos: tuple = None, size: tuple = (1000, 1000)):
        """
        Make cutout from your image with this method.
        pos (tuple) -> position in pixels
        size (tuple) -> size of your image in pixel size, default=(1000,1000)
        """
        out = Cutout2D(
            data=self.image_data,
            position=pos,
            size=size,
            wcs=self.wcs,
            mode='partial'
        )

        return out.data, out.wcs

vmin, vmax = 0.0005177041, 0.01294260291615501

image_1 = Imaging('image_000.app.restored.fits', vmin=vmin, vmax=vmax)
image_2 = Imaging('image_003.app.restored.fits', vmin=vmin, vmax=vmax)
image_3 = Imaging('image_004.app.restored.fits', vmin=vmin, vmax=vmax)
image_4 = Imaging('image_007.app.restored.fits', vmin=vmin, vmax=vmax)


fig, axes = plt.subplots(figsize=(10, 10), nrows=1, ncols=3, subplot_kw={'projection': image_1.wcs})
axes[0].imshow(image_1.image_data, norm=SymLogNorm(linthresh=image_1.vmin / 10, vmin=image_1.vmin / 10, vmax=image_1.vmax / 2),
           origin='lower',
           cmap='CMRmap')
axes[1].imshow(image_2.image_data, norm=SymLogNorm(linthresh=image_2.vmin / 10, vmin=image_2.vmin / 20, vmax=image_2.vmax),
           origin='lower',
           cmap='CMRmap')
im = axes[2].imshow(image_4.image_data, norm=SymLogNorm(linthresh=image_4.vmin / 10, vmin=image_4.vmin / 20, vmax=image_4.vmax),
           origin='lower',
           cmap='CMRmap')
# divider = make_axes_locatable(axes[2])
# cax2 = divider.append_axes("right", size="5%", pad=0.05)
cbar = fig.colorbar(im, ax=axes, orientation='horizontal', shrink=0.6)
cbar.set_label('Flux density [Jy/beam]')
plt.subplots_adjust(left=0.1, bottom=0.3, right=0.9, top=0.5, wspace=0.5, hspace=0.1)
plt.show()