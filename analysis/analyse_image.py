from past.utils import old_div
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
from astropy.nddata import Cutout2D
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import SymLogNorm, LogNorm, PowerNorm
from reproject import reproject_interp
import string
import sys
from astropy.modeling.models import Gaussian2D
from astropy.convolution import convolve, Gaussian2DKernel
from matplotlib.ticker import LogLocator, LogFormatterSciNotation as LogFormatter
import os
from scipy.ndimage import gaussian_filter, filters
import pyregion
from pyregion.mpl_helper import properties_func_default
from astropy.visualization.wcsaxes import WCSAxes
from matplotlib.patches import ConnectionPatch
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from radio_beam import Beams
import warnings
import scipy.ndimage as sn
from scipy.stats.stats import pearsonr, spearmanr, linregress
from scipy.optimize import curve_fit
from scipy import stats
from shapely.geometry import Polygon, Point
from shapely.ops import cascaded_union
from matplotlib.path import Path
from shapely.ops import cascaded_union

warnings.filterwarnings('ignore')
plt.style.use('seaborn-deep')

def flatten(f):
    """ Flatten a fits file so that it becomes a 2D image. Return new header and data """

    naxis = f[0].header['NAXIS']
    if naxis < 2:
        sys.exit('Can\'t make map from this')
    if naxis == 2:
        return fits.PrimaryHDU(header=f[0].header, data=f[0].data)

    w = WCS(f[0].header)
    wn = WCS(naxis=2)

    wn.wcs.crpix[0] = w.wcs.crpix[0]
    wn.wcs.crpix[1] = w.wcs.crpix[1]
    wn.wcs.cdelt = w.wcs.cdelt[0:2]
    wn.wcs.crval = w.wcs.crval[0:2]
    wn.wcs.ctype[0] = w.wcs.ctype[0]
    wn.wcs.ctype[1] = w.wcs.ctype[1]

    header = wn.to_header()
    header["NAXIS"] = 2
    copy = ('EQUINOX', 'EPOCH', 'BMAJ', 'BMIN', 'BPA', 'RESTFRQ', 'TELESCOP', 'OBSERVER')
    for k in copy:
        r = f[0].header.get(k)
        if r is not None:
            header[k] = r

    slice = []
    for i in range(naxis, 0, -1):
        if i <= 2:
            slice.append(np.s_[:], )
        else:
            slice.append(0)

    hdu = fits.PrimaryHDU(header=header, data=f[0].data[slice])
    return hdu

def reproject_hdu(hdu1, hdu2, outputname):
    array, _ = reproject_interp(hdu2, hdu1.header)
    fits.writeto(outputname, array, hdu1.header, overwrite=True)

    return

class Imaging:

    def __init__(self, fits_file: str = None, resolution: int = None):
        self.fitsfile = fits_file
        self.hdu = fits.open(fits_file)
        self.image_data = self.hdu[0].data
        while len(self.image_data.shape) != 2:
            self.image_data = self.image_data[0]
        self.wcs = WCS(self.hdu[0].header, naxis=2)
        self.header = self.wcs.to_header()
        self.rms = self.noise
        self.rms_full = self.rms.copy()
        self.cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)
        self.resolution = resolution
        self.beamarea_copy = self.beamarea

    @property
    def beamarea(self):
        try:
            # Given a fitsfile this calculates the beamarea in pixels

            if 'median' in self.fitsfile:
                hdu = fits.open(self.fitsfile.replace('median.fits', 'all.fits'))
                bmaj = hdu[0].header['BMAJ']
                bmin = hdu[0].header['BMIN']

                beammaj = bmaj / (2.0 * (2 * np.log(2)) ** 0.5)  # Convert to sigma
                beammin = bmin / (2.0 * (2 * np.log(2)) ** 0.5)  # Convert to sigma
                pixarea = abs(hdu[0].header['CDELT1']/2 * hdu[0].header['CDELT2']/2)

                beamarea = 2 * np.pi * 1.0 * beammaj * beammin  # Note that the volume of a two dimensional gaus$
                beamarea_pix = beamarea / pixarea

                return beamarea_pix

            bmaj = self.hdu[0].header['BMAJ']
            bmin = self.hdu[0].header['BMIN']

            beammaj = bmaj / (2.0 * (2 * np.log(2)) ** 0.5)  # Convert to sigma
            beammin = bmin / (2.0 * (2 * np.log(2)) ** 0.5)  # Convert to sigma
            pixarea = abs(self.hdu[0].header['CDELT1'] * self.hdu[0].header['CDELT2'])

            beamarea = 2 * np.pi * 1.0 * beammaj * beammin  # Note that the volume of a two dimensional gaus$
            beamarea_pix = beamarea / pixarea

            return beamarea_pix
        except:
            return self.beamarea_copy

    @staticmethod
    def get_beamarea(hdu):

        bmaj = hdu[0].header['BMAJ']
        bmin = hdu[0].header['BMIN']

        beammaj = bmaj / (2.0 * (2 * np.log(2)) ** 0.5)  # Convert to sigma
        beammin = bmin / (2.0 * (2 * np.log(2)) ** 0.5)  # Convert to sigma
        pixarea = abs(hdu[0].header['CDELT1'] * hdu[0].header['CDELT2'])

        beamarea = 2 * np.pi * 1.0 * beammaj * beammin  # Note that the volume of a two dimensional gaus$
        beamarea_pix = beamarea / pixarea

        return beamarea_pix


    def medianfilter(self, kpc_scale=None, write=None):
        kernelsize = int(kpc_scale / (1000 * self.pix_to_size(0.072).value))
        print(f"Apply median kernel with size {kernelsize}x{kernelsize}")
        self.image_data = sn.median_filter(self.image_data, kernelsize)
        if write:
            self.hdu[0].data = np.expand_dims(np.expand_dims(self.image_data, axis=0), axis=0)
            try:
                self.hdu[0].writeto(write, overwrite=True)
            except:
                self.hdu.writeto(write, overwrite=True)
            print('Saved: '+write)
        self.rms = self.noise
        return self

    def rudnickfilter(self, kpc_scale=None, write=None, open=None):
        """Multi-resolution filtering of radio images
           technique described in Rudnick, 2002 https://iopscience.iop.org/article/10.1086/342499/pdf
           larry@umn.edu -- please contact for assistance, as needed
           code below courtesy of Viral Parekh, SARAO , vparekh@ska.ac.za

        technique creates a diffuse emission map, called “open”;
             for small scale features,  filtered = original_map - open

         pick a box size 3x the beam-size or 3x size of features you want to remove
         open map has an offset zero level - determine it and correct
         open map is at original resolution - units are the same in Jy/beam
         open map will show sharp edges to diffuse regions, but is boxy; convolve for aesthetics
         """

        kernelsize = kpc_scale/(1000*self.pix_to_size(0.072).value)
        print(f'Kernel size is {int(kernelsize)} pixels ({self.header["CDELT2"]*int(kernelsize)*3600})')
        mins = filters.minimum_filter(self.image_data, size=(kernelsize, kernelsize))
        openmp = filters.maximum_filter(mins, size=(kernelsize, kernelsize))
        if open:
            self.image_data = openmp
        else:
            self.image_data -= openmp
        if write:
            self.hdu[0].data = np.expand_dims(np.expand_dims(self.image_data, axis=0), axis=0)
            try:
                self.hdu[0].writeto(write, overwrite=True)
            except:
                self.hdu.writeto(write, overwrite=True)
            print('Saved: '+write)
        self.rms = self.noise
        return self

    def make_image(self, image_data=None, cmap: str = 'CMRmap', vmin=None, vmax=None, show_regions=None, wcs=None,
                   colorbar=True, save=None, text=None, subim=None, beam=True, give_scale=True, convolve=None, show_grid=None, ticks=None,
                   savefits=None, sub='', show_clustername=None, bigim=None):
        plt.style.use('ggplot')
        """
        Image your data with this method.
        image_data -> insert your image_data or plot full image
        cmap -> choose your preferred cmap
        """

        if image_data is None:
            image_data = self.image_data
        if vmin is None:
            vmin = self.rms
        else:
            vmin = vmin
        if vmax is None:
            vmax = self.rms*25

        if wcs is None:
            wcs = self.wcs

        if convolve:
            image_data = self.convolve_image(image_data, convolve)

        if show_regions is not None:
            fig = plt.figure(figsize=(7, 10), dpi=200)
            plt.subplot(projection=wcs)
            WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=wcs)

            if self.resolution==6 or cmap=='Blues':
                color='black'
            else:
                color='darkblue'

            def fixed_color(shape, saved_attrs):
                attr_list, attr_dict = saved_attrs
                attr_dict["color"] = color
                kwargs = properties_func_default(shape, (attr_list, attr_dict))

                return kwargs

            r = pyregion.open(show_regions).as_imagecoord(header=self.hdu[0].header)
            patch_list, artist_list = r.get_mpl_patches_texts(fixed_color)

            # fig.add_axes(ax)
            for patch in patch_list:
                print(patch)
                plt.gcf().gca().add_patch(patch)
            for artist in artist_list:
                if not subim:
                    plt.gca().add_artist(artist)
            if subim:
                plt.text(3000/3400*(1500+120)-30, 3000/3400*(2880-500), 'C', color='black', fontsize=10)
                plt.text(3000/3400*(2016.59+105), 3000/3400*(2821.64-130)+20, 'B', color='black', fontsize=10)
                plt.text(3000/3400*(1539.66+78)-60, 3000/3400*(2049.08-142)-50, 'D', color='black', fontsize=10)
                plt.text(3000/3400*(2309.1+175)+20, 3000/3400*(1808.42-150)-70, 'E', color='black', fontsize=10)
                plt.text(3000/3400*(2100.69+95)-30, 3000/3400*(1088.61-330)-180, 'F', color='black', fontsize=10)
                plt.text(3000/3400*(2114.85+100), 3000/3400*(3500-145)+130, 'A', color='black', fontsize=10)

            if show_clustername is not None:

                def white_color(shape, saved_attrs):
                    attr_list, attr_dict = saved_attrs
                    attr_dict["color"] = 'darkblue'
                    kwargs = properties_func_default(shape, (attr_list, attr_dict))

                    return kwargs

                r = pyregion.open('../regions/clusters.reg').as_imagecoord(header=self.hdu[0].header)
                patch_list, artist_list = r.get_mpl_patches_texts(white_color)

                # fig.add_axes(ax)
                for patch in patch_list:
                    plt.gcf().gca().add_patch(patch)
                for artist in artist_list:
                    plt.gca().add_artist(artist)


        elif show_grid is not None:
            fig = plt.figure(figsize=(7, 10), dpi=200)
            plt.subplot(projection=wcs)
            WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=wcs)
            if cmap=='Blues':
                objts = [('A399', 'orangered'), ('A401', 'orangered'), ('bridge', 'firebrick'), ('A399trail', 'yellow')]
            else:
                objts = [('A399', 'lightcyan'), ('A401', 'lightcyan'), ('bridge', 'green'), ('A399trail', 'yellow')]
            for area in objts:

                # try:

                def colorr(shape, saved_attrs):
                    attr_list, attr_dict = saved_attrs
                    attr_dict["color"] = area[1]
                    kwargs = properties_func_default(shape, (attr_list, attr_dict))

                    return kwargs

                # if area[0]=='trail':
                #     r = pyregion.open(f'../ptp_results_a399_trail/gridA399_cb_.reg').as_imagecoord(
                #         header=self.hdu[0].header)
                #     patch_list, artist_list = r.get_mpl_patches_texts(colorr)
                #     for n, patch in enumerate(patch_list):
                #         f = fits.open(f'../ptp_results_a399_trail/A399_results_cb_.fits')
                #         t = f[1].data
                #         t = t[(t['xray_sb'] > 0) & (t['radio1_sb'] > 0) & (t['xray_sb_err'] > 0) & (t['radio1_sb_err'] > 0)]
                #         if n in [int(l.replace('xaf_','')) for l in list(t['region_name'])]:
                #             plt.gca().add_patch(patch)
                # else:
                if 'rudnick' in self.fitsfile:
                    ext = 'rudnick'
                elif 'cleanbridge' in self.fitsfile:
                    ext = 'cb'
                else:
                    ext = 'cb'
                # try:

                r = pyregion.open(f'../regions2/grid{area[0]}_{ext}{sub}.reg').as_imagecoord(header=self.hdu[0].header)
                # except:
                #     r = pyregion.open(f'../regions/grid{area[0]}_{ext}_1.reg').as_imagecoord(
                #         header=self.hdu[0].header)
                #     sub="_1"
                patch_list, artist_list = r.get_mpl_patches_texts(colorr)

                for n, patch in enumerate(patch_list):
                    if 'trail' in area[0]:
                        f = fits.open(f'../ptp_results_27/{area[0]}_results_{ext}_60.fits')
                    else:
                        f = fits.open(f'../ptp_results_27/{area[0]}_results_{ext}{sub}.fits')
                    t = f[1].data
                    t = t[(t['xray_sb'] > 0) & (t['radio1_sb'] > 0) & (t['xray_sb_err'] > 0) & (t['radio1_sb_err'] > 0)]
                    if n in [int(l.replace('xaf_','')) for l in list(t['region_name'])]:
                        plt.gca().add_patch(patch)
                # except:
                #     print(f'{area[0]} does not exist.')

        else:
            figure, ax = plt.subplots(figsize=(7, 10), dpi=200)
            plt.subplot(projection=wcs)
        im = plt.imshow(image_data, origin='lower', cmap=cmap)
        # if self.resolution>20:
        #     # im.set_norm(SymLogNorm(linthresh=vmin*10, vmin=vmin, vmax=vmax, base=10))
        #     if cmap in ['Blues']:
        #         im.set_norm(SymLogNorm(linthresh=self.rms*2, vmin=self.rms, vmax=self.rms*15, base=10))
        #     else:
        #         im.set_norm(PowerNorm(vmin=self.rms*1.8, vmax=vmax, gamma=1 / 2))
        #
        # else:
        #     # im.set_norm(SymLogNorm(linthresh=vmin*10, vmin=vmin, vmax=vmax, base=10))
        #     # im.set_norm(PowerNorm(vmin=vmin, vmax=vmax, gamma=1/2))
        if cmap=='Blues':
            im.set_norm(SymLogNorm(linthresh = self.rms * 2, vmin=self.rms, vmax = self.rms * 15, base = 10))
        else:
            im.set_norm(PowerNorm(vmin=0, vmax=vmax, gamma=1 / 2))

        if bigim:
            plt.xlabel('Right Ascension (J2000)', size=11)
            plt.ylabel('Declination (J2000)', size=11)
            plt.tick_params(axis='both', which='major', labelsize=10)
        else:
            plt.xlabel('Right Ascension (J2000)', size=14)
            plt.ylabel('Declination (J2000)', size=14)
            plt.tick_params(axis='both', which='major', labelsize=12)
        if colorbar:
            if ticks:
                cbar = plt.colorbar(orientation='horizontal', shrink=1, ticks=ticks)
            else:
                cbar = plt.colorbar(orientation='horizontal', shrink=1)
            # cbar.ax.set_xscale('log')
            # cbar.locator = LogLocator()
            if bigim:
                cbar.set_label('Surface brightness [Jy/beam]', size=11)

            else:
                cbar.set_label('Surface brightness [Jy/beam]', size=14)
            # cbar.formatter = LogFormatter()
            cbar.update_normal(im)
        if text:
            if self.resolution==60:
                plt.text(3450/4 - (6000/4-image_data.shape[0])/2, 2450/4 - (6000/4-image_data.shape[1])/2, 'A399', color='pink', fontsize=14)
                plt.text(2200/4 - (6000/4-image_data.shape[0])/2, 3750/4 - (6000/4-image_data.shape[1])/2, 'A401', color='pink', fontsize=14)
            if self.resolution==20:
                plt.text(3450/2 - (6000/2-image_data.shape[0])/2, 2450/2 - (6000/2-image_data.shape[1])/2, 'A399', color='pink', fontsize=14)
                plt.text(2200/2 - (6000/2-image_data.shape[0])/2, 3750/2 - (6000/2-image_data.shape[1])/2, 'A401', color='pink', fontsize=14)


        if give_scale:
            if cmap=='Blues':
                c = 'brown'
            else:
                c = 'black'

            Mpc_pixels = 1 / (abs((self.header['CDELT2'] * u.deg).to(u.rad).value) * self.cosmo.angular_diameter_distance(0.072)).value
            l1, l2 = [int(image_data.shape[0]*0.1), int(image_data.shape[0]*0.1+Mpc_pixels)], [int(image_data.shape[1]*0.9), int(image_data.shape[1]*0.9)]
            plt.plot(l1, l2, color=c, linewidth=1)
            plt.text((l1[0]+l1[1])/2, 1.02*(l2[0]+l2[1])/2, '1 Mpc', color=c, fontsize=10, horizontalalignment='center')


        plt.xlim(0, image_data.shape[0])
        plt.ylim(0, image_data.shape[1])


        if type(self.resolution) == int and beam:
            beampix = self.resolution / (self.header['CDELT2'] * u.deg).to(u.arcsec).value/2
            x, y = beampix * 1.5 + image_data.shape[0] * 0.03, beampix * 1.5 + image_data.shape[1] * 0.03
            circle = plt.Circle((x, y), beampix, color='g',
                                fill=True)
            if self.resolution<=6:
                rectanglefill = plt.Rectangle(
                    (x - beampix * 5 / 2, y - beampix * 5 / 2), beampix * 5,
                                                            beampix * 5, fill=True, color='white')
                rectangle = plt.Rectangle(
                    (x - beampix * 5 / 2, y - beampix * 5 / 2), beampix * 5,
                                                                beampix * 5, fill=False, color='black', linewidth=1)
            else:
                rectanglefill = plt.Rectangle(
                    (x - beampix * 3 / 2, y - beampix * 3 / 2), beampix * 3,
                                                            beampix * 3, fill=True, color='white')
                rectangle = plt.Rectangle(
                    (x - beampix * 3 / 2, y - beampix * 3 / 2), beampix * 3,
                                                                beampix * 3, fill=False, color='black', linewidth=1)
            plt.gcf().gca().add_artist(rectangle)
            plt.gcf().gca().add_artist(rectanglefill)
            plt.gcf().gca().add_artist(circle)

        plt.grid(False)

        if bigim:
            plt.rc('font', size=12)

        if savefits:
            self.header['CRVAL3'] = 143650817.871094
            image_data = np.expand_dims(np.expand_dims(image_data, axis=0), axis=0)
            fits.writeto(savefits, image_data, self.header, overwrite=True)

        if save:
            plt.savefig(save, dpi=250, bbox_inches="tight")
            plt.close()

        else:
            plt.show()

        return self

    def make_subimages(self, regionfile, cmap='CMRmap', save=None, beamsize=None, convolve=None):

        r = pyregion.open(regionfile).as_imagecoord(header=self.hdu[0].header)

        fig = plt.figure(figsize=(9, 15))
        fig.subplots_adjust(hspace=0.2, wspace=0.4)

        rows, cols = len(r)//2, 2

        for k, shape in enumerate(r):

            if k==0:
                k=1
            elif k==1:
                k=3
            elif k==2:
                k=4
            elif k==3:
                k=5
            elif k==4:
                k=0
            elif k==5:
                k=2

            s = np.array(shape.coord_list)

            out = Cutout2D(
                data=self.image_data,
                position=(s[0], s[1]),
                size=(s[3], s[2]),
                wcs=self.wcs,
                mode='partial'
            )
            norm = SymLogNorm(linthresh=self.rms * 5, vmin=self.rms, vmax=self.rms*30, base=10)

            if convolve:
                image_data = self.convolve_image(out.data, convolve)
            else:
                image_data = out.data

            plt.subplot(rows, cols, k+1, projection=out.wcs)
            im = plt.imshow(image_data, origin='lower', cmap=cmap, norm=norm)
            if k%2==0 and self.resolution==6:
                plt.ylabel('Declination (J2000)', size=14)
            else:
                plt.ylabel(' ')

            if k>=4:
                plt.xlabel('Right Ascension (J2000)', size=14)
            else:
                plt.xlabel(' ')
            plt.tick_params(axis='both', which='major', labelsize=12)

            if type(self.resolution) == int and beamsize:
                beampix = self.resolution / (self.header['CDELT2'] * u.deg).to(u.arcsec).value/2
                x, y = beampix*1.5+out.data.shape[0]*0.03, beampix*1.5+out.data.shape[1]*0.03
                circle = plt.Circle((x, y), beampix, color='g',
                                    fill=True)
                rectanglefill = plt.Rectangle(
                    (x - beampix*3/2, y - beampix*3/2), beampix * 3,
                    beampix * 3, fill=True, color='white')
                rectangle = plt.Rectangle(
                    (x - beampix*3/2, y - beampix*3/2), beampix * 3,
                    beampix * 3, fill=False, color='black', linewidth=2)
                plt.gcf().gca().add_artist(rectangle)
                plt.gcf().gca().add_artist(rectanglefill)
                plt.gcf().gca().add_artist(circle)
                plt.grid(False)

        fig.subplots_adjust(top=0.8)
        cbar_ax = fig.add_axes([0.22, 0.88, 0.6, 0.03]) # l, b, w, h
        cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
        cbar.ax.set_xscale('log')
        cbar.locator = LogLocator()
        cbar.formatter = LogFormatter()
        cbar.update_normal(im)
        cbar.set_label('Surface brightness [Jy/beam]')

        plt.grid(False)
        if save:
            plt.savefig(save, dpi=250, bbox_inches="tight")
            plt.close()

        else:
            plt.show()

    def make_multip_cuts(self, regionfile, cmap='CMRmap', save=None, beamsize=None, vmin=None, vmax=None):

        if vmin is None:
            vmin = self.rms
        else:
            vmin = vmin
        if vmax is None:
            vmax = self.rms*25

        r = pyregion.open(regionfile).as_imagecoord(header=self.hdu[0].header)

        fig = plt.figure(figsize=(9, 15))
        fig.subplots_adjust(hspace=0.2, wspace=0.4)

        for k, shape in enumerate(r):

            s = np.array(shape.coord_list)

            out = Cutout2D(
                data=self.image_data,
                position=(s[0], s[1]),
                size=(s[3], s[2]),
                wcs=self.wcs,
                mode='partial'
            )
            figure, ax = plt.subplots(figsize=(7, 10), dpi=200)
            plt.subplot(projection=out.wcs)
            im = plt.imshow(out.data, origin='lower', cmap=cmap)
            im.set_norm(PowerNorm(vmin=0, vmax=vmax, gamma=1/2))
            plt.xlabel('Right Ascension (J2000)', size=14)
            plt.ylabel('Declination (J2000)', size=14)
            plt.tick_params(axis='both', which='major', labelsize=12)
            plt.grid(False)
            cbar = plt.colorbar(orientation='horizontal', shrink=1, ticks=[7e-4, 1e-3])
            # cbar.ax.set_xscale('log')
            cbar.locator = LogLocator()
            cbar.formatter = LogFormatter()
            cbar.update_normal(im)
            cbar.set_label('Surface brightness [Jy/beam]', size=14)

            if save:
                plt.savefig('subims/'+str(k)+save, dpi=250, bbox_inches="tight")
                plt.close()

            else:
                plt.show()

    def compare_lotss(self, regionfile, cmap='CMRmap', save=None, beamsize=None, vmin=None, vmax=None, colorbar=None):

        hdu2 = fits.open('../fits/lotss.fits')
        image_data2 = hdu2[0].data[0][0]
        wcs2 = WCS(hdu2[0].header, naxis=2)

        if vmin is None:
            vmin = self.rms
        else:
            vmin = vmin
        if vmax is None:
            vmax = self.rms*25

        r = pyregion.open(regionfile).as_imagecoord(header=self.hdu[0].header)

        r2 = pyregion.open(regionfile).as_imagecoord(header=hdu2[0].header)

        fig = plt.figure(figsize=(9, 15))
        fig.subplots_adjust(hspace=0.2, wspace=0.4)

        for k, shape in enumerate(r):

            s = np.array(shape.coord_list)
            s2 = np.array(r2[k].coord_list)

            out = Cutout2D(
                data=self.image_data,
                position=(s[0], s[1]),
                size=(s[3], s[2]),
                wcs=self.wcs,
                mode='partial'
            )

            out2 = Cutout2D(
                data=image_data2,
                position=(s2[0], s2[1]),
                size=(s2[3], s2[2]),
                wcs=wcs2,
                mode='partial'
            )

            norm = PowerNorm(vmin=0, vmax=vmax, gamma=1/2)

            fig, axes = plt.subplots(figsize=(10, 10), nrows=1, ncols=2, subplot_kw={'projection': out.wcs},
                                     sharey='all')
            fig.subplots_adjust(hspace=0.2, wspace=0.4)
            axes[0].imshow(out.data,
                           norm=norm,
                           origin='lower',
                           cmap='cubehelix_r')
            axes[0].set_xlabel('Right Ascension (J2000)', size=14)
            axes[0].set_ylabel('Declination (J2000)', size=14)
            axes[0].tick_params(axis='both', which='major', labelsize=12)
            axes[0].grid(False)
            im = axes[1].imshow(out2.data,
                                norm=norm,
                                origin='lower',
                                cmap='cubehelix_r')
            axes[1].set_xlabel('Right Ascension (J2000)', size=14)
            axes[1].tick_params(axis='both', which='major', labelsize=12)
            axes[1].set_yticks([])
            axes[1].set_ylabel(' ')
            axes[1].yaxis.set_visible(False)
            axes[1].grid(False)

            if colorbar:

                cbar = fig.colorbar(im, ax=axes, orientation='horizontal', shrink=1, ticks=[0., 1e-4, 0.0005, 1e-3])
                cbar.ax.tick_params(labelsize=14)
                cbar.set_label('Surface brightness [Jy/beam]', size=14)
                # cbar.locator = LogLocator()
                # cbar.formatter = LogFormatter()
                cbar.update_normal(im)

            # figure, ax = plt.subplots(figsize=(7, 10), dpi=200)
            # plt.subplot(projection=out.wcs)
            # im = plt.imshow(out.data, origin='lower', cmap=cmap)
            # im.set_norm()
            # plt.xlabel('Right Ascension (J2000)', size=14)
            # plt.ylabel('Declination (J2000)', size=14)
            # plt.tick_params(axis='both', which='major', labelsize=12)
            # plt.grid(False)
            # cbar = plt.colorbar(orientation='horizontal', shrink=1, ticks=[7e-4, 1e-3])
            # # cbar.ax.set_xscale('log')

            # cbar.set_label('Surface brightness [Jy/beam]', size=14)

            if save:
                plt.savefig('subims/'+str(k)+save, dpi=250, bbox_inches="tight")
                plt.close()

            else:
                plt.show()


    def make_subcontour(self, regionfile, save=None, fits_lowres=None, beamsize=None):

        r = pyregion.open(regionfile).as_imagecoord(header=self.hdu[0].header)

        maxlevel=200*self.rms_full
        minlevel=self.rms_full*5

        if fits_lowres:
            self.reproject_map(fits_lowres, 'temp.fits')
            hdu2 = fits.open('temp.fits')[0]

            image_data2 = hdu2.data
            while len(image_data2.shape) != 2:
                image_data2 = image_data2[0]
            wcs2 = WCS(hdu2.header, naxis=2)
            rms2 = self.get_noise(image_data2)
            levelslowres = [rms2 * 5]
            for _ in range(10):
                levelslowres.append(levelslowres[-1] * 2)
            r2 = pyregion.open(regionfile).as_imagecoord(header=hdu2.header)

        levels = [minlevel]
        for _ in range(10):
            levels.append(levels[-1]*2)

        levels2 = np.linspace(0.5*self.rms_full, maxlevel, 1000)

        fig = plt.figure(figsize=(9, 15))
        fig.subplots_adjust(hspace=0.2, wspace=0.4)

        rows, cols = len(r)//2, 2
        for k, shape in enumerate(r):

            if fits_lowres:
                s2 = np.array(r2[k].coord_list)

                out2 = Cutout2D(
                    data=image_data2,
                    position=(s2[0], s2[1]),
                    size=(s2[2], s2[3]),
                    wcs=wcs2,
                    mode='partial'
                )

            #REORDER see govoni et al. 2019
            if k==0:
                k=1
            elif k==1:
                k=3
            elif k==2:
                k=4
            elif k==3:
                k=5
            elif k==4:
                k=0
            elif k==5:
                k=2

            s = np.array(shape.coord_list)

            out = Cutout2D(
                data=self.image_data,
                position=(s[0], s[1]),
                size=(s[3], s[2]),
                wcs=self.wcs,
                mode='partial'
            )

            image_data = np.clip(out.data, a_min=levels2[0], a_max=maxlevel * 0.99)
            plt.subplot(rows, cols, k+1, projection=out.wcs)
            norm = LogNorm(vmin=minlevel, vmax=maxlevel)

            plt.contour(out.data, levels, colors=('k'), linestyles=('-'), linewidths=(0.3,))
            im = plt.contourf(image_data, levels2, cmap='Blues', norm=norm)
            # plt.xlabel('Right Ascension (J2000)')
            plt.text(out.data.shape[1]*0.9, out.data.shape[0]*0.9, string.ascii_uppercase[k],
                     color='black', fontsize=14)

            if fits_lowres:
                plt.contour(out2.data, levelslowres, linestyles=('-'), linewidths=(0.3,), colors=['#A52A2A'])


            if k%2==0 and self.resolution==6:
                plt.ylabel('Declination (J2000)', size=14)
            else:
                plt.ylabel(' ')

            if k>=4:
                plt.xlabel('Right Ascension (J2000)', size=14)
            else:
                plt.xlabel(' ')

            if type(self.resolution) == int and beamsize:
                beampix = self.resolution / (self.header['CDELT2'] * u.deg).to(u.arcsec).value/2
                x, y = beampix*1.5+out.data.shape[0]*0.03, beampix*1.5+out.data.shape[1]*0.03
                circle = plt.Circle((x, y), beampix, color='g',
                                    fill=True)
                rectanglefill = plt.Rectangle(
                    (x - beampix*3/2, y - beampix*3/2), beampix * 3,
                    beampix * 3, fill=True, color='white')
                rectangle = plt.Rectangle(
                    (x - beampix*3/2, y - beampix*3/2), beampix * 3,
                    beampix * 3, fill=False, color='black', linewidth=2)
                plt.gcf().gca().add_artist(rectangle)
                plt.gcf().gca().add_artist(rectanglefill)
                plt.gcf().gca().add_artist(circle)

            def fixed_color(shape, saved_attrs):
                attr_list, attr_dict = saved_attrs
                attr_dict["color"] = "darkgreen"
                kwargs = properties_func_default(shape, (attr_list, attr_dict))

                return kwargs

            r3 = pyregion.open('../regions/optical.reg').as_imagecoord(header=out.wcs.to_header())

            optical_sources = [np.array(r3[i].coord_list) for i in range(len(r3))]
            optical_sources = [i for i in optical_sources if i[0]<s[3] and i[1]<s[2] and i[0]>0 and i[1]>0]
            plt.scatter([i[0] for i in optical_sources], [i[1] for i in optical_sources], color='red', marker='x', s=80)

            r4 = pyregion.open('../regions/sourcelabels.reg').as_imagecoord(header=out.wcs.to_header())

            patch_list, artist_list = r4.get_mpl_patches_texts(fixed_color)

            # fig.add_axes(ax)
            for patch in patch_list:
                print(patch)
                plt.gcf().gca().add_patch(patch)
            for artist in artist_list:
                plt.gca().add_artist(artist)


        fig.subplots_adjust(top=0.8)
        cbar_ax = fig.add_axes([0.22, 0.88, 0.6, 0.03]) # l, b, w, h
        cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal', ticks=[round(a, 5) for a in
                                   np.logspace(np.log(minlevel), np.log(maxlevel), 4, endpoint=True)])

        cbar.set_label('Surface brightness [Jy/beam]', size=14)
        cbar.ax.set_xscale('log')

        plt.grid(False)
        if save:
            plt.savefig(save, dpi=250, bbox_inches="tight")
            plt.close()

        else:
            plt.show()

    def make_3d_map(self):
        pass

    def pix_to_size(self, z):
        return abs((self.header['CDELT2'] * u.deg).to(u.rad).value) * self.cosmo.angular_diameter_distance(z)

    def make_polygon(self, points=None, subtract_points=None, sigma=None, sigma2=None, do_science=None, make_image=None, spectralindex=None, save=None, regionmask=None, size=None, make_cutout=True, savefits=None):
        if not points:
            sys.exit('Error: no point given in pix')

        if not sigma:
            sigma=3
            # 20*self.rms

        mask = self.image_data*0

        for point in points:
            point = WCS(self.header).wcs_world2pix(point[0], point[1], 1)

            cs = plt.contour(self.image_data, [self.rms*sigma], colors='white', linewidths=0.7)
            plt.close()
            cs_list = cs.collections[0].get_paths()
            correct_cs = []
            for cs in cs_list:
                if len(cs) > 2:
                    cs = cs.vertices
                    if Polygon(cs).contains(Point(point[0], point[1])):
                        correct_cs = cs

            x, y = np.meshgrid(np.arange(self.image_data.shape[0]), np.arange(self.image_data.shape[1]))
            x, y = x.flatten(), y.flatten()

            points = np.vstack((x, y)).T

            path = Path(correct_cs)
            grid = path.contains_points(points).astype(int)
            mask += grid.reshape(self.image_data.shape[0], self.image_data.shape[1])

            if subtract_points:
                for subtract_point in subtract_points:
                    if not sigma2:
                        sigma2=10
                    subtract_point = WCS(self.header).wcs_world2pix(subtract_point[0], subtract_point[1], 1)
                    cs = plt.contour(self.image_data, [self.rms * sigma2], colors='white', linewidths=0.7)
                    plt.close()
                    cs_list = cs.collections[0].get_paths()
                    correct_cs = []
                    for cs in cs_list:
                        if len(cs) > 2:
                            cs = cs.vertices
                            if Polygon(cs).contains(Point(subtract_point[0], subtract_point[1])):
                                correct_cs = cs

                    x, y = np.meshgrid(np.arange(self.image_data.shape[0]), np.arange(self.image_data.shape[1]))
                    x, y = x.flatten(), y.flatten()

                    points = np.vstack((x, y)).T
                    try:
                        path = Path(correct_cs)
                        grid = path.contains_points(points).reshape(self.image_data.shape[0], self.image_data.shape[1])
                        mask -= grid.astype(int)
                    except:
                        pass

        if regionmask:
            r = pyregion.open(regionmask).as_imagecoord(header=self.hdu[0].header)
            mask -= r.get_mask(hdu=self.hdu[0], shape=self.image_data.shape).astype(int)

        image_data = np.where(mask>0, 1, 0)*self.image_data

        if make_image:
            rms = self.rms
            if size is None:
                size = (180, 180)
            if make_cutout:
                self.make_cutout(image_data=image_data, pos=(int(point[0]), int(point[1])), size=size)
                # if 'Bridge' in save:
                #     self.make_image(vmin=rms, vmax=rms*25, save=save, beam=False, give_scale=False, savefits=savefits, show_regions='../regions/lines.reg')
                # else:
                self.make_image(vmin=rms, vmax=rms * 25, save=save, beam=False, give_scale=False, savefits=savefits)
            else:
                self.make_image(image_data=image_data, vmin=rms, vmax=rms*25, save=save, beam=False, give_scale=False, savefits=savefits)

        if do_science:
            if spectralindex:
                self.do_science(image_data=image_data, spectralindex=spectralindex, field='halo')
            else:
                self.do_science(image_data=image_data, field='bridge')


    def do_science(self, image_data=None, region='', results=None, spectralindex=1.5, field=None):

        z=0.072

        def radiopower(flux):
            d_L = self.cosmo.luminosity_distance(z)
            return (4 * np.pi * d_L ** 2. * ((1. + z) ** ((spectralindex) - 1.)) * flux).to(u.W / u.Hz)

        # def radiopower(obs_brightness):
        #     return (obs_brightness*4*np.pi*self.pix_to_size(z) ** 2*(1+z)**(0.7-1)* (1 + z) ** 4/ (
        #             (abs(self.header['CDELT1']*3600) * u.arcsec.to(u.rad)) ** 2)).to(u.W / u.Hz)

        if not image_data is None:
            image_data = image_data[image_data!=0]
            N_pixels_bridge = len(image_data)
            av_sb = np.median(image_data) * u.Jy/u.beam # average surface brigthness
            av_sb_arcsec = av_sb*u.beam / self.beamarea/(abs(self.hdu[0].header['CDELT1'] * 3600.)*u.arcsec)**2
            err_sb = np.std(image_data)/np.sqrt(N_pixels_bridge/self.beamarea)*u.Jy/u.beam # 1*sigma of every area element
            err_sb_arcsec = np.std(image_data)/np.sqrt(N_pixels_bridge/self.beamarea)*u.Jy/self.beamarea/(abs(self.hdu[0].header['CDELT1'] * 3600.)*u.arcsec * abs(self.hdu[0].header['CDELT1'] * 3600.)*u.arcsec)
            bridge_area = N_pixels_bridge * self.pix_to_size(0.072)**2 # from govoni et al
            integr_sb = N_pixels_bridge*av_sb_arcsec*(abs((self.header['CDELT2'] * u.deg).to(u.arcsec))**2)

            L=radiopower(integr_sb)

            print(f'\nArea: {bridge_area}, {len(image_data)*(self.header["CDELT1"]*u.deg).to(u.arcmin)**2}')
            print(f'# of pixels: {N_pixels_bridge}')
            print(f'Average surface brightness: {av_sb} $\pm$ {err_sb}')
            print(f'Average surface brightness: {av_sb_arcsec} $\pm$ {err_sb_arcsec}')
            flux_density_err = N_pixels_bridge * err_sb_arcsec * abs((self.header["CDELT2"] * u.deg).to(u.arcsec)) ** 2
            print(f'Total flux density is: {integr_sb} $\pm$ {flux_density_err}')
            radiopower_err = radiopower(flux_density_err)

            print(f'Radio power {L} $\pm$ {radiopower_err}')
            # volume = (1.3/2)**2*np.pi*3*(u.Mpc**3)
            if 'bridge' in field:
                volume = (bridge_area * 12.1 * u.Mpc)
                print(f'Emissivity is {(L / (volume)).to(u.erg / (u.s * u.cm * u.cm * u.cm * u.Hz))} $\pm$ {(radiopower_err / (volume)).to(u.erg / (u.s * u.cm * u.cm * u.cm * u.Hz))}')
            elif 'a399' in field.lower():
                volume = 4/3*np.pi*(3*0.186)*u.Mpc
                print(f'Emissivity is {(L / (volume)).to(u.erg / (u.s * u.cm * u.cm * u.cm * u.Hz))} $\pm$ {(radiopower_err / (volume)).to(u.erg / (u.s * u.cm * u.cm * u.cm * u.Hz))}')
            elif 'a401' in field.lower():
                volume = 4/3*np.pi*(3*0.109)*u.Mpc
                print(f'Emissivity is {(L / (volume)).to(u.erg / (u.s * u.cm * u.cm * u.cm * u.Hz))} $\pm$ {(radiopower_err / (volume)).to(u.erg / (u.s * u.cm * u.cm * u.cm * u.Hz))}')


        if results:
            print('\n' + results.split('_')[0])
            t = fits.open(results)[1].data
            t = t[(t['xray_sb'] > 0)
                  & (t['radio1_sb']>0)
                  & (t['xray_sb_err'] > 0)
                  & (t['radio1_sb_err'] > 0)]

            num_areas = len(t)

            t = t[(t['radio1_fluxdensity']<5*self.rms)]


            av_sb = np.mean(t['radio1_fluxdensity'])*u.Jy/u.beam # average surface brigthness
            av_sb_arcsec = np.mean(t['radio1_sb'])*u.Jy/(u.arcsec)**2
            err_sb = np.mean(t['radio1_fluxdensity_err'])*u.Jy/u.beam # 1*sigma of every area element
            err_sb_arcsec = np.mean(t['radio1_sb_err'])*u.Jy/(u.arcsec)**2
            # bridge_area = 1.3*3*(u.Mpc**2) # from govoni et al
            bridge_area = num_areas*(u.Mpc)**2*(105/(self.header['CDELT1']*u.deg).to(u.arcsec).value * (abs((self.header['CDELT2'] * u.deg).to(u.rad).value) * self.cosmo.angular_diameter_distance(0.072)).value)**2
            N_pixels_bridge = int(bridge_area/(self.pix_to_size(0.072)**2))
            N_beams = N_pixels_bridge/self.beamarea
            # integr_sb = N_pixels_bridge*av_sb_arcsec*(abs((self.header['CDELT2'] * u.deg).to(u.arcsec))**2)
            integr_sb = N_beams * av_sb.value * u.Jy

            L = radiopower(integr_sb)

            print(f'\n# of pixels: {N_pixels_bridge}')
            print(f'Bridge area: {bridge_area} or {(105*u.arcsec).to(u.arcmin)**2*num_areas}')
            print(f'Average surface brightness: {av_sb} $\pm$ {err_sb}')
            print(f'Average surface brightness: {av_sb_arcsec} $\pm$ {err_sb_arcsec}')
            flux_density_err = np.sqrt(N_pixels_bridge)*err_sb_arcsec*abs((self.header["CDELT2"] * u.deg).to(u.arcsec))**2
            print(f'Total flux density is: {integr_sb} $\pm$ {flux_density_err}')
            radiopower_err = radiopower(flux_density_err)
            print(f'Radio power {L} $\pm$ {radiopower_err}')
            # volume = (1.3/2)**2*np.pi*3*(u.Mpc**3)
            volume = (bridge_area*12.1*u.Mpc)
            print(f'Emissivity is {(L/(volume)).to(u.erg/(u.s*u.cm*u.cm*u.cm*u.Hz))} $\pm$ {(radiopower_err/(volume)).to(u.erg/(u.s*u.cm*u.cm*u.cm*u.Hz))}')

        if region:
            r = pyregion.open(region).as_imagecoord(header=self.hdu[0].header)
            mask = r.get_mask(hdu=self.hdu[0], shape=self.image_data.shape)
            if 'bridge' in region:
                image_data = np.where(self.image_data<self.rms*5, self.image_data, 0)
            else:
                image_data = np.where((self.image_data < self.rms * 5), self.image_data, 0)
            image_data = np.where(mask, image_data, 0)
            # self.make_image(image_data=image_data, colorbar=False, show_regions=region)
            image_data = image_data[image_data!=0]
            N_pixels_bridge = len(self.image_data[mask])
            av_sb = np.mean(image_data) * u.Jy/u.beam # average surface brigthness
            av_sb_arcsec = av_sb*u.beam / self.beamarea/(abs(self.hdu[0].header['CDELT1'] * 3600.)*u.arcsec * abs(self.hdu[0].header['CDELT1'] * 3600.)*u.arcsec)
            err_sb = np.std(image_data)/np.sqrt(N_pixels_bridge/self.beamarea)*u.Jy/u.beam # 1*sigma of every area element
            err_sb_arcsec = np.std(image_data)/np.sqrt(N_pixels_bridge/self.beamarea)*u.Jy/self.beamarea/(abs(self.hdu[0].header['CDELT1'] * 3600.)*u.arcsec * abs(self.hdu[0].header['CDELT1'] * 3600.)*u.arcsec)
            bridge_area = N_pixels_bridge * self.pix_to_size(0.072)**2 # from govoni et al
            integr_sb = N_pixels_bridge*av_sb_arcsec*(abs((self.header['CDELT2'] * u.deg).to(u.arcsec))**2)

            L=radiopower(integr_sb)
            print(f'\nArea: {bridge_area}, {N_pixels_bridge*(self.header["CDELT1"]*u.deg).to(u.arcmin)**2}')
            print(f'# of pixels: {N_pixels_bridge}')
            print(f'Average surface brightness: {av_sb} $\pm$ {err_sb}')
            print(f'Average surface brightness: {av_sb_arcsec} $\pm$ {err_sb_arcsec}')
            flux_density_err = N_pixels_bridge * err_sb_arcsec * abs((self.header["CDELT2"] * u.deg).to(u.arcsec)) ** 2
            print(f'Total flux density is: {integr_sb} $\pm$ {flux_density_err}')
            radiopower_err = radiopower(flux_density_err)

            print(f'Radio power {L} $\pm$ {radiopower_err}')
            # volume = (1.3/2)**2*np.pi*3*(u.Mpc**3)
            volume = (bridge_area * 12.1 * u.Mpc)
            print(f'Emissivity is {(L / (volume)).to(u.erg / (u.s * u.cm * u.cm * u.cm * u.Hz))} $\pm$ {(radiopower_err / (volume)).to(u.erg / (u.s * u.cm * u.cm * u.cm * u.Hz))}')

    def convolve_image(self, image_data=None, sigma=None):
        # gauss_kernel = Gaussian2DKernel(sigma)
        # self.image_data = convolve(self.image_data, gauss_kernel)
        if image_data is None:
            image_data = self.image_data
        if sigma:
            image_data = gaussian_filter(image_data, sigma=sigma)
            self.rms = self.noise
        else:
            print('No tapering because no value given.')
        return image_data

    def reproject_map(self, input, output):

        hduflat1 = flatten(self.hdu)

        hdu2 = fits.open(input)
        hduflat2 = flatten(hdu2)

        reproject_hdu(hduflat1, hduflat2, output)

        return self

    def make_contourplot(self, image_data=None, wcs=None, title=None,
                         regions=None, scale=None, save=None, fits_lowres=None):
        if fits_lowres:
            self.reproject_map(fits_lowres, 'temp.fits')
            hdu2 = fits.open('temp.fits')[0]

            image_data2 = hdu2.data
            while len(image_data2.shape) != 2:
                image_data2 = image_data2[0]
            wcs2 = WCS(hdu2.header, naxis=2)
            rms2 = self.get_noise(image_data2)
            levelslowres = [rms2]
            for _ in range(10):
                levelslowres.append(levelslowres[-1] * 2)

        if image_data is None:
            image_data = self.image_data

        if wcs is None:
            wcs = self.wcs


        if regions:
            r = pyregion.open(regions).as_imagecoord(header=self.hdu[0].header)
            patch_list, artist_list = r.get_mpl_patches_texts()
            fig = plt.figure(figsize=(7, 10))

            ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=wcs)
            fig.add_axes(ax)
            for patch in patch_list:
                ax.add_patch(patch)
        else:
            plt.figure(figsize=(7, 10))
            plt.subplot(projection=wcs)

        maxlevel=200*self.rms_full
        minlevel=3*self.rms_full

        image_data = np.clip(image_data, a_min=np.min(image_data), a_max=maxlevel*0.99)

        levels = [minlevel]
        for _ in range(10):
            levels.append(levels[-1]*2)

        levels2 = np.linspace(minlevel/3, maxlevel, 1000)

        norm = LogNorm(vmin=minlevel, vmax=maxlevel)

        if fits_lowres:
            plt.contour(image_data2, levelslowres, linestyles=('-'), linewidths=(0.3,), colors=['#A52A2A'])

        plt.contour(image_data, levels, colors=('k'), linestyles=('-'), linewidths=(0.3,))
        plt.contourf(image_data, levels2, cmap='Blues', norm=norm)
        cbar = plt.colorbar(orientation='horizontal', shrink=1)
        cbar.set_label('Surface brightness [Jy/beam]', size=14)
        cbar.ax.set_xscale('log')
        # cbar.set_ticks([10e-4, 10e-3, 10e-2, 10e-1])
        if scale:
            Mpc_pixels = 1 / (abs((self.header['CDELT2'] * u.deg).to(u.rad).value) * self.cosmo.angular_diameter_distance(0.072)).value
            l1, l2 = [int(image_data.shape[0]*0.1), int(image_data.shape[0]*0.1+Mpc_pixels)], [int(image_data.shape[1]*0.9), int(image_data.shape[1]*0.9)]
            plt.plot(l1, l2, color='pink', linewidth=1)
            plt.text((l1[0] + l1[1]) / 2, 1.02 * (l2[0] + l2[1]) / 2, '1 Mpc', color='cyan', fontsize=10,
                     horizontalalignment='center')
        if type(self.resolution) == int:
            beampix = self.resolution / (self.header['CDELT2'] * u.deg).to(u.arcsec).value/2
            x, y = beampix * 1.5 + image_data.shape[0] * 0.03, beampix * 1.5 + image_data.shape[1] * 0.03
            circle = plt.Circle((x, y), beampix, color='g',
                                fill=True)
            rectangle = plt.Rectangle(
                (x - beampix * 3 / 2, y - beampix * 3 / 2), beampix * 3,
                                                            beampix * 3, fill=True, color='white')
            plt.gcf().gca().add_artist(rectangle)
            plt.gcf().gca().add_artist(circle)

        plt.xlabel('Right Ascension (J2000)', size=14)
        plt.ylabel('Declination (J2000)', size=14)
        plt.title(title)

        plt.grid(False)
        if save:
            plt.savefig(save, dpi=250, bbox_inches="tight")
            plt.close()

        else:
            plt.show()
        return self

    def get_xray(self, fitsfile, reproject=True, plot3d=None):
        if reproject:
            for f in [fitsfile, fitsfile.replace('.fits','_bkg.fits'), fitsfile.replace('.fits', '_exp.fits')]:
                fnew = f.replace('.fits','_reproj.fits')
                self.reproject_map(f, fnew)
                if 'bkg' in f:
                    hdu_bkg = fits.open(fnew)
                    data_bkg = hdu_bkg[0].data
                elif 'exp' in f:
                    hdu_exp = fits.open(fnew)
                    data_exp = hdu_exp[0].data
                else:
                    hdu = fits.open(fnew)
                    data = hdu[0].data

        else:
            hdu = fits.open(fitsfile)
            hdu_bkg = fits.open(fitsfile.replace('.fits', '_bkg.fits'))
            hdu_exp = fits.open(fitsfile.replace('.fits', '_exp.fits'))

            data = hdu[0].data
            data_bkg = hdu_bkg[0].data
            data_exp = hdu_exp[0].data

        data[np.isnan(data)] = 0
        data_bkg[np.isnan(data_bkg)] = 0
        data_exp[np.isnan(data_exp)] = 1
        data_exp[data_exp==0] = 1

        if plot3d:
            return data, data_bkg, data_exp

        return (data - data_bkg)/data_exp

    def make_bridge_overlay_yxr_contourplot(self, fits1, fits2=None, save=None, beamsize=True, show_regions=None):
        plt.style.use('ggplot')

        self.reproject_map(fits1, 'fits1.fits')
        # self.reproject_map(fits2, 'fits2.fits')


        hdu1 = fits.open('fits1.fits')
        image_data1 = self.get_xray(fits1, reproject=True)
        # hdu2 = fits.open('fits2.fits')
        # image_data2 = hdu2[0].data

        while len(image_data1.shape) != 2:
            image_data1 = image_data1[0]
        image_data1 = self.convolve_image(image_data1, 15)
        # plt.figure(figsize=(7, 10))
        # plt.subplot(projection=self.wcs)
        # plt.imshow(image_data1, norm = LogNorm(vmin=self.get_noise(image_data1), vmax=self.get_noise(image_data1)*10))
        # plt.savefig('xray.png')
        wcs_1 = WCS(hdu1[0].header)
        # while len(image_data2.shape) != 2:
        #     image_data2 = image_data2[0]
        # image_data2 = self.convolve_image(image_data2, 30)
        # wcs_2 = WCS(hdu2[0].header)


        plt.figure(figsize=(7, 10))
        plt.subplot(projection=self.wcs)

        if show_regions is not None:
            fig = plt.figure(figsize=(7, 10), dpi=200)
            plt.subplot(projection=self.wcs)
            WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=self.wcs)

            def fixed_color(shape, saved_attrs):
                attr_list, attr_dict = saved_attrs
                attr_dict["color"] = "red"
                kwargs = properties_func_default(shape, (attr_list, attr_dict))

                return kwargs

            r = pyregion.open(show_regions).as_imagecoord(header=self.hdu[0].header)
            patch_list, artist_list = r.get_mpl_patches_texts(fixed_color)

            # fig.add_axes(ax)
            for patch in patch_list:
                print(patch)
                plt.gcf().gca().add_patch(patch)
            for artist in artist_list:
                plt.gca().add_artist(artist)


            # def fixed_color(shape, saved_attrs):
            #     attr_list, attr_dict = saved_attrs
            #     attr_dict["color"] = "red"
            #     kwargs = properties_func_default(shape, (attr_list, attr_dict))
            #
            #     return kwargs

            # r = pyregion.open('../regions/slices_inverse.reg').as_imagecoord(header=self.hdu[0].header)
            # patch_list, artist_list = r.get_mpl_patches_texts(fixed_color)
            #
            # # fig.add_axes(ax)
            # for patch in patch_list:
            #     print(patch)
            #     plt.gcf().gca().add_patch(patch)
            # for artist in artist_list:
            #     plt.gca().add_artist(artist)
            #
            # def fixed_color(shape, saved_attrs):
            #     attr_list, attr_dict = saved_attrs
            #     attr_dict["color"] = "orange"
            #     kwargs = properties_func_default(shape, (attr_list, attr_dict))
            #
            #     return kwargs

            # r = pyregion.open('../regions/a399extensionslices.reg').as_imagecoord(header=self.hdu[0].header)
            # patch_list, artist_list = r.get_mpl_patches_texts(fixed_color)
            #
            # # fig.add_axes(ax)
            # for patch in patch_list:
            #     print(patch)
            #     plt.gcf().gca().add_patch(patch)
            # for artist in artist_list:
            #     plt.gca().add_artist(artist)

        levels_1 = [self.get_noise(image_data1)]
        for _ in range(10):
            levels_1.append(levels_1[-1]*2)

        levels_2 = [(10**(-5))/2]
        for _ in range(10):
            levels_2.append(levels_2[-1]*np.sqrt(2))

        plt.contour(image_data1, levels_1, colors=('darkgreen'), linestyles=('-'), linewidths=(1,))
        # plt.contour(image_data2, levels_2, colors=('orange'), linestyles=('-'), linewidths=(2,))


        norm = SymLogNorm(linthresh = self.rms * 2, vmin=self.rms, vmax = self.rms * 15, base = 10)
        # norm = PowerNorm(vmin=0, vmax=self.rms*25, gamma=1/2)
        plt.imshow(self.image_data, cmap='Blues', norm=norm)
        # plt.contourf(self.image_data, levels, cmap='Blues', norm=norm)
        ticks = [0, 1e-3, 1e-2]
        cbar = plt.colorbar(orientation='horizontal', shrink=1, ticks=ticks)
        cbar.set_label('Surface brightness  [Jy/beam]', size=14)
        cbar.formatter = LogFormatter()
        # cbar.ax.set_xscale('log')
        plt.xlabel('Right Ascension (J2000)')
        plt.ylabel('Declination (J2000)')

        Mpc_pixels = 1 / (abs((self.header['CDELT2'] * u.deg).to(u.rad).value) * self.cosmo.angular_diameter_distance(0.072)).value
        l1, l2 = [int(self.image_data.shape[0]*0.05), int(self.image_data.shape[0]*0.05+Mpc_pixels)], [int(self.image_data.shape[1]*0.9), int(self.image_data.shape[1]*0.9)]
        plt.plot(l1, l2, color='brown', linewidth=1)
        plt.text((l1[0]+l1[1])/2, 1.02*(l2[0]+l2[1])/2, '1 Mpc', color='brown', fontsize=10, horizontalalignment='center')
        if type(self.resolution) == int and beamsize:
            beampix = self.resolution / (self.header['CDELT2'] * u.deg).to(u.arcsec).value/2
            x, y = beampix * 1.5 + self.image_data.shape[0] * 0.03, beampix * 1.5 + self.image_data.shape[1] * 0.03
            circle = plt.Circle((x, y), beampix, color='g',
                                fill=True)
            rectanglefill = plt.Rectangle(
                (x - beampix * 3 / 2, y - beampix * 3 / 2), beampix * 3,
                                                            beampix * 3, fill=True, color='white')
            rectangle = plt.Rectangle(
                (x - beampix * 3 / 2, y - beampix * 3 / 2), beampix * 3,
                                                            beampix * 3, fill=False, color='black', linewidth=2)
            plt.gcf().gca().add_artist(rectangle)
            plt.gcf().gca().add_artist(rectanglefill)
            plt.gcf().gca().add_artist(circle)

        plt.grid(False)
        if save:
            plt.savefig(save, dpi=250, bbox_inches="tight")
            plt.close()
        else:
            plt.show()

    def make_bridge_overlay_contourplot(self, fitsfile=None, minlevel_1=None, maxlevel_1=None, maxlevel_2=None,
                                 minlevel_2=None, steps_1=1000, title=None, convolve_2=None, xray=False,
                                        save=None, steps_2=None):

        self.reproject_map(fitsfile, 'test.fits')

        hdu = fits.open('test.fits')
        if xray:
            image_data_2 = self.get_xray(fitsfile, reproject=True)
        else:
            image_data_2 = hdu[0].data

        while len(image_data_2.shape) != 2:
            image_data_2 = image_data_2[0]
        wcs_2 = WCS(hdu[0].header)

        if convolve_2:
            image_data_2 = self.convolve_image(image_data_2, convolve_2)

        if maxlevel_2 is None:
            maxlevel_2 = np.max(image_data_2)

        if minlevel_2 is None:
            minlevel_2 = self.get_noise(image_data_2)
            if xray:
                minlevel_2*=10

        if maxlevel_1 is None:
            maxlevel_1 = np.max(self.image_data)

        if minlevel_1 is None:
            minlevel_1 = self.rms*3


        # image_data_2 = np.clip(image_data_2, a_min=0, a_max=maxlevel_2)
        plt.figure(figsize=(7, 10))
        plt.subplot(projection=wcs_2)

        levels_2 = [minlevel_2]
        if steps_2 is None:
            steps_2=2
        for _ in range(10):
            levels_2.append(levels_2[-1]*steps_2)

        levels_1 = np.linspace(minlevel_1, maxlevel_1, steps_1)
        # print(self.rms*3)
        plt.contour(image_data_2, levels_2, colors=('r'), linestyles=('-'), linewidths=(1,))
        plt.imshow(np.where(self.image_data<self.rms*3, 0, 1), cmap='Greys')
        plt.contourf(np.clip(self.image_data, a_min=levels_1[0], a_max=np.max(self.image_data)), levels_1, cmap='Blues')
        cbar = plt.colorbar(orientation='horizontal', shrink=1, ticks=[0.0, 0.002, 0.004])
        cbar.set_label('Surface brightness  [Jy/beam]', size=14)
        # cbar.ax.set_xscale('log')
        plt.xlabel('Right Ascension (J2000)')
        plt.ylabel('Declination (J2000)')

        Mpc_pixels = 1 / (abs((self.header['CDELT2'] * u.deg).to(u.rad).value) * self.cosmo.angular_diameter_distance(0.072)).value
        l1, l2 = [int(self.image_data.shape[0]*0.05), int(self.image_data.shape[0]*0.05+Mpc_pixels)], [int(self.image_data.shape[1]*0.9), int(self.image_data.shape[1]*0.9)]
        plt.plot(l1, l2, color='brown', linewidth=1)
        plt.text((l1[0]+l1[1])/2, 1.02*(l2[0]+l2[1])/2, '1 Mpc', color='brown', fontsize=10, horizontalalignment='center')
        # if type(resolution) == int:
        #     beampix = resolution / (self.header['CDELT2'] * u.deg).to(u.arcsec).value/2
        #     circle = plt.Circle((image_data.shape[0] * 0.03, image_data.shape[1] * 0.03), beampix, color='g',
        #                         fill=True)
        #     rectangle = plt.Rectangle(
        #         (image_data.shape[0] * 0.03 - beampix * 4, image_data.shape[1] * 0.03 - beampix * 4), beampix * 8,
        #                                                                                           beampix * 8,
        #         fill=True, color='white')
        #     plt.gcf().gca().add_artist(rectangle)
        #     plt.gcf().gca().add_artist(circle)
        plt.title(title)

        plt.grid(False)
        if save:
            plt.savefig(save, dpi=250, bbox_inches="tight")
            plt.close()
        else:
            plt.show()

    def analyse_corr(self, savefig=None, grid=None, A1758=False):

        def pearsonr_ci(x, y, alpha=0.05):
            ''' calculate Pearson correlation along with the confidence interval using scipy and numpy
            Parameters
            ----------
            x, y : iterable object such as a list or np.array
              Input for correlation calculation
            alpha : float
              Significance level. 0.05 by default
            Returns
            -------
            r : float
              Pearson's correlation coefficient
            pval : float
              The corresponding p value
            lo, hi : float
              The lower and upper bound of confidence intervals
            '''

            r, p = pearsonr(x, y)
            r_z = np.arctanh(r)
            se = 1 / np.sqrt(x.size - 3)
            z = stats.norm.ppf(1 - alpha / 2)
            lo_z, hi_z = r_z - z * se, r_z + z * se
            lo, hi = np.tanh((lo_z, hi_z))
            return r, p, lo, hi

        def spearmanr_ci(x, y, alpha=0.05):
            ''' calculate Spearman correlation along with the confidence interval using scipy and numpy
            Parameters
            ----------
            x, y : iterable object such as a list or np.array
              Input for correlation calculation
            alpha : float
              Significance level. 0.05 by default
            Returns
            -------
            r : float
              Spearman's correlation coefficient
            pval : float
              The corresponding p value
            lo, hi : float
              The lower and upper bound of confidence intervals
            '''

            r, p = spearmanr(x, y)
            r_z = np.arctanh(r)
            se = 1 / np.sqrt(x.size - 3)
            z = stats.norm.ppf(1 - alpha / 2)
            lo_z, hi_z = r_z - z * se, r_z + z * se
            lo, hi = np.tanh((lo_z, hi_z))
            return r, p, lo, hi

        def objective(x, a, b):
            return a * x + b

        def fit(x, y):
            popt, _ = curve_fit(objective, x, y)
            a, b = popt
            x_line = np.arange(min(x) - 10, max(x) + 10, 0.01)
            y_line = objective(x_line, a, b)
            print('y = %.5f * x + %.5f' % (a, b))
            return x_line, y_line

        def linreg(x, y):
            res = linregress(x, y)
            print(f'Slope is {res.slope} +- {res.stderr}')
            print(f'Based on {len(x)} elements')
            return res.slope, res.stderr

        if grid:
            radio = np.load(grid.lower()+'radio.npy')
            radio2 = np.load(grid.lower()+'radio2.npy')
            y = np.power(np.load(grid.lower()+'y.npy'), 2)
            xray = np.load(grid.lower()+'xray.npy')

            radio_err = np.load(grid.lower()+'radio_err.npy')
            radio2_err = np.load(grid.lower()+'radio2_err.npy')
            y_err = 2*np.sqrt(y)*np.load(grid.lower()+'y_err.npy')
            xray_err = np.load(grid.lower()+'xray_err.npy')

        elif A1758:
            import pandas as pd
            df = pd.read_csv('A1758.csv')
            radio = df.radio_sb
            radio_err = df.radio_sb_err
            xray = df.xray_sb
            xray_err = df.xray_sb_err

            print('radio vs. xray')
            slopex, errx = linreg(np.log10(xray), np.log10(radio))
            fitxray = fit(np.log10(xray), np.log10(radio))
            print('Pearson R (x-ray vs radio): ' + str(pearsonr_ci(np.log(xray), np.log(radio))))
            print('Spearman R (x-ray vs radio): ' + str(spearmanr_ci(np.log(xray), np.log(radio))))

            fig, ax = plt.subplots(constrained_layout=True)
            ax.errorbar(np.log10(xray), np.log10(radio), xerr=(0.434 * xray_err / xray),
                        yerr=0.434 * radio_err / radio, fmt='.', ecolor='red', elinewidth=0.4,
                        color='darkred', capsize=2, markersize=4)

            plt.grid(False)
            # ax2 = ax.twiny()
            ax.set_ylim(np.min([np.min(np.log10(radio) - (0.434 * radio_err / radio)),
                                np.min(np.log10(radio) - (0.434 * radio_err / radio))]) - 0.05,
                        np.max([np.max(np.log10(radio) + (0.434 * radio_err / radio)),
                                np.max(np.log10(radio) + (0.434 * radio_err / radio))]) + 0.05)
            ax.set_xlim(np.min(np.min(np.log10(xray) - (xray_err / xray))) - 0.05,
                        np.max(np.max(np.log10(xray) + (xray_err / xray))) + 0.05)
            # ax.plot(fitline[0], fitline[1], color='darkslateblue', linestyle='--')
            # ax.set_xlim(-1, 1)
            # ax.set_ylim(-0.3, 0.45)
            ax.set_ylabel('log($I_{R}$) [SB/mean(SB)]')
            # ax.set_xlabel('X-ray [SB/mean(SB)]')
            ax.set_xlabel('log($I_{X}$) [SB/mean(SB)]')
            ax.legend(['Radio vs. X-ray'], loc='upper left')
            ax.plot(fitxray[0], fitxray[1], color='darkred', linestyle='--')

            plt.tight_layout()
            plt.grid(False)
            plt.savefig('A1758corr.png', bbox_inches='tight')

            return

        else:
            radio = np.load('radio3d.npy')
            y = np.power(np.load('y.npy'), 2)
            xray = np.load('xray.npy')

            radio_err = np.load('radio3d_err.npy').flatten()
            y_err = (2*np.sqrt(y)*np.load('y_err.npy')).flatten()
            xray_err = np.load('xray_err.npy').flatten()

            # X, Y = np.meshgrid(range(radio.shape[1]), range(radio.shape[0]))
            # hf = plt.figure()
            # ha = hf.add_subplot(111, projection='3d')
            # ha.plot_surface(X, Y, radio-1, color='red', alpha=0.5)
            # ha.plot_surface(X, Y, y-1, color='green', alpha=0.5)
            # ha.plot_surface(X, Y, np.clip(xray-1, a_max=1.5, a_min=-1.5), color='blue', alpha=0.5)
            # plt.xticks([])
            # plt.yticks([])
            # plt.grid('off')
            # if savefig:
            #     plt.savefig(savefig, bbox_inches="tight")
            # else:
            #     plt.show()

            radio = radio.flatten()
            xray = xray.flatten()
            y = y.flatten()

        msk = (radio>self.rms*3)

        radio=radio[msk]
        radio_err=radio_err[msk]
        # y=y[msk]
        # y_err=y_err[msk]
        xray=xray[msk]
        xray_err=xray_err[msk]


        radio_err[np.isnan(radio_err)] = 0
        radio2_err[np.isnan(radio2_err)] = 0
        xray_err[np.isnan(xray_err)] = 0
        y_err[np.isnan(y_err)] = 0

        print('Lin reg before dividing by mean:')

        linreg(np.log10(xray), np.log10(radio))
        linreg(np.log10(y), np.log10(radio2))

        radio_err/=np.mean(radio)
        radio/=np.mean(radio)
        radio2_err/=np.mean(radio2)
        radio2/=np.mean(radio2)
        y_err/=np.mean(y)
        y/=np.mean(y)
        xray_err/=np.mean(xray)
        xray/=np.mean(xray)

        mask = 0.434 * xray_err / xray < 0.5
        radio = radio[mask]
        xray = xray[mask]
        # y = y[mask]
        radio_err = radio_err[mask]
        xray_err = xray_err[mask]
        # y_err = y_err[mask]

        print('radio vs. xray')
        slopex, errx = linreg(np.log10(xray), np.log10(radio))
        fitxray = fit(np.log10(xray), np.log10(radio))
        print('radio vs y')
        slopey, erry =linreg(np.log10(y), np.log10(radio2))
        fity = fit(np.log10(y), np.log10(radio2))


        print('Pearson R (x-ray vs radio): ' + str(pearsonr(np.log(xray), np.log(radio))))
        print('Pearson R (ymap vs radio): ' + str(pearsonr(np.log(radio2), np.log(y))))

        print('Spearman R (x-ray vs radio): ' + str(spearmanr(np.log(xray), np.log(radio))))
        print('Spearman R (ymap vs radio): ' + str(spearmanr(np.log(radio2), np.log(y))))

        fig, ax = plt.subplots(constrained_layout=True)
        ax.errorbar(np.log10(xray), np.log10(radio), xerr=(0.434 * xray_err / xray),
                    yerr=0.434 * radio_err / radio, fmt='.', ecolor='red', elinewidth=0.4,
                    color='darkred', capsize=2, markersize=4)

        plt.grid(False)
        # ax2 = ax.twiny()
        ax.errorbar(np.log10(y), np.log10(radio2), xerr=(0.434 * y_err / y),
                     yerr=0.434 * radio2_err / radio2, fmt='.', ecolor='blue', elinewidth=0.4,
                     color='darkblue', capsize=2, markersize=4)
        ax.set_ylim(np.min([np.min(np.log10(radio) - (0.434 * radio_err / radio)),
                            np.min(np.log10(radio) - (0.434 * radio_err / radio))]) - 0.1,
                    np.max([np.max(np.log10(radio) + (0.434 * radio_err / radio)),
                            np.max(np.log10(radio) + (0.434 * radio_err / radio))]) + 0.1)
        ax.set_xlim(np.min([np.min(np.log10(y) - (0.434 * y_err / y)),
                            np.min(np.log10(xray) - (0.434 * xray_err / xray))]) - 0.1,
                    np.max([np.max(np.log10(y) + (0.434 * y_err / y)),
                            np.max(np.log10(xray) + (0.434 * xray_err / xray))]) + 0.1)
        # ax.plot(fitline[0], fitline[1], color='darkslateblue', linestyle='--')
        # ax.set_xlim(-1, 1)
        # ax.set_ylim(-0.3, 0.45)
        ax.set_ylabel('log($I_{R}$) [SB/mean(SB)]')
        # ax.set_xlabel('X-ray [SB/mean(SB)]')
        ax.set_xlabel('log($I_{X}$) [SB/mean(SB)] and log($I_{SZ}^{2}$) [SZ/mean(SZ)]')
        ax.legend(['Radio vs. X-ray', 'Radio vs. SZ'], loc='upper left')
        ax.plot(fity[0], fity[1], color='darkblue', linestyle='--')
        ax.plot(fitxray[0], fitxray[1], color='darkred', linestyle='--')

        plt.tight_layout()
        plt.grid(False)
        plt.savefig(savefig, bbox_inches='tight')

    def ps(self, xray=None, region=None, save=None, rudnick=None):
        if region is None:
            region = f'../regions/lines.reg'
        if save is None:
            save = 'slicesbridge.png'
        region = open(region).read()
        region = region.split('\n')
        region_head = region[0:2]
        structures = region[3:]
        data_radio_cb = []
        data_radio_err_cb = []
        if rudnick:
            hdu_rud = fits.open('../fits/60rudnick.fits')
            image_data_rud = hdu_rud[0].data
            while len(image_data_rud.shape) != 2:
                image_data_rud = image_data_rud[0]
            data_radio_rudnick = []
            data_radio_err_rudnick = []
            rud_noise = self.get_noise(image_data_rud)
            rud_beam = self.get_beamarea(hdu_rud)
        data_xray = []
        data_xray_err = []
        for n, structure in enumerate(structures):
            if (not 'box' in structure) and (not 'annulus' in structure):
                continue
            r = pyregion.parse('\n'.join(region_head + [structure])).as_imagecoord(header=self.hdu[0].header)
            mask_radio = r.get_mask(hdu=self.hdu[0], shape=self.image_data.shape)
            # im_mask_radio = np.clip(self.image_data[mask_radio], a_max=self.rms*5, a_min=None)
            if 'bridge' in save:
                im_mask_radio_sub = np.where((self.image_data<self.rms*5), self.image_data, 0)[mask_radio]
            else:
                im_mask_radio_sub = self.image_data[mask_radio]

            im_mask_radio_sub = im_mask_radio_sub[im_mask_radio_sub!=0]
            fluxdens = len(self.image_data[mask_radio])*np.mean(im_mask_radio_sub)/self.beamarea
            data_radio_cb.append(fluxdens)
            err = np.sqrt((self.rms * np.sqrt(len(im_mask_radio_sub) / self.beamarea)) ** 2 +
                          (0.1 * fluxdens) ** 2)
            # err = np.std(im_mask_radio_sub)/np.sqrt(len(im_mask_radio_sub) / self.beamarea)
            data_radio_err_cb.append(err)

            if rudnick:
                r = pyregion.parse('\n'.join(region_head + [structure])).as_imagecoord(header=hdu_rud[0].header)
                mask_radio = r.get_mask(hdu=hdu_rud[0], shape=image_data_rud.shape)
                # im_mask_radio = np.clip(self.image_data[mask_radio], a_max=self.rms*5, a_min=None)
                image_data_rud_sub = np.where((image_data_rud < rud_noise * 4), image_data_rud, 0)[mask_radio]
                image_data_rud_sub = image_data_rud_sub[image_data_rud_sub != 0]
                fluxdens_rud = len(image_data_rud[mask_radio])*np.mean(image_data_rud_sub)/rud_beam
                data_radio_rudnick.append(fluxdens_rud)
                # err = np.sqrt((rud_noise * np.sqrt(len(image_data_rud_sub) / self.get_beamarea(hdu_rud))) ** 2 +
                #               (0.15 * np.sum(image_data_rud_sub) / self.get_beamarea(hdu_rud)) ** 2)/(len(im_mask_radio_sub) / self.beamarea)
                err = np.sqrt((rud_noise * np.sqrt(len(image_data_rud_sub) / rud_beam)) ** 2 +
                          (0.1 * fluxdens_rud) ** 2)

                data_radio_err_rudnick.append(err)


            if not xray is None:

                self.reproject_map(xray, 'test.fits')
                hdu = fits.open('test.fits')
                # wcs = WCS(hdu[0].header, naxis=2)
                image_data, xray_back, xray_exp = self.get_xray(xray, reproject=True, plot3d=True)
                r = pyregion.parse('\n'.join(region_head + [structure])).as_imagecoord(header=hdu[0].header)
                mask_xray = r.get_mask(hdu=hdu[0], shape=image_data.shape)
                im_mask_xray = image_data[mask_xray]
                xray_exp[xray_exp == 0.] = 1.
                xray_back[np.isnan(xray_back)] = 0.
                im_mask_xray[np.isnan(im_mask_xray)] = 0
                xray_exp_mask = xray_exp[mask_xray]
                xray_back_mask = xray_back[mask_xray]
                xray_exp_mask[np.isnan(xray_exp_mask)] = 1
                xray_back_mask[np.isnan(xray_back_mask)] = 0
                xr = (np.sum(im_mask_xray) - np.sum(xray_back_mask)) / np.mean(xray_exp_mask)
                xr_err = 4 * np.sqrt(np.sum(im_mask_xray) + np.sum(xray_back_mask)) / np.mean(xray_exp_mask)
                data_xray.append(xr)
                data_xray_err.append(xr_err)

        if not xray is None:
            xraydat = [[n, d / np.mean(data_xray) - 1] for n, d in enumerate(data_xray)]
            xraydat_err = [d / np.mean(data_xray) for d in data_xray_err]

        radiodat = [[n, d / np.mean(data_radio_cb) - 1] for n, d in enumerate(data_radio_cb)]
        radiodat_err = [d / np.mean(data_radio_cb) for d in data_radio_err_cb]

        lgnd = ['UV-subtract']

        plt.errorbar([i[0] for i in radiodat], [i[1] for i in radiodat],
                     yerr=radiodat_err, fmt='.', ecolor='darkblue', elinewidth=0.4,
                     color='darkblue', capsize=2, markersize=4)

        if rudnick:
            radiodat2 = [[n, d / np.mean(data_radio_rudnick) - 1] for n, d in enumerate(data_radio_rudnick)]
            radiodat_err2 = [d / np.mean(data_radio_rudnick) for d in data_radio_err_rudnick]
            lgnd.append('Rudnick')
            plt.errorbar([i[0] for i in radiodat], [i[1] for i in radiodat2],
                         yerr=radiodat_err2, fmt='.', ecolor='darkgreen', elinewidth=0.4,
                         color='darkgreen', capsize=2, markersize=4)

        if not xray is None:
            plt.errorbar([i[0] for i in xraydat], [i[1] for i in xraydat],
                        yerr=xraydat_err, fmt='.', ecolor='darkred', elinewidth=0.4,
                        color='darkred', capsize=2, markersize=4)
            lgnd.append('X-ray')

        if len(lgnd)>1:
            plt.legend(lgnd)
        plt.xticks(list(range(0,len(radiodat))), list(range(1,len(radiodat)+1)))
        plt.xlabel('Slice')
        if not xray is None:
            plt.ylabel('SB/mean(SB)-1')
        else:
            plt.ylabel('$I_{R}$ (Jy/arcsec$^{2}$)')
        plt.grid(False)
        plt.savefig(save, bbox_inches="tight")
        plt.close()

        # hdu = fits.open(f'../fits/Bridge.fits')
        # image_data = hdu[0].data
        # while len(image_data.shape) != 2:
        #     image_data = image_data[0]
        # wcs = WCS(hdu[0].header, naxis=2)



    def ptp(self, pixelsize=None, dYpix=300, dXpix=210, start=(1495,1275), fitsfile=None, y=None, xray=None, savenumpy=None, savefig=None, maskregion=None, grid=None):

        if pixelsize is None:
            pixelsize = 4*int(np.sqrt(self.beamarea/np.pi))

        stepsize = int(pixelsize*np.sin(40))

        if not grid:
            print(f'Cell size is {self.pix_to_size(0.072)*pixelsize} X {self.pix_to_size(0.072)*pixelsize}')
            print(f'Cell size is {pixelsize*abs((self.header["CDELT2"] * u.deg).to(u.arcsec))} X {pixelsize*abs((self.header["CDELT2"] * u.deg).to(u.arcsec))}')

        if self.resolution >= 60:
            start = (start[0]//2, start[1]//2)
            dYpix//=2
            dXpix//=2

        if self.resolution < 20:
            start = (int(start[0]*2), int(start[1]*2))
            dYpix *= 2
            dXpix *= 2

        if fitsfile:
            self.reproject_map(fitsfile, 'test.fits')
            hdu = fits.open('test.fits')
            wcs = WCS(hdu[0].header, naxis=2)
            header = wcs.to_header()
            if xray:
                image_data, xray_back, xray_exp = self.get_xray(fitsfile, reproject=True, plot3d=True)
                xray_exp[xray_exp == 0.] = 1.
                xray_back[np.isnan(xray_back)] = 0.
                # image_data = self.convolve_image(image_data, 3)
            else:
                image_data = hdu[0].data
        else:
            image_data = self.image_data
            hdu = self.hdu
            header = self.wcs.to_header()

        image_data[np.isnan(image_data)] = 0

        if maskregion:
            r = pyregion.open(maskregion)
            manualmask = r.get_mask(hdu=hdu[0], shape=image_data.shape)
            image_data[manualmask] = 0

        def next_line(d):
            original = d.replace('box(', '').split(',')[0:2]
            return f'box({int(original[0]) - stepsize},{int(original[1]) + stepsize},{pixelsize},{pixelsize},45)'

        if grid:

            if y:
                region = open(f'../regions/{grid.lower()}5arcmin.reg').read()
                region = region.split('\n')
                region_head = region[0:2]
                structures = region[3:]

            else:

                gridsize=11

                region = open('../regions/'+grid+f'6arc.reg', 'r').read()
                print(f'Cell size is {self.pix_to_size(0.072)*gridsize} X {self.pix_to_size(0.072)*gridsize}')
                print(f'Cell size is {gridsize*abs((self.header["CDELT2"] * u.deg).to(u.arcsec))} X {gridsize*abs((self.header["CDELT2"] * u.deg).to(u.arcsec))}')
                region = region.split('\n')
                region_head = region[0:2]
                structures = region[3:]

                f = fits.open(f'../ptp_dir/'+grid+f'/grid_{gridsize}x{gridsize}_results.fits')
                t = f[1].data
                t = t[(t['xray_sb'] > 0) & (t['radio1_sb'] > 0) & (t['xray_sb_err'] > 0) & (t['radio1_sb_err'] > 0)]

                indices = [int(i[0].replace('xaf_','')) for i in t]

            data = []
            data_error = []
            for n, structure in enumerate(structures):
                if (not 'box' in structure):
                    continue
                if not y:
                     if n not in indices:
                        continue
                r = pyregion.parse('\n'.join(region_head + [structure])).as_imagecoord(header=hdu[0].header)
                mask = r.get_mask(hdu=hdu[0], shape=image_data.shape)
                im_mask = image_data[
                    mask]  # /(np.sum(mask) * (abs((header['CDELT1'] * u.deg).to(u.arcsec)).value**2))

                if xray:
                    im_mask[np.isnan(im_mask)] = 0
                    xray_exp_mask = xray_exp[mask]
                    xray_back_mask = xray_back[mask]
                    xray_exp_mask[np.isnan(xray_exp_mask)] = 1
                    xray_back_mask[np.isnan(xray_back_mask)] = 0
                    xr = (np.sum(im_mask) - np.sum(xray_back_mask)) / np.mean(xray_exp_mask)
                    xr_err = 4 * np.sqrt(np.sum(im_mask) + np.sum(xray_back_mask)) / np.mean(xray_exp_mask)
                    data.append(xr)
                    data_error.append(xr_err)
                else:
                    im_mask = im_mask[im_mask != 0]
                    data.append(np.mean(im_mask))
                    data_error.append(np.std(im_mask))
            arr = np.array(data)
            arr_err = np.array(data_error)
            if savenumpy:
                np.save(savenumpy, arr)
                np.save(savenumpy.replace('.npy', '_err.npy'), arr_err)

        else:
            region = \
            """
            # Region file format: DS9 version 4.1
            global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
            fk5"""
            for n in range(0, dXpix, stepsize):
                region+= f"\nbox({start[0]+n},{start[1]+n},{pixelsize},{pixelsize},45)"

            region = region.replace('        ','').split('\n')[1:]
            regionoutput = '\n'.join(region)

            region_head = region[0:2]
            structures = region[3:]

            datatotal= [[] for _ in range(0, dYpix, stepsize)]
            for structure in structures:
                if 'box' in structure:
                    for N in range(len(datatotal)):
                        if N>0:
                            structure = next_line(structure)
                            # im_mask_old = im_mask
                        r = pyregion.parse('\n'.join(region_head + [structure])).as_imagecoord(header=hdu[0].header)
                        mask = r.get_mask(hdu=hdu[0], shape=image_data.shape)
                        im_mask = image_data[mask]#/(np.sum(mask) * (abs((header['CDELT1'] * u.deg).to(u.arcsec)).value**2))

                        if xray:
                            im_mask[np.isnan(im_mask)] = 0
                            xray_exp_mask = xray_exp[mask]
                            xray_back_mask = xray_back[mask]
                            xray_exp_mask[np.isnan(xray_exp_mask)] = 1
                            xray_back_mask[np.isnan(xray_back_mask)] = 0
                            xr = (np.sum(im_mask)-np.sum(xray_back_mask)) / np.mean(xray_exp_mask)
                            xr_err = 4*np.sqrt(np.sum(im_mask)+np.sum(xray_back_mask)) / np.mean(xray_exp_mask)
                            datatotal[N].append(
                                (structure.replace('box(', '').split(',')[0:2], xr, xr_err))
                        # if np.sum(im_mask)==0:
                        #     im_mask = im_mask_old
                        else:
                            im_mask = im_mask[im_mask != 0]
                            datatotal[N].append((structure.replace('box(', '').split(',')[0:2], np.mean(im_mask), np.std(im_mask)))
                        regionoutput+='\n'+structure

            arr = np.array([[d[1] for d in data] for data in datatotal])
            arr_err = np.array([[d[2] for d in data] for data in datatotal])
            arr[np.isnan(arr)]=0
            arr[arr==0]=np.mean(arr)
            arr_err[arr_err==0]=np.mean(arr_err)
            arr_err[np.isnan(arr)]=np.mean(arr_err)
            if savenumpy:
                np.save(savenumpy, arr)
                np.save(savenumpy.replace('.npy', '_err.npy'), arr_err)
            # ticks = np.array([[d[0] for d in data] for data in datatotal])
            X, Y = np.meshgrid(range(arr.shape[1]), range(arr.shape[0]))
            hf = plt.figure()
            ha = hf.add_subplot(111, projection='3d')
            ha.plot_surface(X, Y, arr/np.mean(arr[arr!=0]), cmap ='viridis', edgecolor ='green')
            plt.xticks([])
            plt.yticks([])
            ha.text(len(X)//8, len(Y), 1, 'A401', color='red', fontsize=14)
            ha.text(len(X)//4, 0, 1, 'A399', color='red', fontsize=14)

            with open('corr_area.reg', 'w') as f:
                f.write('\n'.join(regionoutput.split('\n')[3:]))

            # print(ticks[0, :, 1])
            # plt.xticks(list(range(arr.shape[1])), [self.wcs.wcs_pix2world(0, int(i), 0)[1] for i in ticks[0, :, 1]], rotation=70)
            # plt.yticks(list(range(arr.shape[0])), [self.wcs.wcs_pix2world(int(i), 0, 0)[0] for i in ticks[:, 0, 0]], rotation=70)

            # plt.ylabel('DEC')
            # plt.xlabel('RA')
            plt.grid(False)
            if savefig:
                plt.savefig(savefig, bbox_inches="tight")
                plt.close()
            else:
                plt.show()

            # for data in datatotal:
            #     y = [i[1] for i in data]
            #     x = [i[0][0] for i in data]
            #     plt.plot(y)
            #     # plt.xticks(list(range(len(y))), x, rotation=90)
            #     # plt.ylabel('surface brightness [Jy/beam]')
            #     # plt.show()
            # plt.legend(['line '+str(i) for i in list(range(Nlines))])
            # plt.show()

    def make_cutout(self, image_data=None, pos: tuple = None, size: tuple = (1000, 1000), savefits=None):
        """
        Make cutout from your image with this method.
        pos (tuple) -> position in pixels
        size (tuple) -> size of your image in pixel size, default=(1000,1000)
        """
        if image_data is None:
            image_data = self.image_data

        out = Cutout2D(
            data=image_data,
            position=pos,
            size=size,
            wcs=self.wcs,
            mode='partial'
        )

        self.wcs = out.wcs
        header = out.wcs.to_header()
        for key in list(header.keys()):
            self.header[key] = header[key]
        self.header["BUNIT"] = "JY/BEAM"
        if 'BMAJ' not in list(header.keys()):
            if self.resolution==60:
                self.header['BMAJ'] = 0.0210939002735131
                self.header['BMIN'] = 0.0202202375560568
            else:
                self.header['BMAJ'] = .00292636961915446
                self.header['BMIN'] = .00164009373956384
            # self.header['CRVAL3'] = 143650817.871094
        self.image_data = out.data
        self.rms = self.noise
        self.hdu = [fits.PrimaryHDU(self.image_data, header=self.header)]

        if savefits:
            image_data = out.data
            self.header['CRVAL3'] = 143650817.871094
            image_data = np.expand_dims(np.expand_dims(image_data, axis=0), axis=0)
            fits.writeto(savefits, image_data, self.header, overwrite=True)


        return self

    @property
    def noise(self):
        """
        from Cyril Tasse/kMS
        """
        maskSup = 1e-7
        m = self.image_data[np.abs(self.image_data)>maskSup]
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
        self.rms = rms
        return rms

    def get_noise(self, image_data):
        """
        from Cyril Tasse/kMS
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
        print('Noise : ' + str(round(rms * 1000, 2)) + f'{u.mJy/u.beam}')
        return rms

if __name__ == '__main__':

    #lotss
    # Image = Imaging('../fits/lotss.fits', resolution=6)
    # Image.make_cutout(pos=(int(Image.image_data.shape[0] / 2), int(Image.image_data.shape[0] / 2)), size=(3400, 3400))
    # Image.make_image(convolve=True)
    # Image.make_multip_cuts(regionfile='../regions/lotsscomp.reg', save='lotsscut.png', beamsize=True, vmin=2*7.8e-05, vmax=20*7.8e-05)
    #, vmin=7.8e-05, vmax=20*7.8e-05
    # Image.make_image(vmin=0.00005, show_regions='../regions.reg', save='lotsscut.png', subim=True, colorbar=True)

    #6"
    # Image = Imaging('../fits/6all.fits', resolution=6)
    # Image.make_cutout(pos=(int(Image.image_data.shape[0] / 2), int(Image.image_data.shape[0] / 2)), size=(3000, 3000))
    # Image.make_polygon(points=[[44.497, 13.03]], do_science=True, make_image=True, spectralindex=1.75, save='A399.png',
    #                    savefits='../fits/A399.fits', sigma=5)
    # Image.do_science(region='../regions/A399.reg')
    # Image.make_image(convolve=5)

    # Image.compare_lotss(regionfile='../regions/lotsscomp.reg', save='comparebar.png', beamsize=True, vmin=7.8e-05, vmax=25*7.8e-05, colorbar=True)
    # Image.compare_lotss(regionfile='../regions/lotsscomp.reg', save='compare.png', beamsize=True, vmin=7.8e-05, vmax=25 * 7.8e-05)

    # Image.make_image()
    # Image.ptp(pixelsize=70, savenumpy='radio3d_6.npy', savefig='radio3d_6.png')
    # Image.do_science(region='../regions/bridge.reg')
    # Image.make_image(show_regions='../boxlayout.reg', vmin=0.00005, save='layout.png', colorbar=False, beam=False, give_scale=True)
    # Image.make_image(show_regions='../tessupdate.reg', vmin=0.00005, save='tess.png', colorbar=False, beam=False, give_scale=False)
    # Image.make_cutout(pos=(int(Image.image_data.shape[0] / 2), int(Image.image_data.shape[0] / 2)), size=(3400, 3400))
    # Image.make_subimages(regionfile='../regions.reg', save='6cutouts.png', beamsize=True)
    # Image.make_image(show_regions='../regions.reg', save='6image.png', subim=True, convolve=2, cmap='cubehelix_r', ticks=[0., 0.0001, 0.0005, 0.007, 0.001, 0.002], show_clustername=True, bigim=True)
    # Image.make_subcontour('../regions.reg', save='6subimages.png', fits_lowres='../fits/60all.fits', beamsize=True, cmap='cubehelix')
    # Image = Imaging('../fits/6all.fits', resolution=6)
    # Image.make_cutout(pos=(int(Image.image_data.shape[0] / 2), int(Image.image_data.shape[0] / 2)),
    #                   size=(int(Image.image_data.shape[0] / 2), int(Image.image_data.shape[0] / 2)))
    # Image.make_image(vmin=0.00005, text=True, save='a399a401.png')
    # Image.make_image()


    #20"
    # Image = Imaging('../fits/20all.fits', resolution=20)
    # Image.make_cutout(pos=(int(Image.image_data.shape[0] / 2), int(Image.image_data.shape[0] / 2)), size=(1700, 1700))
    # Image.make_image(save='20image.png', cmap='cubehelix_r', show_regions='../regions/20regions.reg')
    # Image.rudnickfilter(300, open=True)
    # Image.make_image()
    # Image.make_image(show_grid=True, save='grids.png')
    # Image.make_image(convolve=True, save='test20.png', text=True)
    # Image.ptp(pixelsize=35, savenumpy='radio3d_20.npy', savefig='radio3d_20.png')
    # Image.make_cutout(pos=(int(Image.image_data.shape[0] / 2), int(Image.image_data.shape[0] / 2)),
    #                   size=(int(Image.image_data.shape[0] / 2), int(Image.image_data.shape[0] / 2)))
    # Image.make_contourplot(regions='../regions.reg')
    # Image.make_subcontour('../regions.reg', save='20subimages.png', fits_lowres='../fits/80all.fits', beamsize=False)

    # Image.make_image(save='20image.png', vmin=0.0001)
    # Image.rudnickfilter(70, open=True, write='../fits/20rudnick.fits')
    # Image.medianfilter(kpc_scale=70, write='../fits/20median.fits')
    # Image.make_image()

    #BRIDGE ANALYSIS
    # Image = Imaging('../fits/20all.fits', resolution=20)
    # Image.rudnickfilter(kpc_scale=60, open=True, write='../fits/20rudnick.fits')
    # Image.make_image()



    # Image = Imaging('../fits/60cleanbridge_300kpc.fits', resolution=60)
    # Image.make_cutout(pos=(int(Image.image_data.shape[0] / 2), int(Image.image_data.shape[0] / 2)), size=(750, 750))
    # Image.make_image(vmin=0., show_regions='../regions/slicegrid.reg', save='cleanbridgesliced.png', cmap='Blues', ticks=[2e-3, 4e-3, 1e-2])
    # Image.ps(xray='../fits/mosaic_a399_a401.fits', region='../regions/slices_inverse.reg', save='slicesbridge_inv.png')
    # Image.ps(region='../regions/slices.reg', save='slicesbridge.png', rudnick=True)

    # Image.ps(xray='../fits/mosaic_a399_a401.fits', region='../regions/a399extensionslices.reg', save='slicesa399.png')


    # Image.ps(region='../regions/circles_A399.reg', save='A399annulus.png', xray='../fits/mosaic_a399_a401.fits')
    # Image.ps(region='../regions/circles_A401.reg', save='A401annulus.png', xray='../fits/mosaic_a399_a401.fits')

    # Image.make_bridge_overlay_yxr_contourplot(fits1='../fits/mosaic_a399_a401.fits',
    #                                           save='radioxray.png', show_regions='../regions/slices.reg')

    # Image.make_polygon(points=[[44.75, 13.56]], subtract_points=[[44.82, 13.49]],
    #                    do_science=True, make_image=True, spectralindex=1.63, sigma2=8.53, save='A401.png', regionmask='../regions/maska401.reg',
    #                    savefits='../fits/A401.fits', sigma=3) #A401
    # # #
    # Image = Imaging('../fits/60cleanbridge_300kpc.fits', resolution=60)
    # Image.make_cutout(pos=(int(Image.image_data.shape[0] / 2), int(Image.image_data.shape[0] / 2)), size=(750, 750))
    # Image.make_polygon(points=[[44.497, 13.03]], sigma=3, do_science=True, make_image=True, spectralindex=1.75, save='A399.png', savefits='../fits/A399.fits')  # A399
    # # #
    # Image = Imaging('../fits/60cleanbridge_300kpc.fits', resolution=60)
    # Image.make_cutout(pos=(int(Image.image_data.shape[0] / 2), int(Image.image_data.shape[0] / 2)), size=(750, 750))
    # Image.make_polygon(points=[[44.587, 13.27598]],
    #                    subtract_points=[[44.5377, 13.2926], [44.497, 13.013], [44.82, 13.49],
    #                                     [44.7857, 13.1997], [44.597, 13.507], [44.587, 12.856], [44.441, 13.4262]],
    #                    sigma=2, sigma2=5, do_science=True, save='Bridge.png', make_image=True, make_cutout=True,
    #                    regionmask='../regions/maskbridge.reg', savefits='../fits/Bridge.fits', size=(400, 400)) #bridge

    # Image.do_science(region='../regions/bridge.reg')
    # Image = Imaging('../fits/60rudnick.fits', resolution=60)
    # Image.make_cutout(pos=(int(Image.image_data.shape[0] / 2), int(Image.image_data.shape[0] / 2)), size=(1500, 1500))
    # Image.medianfilter(kpc_scale=200)
    # Image.rudnickfilter(kpc_scale=100, open=True, write='../fits/20rudnick.fits')
    # Image.make_image(show_grid=True, save='60rudnick_grid.png', ticks=[1e-3, 2e-2, 1e-1])
    # Image.analyse_corr(A1758=True)
    # Image.make_image(save='justbridge.png', text=True)

    Image = Imaging('../fits/60rudnick.fits', resolution=60)
    # Image.make_cutout(pos=(int(Image.image_data.shape[0] / 2), int(Image.image_data.shape[0] / 2)), size=(1500, 1500))
    # Image.make_image(show_regions='../regions/slices_inverse.reg', save='rudnicksliced.png', cmap='Blues', ticks=[2e-3, 4e-3, 1e-2])
    # Image.ps(xray='../fits/mosaic_a399_a401.fits', region='../regions/slices_inverse.reg', save='slicesbridge_inv_rudnick.png')
    # Image.ps(xray='../fits/mosaic_a399_a401.fits', region='../regions/slices.reg', save='slicesbridge_rudnick.png')

    # Image.ps(region='../regions/circles_A399.reg', save='A399annulus.png', xray='../fits/mosaic_a399_a401.fits')
    # Image.ps(region='../regions/circles_A401.reg', save='A401annulus.png', xray='../fits/mosaic_a399_a401.fits')
    # Image.do_science(results='A401_results_rudnick.fits', spectralindex=1.63)
    # Image.do_science(results='A399_results_rudnick.fits', spectralindex=1.75)
    # Image.do_science(region='../regions/bridge.reg', results='bridge_results_rudnick.fits')
    Image.make_cutout(pos=(int(Image.image_data.shape[0] / 2), int(Image.image_data.shape[0] / 2)), size=(1500, 1500))
    # Image.make_image(show_grid=True, save=f'../ptp_results/60rudnick_grid.png', vmin=1e-3, ticks=[1e-3, 2e-2, 1e-1])

    for i in range(11,12):
        try:
            Image.make_image(show_grid=True, save=f'../analysis/60rudnick_grid.png', ticks=[2e-3, 1e-2], sub=f'_11', cmap='Blues')
        except:
            pass

    #BRIDGE
    # Image.ptp(savenumpy='bridgeradio.npy', grid='bridge')
    # Image.ptp(savenumpy='bridgeradio2.npy', grid='bridge', y=True)
    # Image.ptp(savenumpy='bridgey.npy', fitsfile='../fits/a401_curdecmaps_0.2_1.5s_sz.fits', grid='bridge', y=True)
    # Image.ptp(savenumpy='bridgexray.npy', fitsfile='../fits/mosaic_a399_a401.fits', xray=True, grid='bridge')
    # Image.analyse_corr(grid='bridge', savefig='bridgecorr.png')

    # Image = Imaging('../fits/60cleanbridge_300kpc.fits', resolution=60)
    # Image.ps(xray='../fits/mosaic_a399_a401.fits')

    #HALO A399
    # Image.ptp(savenumpy='a399radio.npy', grid='A399')
    # Image.ptp(savenumpy='a399radio2.npy', grid='A399', y=True)
    # Image.ptp(savenumpy='a399y.npy', fitsfile='../fits/a401_curdecmaps_0.2_1.5s_sz.fits', grid='A399', y=True)
    # Image.ptp(savenumpy='a399xray.npy', fitsfile='../fits/mosaic_a399_a401.fits', xray=True, grid='A399')
    # Image.analyse_corr(grid='A399', savefig='A399corr.png')

    #HALO A401
    # Image.ptp(savenumpy='a401radio.npy', grid='A401')
    # Image.ptp(savenumpy='a401radio2.npy', grid='A401', y=True)
    # Image.ptp(savenumpy='a401y.npy', fitsfile='../fits/a401_curdecmaps_0.2_1.5s_sz.fits', grid='A401', y=True)
    # Image.ptp(savenumpy='a401xray.npy', fitsfile='../fits/mosaic_a399_a401.fits', xray=True, grid='A401')
    # Image.analyse_corr(grid='A401', savefig='A401corr.png')

    # Image.make_cutout(pos=(650, 960), size=(300, 300), savefits='../fits/A401.fits')
    # Image.make_cutout(pos=(800, 600), size=(300, 300), savefits='../fits/A399.fits')

    # Image.make_image()
    # Image.make_bridge_overlay_yxr_contourplot(fits1='../fits/mosaic_a399_a401.fits',
    #                                           save='60cleanbridge.png', show_regions='../regions/lines_2.reg')

    # Image.do_science(region='../regions/bridge.reg')
    # Image.make_cutout(pos=(int(Image.image_data.shape[0]/2), int(Image.image_data.shape[0]/2)), size=(1500, 1500))
    # Image.make_bridge_overlay_contourplot('../fits/a401_curdecmaps_0.2_1.5s_sz.fits', title=' ',
    #                                       minlevel_1=0, maxlevel_1=0.005, steps_1=100, steps_2=np.sqrt(2), minlevel_2=(10**(-5))/2, convolve_2=30, save='ymap.png')
    # Image.make_bridge_overlay_contourplot('../fits/mosaic_a399_a401.fits', title=' ',
    #                                       minlevel_1=0, maxlevel_1=0.005, steps_1=100, convolve_2=2, steps_2=2, maxlevel_2=50, xray=True)
    # Image.make_bridge_overlay_yxr_contourplot(fits2='../fits/a401_curdecmaps_0.2_1.5s_sz.fits', fits1='../fits/mosaic_a399_a401.fits',
    #                                           show_regions='../ptp_dir/grid_35x35/grid_35x35_ds9_image.reg', save='ymapxray.png')

    #60"
    # Image = Imaging('../fits/60all.fits', resolution=60)
    # Image.make_cutout(pos=(int(Image.image_data.shape[0] / 2), int(Image.image_data.shape[0] / 2)), size=(850, 850))
    # Image.rudnickfilter(100, open=True)
    # Image.ptp(savenumpy='bridgegridradio.npy', grid='bridge')
    # Image.ptp(savenumpy='bridgegridy.npy', fitsfile='../fits/a401_curdecmaps_0.2_1.5s_sz.fits', grid='bridge')
    # Image.ptp(savenumpy='bridgegridxray.npy', fitsfile='../fits/mosaic_a399_a401.fits', xray=True, grid='bridge')
    # Image.analyse_corr(grid='bridgegrid', savefig='bridgecorr.png')
    # Image.ptp(savenumpy='radio3d.npy', savefig='radio3d.png', maskregion='../regions/excluderegions60.reg')
    # Image.ptp(savenumpy='y.npy', savefig='y3d.png', fitsfile='../fits/a401_curdecmaps_0.2_1.5s_sz.fits', maskregion='../regions/excluderegions60.reg')
    # Image.ptp(savenumpy='xray.npy', savefig='xray3d.png', fitsfile='../fits/mosaic_a399_a401.fits', xray=True, maskregion='../regions/excluderegions60.reg')
    # Image.make_cutout(pos=(int(Image.image_data.shape[0]/2), int(Image.image_data.shape[0]/2)), size=(850, 850))
    # Image.make_image(show_regions='../regions/60regions.reg', save='60image.png', cmap='cubehelix_r')
    # Image.make_contourplot(regions='../regions.reg')
    # Image.medianfilter(kernelsize=51, write='../fits/60median.fits')
    # Image.make_image()
