from past.utils import old_div
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
from astropy.nddata import Cutout2D
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm, LogNorm
from reproject import reproject_interp
import string
import sys
from astropy.modeling.models import Gaussian2D
from astropy.convolution import convolve, Gaussian2DKernel
from matplotlib.ticker import LogLocator, LogFormatterSciNotation as LogFormatter
import os
from scipy.ndimage import gaussian_filter
import pyregion
from pyregion.mpl_helper import properties_func_default
from astropy.visualization.wcsaxes import WCSAxes
from matplotlib.patches import ConnectionPatch
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from radioflux import Radiomap
from radio_beam import Beams
import warnings

warnings.filterwarnings('ignore')
plt.style.use('ggplot')


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

    def make_image(self, image_data=None, cmap: str = 'CMRmap', vmin=None, vmax=None, show_regions=None, wcs=None,
                   colorbar=True, save=None, text=None, subim=None, beam=True, give_scale=True):
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

        print(vmin, vmax)

        if show_regions is not None:

            def fixed_color(shape, saved_attrs):
                attr_list, attr_dict = saved_attrs
                attr_dict["color"] = "white"
                kwargs = properties_func_default(shape, (attr_list, attr_dict))

                return kwargs

            r = pyregion.open(show_regions).as_imagecoord(header=self.hdu[0].header)
            patch_list, artist_list = r.get_mpl_patches_texts(fixed_color)
            fig = plt.figure(figsize=(7, 10), dpi=200)

            ax = WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=wcs)

            fig.add_axes(ax)
            for patch in patch_list:
                ax.add_patch(patch)
            # for artist in artist_list:
            #     ax.add_artist(artist)
            if subim:
                plt.text(1500+120, 2880-500, 'C', color='white', fontsize=10)
                plt.text(2016.59+105, 2821.64-130, 'B', color='white', fontsize=10)
                plt.text(1539.66+78, 2049.08-142, 'D', color='white', fontsize=10)
                plt.text(2309.1+175, 1808.42-150, 'E', color='white', fontsize=10)
                plt.text(2100.69+95, 1138.61-330, 'F', color='white', fontsize=10)
                plt.text(2114.85+100, 3500-145, 'A', color='white', fontsize=10)

        else:
            figure, ax = plt.subplots(figsize=(7, 10), dpi=200)
            plt.subplot(projection=wcs)
        im = plt.imshow(image_data, origin='lower', cmap=cmap)
        im.set_norm(SymLogNorm(linthresh=vmin*10, vmin=vmin, vmax=vmax, base=10))
        plt.xlabel('Galactic Longitude', size=15)
        plt.ylabel('Galactic Latitude', size=15)
        plt.tick_params(axis='both', which='major', labelsize=12)
        if colorbar:
            cbar = plt.colorbar(orientation='horizontal', shrink=1)
            # cbar.ax.set_xscale('log')
            cbar.locator = LogLocator()
            cbar.formatter = LogFormatter()
            cbar.update_normal(im)
            cbar.set_label('Surface brightness [Jy/beam]', size=15)
        if text:
            plt.text(3350 - (6000-image_data.shape[0])/2, 1950 - (6000-image_data.shape[1])/2, 'A399', color='pink', fontsize=14)
            plt.text(2200 - (6000-image_data.shape[0])/2, 4000 - (6000-image_data.shape[1])/2, 'A401', color='pink', fontsize=14)

        if give_scale:
            Mpc_pixels = 1 / (abs((self.header['CDELT2'] * u.deg).to(u.rad).value) * self.cosmo.angular_diameter_distance(0.072)).value
            l1, l2 = [int(image_data.shape[0]*0.1), int(image_data.shape[0]*0.1+Mpc_pixels)], [int(image_data.shape[1]*0.9), int(image_data.shape[1]*0.9)]
            plt.plot(l1, l2, color='cyan', linewidth=1)
            plt.text((l1[0]+l1[1])/2, 1.02*(l2[0]+l2[1])/2, 'Mpc', color='cyan', fontsize=12, horizontalalignment='center')


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
        if save:
            plt.savefig(save, dpi=250, bbox_inches="tight")
            plt.close()

        else:
            plt.show()

        return self

    def make_subimages(self, regionfile, cmap='CMRmap', save=None):

        r = pyregion.open(regionfile).as_imagecoord(header=self.hdu[0].header)

        fig = plt.figure(figsize=(10, 10))

        rows, cols = len(r)//2, 2

        for k, shape in enumerate(r):
            s = np.array(shape.coord_list)

            out = Cutout2D(
                data=self.image_data,
                position=(s[0], s[1]),
                size=(s[2], s[3]),
                wcs=self.wcs,
                mode='partial'
            )
            norm = SymLogNorm(linthresh=self.rms * 10, vmin=self.rms, vmax=self.rms*25, base=10)

            plt.subplot(rows, cols, k+1, projection=self.wcs)
            im = plt.imshow(out.data, origin='lower', cmap=cmap, norm=norm)
            plt.xlabel('Right Ascension (J2000)')
            plt.ylabel('Declination (J2000)')
            plt.tick_params(axis='both', which='major', labelsize=12)

            if type(self.resolution) == int:
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

    def make_subcontour(self, regionfile, save=None, fits_lowres=None, beamsize=None):

        r = pyregion.open(regionfile).as_imagecoord(header=self.hdu[0].header)

        maxlevel=200*self.rms_full
        minlevel=self.rms_full*5

        if fits_lowres:
            # self.reproject_map(fits_lowres, 'temp.fits')
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
                size=(s[2], s[3]),
                wcs=self.wcs,
                mode='partial'
            )

            image_data = np.clip(out.data, a_min=levels2[0], a_max=maxlevel * 0.99)
            plt.subplot(rows, cols, k+1, projection=self.wcs)
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

        fig.subplots_adjust(top=0.8)
        cbar_ax = fig.add_axes([0.22, 0.88, 0.6, 0.03]) # l, b, w, h
        cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal', ticks=[round(a, 5) for a in
                                   np.logspace(np.log(minlevel), np.log(maxlevel), 4, endpoint=True)])

        cbar.set_label('Surface brightness [Jy/beam]')
        cbar.ax.set_xscale('log')

        plt.grid(False)
        if save:
            plt.savefig(save, dpi=250, bbox_inches="tight")
            plt.close()

        else:
            plt.show()
    def pix_to_size(self, z):
        return abs((self.header['CDELT2'] * u.deg).to(u.rad).value) * self.cosmo.angular_diameter_distance(z)



    def do_science(self, region, objects='bridge'):

        rm = Radiomap(fits.open(self.fitsfile))

        if objects=='bridge':
            # mask1 = np.where(self.image_data<self.rms, True, False)
            r = pyregion.open(region).as_imagecoord(header=self.hdu[0].header)
            extrapolatedregion = pyregion.open('../extrapolatedregion.reg').as_imagecoord(header=self.hdu[0].header).\
                get_mask(hdu=self.hdu[0], shape=self.image_data.shape)
            pix_extrap = np.sum(extrapolatedregion)
            mask = r.get_mask(hdu=self.hdu[0], shape=self.image_data.shape)
            # mask = np.where(mask2, mask1, False)
            image_data = np.where(mask, self.image_data, 0)
            pixnum=np.sum(np.where(image_data!=0, 1, 0))
            print(pixnum)
            # image_data = image_data[image_data>self.rms*3]
            # self.make_image(image_data, vmin=0.0001, vmax=0.003)
            integrated_surface_brightness = np.nansum(image_data)
            av_sb = round(np.mean(image_data[image_data!=0])*1000,10) * u.mJy/u.beam
            # print(f'Integrated surface brightness is: {round(integrated_surface_brightness*1000,2) * u.mJy/u.beam}')
            print(f'Average surface brightness is: {av_sb}')
            flux_density = round(integrated_surface_brightness*1000/rm.area,2) * u.mJy
            flux_density_err = np.std(image_data)/rm.area
            print(f'Total extrapolated flux density is: {flux_density} $\pm$ {flux_density_err}')
            area = self.pix_to_size(0.072)**2*pixnum

            z=0.072
            L = np.nansum(image_data) * self.pix_to_size(z).value**2 * (1 / rm.area) * 10 ** (-26) * 4 * np.pi * (1 + z) ** (0.7 - 1) * (1 + z) ** 4 / (
                    ((1.5 * u.arcsec.to(u.rad)) * (1 * u.m).to(u.Mpc).value) ** 2)*u.W/u.Hz

            print(f'Radio power: {L}')

            # area_govoni=3.9*u.Mpc**2
            # print(f'Calculated area {area} and Govoni area {area_govoni}')
            # print(f'Radio power: {(flux_density*area).to(u.W/u.Hz)}')


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
                         regions=None, scale=None, save=None):

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
        minlevel=self.rms_full*3

        image_data = np.clip(image_data, a_min=np.min(image_data), a_max=maxlevel*0.99)


        levels = [minlevel]
        for _ in range(10):
            levels.append(levels[-1]*2)

        levels2 = np.linspace(minlevel/3, maxlevel, 1000)

        norm = LogNorm(vmin=minlevel, vmax=maxlevel)

        plt.contour(image_data, levels, colors=('k'), linestyles=('-'), linewidths=(0.3,))
        plt.contourf(image_data, levels2, cmap='Blues', norm=norm)
        cbar = plt.colorbar(orientation='horizontal', shrink=1)
        cbar.set_label('Surface brightness [Jy/beam]')
        cbar.ax.set_xscale('log')
        # cbar.set_ticks([10e-4, 10e-3, 10e-2, 10e-1])
        if scale:
            Mpc_pixels = 1 / (abs((self.header['CDELT2'] * u.deg).to(u.rad).value) * self.cosmo.angular_diameter_distance(0.072)).value
            l1, l2 = [int(image_data.shape[0]*0.1), int(image_data.shape[0]*0.1+Mpc_pixels)], [int(image_data.shape[1]*0.9), int(image_data.shape[1]*0.9)]
            plt.plot(l1, l2, color='pink', linewidth=1)
            plt.text((l1[0] + l1[1]) / 2, 1.02 * (l2[0] + l2[1]) / 2, 'Mpc', color='cyan', fontsize=12,
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

        plt.xlabel('Right Ascension (J2000)')
        plt.ylabel('Declination (J2000)')
        plt.title(title)

        plt.grid(False)
        if save:
            plt.savefig(save, dpi=250, bbox_inches="tight")
            plt.close()

        else:
            plt.show()
        return self

    def get_xray(self, fitsfile, reproject=True):
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

        return (data - data_bkg)/data_exp

    def make_bridge_overlay_yxr_contourplot(self, fits1, fits2, save=None, beamsize=True, show_regions=None):

        self.reproject_map(fits1, 'fits1.fits')
        self.reproject_map(fits2, 'fits2.fits')


        hdu1 = fits.open('fits1.fits')
        image_data1 = self.get_xray(fits1, reproject=True)
        hdu2 = fits.open('fits2.fits')
        image_data2 = hdu2[0].data

        while len(image_data1.shape) != 2:
            image_data1 = image_data1[0]
        image_data1 = self.convolve_image(image_data1, 2)
        wcs_1 = WCS(hdu1[0].header)
        while len(image_data2.shape) != 2:
            image_data2 = image_data2[0]
        image_data2 = self.convolve_image(image_data2, 30)
        wcs_2 = WCS(hdu2[0].header)


        plt.figure(figsize=(7, 10))
        plt.subplot(projection=self.wcs)

        if show_regions is not None:

            def fixed_color(shape, saved_attrs):
                attr_list, attr_dict = saved_attrs
                attr_dict["color"] = "indianred"
                kwargs = properties_func_default(shape, (attr_list, attr_dict))

                return kwargs

            r = pyregion.open(show_regions).as_imagecoord(header=self.hdu[0].header)
            patch_list, artist_list = r.get_mpl_patches_texts(fixed_color)

            for n, patch in enumerate(patch_list):
                f = fits.open(show_regions.split('_ds9_image.reg')[0]+'_results.fits')
                t = f[1].data
                t = t[(t['xray_sb'] > 0) & (t['radio1_sb'] > 0) & (t['xray_sb_err'] > 0) & (t['radio1_sb_err'] > 0)]
                if n in [int(l.replace('xaf_','')) for l in list(t['region_name'])]:
                    plt.gca().add_patch(patch)

        levels_1 = [self.get_noise(image_data1)*10]
        for _ in range(10):
            levels_1.append(levels_1[-1]*2)

        levels_2 = [(10**(-5))/2]
        for _ in range(10):
            levels_2.append(levels_2[-1]*np.sqrt(2))

        plt.contour(image_data1, levels_1, colors=('red'), linestyles=('-'), linewidths=(2,))
        plt.contour(image_data2, levels_2, colors=('orange'), linestyles=('-'), linewidths=(2,))

        plt.imshow(np.where(self.image_data<self.rms*3, 0, 1), cmap='Greys')

        levels = np.linspace(0, 0.005, 100)
        plt.contourf(np.clip(self.image_data, a_min=levels[0], a_max=np.max(self.image_data)), levels, cmap='Blues')
        cbar = plt.colorbar(orientation='horizontal', shrink=1, ticks=[0.0, 0.002, 0.004])
        cbar.set_label('Surface brightness  [Jy/beam]')
        # cbar.ax.set_xscale('log')
        plt.xlabel('Right Ascension (J2000)')
        plt.ylabel('Declination (J2000)')

        Mpc_pixels = 1 / (abs((self.header['CDELT2'] * u.deg).to(u.rad).value) * self.cosmo.angular_diameter_distance(0.072)).value
        l1, l2 = [int(self.image_data.shape[0]*0.05), int(self.image_data.shape[0]*0.05+Mpc_pixels)], [int(self.image_data.shape[1]*0.9), int(self.image_data.shape[1]*0.9)]
        plt.plot(l1, l2, color='brown', linewidth=1)
        plt.text((l1[0]+l1[1])/2, 1.02*(l2[0]+l2[1])/2, 'Mpc', color='brown', fontsize=12, horizontalalignment='center')
        if type(self.resolution) == int and beamsize:
            beampix = self.resolution / (self.header['CDELT2'] * u.deg).to(u.arcsec).value/2
            print(beampix)
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
        cbar.set_label('Surface brightness  [Jy/beam]')
        # cbar.ax.set_xscale('log')
        plt.xlabel('Right Ascension (J2000)')
        plt.ylabel('Declination (J2000)')

        Mpc_pixels = 1 / (abs((self.header['CDELT2'] * u.deg).to(u.rad).value) * self.cosmo.angular_diameter_distance(0.072)).value
        l1, l2 = [int(self.image_data.shape[0]*0.05), int(self.image_data.shape[0]*0.05+Mpc_pixels)], [int(self.image_data.shape[1]*0.9), int(self.image_data.shape[1]*0.9)]
        plt.plot(l1, l2, color='brown', linewidth=1)
        plt.text((l1[0]+l1[1])/2, 1.02*(l2[0]+l2[1])/2, 'Mpc', color='brown', fontsize=12, horizontalalignment='center')
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

        self.wcs = out.wcs
        self.header = self.wcs.to_header()
        self.image_data = out.data
        self.rms = self.noise
        self.hdu = [fits.PrimaryHDU(self.image_data, header=self.header)]


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
        print('Noise : ' + str(round(rms * 1000, 2)) + f'{u.mJy/u.beam}')
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

    #6"
    # Image = Imaging('../fits/6all.fits', resolution=6)
    # Image.make_image(show_regions='../boxlayout.reg', vmin=0.00005, save='layout.png', colorbar=False, beam=False, give_scale=True)
    # Image.make_image(show_regions='../tessupdate.reg', vmin=0.00005, save='tess.png', colorbar=False, beam=False, give_scale=False)
    # Image.make_cutout(pos=(int(Image.image_data.shape[0] / 2), int(Image.image_data.shape[0] / 2)),
    #                   size=(int(Image.image_data.shape[0] / 1.5), int(Image.image_data.shape[0] / 1.5)))
    # Image.do_science(region='../bridge.reg', objects='bridge')
    # Image.make_image(vmin=0.00005, show_regions='../regions.reg', save='subimagelayout.png', subim=True, colorbar=True, text=True)
    # Image.make_subcontour('../regions.reg', save='6subimages.png', fits_lowres='../fits/80all.fits', beamsize=False)
    # Image = Imaging('../fits/6all.fits', resolution=6)
    # Image.make_cutout(pos=(int(Image.image_data.shape[0] / 2), int(Image.image_data.shape[0] / 2)),
    #                   size=(int(Image.image_data.shape[0] / 2), int(Image.image_data.shape[0] / 2)))
    # Image.make_image(vmin=0.00005, text=True, save='a399a401.png')

    #10"

    #20"
    # Image = Imaging('../fits/20all.fits', resolution=20)
    # Image.make_subcontour('../regions.reg', save='20subimages.png')

    #60"
    # Image = Imaging('../fits/60all.fits', resolution=60)
    # Image.make_cutout(pos=(int(Image.image_data.shape[0]/2), int(Image.image_data.shape[0]/2)), size=(850, 850))
    # Image.make_bridge_overlay_contourplot('../fits/a401_curdecmaps_0.2_1.5s_sz.fits', title='y-map contour lines and radio filled contour',
    #                                       minlevel_1=0, maxlevel_1=0.005, steps_1=100, steps_2=6, minlevel_2=(10**(-5)/2), convolve_2=True)
    # Image.make_bridge_overlay_contourplot('../fits/mosaic_a399_a401.fits', title='X-ray contour lines and radio filled contour',
    #                                       minlevel_1=0, maxlevel_1=0.005, steps_1=100, steps_2=6, convolve_2=True, maxlevel_2=50, xray=True)

    #80"
    Image = Imaging('../fits/80all.fits', resolution=80)
    Image.make_cutout(pos=(int(Image.image_data.shape[0]/2), int(Image.image_data.shape[0]/2)), size=(850, 850))
    # Image.make_bridge_overlay_contourplot('../fits/a401_curdecmaps_0.2_1.5s_sz.fits', title=' ',
    #                                       minlevel_1=0, maxlevel_1=0.005, steps_1=100, steps_2=np.sqrt(2), minlevel_2=(10**(-5))/2, convolve_2=30, save='ymap.png')
    # Image.make_bridge_overlay_contourplot('../fits/mosaic_a399_a401.fits', title=' ',
    #                                       minlevel_1=0, maxlevel_1=0.005, steps_1=100, convolve_2=2, steps_2=2, maxlevel_2=50, xray=True)
    Image.make_bridge_overlay_yxr_contourplot(fits2='../fits/a401_curdecmaps_0.2_1.5s_sz.fits', fits1='../fits/mosaic_a399_a401.fits', show_regions='../regions/ptp_dir/grid_13x13/grid_13x13_ds9_image.reg', save='ymapxray.png')