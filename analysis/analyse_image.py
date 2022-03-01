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

    @property
    def beamarea(self):
        # Given a fitsfile this calculates the beamarea in pixels

        bmaj = self.hdu[0].header['BMAJ']
        bmin = self.hdu[0].header['BMIN']

        beammaj = bmaj / (2.0 * (2 * np.log(2)) ** 0.5)  # Convert to sigma
        beammin = bmin / (2.0 * (2 * np.log(2)) ** 0.5)  # Convert to sigma
        pixarea = abs(self.hdu[0].header['CDELT1'] * self.hdu[0].header['CDELT2'])

        beamarea = 2 * np.pi * 1.0 * beammaj * beammin  # Note that the volume of a two dimensional gaus$
        beamarea_pix = beamarea / pixarea

        return beamarea_pix

    def remove_compactsources(self, kernelsize=None, write=None):
        print("Apply median kernel")
        self.image_data = sn.median_filter(self.image_data, kernelsize)
        if write:
            self.hdu[0].data = np.expand_dims(np.expand_dims(self.image_data, axis=0), axis=0)
            self.hdu.writeto(write, overwrite=True)
            print('Saved: '+write)
        return self

    def remove_compactsources2(self, kernelsize=None, write=None):
        mins = filters.minimum_filter(self.image_data, size=(kernelsize, kernelsize))
        openmp = filters.maximum_filter(mins, size=(kernelsize, kernelsize))
        self.image_data -= openmp
        if write:
            self.hdu[0].data = np.expand_dims(np.expand_dims(self.image_data, axis=0), axis=0)
            self.hdu.writeto(write, overwrite=True)
            print('Saved: '+write)
        return self

    def make_image(self, image_data=None, cmap: str = 'CMRmap', vmin=None, vmax=None, show_regions=None, wcs=None,
                   colorbar=True, save=None, text=None, subim=None, beam=True, give_scale=True, convolve=None):
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
            image_data = self.convolve_image(image_data, 4)

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
                plt.text(3000/3400*(1500+120)-30, 3000/3400*(2880-500), 'C', color='white', fontsize=10)
                plt.text(3000/3400*(2016.59+105), 3000/3400*(2821.64-130)+20, 'B', color='white', fontsize=10)
                plt.text(3000/3400*(1539.66+78)-60, 3000/3400*(2049.08-142)-50, 'D', color='white', fontsize=10)
                plt.text(3000/3400*(2309.1+175)+20, 3000/3400*(1808.42-150)-70, 'E', color='white', fontsize=10)
                plt.text(3000/3400*(2100.69+95)-30, 3000/3400*(1088.61-330)-180, 'F', color='white', fontsize=10)
                plt.text(3000/3400*(2114.85+100), 3000/3400*(3500-145)+130, 'A', color='white', fontsize=10)

        else:
            figure, ax = plt.subplots(figsize=(7, 10), dpi=200)
            plt.subplot(projection=wcs)
        im = plt.imshow(image_data, origin='lower', cmap=cmap)
        if self.resolution>6:
            im.set_norm(PowerNorm(vmin=vmin, vmax=vmax, gamma=2/3))
        else:
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
            plt.text(3450/4 - (6000/4-image_data.shape[0])/2, 2450/4 - (6000/4-image_data.shape[1])/2, 'A399', color='pink', fontsize=14)
            plt.text(2200/4 - (6000/4-image_data.shape[0])/2, 3750/4 - (6000/4-image_data.shape[1])/2, 'A401', color='pink', fontsize=14)

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

    def make_subimages(self, regionfile, cmap='CMRmap', save=None, beamsize=None):

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

            plt.subplot(rows, cols, k+1, projection=out.wcs)
            im = plt.imshow(out.data, origin='lower', cmap=cmap, norm=norm)
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

    def make_3d_map(self):
        pass

    def pix_to_size(self, z):
        return abs((self.header['CDELT2'] * u.deg).to(u.rad).value) * self.cosmo.angular_diameter_distance(z)

    def do_science(self, region=None):

        z=0.072

        t = fits.open('../regions/ptp_dir/grid_13x13/grid_13x13_results.fits')[1].data
        t = t[(t['xray_sb'] > 0)
              & (t['radio1_sb'] > 0)
              & (t['xray_sb_err'] > 0)
              & (t['radio1_sb_err'] > 0)]

        av_sb = np.mean(t['radio1_fluxdensity'])*u.Jy/u.beam # average surface brigthness
        av_sb_arcsec = np.mean(t['radio1_sb'])*u.Jy/(u.arcsec)**2
        err_sb = np.std(t['radio1_fluxdensity'])*u.Jy/u.beam # 1*sigma of every area element
        err_sb_arcsec = np.std(t['radio1_sb'])*u.Jy/(u.arcsec)**2
        bridge_area = 1.3*3*(u.Mpc**2) # from govoni et al
        N_pixels_bridge = int(bridge_area/(self.pix_to_size(0.072)**2))
        integr_sb = N_pixels_bridge*av_sb_arcsec*(abs((self.header['CDELT2'] * u.deg).to(u.arcsec))**2)
        # integr_sb = 822

        L=(integr_sb* self.pix_to_size(z)**2 * 4 * np.pi * (1 + z) ** (0.7 - 1) * (1 + z) ** 4 / (
                (1.5 * u.arcsec.to(u.rad)) ** 2)).to(u.W/u.Hz)

        print(f'# of pixels: {N_pixels_bridge}')
        print(f'Average surface brightness: {av_sb} $\pm$ {err_sb}')
        print(f'Average surface brightness: {av_sb_arcsec} $\pm$ {err_sb_arcsec}')
        flux_density_err = N_pixels_bridge*err_sb_arcsec*abs((self.header["CDELT2"] * u.deg).to(u.arcsec))**2
        print(f'Total flux density is: {integr_sb} $\pm$ {flux_density_err}')
        radiopower_err = (flux_density_err*self.pix_to_size(z)**2 * 4 * np.pi * (1 + z) ** (0.7 - 1) * (1 + z) ** 4 / ((1.5 * u.arcsec.to(u.rad)) ** 2)).to(u.W/u.Hz)
        print(f'Radio power {L} $\pm$ {radiopower_err}')
        # volume = (1.3/2)**2*np.pi*3*(u.Mpc**3)
        volume = (bridge_area*12.1*u.Mpc)
        print(f'Emissivity is {(L/(volume)).to(u.erg/(u.s*u.cm*u.cm*u.cm*u.Hz))} $\pm$ {(radiopower_err/(volume)).to(u.erg/(u.s*u.cm*u.cm*u.cm*u.Hz))}')

        if 'bridge' in region:
            r = pyregion.open(region).as_imagecoord(header=self.hdu[0].header)
            mask = r.get_mask(hdu=self.hdu[0], shape=self.image_data.shape)
            image_data = np.where(mask, self.image_data, 0)
            image_data = image_data[image_data!=0]
            N_pixels_bridge = len(image_data)
            av_sb = np.median(image_data) * u.Jy/u.beam # average surface brigthness
            av_sb_arcsec = av_sb*u.beam / self.beamarea/(abs(self.hdu[0].header['CDELT1'] * 3600.)*u.arcsec * abs(self.hdu[0].header['CDELT1'] * 3600.)*u.arcsec)
            err_sb = np.std(image_data)*u.Jy/u.beam # 1*sigma of every area element
            err_sb_arcsec = np.std(image_data)*u.Jy/self.beamarea/(abs(self.hdu[0].header['CDELT1'] * 3600.)*u.arcsec * abs(self.hdu[0].header['CDELT1'] * 3600.)*u.arcsec)
            bridge_area = N_pixels_bridge * self.pix_to_size(0.072)**2 # from govoni et al
            integr_sb = N_pixels_bridge*av_sb_arcsec*(abs((self.header['CDELT2'] * u.deg).to(u.arcsec))**2)

            L=(integr_sb* self.pix_to_size(z)**2 * 4 * np.pi * (1 + z) ** (0.7 - 1) * (1 + z) ** 4 / (
                    (1.5 * u.arcsec.to(u.rad)) ** 2)).to(u.W/u.Hz)
            print(f'Bridge area: {bridge_area}')
            print(f'# of pixels: {N_pixels_bridge}')
            print(f'Average surface brightness: {av_sb} $\pm$ {err_sb}')
            print(f'Average surface brightness: {av_sb_arcsec} $\pm$ {err_sb_arcsec}')
            flux_density_err = N_pixels_bridge * err_sb_arcsec * abs((self.header["CDELT2"] * u.deg).to(u.arcsec)) ** 2
            print(f'Total flux density is: {integr_sb} $\pm$ {flux_density_err}')
            radiopower_err = (flux_density_err * self.pix_to_size(z) ** 2 * 4 * np.pi * (1 + z) ** (0.7 - 1) * (
                        1 + z) ** 4 / ((1.5 * u.arcsec.to(u.rad)) ** 2)).to(u.W / u.Hz)
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

    def make_bridge_overlay_yxr_contourplot(self, fits1, fits2, save=None, beamsize=True, show_regions=None):

        self.reproject_map(fits1, 'fits1.fits')
        self.reproject_map(fits2, 'fits2.fits')


        hdu1 = fits.open('fits1.fits')
        image_data1 = self.get_xray(fits1, reproject=True)
        hdu2 = fits.open('fits2.fits')
        image_data2 = hdu2[0].data

        while len(image_data1.shape) != 2:
            image_data1 = image_data1[0]
        image_data1 = self.convolve_image(image_data1, 15)
        # plt.figure(figsize=(7, 10))
        # plt.subplot(projection=self.wcs)
        # plt.imshow(image_data1, norm = LogNorm(vmin=self.get_noise(image_data1), vmax=self.get_noise(image_data1)*10))
        # plt.savefig('xray.png')
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
                attr_dict["color"] = "green"
                kwargs = properties_func_default(shape, (attr_list, attr_dict))

                return kwargs

            r = pyregion.open(show_regions).as_imagecoord(header=self.hdu[0].header)
            patch_list, artist_list = r.get_mpl_patches_texts(fixed_color)

            for n, patch in enumerate(patch_list):
                # f = fits.open(show_regions.split('_ds9_image.reg')[0]+'_results.fits')
                # t = f[1].data
                # t = t[(t['xray_sb'] > 0) & (t['radio1_sb'] > 0) & (t['xray_sb_err'] > 0) & (t['radio1_sb_err'] > 0)]
                # if n in [int(l.replace('xaf_','')) for l in list(t['region_name'])]:
                plt.gca().add_patch(patch)

        levels_1 = [self.get_noise(image_data1)]
        for _ in range(10):
            levels_1.append(levels_1[-1]*2)

        levels_2 = [(10**(-5))/2]
        for _ in range(10):
            levels_2.append(levels_2[-1]*np.sqrt(2))

        plt.contour(image_data1, levels_1, colors=('red'), linestyles=('-'), linewidths=(2,))
        plt.contour(image_data2, levels_2, colors=('orange'), linestyles=('-'), linewidths=(2,))


        levels = np.linspace(self.rms, self.rms*5, 100)
        norm = SymLogNorm(linthresh=self.rms*2, vmin=self.rms, vmax=levels[-1], base=10)
        plt.imshow(np.where(self.image_data<levels[-1], 0, 1), cmap='Greys')
        plt.contourf(self.image_data, levels, cmap='Blues', norm=norm)
        ticks = [0.0002, 0.0004]
        cbar = plt.colorbar(orientation='horizontal', shrink=1, ticks=ticks)
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

    def plot_corr(self, savefig=None, halo=None):

        def objective(x, a, b):
            return a * x + b

        def fit(x, y):
            popt, _ = curve_fit(objective, x, y)
            a, b = popt
            x_line = np.arange(min(x) - 1, max(x) + 1, 0.01)
            y_line = objective(x_line, a, b)
            print('y = %.5f * x + %.5f' % (a, b))
            return x_line, y_line

        def linreg(x, y):
            res = linregress(x, y)
            print(f'Slope is {res.slope} +- {res.stderr}')
            return res.slope, res.stderr

        if halo:
            radio = np.load(halo.lower()+'radio.npy')
            y = np.load(halo.lower()+'y.npy')
            xray = np.load(halo.lower()+'xray.npy')

            radio_err = np.load(halo.lower()+'radio_err.npy')
            y_err = np.load(halo.lower()+'y_err.npy')
            xray_err = np.load(halo.lower()+'xray_err.npy')

        else:
            radio = np.load('radio3d.npy')
            y = np.load('y.npy')
            xray = np.load('xray.npy')

            radio_err = np.load('radio3d_err.npy').flatten()
            y_err = np.load('y_err.npy').flatten()
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


        msk = ((xray>0) & (radio>self.rms*3) & (xray_err/xray<1))

        radio=radio[msk]
        radio_err=radio_err[msk]
        y=y[msk]
        y_err=y_err[msk]
        xray=xray[msk]
        xray_err=xray_err[msk]

        radio_err[np.isnan(radio_err)] = 0
        xray_err[np.isnan(xray_err)] = 0
        y_err[np.isnan(y_err)] = 0

        radio_err/=np.mean(radio)
        radio/=np.mean(radio)
        y_err/=np.mean(y)
        y/=np.mean(y)
        xray_err/=np.mean(xray)
        xray/=np.mean(xray)
        print('radio vs. xray')
        slopex, errx = linreg(np.log10(xray), np.log10(radio))
        fitxray = fit(np.log10(xray), np.log10(radio))
        print('radio vs y')
        slopey, erry =linreg(np.log10(y), np.log10(radio))
        fity = fit(np.log10(y), np.log10(radio))


        print('Pearson R (x-ray vs radio): ' + str(pearsonr(np.log(xray), np.log(radio))))
        print('Pearson R (ymap vs radio): ' + str(pearsonr(np.log(radio), np.log(y))))

        print('Spearman R (x-ray vs radio): ' + str(spearmanr(np.log(xray), np.log(radio))))
        print('Spearman R (ymap vs radio): ' + str(spearmanr(np.log(radio), np.log(y))))

        fig, ax = plt.subplots(constrained_layout=True)
        ax.errorbar(np.log10(xray), np.log10(radio), xerr=(0.434 * xray_err / xray),
                    yerr=0.434 * radio_err / radio, fmt='.', ecolor='red', elinewidth=0.4,
                    color='darkred')

        plt.grid(False)
        # ax2 = ax.twiny()
        ax.errorbar(np.log10(y), np.log10(radio), xerr=(y_err / y),
                     yerr=0.434 * radio_err / radio, fmt='.', ecolor='blue', elinewidth=0.4,
                     color='darkblue')
        ax.set_ylim(np.min([np.min(np.log10(radio) - (0.434 * radio_err / radio)), np.min(np.log10(radio) - (0.434 * radio_err / radio))])-0.05,
                     np.max([np.max(np.log10(radio) + (0.434 * radio_err / radio)),np.max(np.log10(radio) + (0.434 * radio_err / radio))])+0.05)
        ax.set_xlim(np.min([np.min(np.log10(y) - (y_err / y)),np.min(np.log10(xray) - (xray_err / xray))])-0.05,
                     np.max([np.max(np.log10(y) + (y_err / y)),np.max(np.log10(xray) + (xray_err / xray))])+0.05)
        # ax.plot(fitline[0], fitline[1], color='darkslateblue', linestyle='--')
        # ax.set_xlim(-1, 1)
        # ax.set_ylim(-0.3, 0.45)
        ax.set_ylabel('log($I_{R}$) [SB/mean(SB)]')
        # ax.set_xlabel('X-ray [SB/mean(SB)]')
        ax.set_xlabel('log($I_{X}$) [SB/mean(SB)] and log(y) [SZ/mean(SZ)]')
        ax.legend(['Radio vs. X-ray', 'Radio vs. SZ'], loc='upper left')
        ax.plot(fity[0], fity[1], color='darkblue', linestyle='--')
        ax.plot(fitxray[0], fitxray[1], color='darkred', linestyle='--')

        plt.tight_layout()
        plt.grid(False)
        plt.savefig(savefig, bbox_inches='tight')

    def plot3d(self, pixelsize=None, dYpix=300, dXpix=210, start=(1495,1275), fitsfile=None, xray=None, savenumpy=None, savefig=None, maskregion=None, halo=None):

        if pixelsize is None:
            pixelsize = 4*int(np.sqrt(self.beamarea/np.pi))

        stepsize = int(pixelsize*np.sin(40))

        if not halo:
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

        if halo:
            region = open('../ptp_dir/'+halo+'/grid_35x35/grid_35x35_ds9_image.reg', 'r').read()
            print(f'Cell size is {self.pix_to_size(0.072)*35} X {self.pix_to_size(0.072)*35}')
            print(f'Cell size is {35*abs((self.header["CDELT2"] * u.deg).to(u.arcsec))} X {35*abs((self.header["CDELT2"] * u.deg).to(u.arcsec))}')
            region = region.split('\n')
            region_head = region[0:2]
            structures = region[3:]
            data = []
            data_error = []
            for structure in structures:
                if not 'box' in structure:
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
                    data.append(np.median(im_mask))
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
                            datatotal[N].append((structure.replace('box(', '').split(',')[0:2], np.median(im_mask), np.std(im_mask)))
                        regionoutput+='\n'+structure

            arr = np.array([[d[1] for d in data] for data in datatotal])
            arr_err = np.array([[d[2] for d in data] for data in datatotal])
            arr[np.isnan(arr)]=0
            arr[arr==0]=np.median(arr)
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
            ha.text(len(X)//8, len(Y), 1, 'A401', color='red', fontsize=15)
            ha.text(len(X)//4, 0, 1, 'A399', color='red', fontsize=15)

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
        print(f'Noise : {str(round(rms * 1000, 2))} {u.mJy/u.beam}')
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
    # Image.plot3d(pixelsize=70, savenumpy='radio3d_6.npy', savefig='radio3d_6.png')
    # Image.do_science(region='../regions/bridge.reg')
    # Image.make_image(show_regions='../boxlayout.reg', vmin=0.00005, save='layout.png', colorbar=False, beam=False, give_scale=True)
    # Image.make_image(show_regions='../tessupdate.reg', vmin=0.00005, save='tess.png', colorbar=False, beam=False, give_scale=False)
    # Image.make_cutout(pos=(int(Image.image_data.shape[0] / 2), int(Image.image_data.shape[0] / 2)), size=(3400, 3400))
    # Image.make_subimages(regionfile='../regions.reg', save='6cutouts.png', beamsize=True)
    # Image.make_image(vmin=0.00005, show_regions='../regions.reg', save='subimagelayout.png', subim=True, colorbar=True)
    # Image.make_subcontour('../regions.reg', save='6subimages.png', fits_lowres='../fits/80all.fits', beamsize=True)
    # Image = Imaging('../fits/6all.fits', resolution=6)
    # Image.make_cutout(pos=(int(Image.image_data.shape[0] / 2), int(Image.image_data.shape[0] / 2)),
    #                   size=(int(Image.image_data.shape[0] / 2), int(Image.image_data.shape[0] / 2)))
    # Image.make_image(vmin=0.00005, text=True, save='a399a401.png')
    # Image.make_image()

    #10"

    #20"
    Image = Imaging('../fits/20all.fits', resolution=20)
    # Image.make_cutout(pos=(int(Image.image_data.shape[0] / 2), int(Image.image_data.shape[0] / 2)), size=(1500, 1500))
    # Image.make_image(convolve=True, save='test20.png', text=True)
    # Image.plot3d(pixelsize=35, savenumpy='radio3d_20.npy', savefig='radio3d_20.png')
    # Image.make_cutout(pos=(int(Image.image_data.shape[0] / 2), int(Image.image_data.shape[0] / 2)),
    #                   size=(int(Image.image_data.shape[0] / 2), int(Image.image_data.shape[0] / 2)))
    # Image.make_contourplot(regions='../regions.reg')
    # Image.make_subcontour('../regions.reg', save='20subimages.png', fits_lowres='../fits/80all.fits', beamsize=False)

    # Image.make_image(save='20image.png', vmin=0.0001)
    Image.remove_compactsources(kernelsize=100, write='../fits/20median.fits')
    # Image.make_image()

    #20" median (will be substitute with bridge? for correlating)
    # Image = Imaging('../fits/20median.fits', resolution=20)
    # Image.make_cutout(pos=(int(Image.image_data.shape[0] / 2), int(Image.image_data.shape[0] / 2)), size=(1500, 1500))
    # Image.make_image(save='justbridge.png', text=True)
    # Image.plot3d(pixelsize=35, savenumpy='radio3d.npy', savefig='radio3d.png')
    # Image.plot3d(pixelsize=35, savenumpy='y.npy', savefig='y3d.png', fitsfile='../fits/a401_curdecmaps_0.2_1.5s_sz.fits')
    # Image.plot3d(pixelsize=35, savenumpy='xray.npy', savefig='xray3d.png', fitsfile='../fits/mosaic_a399_a401.fits', xray=True)
    # Image.plot_corr(savefig='bridgecorr.png')
    # Image.plot3d(savenumpy='a399radio.npy', halo='A399')
    # Image.plot3d(savenumpy='a399y.npy', fitsfile='../fits/a401_curdecmaps_0.2_1.5s_sz.fits', halo='A399')
    # Image.plot3d(savenumpy='a399xray.npy', fitsfile='../fits/mosaic_a399_a401.fits', xray=True, halo='A399')
    # Image.plot_corr(halo='A399', savefig='A399corr.png')
    # Image.plot3d(savenumpy='a401radio.npy', halo='A401')
    # Image.plot3d(savenumpy='a401y.npy', fitsfile='../fits/a401_curdecmaps_0.2_1.5s_sz.fits', halo='A401')
    # Image.plot3d(savenumpy='a401xray.npy', fitsfile='../fits/mosaic_a399_a401.fits', xray=True, halo='A401')
    # Image.plot_corr(halo='A401', savefig='A401corr.png')
    # Image.make_cutout(pos=(int(Image.image_data.shape[0] / 2), int(Image.image_data.shape[0] / 2)), size=(1500, 1500))
    # Image.make_bridge_overlay_yxr_contourplot(fits2='../fits/a401_curdecmaps_0.2_1.5s_sz.fits', fits1='../fits/mosaic_a399_a401.fits',
    #                                           show_regions='corr_area.reg', save='ymapxray.png')

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
    # Image.plot3d(savenumpy='radio3d.npy', savefig='radio3d.png', maskregion='../regions/excluderegions60.reg')
    # Image.plot3d(savenumpy='y.npy', savefig='y3d.png', fitsfile='../fits/a401_curdecmaps_0.2_1.5s_sz.fits', maskregion='../regions/excluderegions60.reg')
    # Image.plot3d(savenumpy='xray.npy', savefig='xray3d.png', fitsfile='../fits/mosaic_a399_a401.fits', xray=True, maskregion='../regions/excluderegions60.reg')
    # Image.make_cutout(pos=(int(Image.image_data.shape[0]/2), int(Image.image_data.shape[0]/2)), size=(850, 850))
    # Image.make_image(vmin=0.002, vmax=0.03, show_regions='../bridgebroken.reg', save='bridge.png', text=True)
    # Image.make_contourplot(regions='../regions.reg')
    # Image.remove_compactsources(kernelsize=51, write='../fits/60median.fits')
    # Image.make_image()

    #60median
    # Image = Imaging('../fits/60median.fits', resolution=60)
    # Image.make_cutout(pos=(int(Image.image_data.shape[0] / 2), int(Image.image_data.shape[0] / 2)), size=(750, 750))
    # Image.plot3d(savenumpy='radio3d.npy', savefig='radio3d.png', maskregion='../regions/excluderegions60.reg')
    # Image.plot3d(savenumpy='y.npy', savefig='y3d.png', fitsfile='../fits/a401_curdecmaps_0.2_1.5s_sz.fits', maskregion='../regions/excluderegions60.reg')
    # Image.plot3d(savenumpy='xray.npy', savefig='xray3d.png', fitsfile='../fits/mosaic_a399_a401.fits', xray=True, maskregion='../regions/excluderegions60.reg')
    # Image.make_bridge_overlay_yxr_contourplot(fits2='../fits/a401_curdecmaps_0.2_1.5s_sz.fits', fits1='../fits/mosaic_a399_a401.fits',
    #                                           show_regions='corr_area.reg', save='ymapxray.png')
    # Image.plot_corr()

    #80"
    # Image = Imaging('../fits/80all.fits', resolution=20)
    # Image.make_cutout(pos=(int(Image.image_data.shape[0]/2), int(Image.image_data.shape[0]/2)), size=(850, 850))
    # Image.make_bridge_overlay_contourplot('../fits/a401_curdecmaps_0.2_1.5s_sz.fits', title=' ',
    #                                       minlevel_1=0, maxlevel_1=0.005, steps_1=100, steps_2=np.sqrt(2), minlevel_2=(10**(-5))/2, convolve_2=30, save='ymap.png')
    # Image.make_bridge_overlay_contourplot('../fits/mosaic_a399_a401.fits', title=' ',
    #                                       minlevel_1=0, maxlevel_1=0.005, steps_1=100, convolve_2=2, steps_2=2, maxlevel_2=50, xray=True)
    # Image.make_bridge_overlay_yxr_contourplot(fits2='../fits/a401_curdecmaps_0.2_1.5s_sz.fits', fits1='../fits/mosaic_a399_a401.fits',
    #                                           show_regions='../ptp_dir/grid_35x35/grid_35x35_ds9_image.reg', save='ymapxray.png')