import warnings
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy import coordinates
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from matplotlib.colors import LogNorm, SymLogNorm, PowerNorm
import pyregion
from pyregion.mpl_helper import properties_func_default
import sys

warnings.filterwarnings("ignore")

__all__ = ["MakeImages"]


def get_rms(image_data):
    """
    from Cyril Tasse/kMS

    :param image_data: image data array
    :return: rms (noise measure)
    """

    from past.utils import old_div

    maskSup = 1e-7
    m = image_data[np.abs(image_data) > maskSup]
    rmsold = np.std(m)
    diff = 1e-1
    cut = 3.
    med = np.median(m)
    for _ in range(10):
        ind = np.where(np.abs(m - med) < rmsold * cut)[0]
        rms = np.std(m[ind])
        if np.abs(old_div((rms - rmsold), rmsold)) < diff: break
        rmsold = rms
    print(f'Noise : {str(round(rms * 1000, 4))} {u.mJy / u.beam}')
    return rms


class MakeImages:
    def __init__(self, fits_file=None):
        """
        Make LOFAR images
        ------------------------------------------------------------
        :param fits_file: Fits file name and path (can also put url)
        """

        self.fits_file = fits_file
        self.hdu = fits.open(fits_file)
        self.image_data = self.hdu[0].data
        self.header = self.hdu[0].header
        while len(self.image_data.shape) != 2:
            self.image_data = self.image_data[0]
        self.wcs = WCS(self.hdu[0].header, naxis=2)

    def imaging(self, image_name: str = None, dpi: int = 50, cmap: str = "CMRmap",
                norm: str = 'SquareRootNorm', vmin: float = None, vmax: float = None,
                cbar: bool = True):
        """
        Imaging of your data.
        ------------------------------------------------------------
        :param image_name: name of your output image
        :param dpi: dots per inch
        :param cmap: cmap of your image
        :param norm: imaging norm [LogNorm, SymLogNorm, SquareRootNorm, SquareNorm, Linear]
        :param vmin: vmin image
        :param vmax: vmax image
        :param cbar: colorbar
        """

        if vmin is None:
            vmin = get_rms(self.image_data)
        if vmax is None:
            vmax = 20 * vmin

        shape_factor = self.image_data.shape[0]/self.image_data.shape[1]
        figsize=12

        plt.figure(figsize=(figsize, int(figsize/shape_factor)))
        plt.subplot(projection=self.wcs)

        if norm == 'LogNorm':
            im = plt.imshow(self.image_data,
                       norm=LogNorm(vmin=vmin, vmax=vmax),
                       origin="lower",
                       cmap=cmap,
                       )
        elif norm == 'SymLogNorm':
            plt.imshow(self.image_data,
                       norm=SymLogNorm(vmin=vmin / 10, vmax=vmax, linthresh=vmin),
                       origin="lower",
                       cmap=cmap,
                       )

        elif norm == 'SquareRootNorm':
            im = plt.imshow(self.image_data,
                       norm=PowerNorm(gamma=1 / 2, vmax=vmax, vmin=vmin),
                       origin="lower",
                       cmap=cmap,
                       )
        elif norm == 'SquareNorm':
            im = plt.imshow(self.image_data,
                       norm=PowerNorm(gamma=2, vmax=vmax, vmin=vmin),
                       origin="lower",
                       cmap=cmap,
                       )
        elif norm == 'Linear':
            im = plt.imshow(self.image_data,
                       vmax=vmax,
                       vmin=vmin,
                       origin="lower",
                       cmap=cmap,
                       )
        else:
            print("Choose from: [LogNorm, SymLogNorm, SquareRootNorm, SquareNorm, Linear]")

        plt.ylabel("Declination", size=16)
        plt.xlabel("Right ascension", size=16)
        plt.xticks(size=16)
        plt.yticks(size=16)

        if cbar:
            cb = plt.colorbar(im, orientation='horizontal')
            cb.set_label('Surface brightness [mJy/beam]', size=16)
            cb.ax.tick_params(labelsize=16)

        plt.grid(False)
        plt.grid('off')
        plt.tight_layout()
        if image_name is None:
            plt.show()
        else:
            plt.savefig(image_name,
                        bbox_inches="tight",
                        dpi=dpi,
                        )
        plt.close()

        return self

    def make_cutout(self, pos: tuple = None, size: tuple = (1000, 1000), region: str = None):
        """
        Make cutout from your image.
        ------------------------------------------------------------
        :param pos: position in pixel size or degrees (RA, DEC)
        :param size: size of your image in pixels
        :param region: pyregion file (if given it ignores pos and size)
        """

        if region is not None:
            r = pyregion.open(region).as_imagecoord(header=self.hdu[0].header)
            mask = r.get_mask(hdu=self.hdu[0],
                              shape=(self.hdu[0].header["NAXIS1"], self.hdu[0].header["NAXIS2"])).astype(np.int16)
            if len(r) > 1:
                sys.exit('Multiple regions in 1 file given, only one allowed')
            else:
                shape = np.array(r[0].coord_list)
                # center = self.wcs.all_pix2world(shape[0], shape[1], 0)
                if len(shape) == 3: # circle
                    cutout = Cutout2D(data=self.image_data * mask,
                                        position=(shape[0], shape[1]),
                                        size=(shape[2], shape[2]),
                                        wcs=self.wcs,
                                        mode='partial')
                elif len(shape) > 3: # square
                    cutout = Cutout2D(data=self.image_data * mask,
                                        position=(shape[0], shape[1]),
                                        size=(shape[3], shape[2]),
                                        wcs=self.wcs,
                                        mode='partial')

        else:
            if type(pos[0]) != int and type(pos[1]) != int:
                pos = self._to_pixel(pos[0], pos[1])
            cutout = Cutout2D(data=self.image_data,
                              position=pos,
                              size=size,
                              wcs=self.wcs,
                              mode="partial")

        self.image_data = cutout.data
        header = cutout.wcs.to_header()
        for key in self.header.keys():
            try:
                self.header[key] = header[key]
            except KeyError:
                pass
        self.hdu[0].data = self.image_data
        self.wcs = WCS(self.header, naxis=2)
        self.hdu[0].header = self.header
        return self

    def _to_pixel(self, ra: float = None, dec: float = None):
        """
        To pixel position from RA and DEC in degrees
        ------------------------------------------------------------
        :param ra: Right ascension (degrees)
        :param dec: Declination (degrees)
        :return: Pixel of position
        """

        from astropy import wcs
        position = coordinates.SkyCoord(
            ra, dec, frame=wcs.utils.wcs_to_celestial_frame(self.wcs).name, unit=(u.degree, u.degree)
        )
        position = np.array(position.to_pixel(self.wcs))
        return position

    def make_fits(self, filename: str = None):
        """
        Make image cutout and make fitsfile
        ------------------------------------------------------------
        :param filename: name of output fits image
        """

        self.hdu[0].writeto(filename, overwrite=True, output_verify='ignore')
        return self


if __name__ == "__main__":
    print("Cannot call script directly.")
