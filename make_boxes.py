"""
LAST UPDATE: 7-6-2021

Use this script to return region files of directions that need a selfcal.

TODO:
- Sources in same figure check flux with sqrt(a**2+b**2)
- The brighter, the less distance between sources in same image
- Initial size depending on brightness of source.
"""

import numpy as np
import os
from astropy.wcs import WCS
from astropy.io import fits
from astropy.nddata import Cutout2D
from argparse import ArgumentParser
import warnings
from pandas import DataFrame, concat, read_csv

__author__ = "Jurjen de Jong"
__all__ = ['']

warnings.filterwarnings("ignore")

parser = ArgumentParser()
parser.add_argument('-f', '--file', type=str, help='fitsfile name')
parser.add_argument('-l', '--location', type=str, help='data location folder name')
parser.add_argument('-i', '--images', type=bool, default=True, help='return images of boxes')
args = parser.parse_args()

if args.location:
    folder = args.location
    if folder[-1]=='/':
        folder=folder[0:-1]
else:
    folder = ''

#check if folder exists and create if not
folders = folder.split('/')
for i, f in enumerate(folders):
    subpath = '/'.join(folder.split('/')[0:i+1])
    if not os.path.isdir(subpath):
        print(f'Create directory: {subpath}')
        os.system(f'mkdir {subpath}')

if args.images:
    import matplotlib.pyplot as plt
    from matplotlib.colors import SymLogNorm

def resample_pixels(image_data, rows, cols):
    """Resample image by summing pixels together"""
    while image_data.shape[0]%rows!=0:
        rows-=1
    while image_data.shape[0]%cols!=0:
        cols-=1
    return image_data.reshape(int(rows), image_data.shape[0]//rows, int(cols), image_data.shape[1]//cols).sum(axis=1).sum(axis=2)


class Imaging:
    def __init__(self, fits_file: str = None, vmin: float = None, vmax: float = None):
        self.hdu = fits.open(fits_file)[0]
        self.image_data = self.hdu.data
        while len(self.image_data.shape) != 2:
            self.image_data = self.image_data[0]
        self.wcs = WCS(self.hdu.header, naxis=2)
        if vmin is None:
            self.vmin = np.nanstd(self.image_data)
        else:
            self.vmin = vmin
        if vmax is None:
            self.vmax = np.nanstd(self.image_data)*25

    def imaging(self, image_data=None, cmap: str = 'CMRmap'):
        """
        Image your data with this method.
        image_data -> insert your image_data or plot full image
        cmap -> choose your preferred cmap
        """

        if image_data is None:
            image_data = self.image_data
        plt.figure(figsize=(10, 10))
        plt.subplot(projection=self.wcs)
        plt.imshow(image_data, norm=SymLogNorm(linthresh=self.vmin/20, vmin=self.vmin/50, vmax=self.vmax), origin='lower', cmap=cmap)
        plt.xlabel('Galactic Longitude')
        plt.ylabel('Galactic Latitude')
        plt.show()

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

class SetBoxes(Imaging):
    def __init__(self, fits_file: str = None, initial_box_size: float = 0.2, peak_flux=0.07):
        self.pix_y = None
        self.pix_x = None
        self.flux = None
        self.image_number = None
        self.initial_box_size = initial_box_size
        super().__init__(fits_file=fits_file)
        self.wcs_cut = self.wcs

        # get the positions of the flux peaks from full image
        self.peak_flux = peak_flux

        df_peaks = DataFrame([{'pix_y': c[0], 'pix_x': c[1], 'flux': self.image_data[c[0], c[1]], 'resampled':False}
                                   for c in np.argwhere(self.image_data > peak_flux)]). \
                            sort_values('flux', ascending=False)

        # find from resampled data more sources, which are bright sources over bigger area
        resampled = resample_pixels(self.image_data,
                                    rows=int(self.image_data.shape[0] / 100),
                                    cols=int(self.image_data.shape[1] / 100))
        # plt.imshow(resampled, norm=SymLogNorm(linthresh=np.nanstd(resampled)/20, vmin=np.nanstd(resampled)/50, vmax=np.nanstd(resampled)*25), origin='lower')
        # plt.show()
        resample_x, resample_y = resampled.shape
        original_x, original_y = self.image_data.shape
        resample_scale_x = original_x//resample_x
        resample_scale_y = original_y//resample_y

        resampled_df_peaks = DataFrame([{'pix_y': c[0]*resample_scale_x, 'pix_x': c[1]*resample_scale_y,
                                              'flux': resampled[c[0], c[1]], 'resampled':True}
                                   for c in np.argwhere(resampled > 15)])

        self.df_peaks = concat([df_peaks, resampled_df_peaks], axis=0).\
            drop_duplicates(subset=['pix_y', 'pix_x']).reset_index(drop=True)

    # euclidean distance
    @staticmethod
    def ed_array(pos=None, lst=None):
        """Euclidean distance between position and list"""
        return np.sqrt(np.sum(np.square(pos - lst), axis=1))

    def reposition(self):
        """Reposition image by looking at the data points near the boundaries."""

        # calculate percentage of high flux points at the boundaries
        outlier_threshold = 4 * np.std(self.image_data)  # use five times the standard deviation of full image as outlier threshold

        def boundary_perc(image_data):
            """Percentage of high flux"""
            outliers = (image_data > outlier_threshold).astype(
                int)  # points are considered outliers when 2 times bigger than the rms/std
            boundary_size = int(image_data.shape[0] / 10)
            left = outliers.T[0:boundary_size]  # left boundary
            lower = outliers[0:boundary_size]  # lower boundary
            right = outliers.T[len(outliers) - boundary_size:len(outliers)]  # right boundary
            upper = outliers[len(outliers) - boundary_size:len(outliers)]  # upper boundary

            left_p = np.sum(left) / len(left.flatten())  # percentage of high pixels left boundary
            lower_p = np.sum(lower) / len(lower.flatten())  # percentage of high pixels lower boundary
            right_p = np.sum(right) / len(right.flatten())  # percentage of high pixels right boundary
            upper_p = np.sum(upper) / len(upper.flatten())  # percentage of high pixels upper boundary

            return left_p, lower_p, right_p, upper_p

        def boundary_sources(image_data, threshold=0.01):
            """Sources within the boundaries of the box"""
            other_source = (image_data > threshold).astype(int)
            boundary_size = int(image_data.shape[0] / 10)
            left = other_source.T[0:boundary_size]  # left boundary
            lower = other_source[0:boundary_size]  # lower boundary
            right = other_source.T[len(other_source) - boundary_size:len(other_source)]  # right boundary
            upper = other_source[len(other_source) - boundary_size:len(other_source)]  # upper boundary
            total = np.sum(left)+np.sum(right)+np.sum(upper)+np.sum(lower)
            return total>0

        self.after = self.before.copy()  # image data before data after repositioning

        #get pixel scale
        pixscale = self.wcs.to_header()['CDELT1']
        max_size = abs(int(0.4//pixscale))
        min_size = abs(int(0.2//pixscale))

        step_size = np.max(self.wcs.pixel_scale_matrix)  # step size in degrees per pixel
        im_size = int(self.initial_box_size / step_size)  # start with square boxes with a size of 0.4 degrees
        threshold_p = 0.000005  # max percentage of boundary elements

        for N in range(3):#looping multiple times
            # STEP 1: Reposition for higher flux around the borders
            n = 0
            left_p, lower_p, right_p, upper_p = boundary_perc(self.after)  # get boundary percentages
            while (left_p > threshold_p or lower_p > threshold_p or right_p > threshold_p or upper_p > threshold_p or
                boundary_sources(image_data=self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))[0], threshold=0.01)) \
                    and n < 100 and max_size > im_size > min_size:  # image cannot be smaller than 0.3 arcsec
                if lower_p > threshold_p and lower_p > upper_p and self.pix_y - int(im_size / 2) > 0 \
                        and not boundary_sources(image_data=self.make_cutout((self.pix_x, self.pix_y-int(im_size/5)), (im_size, im_size))[0], threshold=0.01):
                    self.pix_y -= int(self.after.shape[0] / 200)  # size shifted down
                    self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                    left_p, lower_p, right_p, upper_p = boundary_perc(self.after)
                elif upper_p > threshold_p and int(im_size / 2) + self.pix_y < self.image_data.shape[0] \
                        and not boundary_sources(image_data=self.make_cutout((self.pix_x, self.pix_y+int(im_size/5)), (im_size, im_size))[0], threshold=0.01):
                    self.pix_y += int(self.after.shape[0] / 200)  # size shifted above
                    self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                    left_p, lower_p, right_p, upper_p = boundary_perc(self.after)
                elif boundary_sources(image_data=self.make_cutout((self.pix_x, self.pix_y+int(im_size/5)), (im_size, im_size))[0], threshold=0.01):
                    im_size -= int(self.after.shape[0]) / 200  # size increase
                    self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y),
                                                                (im_size, im_size))  # new image
                    left_p, lower_p, right_p, upper_p = boundary_perc(self.after)
                n += 1
                if left_p > threshold_p and left_p > right_p and int(im_size / 2) + self.pix_x < self.image_data.shape[0] \
                        and not boundary_sources(image_data=self.make_cutout((self.pix_x-int(im_size/5), self.pix_y), (im_size, im_size))[0], threshold=0.01):
                    self.pix_x -= int(self.after.shape[0] / 200)  # size shifted to the left
                    self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                    left_p, lower_p, right_p, upper_p = boundary_perc(self.after)
                elif right_p > threshold_p and self.pix_x - int(im_size / 2) > 0 and not \
                        boundary_sources(image_data=self.make_cutout((self.pix_x+int(im_size/5), self.pix_y), (im_size, im_size))[0], threshold=0.01):
                    self.pix_x += int(self.after.shape[0] / 200)  # size shifted to the right
                    self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                    left_p, lower_p, right_p, upper_p = boundary_perc(self.after)
                elif boundary_sources(image_data=self.make_cutout((self.pix_x, self.pix_y+int(im_size/5)), (im_size, im_size))[0], threshold=0.01):
                    im_size -= int(self.after.shape[0]) / 200  # size increase
                    self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y),
                                                                (im_size, im_size))  # new image
                    left_p, lower_p, right_p, upper_p = boundary_perc(self.after)
                n += 1

            # STEP 2: Resizing
            while ((left_p < threshold_p*7 and lower_p < threshold_p*7 and right_p < threshold_p*7 and upper_p < threshold_p*7) or
                   boundary_sources(image_data=self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))[0],
                       threshold=0.01)) \
                    and max_size > im_size > min_size \
                    and n<100:
                if im_size<=(max_size+min_size)/2:
                    im_size += int(self.after.shape[0] / 200)  # size reduced image
                    self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                    left_p, lower_p, right_p, upper_p = boundary_perc(self.after)
                else:
                    im_size -= int(self.after.shape[0] / 200)  # size reduced image
                    self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                    left_p, lower_p, right_p, upper_p = boundary_perc(self.after)
                n += 1

            # STEP 3: Reposition for high flux near the borders. These are bright sources that we want to have more in the center
            peak_y, peak_x = np.argwhere(self.after > self.peak_flux / 2).T / im_size
            n = 0
            while (np.sum(peak_y > 0.75) or np.sum(peak_y < 0.25) or np.sum(peak_x > 0.75) or np.sum(
                    peak_x < 0.25)) and n < 100 and max_size > im_size > min_size:
                if np.sum(peak_y > 0.75) and not np.sum(peak_y < 0.25) and not \
                        boundary_sources(image_data=self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))[0], threshold=0.02):
                    self.pix_y += int(self.after.shape[0] / 200)
                    self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                    peak_y, peak_x = np.argwhere(self.after > self.peak_flux / 2).T / im_size
                elif np.sum(peak_y < 0.25) and not np.sum(peak_y > 0.75) and not \
                        boundary_sources(image_data=self.make_cutout((self.pix_x, self.pix_y-int(im_size/8)), (im_size, im_size))[0], threshold=0.02):
                    self.pix_y -= int(self.after.shape[0] / 200)
                    self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                    peak_y, peak_x = np.argwhere(self.after > self.peak_flux / 2).T / im_size
                elif np.sum(peak_y > 0.75) and np.sum(peak_y < 0.25) and not \
                        boundary_sources(image_data=self.make_cutout((self.pix_x, self.pix_y), (im_size-int(im_size/8), im_size-int(im_size/8)))[0], threshold=0.02):
                    im_size += int(self.after.shape[0]) / 200  # size increase
                    self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                    peak_y, peak_x = np.argwhere(self.after > self.peak_flux / 2).T / im_size
                if np.sum(peak_x > 0.75) and not np.sum(peak_x < 0.25) and not \
                        boundary_sources(image_data=self.make_cutout((self.pix_x+int(im_size/8), self.pix_y), (im_size, im_size))[0], threshold=0.02):
                    self.pix_x += int(self.after.shape[0] / 200)
                    self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                    peak_y, peak_x = np.argwhere(self.after > self.peak_flux / 2).T / im_size
                elif np.sum(peak_x < 0.25) and not np.sum(peak_x > 0.75) and not \
                        boundary_sources(image_data=self.make_cutout((self.pix_x-int(im_size/8), self.pix_y), (im_size, im_size))[0], threshold=0.02):
                    self.pix_x -= int(self.after.shape[0] / 200)
                    self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                    peak_y, peak_x = np.argwhere(self.after > 0.07).T / im_size
                elif np.sum(peak_x < 0.25) and np.sum(peak_x > 0.75) and not \
                        boundary_sources(image_data=self.make_cutout((self.pix_x, self.pix_y), (im_size-int(im_size/8), im_size-int(im_size/8)))[0], threshold=0.02):
                    im_size += int(self.after.shape[0]) / 200  # size increase
                    self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                    peak_y, peak_x = np.argwhere(self.after > self.peak_flux / 2).T / im_size
                n += 1

            # STEP 4: Resizing based on flux on boundary
            while (left_p * 2 < threshold_p and lower_p * 2 < threshold_p and right_p * 2 < threshold_p and upper_p * 2 < threshold_p) \
                    and max_size > im_size > min_size and n < 200:
                if im_size<=(max_size+min_size)/2:
                    im_size += int(self.after.shape[0] / 400)  # size reduced image
                    self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                    left_p, lower_p, right_p, upper_p = boundary_perc(self.after)
                else:
                    im_size -= int(self.after.shape[0] / 400)  # size reduced image
                    self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                    left_p, lower_p, right_p, upper_p = boundary_perc(self.after)
                n+=1

            # STEP 5: Resizing based on flux outside of image
            while boundary_sources(image_data=self.make_cutout((self.pix_x, self.pix_y), (im_size-int(im_size/5), im_size-int(im_size/5)))[0], threshold=0.01) \
                and n < 400 and max_size > im_size > min_size:
                im_size -= int(self.after.shape[0]) / 400  # size increase
                self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y),
                                                            (im_size, im_size))
                n+=1

            # Step 6: Reposition
            n, t2 = 0, 0
            left_p, lower_p, right_p, upper_p = boundary_perc(self.after)  # get boundary percentages
            while (left_p > threshold_p or lower_p > threshold_p or right_p > threshold_p or upper_p > threshold_p or
                   boundary_sources(image_data=self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))[0],
                       threshold=0.01))\
                    and n < 200 and max_size > im_size > min_size:
                t1 = 0
                if upper_p > threshold_p and int(im_size / 2) + self.pix_y < self.image_data.shape[0] and not \
                        boundary_sources(image_data=self.make_cutout((self.pix_x, self.pix_y+int(im_size/8)), (im_size, im_size))[0], threshold=0.04):
                    self.pix_y += int(self.after.shape[0] / 300)  # size shifted above
                    new_image, _ = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))
                    if np.sum(boundary_perc(new_image)) < np.sum(boundary_perc(self.after)):  # check if improvement
                        self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                        left_p, lower_p, right_p, upper_p = boundary_perc(self.after)
                    else:  # correct back
                        self.pix_y -= int(self.after.shape[0] / 300)
                        t1 += 1
                elif lower_p > threshold_p and self.pix_y - int(im_size / 2) > 0 and not \
                        boundary_sources(image_data=self.make_cutout((self.pix_x, self.pix_y-int(im_size/8)), (im_size, im_size))[0], threshold=0.04):
                    self.pix_y -= int(self.after.shape[0] / 300)  # size shifted down
                    new_image, _ = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))
                    if np.sum(boundary_perc(new_image)) < np.sum(boundary_perc(self.after)):  # check if improvement
                        self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                        left_p, lower_p, right_p, upper_p = boundary_perc(self.after)
                    else:  # correct back
                        self.pix_y += int(self.after.shape[0] / 300)
                        t1 += 1
                n += 1
                if right_p > threshold_p and self.pix_x - int(im_size / 2) > 0 and not \
                        boundary_sources(image_data=self.make_cutout((self.pix_x+int(im_size/8), self.pix_y), (im_size, im_size))[0], threshold=0.04):
                    self.pix_x += int(self.after.shape[0] / 300)  # size shifted to the right
                    new_image, _ = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))
                    if np.sum(boundary_perc(new_image)) < np.sum(boundary_perc(self.after)):  # check if improvement
                        self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                        left_p, lower_p, right_p, upper_p = boundary_perc(self.after)
                    else:  # correct back
                        self.pix_x -= int(self.after.shape[0] / 300)
                        t1 += 1
                elif left_p > threshold_p and left_p > right_p and int(im_size / 2) + self.pix_x < self.image_data.shape[0] and not \
                        boundary_sources(image_data=self.make_cutout((self.pix_x-int(im_size/8), self.pix_y), (im_size, im_size))[0], threshold=0.04):
                    self.pix_x -= int(self.after.shape[0] / 300)  # size shifted to the left
                    new_image, _ = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))
                    if np.sum(boundary_perc(new_image)) < np.sum(boundary_perc(self.after)):  # check if improvement
                        self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                        left_p, lower_p, right_p, upper_p = boundary_perc(self.after)
                    else:  # correct back
                        self.pix_x += int(self.after.shape[0] / 300)
                        t1 += 1
                n += 1

            # STEP 7: Resizing based on flux within and on the borders
            while im_size<max_size and np.sum(self.after>np.max(self.after)/10)/self.after.size*self.before.size>50:
                im_size += int(self.after.shape[0] / 100)  # size reduced image
                self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image

            n=0
            # STEP 8: Resizing and moving if needed
            while boundary_sources(image_data=self.make_cutout((self.pix_x, self.pix_y), (im_size+int(im_size/30), im_size+int(im_size/30)))[0], threshold=0.01)\
                    and n<200:
                if im_size>min_size:
                    im_size -= int(self.after.shape[0]) / 100  # size increase
                    self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y),
                                                                (im_size, im_size))
                else:
                    im_size += int(self.after.shape[0]) / 200  # size increase
                    self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y),
                                                                (im_size, im_size))
                if n>100:#if still not optimized we can move around again.
                    if not boundary_sources(image_data=self.make_cutout((self.pix_x, self.pix_y - int(im_size / 10)),
                                                                      (im_size, im_size))[0], threshold=0.01):
                        self.pix_y -= int(self.after.shape[0] / 300)  # size shifted down
                        self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y),
                                                                    (im_size, im_size))  # new image
                    elif not boundary_sources(image_data=self.make_cutout((self.pix_x, self.pix_y + int(im_size / 10)),
                                                                      (im_size, im_size))[0], threshold=0.01):
                        self.pix_y += int(self.after.shape[0] / 300)  # size shifted above
                        self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y),
                                                                    (im_size, im_size))  # new image
                    if not boundary_sources(image_data=self.make_cutout((self.pix_x - int(im_size / 10), self.pix_y),
                                                                      (im_size, im_size))[0], threshold=0.01):
                        self.pix_x -= int(self.after.shape[0] / 300)  # size shifted to the left
                        self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y),
                                                                    (im_size, im_size))  # new image
                    elif not boundary_sources(image_data=self.make_cutout((self.pix_x + int(im_size / 10), self.pix_y),
                                                              (im_size, im_size))[0], threshold=0.01):
                        self.pix_x += int(self.after.shape[0] / 300)  # size shifted to the right
                        self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y),
                                                                    (im_size, im_size))  # new image
                n+=1
        return self

    def save_box(self, box_name: str = 'box.reg'):
        """
        save the box as an .reg file
        """
        position_h = [v.replace('h', ':').replace('m', ':').replace('s', '').replace('d', ':') for v in
                      self.wcs.pixel_to_world(self.pix_x, self.pix_y).to_string('hmsdms').split()]
        arcsec_size = round(np.max(self.wcs.pixel_scale_matrix) * self.after.shape[0] * 3600, 3)
        f = open(box_name, "a")
        f.write(f'box({position_h[0]},{position_h[1]},{str(arcsec_size)}\",{str(arcsec_size)}\",0) # text={{{box_name.split("/")[-1].split(".reg")[0]}}}')
        f.close()

        return self

    def source_to_csv(self, sources):
        """Write source to a csv file"""
        with open(f'{folder}/source_file.csv', 'a+') as source_file:
            source_file.write(f'{self.image_number},{self.pix_x},{self.pix_y},{self.flux},{str(self.after.shape).replace(",", ";")},{str(sources).replace(",", ";")}\n')
        return self

    @staticmethod
    def intersecting_boxes(position, positions):
        """
        This method checks if the new box will intersect with another box
        position -> consists of the x-position, y-position, image size respectively
        positions -> consists of a list with already accepted positions to compare with
        """
        if len(positions) == 0:
            return False
        positions = np.array(positions).T
        if np.all((position[2] + positions[2]) / 2 < abs(position[0] - positions[0])) and \
                np.all((position[2] + positions[2]) / 2 < abs(position[1] - positions[1])):
            return False
        else:
            return True

    @property
    def other_sources_in_image(self):
        """Check if there are other sources in your image.
        Score of 4 means that it is in the same image
        Score of 2 means that it is nearby (near boundary)"""
        self.df_peaks['nearby'] = (abs(self.pix_x - self.df_peaks.pix_x) <= abs(self.after.shape[0]) / 2).astype(int) + \
                                    (abs(self.pix_y - self.df_peaks.pix_y) <= abs(self.after.shape[0]) / 2).astype(int) + \
                                      (abs(self.pix_x - self.df_peaks.pix_x) <= abs(self.after.shape[0])).astype(int) + \
                                      (abs(self.pix_y - self.df_peaks.pix_y) <= abs(self.after.shape[0])).astype(int)

        other_sources_in_image = self.df_peaks[(self.df_peaks.nearby == 4) & (self.df_peaks.index != self.image_number)]
        sources_nearby = self.df_peaks[(self.df_peaks.nearby == 2) & (self.df_peaks.index != self.image_number)]

        self.number_of_sources(list(other_sources_in_image.index))

        return list(other_sources_in_image.index), list(sources_nearby.index)

    @staticmethod
    def center_cutout(image_data):
        return Cutout2D(data=image_data,
                        position=(int(image_data.shape[0] / 2), int(image_data.shape[0] / 2)),
                        size=(int(image_data.shape[0] / 8), int(image_data.shape[0] / 8))).data

    def make_initial_box(self):
        """Check if source is interesting for selfcal."""

        im_size = int(self.initial_box_size / np.max(self.wcs.pixel_scale_matrix))

        self.initial_pos = (self.pix_x, self.pix_y)
        self.before, self.wcs_cut = self.make_cutout(pos=self.initial_pos,
                                       size=(im_size, im_size))  # image data after repositioning

        return self

    def number_of_sources(self, sources):
        """Count number of sources on a distance from each other, to see if split is possible"""
        peaks_in_image = self.df_peaks.iloc[sources]
        positions = peaks_in_image[['pix_x', 'pix_y']].to_numpy()
        n=0
        #count number of sources in image
        while len(positions)>n+1:
            peaks_in_image['distance'] = self.ed_array(positions[n], positions)
            peaks_in_image = peaks_in_image[peaks_in_image['distance']>200]
            positions = peaks_in_image[['pix_x', 'pix_y']].to_numpy()
            n+=1
        if n>1:
            print(f'There are multiple sources in same box.\n Splitting possible.')
        return self


if __name__ == '__main__':

    image = SetBoxes(fits_file=args.file, initial_box_size=0.2)

    if args.images:
        os.system(f'rm -rf {folder}/box_images; mkdir {folder}/box_images')  # make folder with box images
    os.system(f'rm -rf {folder}/boxes; mkdir {folder}/boxes')  # make folder with the .reg files
    os.system(f'rm source_file.csv') # clean up

    with open(f'{folder}/source_file.csv', 'w') as source_file:
        source_file.write('id,pix_x,pix_y,flux,size,sources\n')

    sources_done = []
    print(f'We found {len(image.df_peaks)} high flux peaks.\n')

    m, r = 0, 0
    for n, p in enumerate(image.df_peaks.to_dict(orient="records")):

        replace=False

        # skip sources that are already displayed in other boxes
        if n in sources_done:
            continue

        # set values for source of interest
        image.pix_x, image.pix_y, image.flux, image.image_number = p['pix_x'], p['pix_y'], p['flux'], n

        image.make_initial_box()

        # reposition box
        image.reposition()

        # we now check if there are in our box sources from our list of peak sources, which we can skip later on
        other_sources, _ = image.other_sources_in_image

        # check if boxes contain multiple times the same source. If so, we replace this box with a better new one.
        found=False
        for source in other_sources:
            if source in sources_done:
                sources = read_csv('source_file.csv')['sources']
                for M, source_list in enumerate(sources):
                    if len(source_list.replace('[','').replace(']',''))>0:
                        source_list = [int(s) for s in source_list.replace('[','').replace(']','').replace(' ','').split(';')]
                        if bool(set(other_sources) & set(source_list)):
                            os.system(f'rm {folder}/box_images/box_{M+1}.png')
                            os.system(f'rm {folder}/boxes/box_{M+1}.reg')
                            replace, found = True, True
                            break
            if found:
                break

        sources_done += other_sources
        image.source_to_csv(other_sources)

        # make image with before and after repositioning of our box
        if not replace:
            m += 1

        if args.images:
            fig = plt.figure(figsize=(10, 10))
            plt.subplot(1, 2, 1, projection = image.wcs_cut)
            plt.title(f'Initial image')
            plt.imshow(image.before, norm=SymLogNorm(linthresh=image.vmin/10, vmin=image.vmin/10, vmax=image.vmax/2), origin='lower',
                          cmap='CMRmap')
            plt.subplot(1, 2, 2, projection = image.wcs_cut)
            plt.title('Repositioned')
            plt.imshow(image.after, norm=SymLogNorm(linthresh=image.vmin/10, vmin=image.vmin/20, vmax=image.vmax), origin='lower',
                          cmap='CMRmap')
            if replace:
                fig.savefig(f'{folder}/box_images/box_{M+1}.png')
            else:
                fig.savefig(f'{folder}/box_images/box_{m}.png')

        if replace:
            print(f'Replace box {M+1}.')
            image.save_box(box_name=f'{folder}/boxes/box_{M+1}.reg')
        else:
            print(f'Create box {m}.')
            image.save_box(box_name=f'{folder}/boxes/box_{m}.reg')

    print('-------------------------------------------------')
    print(f'Made succesfully {m} boxes.')
    if args.images:
        print(f'Images of boxes are in {folder}/box_images.')
    print(f'Region files are in {folder}/boxes.')