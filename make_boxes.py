"""
Use this script to return region files of directions that need a selfcal.
"""

import numpy as np
import os
from astropy.wcs import WCS
from astropy.io import fits
from astropy.nddata import Cutout2D
from scipy.ndimage import gaussian_filter
from argparse import ArgumentParser
import warnings
from pandas import DataFrame

__author__ = "Jurjen de Jong"
__all__ = ['SetBoxes']

warnings.filterwarnings("ignore")

parser = ArgumentParser()
parser.add_argument('-f', '--file', type=str, help='fitsfile name')
parser.add_argument('-l', '--location', type=str, help='data location folder name')
parser.add_argument('-i', '--images', type=bool, default=False, help='return images of boxes')
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
    def __init__(self, fits_file: str = None, initial_box_size: float = 0.35, peak_flux=0.07):
        self.flux = None
        self.image_number = None
        self.initial_box_size = initial_box_size
        super().__init__(fits_file=fits_file)
        self.wcs_cut = self.wcs


        # get the positions of the flux peaks from full image
        self.peak_flux = peak_flux
        self.df_peaks = DataFrame([{'pix_x': c[0], 'pix_y': c[1], 'flux': self.image_data[c[0], c[1]]} for c in
                                   np.argwhere(self.image_data > peak_flux)]). \
            sort_values('flux', ascending=False)

        peaks = self.df_peaks[['pix_x', 'pix_y']].to_numpy()

        # peaks sources close to each other are considered from the same sources, so these can be removed
        for i in range(len(peaks)):
            if i > len(peaks) - 1:
                break
            # points that are too close to each other have to be removed
            peaks = peaks[(self.ed_array(peaks[i], peaks) == 0) |
                          (self.ed_array(peaks[i], peaks) >
                           int(self.initial_box_size / np.max(self.wcs.pixel_scale_matrix)) / 100)]

        self.df_peaks = DataFrame([{'pix_y': c[0], 'pix_x': c[1], 'flux': self.image_data[c[0], c[1]]} for c in peaks]). \
            sort_values('flux', ascending=False)

    # euclidean distance
    @staticmethod
    def ed_array(pos=None, lst=None):
        """Euclidean distance between position and list"""
        return np.sqrt(np.sum(np.square(pos - lst), axis=1))

    def reposition(self):
        """Reposition image by looking at the data points near the boundaries."""

        # calculate percentage of high flux points at the boundaries
        outlier_threshold = 2 * np.std(
            self.image_data)  # use two times the standard deviation of full image as outlier threshold

        def boundary_perc(image_data):

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

        self.after = self.before.copy()  # image data before data after repositioning

        step_size = np.max(self.wcs.pixel_scale_matrix)  # step size in degrees per pixel
        im_size = int(self.initial_box_size / step_size)  # start with square boxes with a size of 0.4 degrees
        threshold_p = 0.0005  # max percentage of boundary elements

        # STEP 1: Reposition for higher flux around the borders
        n = 0
        left_p, lower_p, right_p, upper_p = boundary_perc(self.after)  # get boundary percentages
        while (left_p > threshold_p or lower_p > threshold_p or right_p > threshold_p or upper_p > threshold_p) \
                and n < 100 and im_size > self.before.shape[0] * 3 / 4:  # image cannot be smaller than 0.3 arcsec
            if lower_p > threshold_p and lower_p > upper_p and self.pix_y - int(im_size / 2) > 0:
                self.pix_y -= int(self.after.shape[0] / 200)  # size shifted down
                self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                left_p, lower_p, right_p, upper_p = boundary_perc(self.after)
                n += 1
            elif upper_p > threshold_p and int(im_size / 2) + self.pix_y < self.image_data.shape[0]:
                self.pix_y += int(self.after.shape[0] / 200)  # size shifted above
                self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                left_p, lower_p, right_p, upper_p = boundary_perc(self.after)
                n += 1
            if left_p > threshold_p and left_p > right_p and int(im_size / 2) + self.pix_x < self.image_data.shape[0]:
                self.pix_x -= int(self.after.shape[0] / 200)  # size shifted to the left
                self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                left_p, lower_p, right_p, upper_p = boundary_perc(self.after)
                n += 1
            elif right_p > threshold_p and self.pix_x - int(im_size / 2) > 0:
                self.pix_x += int(self.after.shape[0] / 200)  # size shifted to the right
                self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                left_p, lower_p, right_p, upper_p = boundary_perc(self.after)
                n += 1
            if n > 5:
                im_size_new = im_size - int(min(max(left_p, lower_p, right_p, upper_p) / threshold_p * 10,
                                                self.after.shape[0] / 200))  # size reduced image
                new_image, _ = self.make_cutout((self.pix_x, self.pix_y), (im_size_new, im_size_new))  # new image
                if np.mean(boundary_perc(new_image)) \
                        < np.mean([left_p, lower_p, right_p, upper_p]):
                    im_size = im_size_new  # accept new image size
                    self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size_new, im_size_new))  # accept image
                    left_p, lower_p, right_p, upper_p = boundary_perc(self.after)

        # Don't accept if image is worse
        if np.mean(boundary_perc(self.after)) > np.mean(boundary_perc(self.before)):  # accept image if better
            self.after = self.before.copy()
            self.pix_x, self.pix_y = self.initial_pos

        # STEP 2: Extra resize step because we prefer smaller images (not smaller than 0.2 degrees)
        while (left_p < threshold_p and lower_p < threshold_p and right_p < threshold_p and upper_p < threshold_p) \
                and im_size > self.before.shape[0] * 1 / 2 and n < 200:
            im_size -= int(self.after.shape[0] / 200)  # size reduced image
            self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
            left_p, lower_p, right_p, upper_p = boundary_perc(self.after)
            n += 1

        # STEP 3: Reposition for high flux near the borders. These are bright sources that we want to have more in the center
        peak_y, peak_x = np.argwhere(self.after > self.peak_flux / 2).T / im_size
        n = 0
        while (np.sum(peak_y > 0.75) or np.sum(peak_y < 0.25) or np.sum(peak_x > 0.75) or np.sum(
                peak_x < 0.25)) and n < 100:
            if np.sum(peak_y > 0.75) and not np.sum(peak_y < 0.25):
                self.pix_y += int(self.after.shape[0] / 200)
                self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                peak_y, peak_x = np.argwhere(self.after > self.peak_flux / 2).T / im_size
            elif np.sum(peak_y < 0.25) and not np.sum(peak_y > 0.75):
                self.pix_y -= int(self.after.shape[0] / 200)
                self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                peak_y, peak_x = np.argwhere(self.after > self.peak_flux / 2).T / im_size
            elif np.sum(peak_y > 0.75) and np.sum(peak_y < 0.25) and im_size * step_size < self.initial_box_size * 1.1:
                im_size += int(self.after.shape[0]) / 200  # size increase
                self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                peak_y, peak_x = np.argwhere(self.after > self.peak_flux / 2).T / im_size
            if np.sum(peak_x > 0.75) and not np.sum(peak_x < 0.25):
                self.pix_x += int(self.after.shape[0] / 200)
                self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                peak_y, peak_x = np.argwhere(self.after > self.peak_flux / 2).T / im_size
            elif np.sum(peak_x < 0.25) and not np.sum(peak_x > 0.75):
                self.pix_x -= int(self.after.shape[0] / 200)
                self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                peak_y, peak_x = np.argwhere(self.after > 0.07).T / im_size
            elif np.sum(peak_x < 0.25) and np.sum(peak_x > 0.75) and im_size * step_size < self.initial_box_size * 1.1:
                im_size += int(self.after.shape[0]) / 200  # size increase
                self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                peak_y, peak_x = np.argwhere(self.after > self.peak_flux / 2).T / im_size
            n += 1

        # Step 4: Last border check with resizing
        n, t2 = 0, 0
        left_p, lower_p, right_p, upper_p = boundary_perc(self.after)  # get boundary percentages
        while (left_p > threshold_p or lower_p > threshold_p or right_p > threshold_p or upper_p > threshold_p) \
                and n < 200 and im_size > self.before.shape[0] * 3 / 4 and im_size > self.before.shape[0] * 1 / 2:
            t1 = 0
            if upper_p > threshold_p and int(im_size / 2) + self.pix_y < self.image_data.shape[0]:
                self.pix_y += int(self.after.shape[0] / 300)  # size shifted above
                new_image, _ = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))
                if np.sum(boundary_perc(new_image)) < np.sum(boundary_perc(self.after)):  # check if improvement
                    self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                    left_p, lower_p, right_p, upper_p = boundary_perc(self.after)
                else:  # correct back
                    self.pix_y -= int(self.after.shape[0] / 300)
                    t1 += 1
                n += 1
            elif lower_p > threshold_p and self.pix_y - int(im_size / 2) > 0:
                self.pix_y -= int(self.after.shape[0] / 300)  # size shifted down
                new_image, _ = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))
                if np.sum(boundary_perc(new_image)) < np.sum(boundary_perc(self.after)):  # check if improvement
                    self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                    left_p, lower_p, right_p, upper_p = boundary_perc(self.after)
                else:  # correct back
                    self.pix_y += int(self.after.shape[0] / 300)
                    t1 += 1
                n += 1
            if right_p > threshold_p and self.pix_x - int(im_size / 2) > 0:
                self.pix_x += int(self.after.shape[0] / 300)  # size shifted to the right
                new_image, _ = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))
                if np.sum(boundary_perc(new_image)) < np.sum(boundary_perc(self.after)):  # check if improvement
                    self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                    left_p, lower_p, right_p, upper_p = boundary_perc(self.after)
                else:  # correct back
                    self.pix_x -= int(self.after.shape[0] / 300)
                    t1 += 1
                n += 1
            elif left_p > threshold_p and left_p > right_p and int(im_size / 2) + self.pix_x < self.image_data.shape[0]:
                self.pix_x -= int(self.after.shape[0] / 300)  # size shifted to the left
                new_image, _ = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))
                if np.sum(boundary_perc(new_image)) < np.sum(boundary_perc(self.after)):  # check if improvement
                    self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
                    left_p, lower_p, right_p, upper_p = boundary_perc(self.after)
                else:  # correct back
                    self.pix_x += int(self.after.shape[0] / 300)
                    t1 += 1
                n += 1

        # STEP 5: Last resize step because we prefer smaller images (not smaller than half of the initial size)
        while (
                left_p * 2 < threshold_p and lower_p * 2 < threshold_p and right_p * 2 < threshold_p and upper_p * 2 < threshold_p) \
                and im_size > self.before.shape[0] * 2 / 3 and n < 200:
            im_size -= int(self.after.shape[0] / 400)  # size reduced image
            self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y), (im_size, im_size))  # new image
            left_p, lower_p, right_p, upper_p = boundary_perc(self.after)
            n += 1

        return self

    def save_box(self, box_name: str = 'box.reg'):
        """
        save the box as an .reg file
        """
        position_h = [v.replace('h', ':').replace('m', ':').replace('s', '').replace('d', ':') for v in
                      self.wcs.pixel_to_world(self.pix_x, self.pix_y).to_string('hmsdms').split()]
        arcsec_size = round(np.max(self.wcs.pixel_scale_matrix) * self.after.shape[0] * 3600, 3)
        f = open(f"{folder}/boxes/{box_name}", "a")
        f.write(f'box({position_h[0]},{position_h[1]},{str(arcsec_size)}\",{str(arcsec_size)}\",0)')
        f.close()

        return self

    def source_to_csv(self):
        """Write source to a csv file"""
        with open(f'{folder}/source_file.csv', 'a+') as source_file:
            source_file.write(f'{self.image_number},{self.pix_x},{self.pix_y},{self.flux},{self.after.shape}\n')
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
        """Check if there are other sources in your image."""
        self.df_peaks['same_image'] = (abs(self.pix_x - self.df_peaks.pix_x) <= abs(self.after.shape[0]) / 2).astype(
            int) + \
                                      (abs(self.pix_y - self.df_peaks.pix_y) <= abs(self.after.shape[0]) / 2).astype(
                                          int)
        other_sources_in_image = list(
            self.df_peaks[(self.df_peaks.same_image == 2) & (self.df_peaks.index != self.image_number)].index)

        #TODO: add split boxes if too many high fluxes together

        # print(other_sources_in_image)
        # sources_in_images = self.df_peaks[self.df_peaks.index.isin(other_sources_in_image)].reset_index(keep=False)
        # source_positions = sources_in_images[['pix_y', 'pix_x']].to_numpy()
        # for pos in source_positions:
        #     print([p>100 for p in self.ed_array(pos, source_positions)])

        return other_sources_in_image

    @staticmethod
    def center_cutout(image_data):
        return Cutout2D(data=image_data,
                        position=(int(image_data.shape[0] / 2), int(image_data.shape[0] / 2)),
                        size=(int(image_data.shape[0] / 8), int(image_data.shape[0] / 8))).data

    @property
    def interesting_source(self):
        """Check if source is interesting for selfcal."""

        im_size = int(self.initial_box_size / np.max(self.wcs.pixel_scale_matrix))

        self.initial_pos = (self.pix_x, self.pix_y)
        self.before, self.wcs_cut = self.make_cutout(pos=self.initial_pos,
                                       size=(im_size,
                                             im_size))  # image data after repositioning

        def central_cutout(imdat):
            return Cutout2D(data=imdat,
                            position=(int(imdat.shape[0] / 2), int(imdat.shape[0] / 2)),
                            size=(int(imdat.shape[0] / 8), int(imdat.shape[0] / 8))).data

        outliers = self.before > np.std(self.before)
        outlier_cutout = central_cutout(outliers.astype(int))
        gaussian_image = gaussian_filter(self.before, sigma=3)
        gaussian_outliers = gaussian_image > np.std(gaussian_image)
        gaussian_outliers_cutout = central_cutout(gaussian_outliers.astype(int))
        res_gauss = outlier_cutout - gaussian_outliers_cutout
        res_gauss = np.where(res_gauss < 0, 0, res_gauss)

        perc_outliers_cutout = np.sum(outlier_cutout) / len(outlier_cutout.flatten())
        perc_gaussian_cutout = np.sum(res_gauss) / len(res_gauss.flatten())

        # these threshold values are based on tests
        if perc_outliers_cutout > 0.014 or perc_gaussian_cutout > 0:
            return True
        else:
            return False


if __name__ == '__main__':

    image = SetBoxes(fits_file=args.file, initial_box_size=0.4)

    if args.images:
        os.system(f'rm -rf {folder}/box_images; mkdir {folder}/box_images')  # make folder with box images
    os.system(f'rm -rf {folder}/boxes; mkdir {folder}/boxes')  # make folder with the .reg files

    with open(f'{folder}/source_file.csv', 'w') as source_file:
        source_file.write('id,pix_x,pix_y,flux,size\n')

    sources_done = []
    print(f'We found {len(image.df_peaks)} interesting directions.\n')

    m, r = 0, 0
    for n, p in enumerate(image.df_peaks.to_dict(orient="records")):

        # skip sources that are already displayed in other boxes
        if n in sources_done:
            print(f'Direction {n} already included in box.')
            continue

        # set values for source of interest
        image.pix_x, image.pix_y, image.flux, image.image_number = p['pix_x'], p['pix_y'], p['flux'], n

        # skip source if not interesting to calibrate
        if not image.interesting_source:
            r += 1
            print(f'Direction {n} removed.')
            if args.images:
                fig = plt.figure(figsize=(10, 10))
                plt.subplot(1, 1, 1, projection=image.wcs)
                plt.title(f'REMOVED')
                plt.imshow(image.before, norm=SymLogNorm(linthresh=image.vmin / 10, vmin=image.vmin / 20, vmax=image.vmax),
                           origin='lower',
                           cmap='CMRmap')
                fig.savefig(f'{folder}/box_images/removed_{r}.png')
                continue
        else:
            print(f'Direction {n} in box {m+1}.')

        # reposition box
        image.reposition()

        # we now check if there are in our box sources from our list of peak sources, which we can skip later on
        other_sources = image.other_sources_in_image
        if len(other_sources) == 0 and image.flux < 0.07:  # if flux <0.07 it needs to have other sources in the image, otherwise skip
            continue
        else:
            sources_done += other_sources

        image.source_to_csv()

        # make image with before and after repositioning of our box
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
            fig.savefig(f'{folder}/box_images/box_{m}.png')

        image.save_box(box_name=f'box_{m}.reg')

    print('-------------------------------------------------')
    print(f'Made succesfully {m} boxes.')
    if args.images:
        print(f'See {folder}/box_images for the images of the directions.')
    print(f'See {folder}/boxes for the region files of the boxes that are made.')