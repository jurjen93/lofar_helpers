"""
Use this script to return region files of directions that need a selfcal.

Use make_boxes.py as a standalone script by running on the command line:
python make_boxes.py <FLAGS>
You can use the following flags:

    -f --> followed by the fits file name (and path)
    -l --> followed by the location (path) to store the data
    --no_images --> don't save the images locally
    --ds9 --> interactive mode to validate the box selection in ds9
    -mb --> max number of boxes


The script returns the following:

    directory with .reg region boxes
    directory with box images, to check the quality of the boxes.
"""

__author__ = "Jurjen de Jong (jurjendejong@strw.leidenuniv.nl)"

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

__all__ = ['']

warnings.filterwarnings("ignore")

parser = ArgumentParser()
parser.add_argument('-f', '--file', type=str, help='fitsfile name')
parser.add_argument('-l', '--location', type=str, help='data location folder name', default='.')
parser.add_argument('--no_images', action='store_true', help='store images')
parser.add_argument('--make_DL_food', action='store_true', help='store images for creating the DL model')
parser.add_argument('--ds9', action='store_true', help='open ds9 to interactively check and change the box selection')
parser.add_argument('-ac', '--angular_cutoff', type=float, default=None, help='angular distances higher than this value from the center will be excluded from the box selection')
parser.add_argument('-mb', '--max_boxes', type=int, default=999, help='Set max number of boxes that can be made')
args = parser.parse_args()
print(args)

if sys.version_info.major == 2:
    print('ERROR: This code only works for Python 3. Please switch.\n.....ENDED.....')
    sys.exit()

folder = args.location
if folder[-1] == '/':
    folder = folder[0:-1]

if args.make_DL_food:
    os.system('mkdir -p {LOCATION}/DL_data/numpy'.format(LOCATION=folder))
    os.system('mkdir -p {LOCATION}/DL_data/png'.format(LOCATION=folder))
    if not os.path.isfile("{LOCATION}/DL_data/label_data.csv".format(LOCATION=folder)):
        with open("{LOCATION}/DL_data/label_data.csv".format(LOCATION=folder), 'w') as file:
            writer = csv.writer(file)
            writer.writerow(["Name", "Recalibrate"])
    # try: # Install tkinter/tk for using interactive window mode
    #     import importlib
    #     importlib_found = importlib.util.find_spec("tk") is not None
    #     if not importlib_found:
    #         os.system('pip install --user tk')
    #     try:
    #         import tk
    #         use_tk = True
    #     except ModuleNotFoundError:
    #         use_tk = False
    # except:
    #     print('ERROR: tk cannot be installed')
    #     use_tk = False


#check if folder exists and create if not
folders = folder.split('/')
for i, f in enumerate(folders):
    subpath = '/'.join(folder.split('/')[0:i+1])
    if not os.path.isdir(subpath):
        print(f'Create directory: {subpath}')
        os.system(f'mkdir {subpath}')

if not args.no_images:
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.colors import SymLogNorm
    except ImportError:
        print('Failed to import matplotlib. Check your version.\nNo images will be made.')
        args.no_images = True

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
        self.header = self.wcs.to_header()

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
        try:
            plt.figure(figsize=(10, 10))
            plt.subplot(projection=self.wcs)
            plt.imshow(image_data, norm=SymLogNorm(linthresh=self.vmin/20, vmin=self.vmin/50, vmax=self.vmax), origin='lower', cmap=cmap)
            plt.xlabel('Galactic Longitude')
            plt.ylabel('Galactic Latitude')
            plt.show()
        except:
            print('Error making images with matplotlib. Images will not be made.')
            args.no_images = True

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
    def __init__(self, fits_file: str = None, initial_box_size: float = 0.4, peak_flux=0.07):
        """
        :param fits_file: fits file to find boxes
        :param initial_box_size: initial box size to start with (making this lower, makes it more likely you get a smaller box)
        :param peak_flux: peak fluxes considered for boxes
        """
        self.pix_y = None
        self.pix_x = None
        self.flux = None
        self.image_number = None
        self.initial_box_size = initial_box_size
        super().__init__(fits_file=fits_file)
        self.wcs_cut = self.wcs

        # get the positions of the flux peaks from full image
        self.peak_flux = peak_flux

        # peaks of non-resampled sources
        df_peaks = DataFrame([{'pix_y': c[0], 'pix_x': c[1], 'flux': self.image_data[c[0], c[1]], 'resampled':False}
                                   for c in np.argwhere(self.image_data > peak_flux)]). \
                            sort_values('flux', ascending=False)

        # find from resampled data more sources, which are bright sources over bigger area
        resampled = resample_pixels(self.image_data,
                                    rows=int(self.image_data.shape[0] / 100),
                                    cols=int(self.image_data.shape[1] / 100))

        resample_x, resample_y = resampled.shape
        original_x, original_y = self.image_data.shape
        resample_scale_x = original_x//resample_x
        resample_scale_y = original_y//resample_y

        # peaked resampled sources [flux is resampled values so much higher than non-resampled]
        resampled_df_peaks = DataFrame([{'pix_y': c[0]*resample_scale_x,
                                        'pix_x': c[1]*resample_scale_y,
                                         'flux': resampled[c[0], c[1]]/np.mean([resample_scale_y, resample_scale_x]),
                                         'resampled':True}
                                        for c in np.argwhere(resampled > 15)])

        self.df_peaks = concat([df_peaks, resampled_df_peaks], axis=0).\
            drop_duplicates(subset=['pix_y', 'pix_x'])

        if args.angular_cutoff:
            self.df_peaks['distance_from_center_deg'] = self.df_peaks.apply(lambda x: self.angular_distance((self.header['CRPIX1'], self.header['CRPIX2']),
                               (x['pix_x'], x['pix_y'])).value, axis=1)
            excluded_sources = self.df_peaks[self.df_peaks.distance_from_center_deg > args.angular_cutoff]
            excluded_sources['angular_position'] = excluded_sources.apply(
                lambda x: ';'.join([str(self.degree_to_radian(i)) for i in self.wcs.pixel_to_world(x['pix_x'], x['pix_y']).to_string().split()]), axis=1)
            excluded_sources.to_csv('excluded_sources.csv', index=False)
            self.df_peaks = self.df_peaks[self.df_peaks.distance_from_center_deg<=args.angular_cutoff]

        self.df_peaks = self.df_peaks.reset_index(drop=True)

    # euclidean distance
    @staticmethod
    def ed_array(pos=None, lst=None):
        """Euclidean distance between position and list"""


    def angular_distance(self, p1, p2):
        """Angular distance between two points"""
        c1 = self.wcs.pixel_to_world(p1[0], p1[1])
        c2 = self.wcs.pixel_to_world(p2[0], p2[1])
        return c1.separation(c2)*u.degree/u.degree

    @staticmethod
    def degree_to_radian(inp):
        """degree to radian"""
        return float(inp)/360*np.pi*2

    def reposition(self):
        """Reposition image by looking at the data points near the boundaries."""

        noise_threshold = np.std(self.image_data)*3

        def boundary_sources(image_data, threshold=noise_threshold):
            """Sources within the boundaries of the box"""
            other_source = (image_data > threshold).astype(int)
            boundary_size = int(image_data.shape[0] / 6)
            left = other_source.T[0:boundary_size]  # left boundary
            lower = other_source[0:boundary_size]  # lower boundary
            right = other_source.T[len(other_source) - boundary_size:len(other_source)]  # right boundary
            upper = other_source[len(other_source) - boundary_size:len(other_source)]  # upper boundary
            total = np.sum(left)+np.sum(right)+np.sum(upper)+np.sum(lower)
            return total>0

        def check_position(position1, position2, data):
            """Check if position1 is not more than half the diameter away from position2"""
            dist = np.sqrt(np.square(position1[0]-position2[0])+np.square(position1[1]-position2[1]))
            shape = data.shape
            max_len = np.sqrt(shape[0]**2+shape[1]**2)/4
            return dist < max_len

        self.after = self.before.copy()  # image data before data after repositioning

        #get pixel scale
        pixscale = self.header['CDELT1']
        max_size = abs(int(0.4//pixscale))
        #min size depends on the flux in the image, as low flux needs a bigger box
        if np.sum(self.after) * self.flux < 0.4:
            min_size = abs(int(0.35 // pixscale))
        elif np.sum(self.after) * self.flux < 0.65:
            min_size = abs(int(0.3 // pixscale))
        else:
            min_size = abs(int(0.25 // pixscale))

        im_size = max(self.initial_box_size, min_size*1.05) # starting image size

        start_pos = (self.pix_x, self.pix_y)

        for N in range(3):#looping multiple times

            m = 0
            # Shrink size when noise at boundary
            while boundary_sources(image_data=self.make_cutout((self.pix_x, self.pix_y), (im_size+int(im_size/10), im_size+int(im_size/10)))[0])\
                    and m<100 \
                    and im_size > min_size \
                    and check_position(start_pos, (self.pix_x, self.pix_y), self.make_cutout((self.pix_x, self.pix_y),
                        (im_size-3, im_size-3))[0]):
                im_size -= 3 # size decrease
                self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y),
                                                            (im_size, im_size))
                m+=1

            m = 0
            # Shift to north when noise outside northern boundary
            while not boundary_sources(image_data=self.make_cutout((self.pix_x, self.pix_y + int(im_size / 5)),
                                                              (im_size, im_size))[0]) \
                    and m<100 \
                    and check_position(start_pos, (self.pix_x, self.pix_y+1), self.after):
                self.pix_y += 1 # shift north
                self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y),
                                                            (im_size, im_size))  # new image
                m+=1

            m = 0
            # Shift to south when noise outside southern boundary
            while not boundary_sources(image_data=self.make_cutout((self.pix_x, self.pix_y - int(im_size / 5)),
                                                              (im_size, im_size))[0]) \
                    and m<100\
                    and check_position(start_pos, (self.pix_x, self.pix_y-1), self.after):
                self.pix_y -= 1 # shift south
                self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y),
                                                            (im_size, im_size))  # new image
                m+=1

            m = 0
            # Shift to left when noise outside left boundary
            while not boundary_sources(image_data=self.make_cutout((self.pix_x - int(im_size / 5), self.pix_y),
                                                              (im_size, im_size))[0]) \
                    and m<100\
                    and check_position(start_pos, (self.pix_x-1, self.pix_y), self.after):
                self.pix_x -= 1 # shift to the left
                self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y),
                                                            (im_size, im_size))  # new image
                m+=1

            m = 0
            # Shift to right when noise outside right boundary
            while not boundary_sources(image_data=self.make_cutout((self.pix_x + int(im_size / 5), self.pix_y),
                                                      (im_size, im_size))[0]) \
                    and m<100\
                    and check_position(start_pos, (self.pix_x+1, self.pix_y), self.after):
                self.pix_x += 1 # shift to the right
                self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y),
                                                            (im_size, im_size))  # new image
                m+=1

            m = 0
            # Increase size when noise outside boundary
            while boundary_sources(image_data=self.make_cutout((self.pix_x, self.pix_y), (im_size+int(im_size/10), im_size+int(im_size/10)))[0])\
                    and m<100 \
                    and im_size < max_size \
                    and check_position(start_pos, (self.pix_x, self.pix_y), self.make_cutout((self.pix_x, self.pix_y),
                                                            (im_size+2, im_size+2))[0]):
                im_size += 2 # size increase
                self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y),
                                                            (im_size, im_size))
                m+=1

            m = 0
            # Shrink size when peak flux outside boundary
            while boundary_sources(image_data=self.make_cutout((self.pix_x, self.pix_y),
                                                               (im_size+int(im_size/10), im_size+int(im_size/10)))[0], threshold=0.01)\
                    and m<100 \
                    and im_size > min_size \
                    and check_position(start_pos, (self.pix_x, self.pix_y), self.make_cutout((self.pix_x, self.pix_y),
                        (im_size-2, im_size-2))[0]):
                im_size -= 2 # size decrease
                self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y),
                                                            (im_size, im_size))
                m+=1

            m = 0
            # Shift to north when peak flux outside northern boundary
            while not boundary_sources(image_data=self.make_cutout((self.pix_x, self.pix_y + int(im_size / 5)),
                                                              (im_size, im_size))[0], threshold=0.01) \
                    and m<100 \
                    and check_position(start_pos, (self.pix_x, self.pix_y+1), self.after):
                self.pix_y += 1 # shift north
                self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y),
                                                            (im_size, im_size))  # new image
                m+=1

            m = 0
            # Shift to south when peak flux outside southern boundary
            while not boundary_sources(image_data=self.make_cutout((self.pix_x, self.pix_y - int(im_size / 5)),
                                                              (im_size, im_size))[0], threshold=0.01) \
                    and m<100\
                    and check_position(start_pos, (self.pix_x, self.pix_y-1), self.after):
                self.pix_y -= 1 # shift south
                self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y),
                                                            (im_size, im_size))  # new image
                m+=1

            m = 0
            # Shift to left when noise outside left boundary
            while not boundary_sources(image_data=self.make_cutout((self.pix_x - int(im_size / 5), self.pix_y),
                                                              (im_size, im_size))[0], threshold=0.01) \
                    and m<100\
                    and check_position(start_pos, (self.pix_x-1, self.pix_y), self.after):
                self.pix_x -= 1 # shift to the left
                self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y),
                                                            (im_size, im_size))  # new image
                m+=1

            m = 0
            # Shift to right when noise outside rightern boundary
            while not boundary_sources(image_data=self.make_cutout((self.pix_x + int(im_size / 5), self.pix_y),
                                                      (im_size, im_size))[0], threshold=0.01) \
                    and m<100\
                    and check_position(start_pos, (self.pix_x+1, self.pix_y), self.after):
                self.pix_x += 1 # shift to the right
                self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y),
                                                            (im_size, im_size))  # new image
                m+=1

            m = 0
            # Increase size when noise outside boundary
            while boundary_sources(image_data=self.make_cutout((self.pix_x, self.pix_y),
                                                               (im_size+int(im_size/10), im_size+int(im_size/10)))[0], threshold=0.01)\
                    and m<100 \
                    and im_size < max_size \
                    and check_position(start_pos, (self.pix_x, self.pix_y), self.make_cutout((self.pix_x, self.pix_y),
                                                            (im_size+1, im_size+1))[0]):
                im_size += 1 # size increase
                self.after, self.wcs_cut = self.make_cutout((self.pix_x, self.pix_y),
                                                            (im_size, im_size))
                m+=1

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

    image = SetBoxes(fits_file=args.file)

    if not args.no_images:
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

        # make initial cutout
        image.make_initial_box()

        # make DL data
        if args.make_DL_food:
            np.save("{LOCATION}/DL_data/numpy/{NAME}_box_".format(LOCATION=folder, NAME=args.file.split('/')[-1].split('.fits')[0])+str(m), image.before)
            matplotlib.use("TkAgg")
            plt.ion()
            plt.imshow(image.before,
                       norm=SymLogNorm(linthresh=image.vmin / 10, vmin=image.vmin / 10, vmax=image.vmax / 2),
                       origin='lower',
                       cmap='CMRmap')
            plt.axis('off')
            plt.show()
            plt.savefig("{LOCATION}/DL_data/png/{NAME}_box_".format(LOCATION=folder, NAME=args.file.split('/')[-1].split('.fits')[0])+str(m)+'.png',)

            with open("{LOCATION}/DL_data/label_data.csv".format(LOCATION=folder), 'a+') as file:
                writer = csv.writer(file)
                writer.writerow([args.file.split('/')[-1].split('.fits')[0]+"_box_"+str(m),
                                 int(input("Recalibration? (y/n):")=='y')])
            matplotlib.use('Agg')

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
                            if not args.no_images:
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

        if not args.no_images:
            try:
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
            except:
                print('Error making images with matplotlib. Images will not be made.')
                args.no_images = True

        if replace:
            print(f'Replace box {M+1}.')
            image.save_box(box_name=f'{folder}/boxes/box_{M+1}.reg')
        else:
            print(f'Create box {m}.')
            image.save_box(box_name=f'{folder}/boxes/box_{m}.reg')

        if m == args.max_boxes: # finish if max boxes reached
            break # break for loop

    print('-------------------------------------------------')
    print(f'Made succesfully {m} boxes.')
    if not args.no_images:
        print(f'Images of boxes are in {folder}/box_images.')
    print(f'Region files are in {folder}/boxes.')

    os.system('rm {DATALOC}/source_file.csv && rm {DATALOC}/excluded_sources.csv'.format(DATALOC=folder))


    if args.ds9:
        """
        With the following part you can move ds9 region files
        """
        try:
            from glob import glob
            print('Opening ds9 to verify box selections and make manual changes if needed.'
                  '\nIf you wish to make changes, please save the new full regions file under a new name in '+folder+'/boxes.')
            current_boxes = glob(folder+'/boxes/*')
            os.system("ds9 {FILE} -regions load all '{DATALOC}/boxes/*.reg'".format(FILE=args.file, DATALOC=folder))
            new_box = [b for b in glob(folder+'/boxes/*') if b not in current_boxes]
            if len(new_box)==1:
                os.system('mkdir '+folder+'/boxestemp')
                for nb in new_box:
                    with open(nb, 'r') as f:
                        for line in f:
                            if '{box' in line:
                                g = open(folder+'/boxestemp/'+line.strip().split()[-1].replace('text={','').replace('}','')+'.reg', "a")
                                g.write(line)
                                g.close()
                os.system('rm -rf '+folder+'/boxes && mv '+folder+'/boxestemp '+folder+'/boxes')

            print('Closed ds9.')
        except:
            print("Failing to open ds9 to verify box selection, check if installed and try to run on the commandline"
                  "\nds9 {FILE} -regions load all '{DATALOC}/boxes/*.reg'".format(FILE=args.file, DATALOC=folder))