import numpy as np
from shapely.geometry import Polygon, Point, MultiPolygon
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
import sys
from shapely.ops import cascaded_union
from matplotlib.path import Path
from itertools import combinations
from astropy.table import Table
from past.utils import old_div
import argparse
from astropy.visualization.wcsaxes import WCSAxes
import os

warnings.filterwarnings("ignore")


def get_rms(image_data):
    """
    from Cyril Tasse/kMS

    :param image_data: image data array
    :return: rms (noise measure)
    """

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
    return rms


class MeasureSource:

    def __init__(self, fitsfile=None, region_mask=None, rms=None, rms_peak_threshold=3,
                 rms_island_threshold=1, max_n_polygons=None):
        """
        Measure properties from a source
        ----------------------------------------
        :param fitsfile: Fits file
        :param region_mask: ds9 region mask
        :param rms: RMS noise (default: calculated by code)
        :param rms_max_threshold: RMS threshold above which flux is fitted for islands
        :param max_n_polygons: set a maximum number of polygons. For example 1 for compact sources
        """

        self.hdu = fits.open(fitsfile)
        self.image_data = self.hdu[0].data
        while self.image_data.ndim > 2:
            self.image_data = self.image_data[0]
        self.header = self.hdu[0].header
        self.wcs = WCS(self.header, naxis=2)

        self.poly_list = []
        self.max_n_polygons = max_n_polygons

        self.region_mask = region_mask

        if rms is None:
            self.rms = get_rms(self.image_data)
            self.rms_peak_threshold = rms_peak_threshold * self.rms
            self.rms_island_threshold = rms_island_threshold * self.rms
        else:
            self.rms = rms
            self.rms_peak_threshold = rms_peak_threshold * rms
            self.rms_island_threshold = rms_island_threshold * rms

        self.beam_pixels = self.beamarea
        self.peak = 0
        self.peak_flux_pos = [0, 0]
        self.cut_image_data = None
        self.cut_header = None

    @property
    def beamarea(self):
        """
        Calculate beam area in pixels
        -----------------------------
        :param bmaj: beam major axis
        :param bmin: beam minor axis
        :return:
        """

        bmaj = self.hdu[0].header['BMAJ']
        bmin = self.hdu[0].header['BMIN']

        beammaj = bmaj / (2.0 * (2 * np.log(2)) ** 0.5)  # Convert to sigma
        beammin = bmin / (2.0 * (2 * np.log(2)) ** 0.5)  # Convert to sigma
        pixarea = abs(self.hdu[0].header['CDELT1'] * self.hdu[0].header['CDELT2'])

        beamarea = 2 * np.pi * 1.0 * beammaj * beammin  # Note that the volume of a two dimensional gaus$
        beamarea_pix = beamarea / pixarea

        return beamarea_pix

    def _to_pixel(self, ra: float = None, dec: float = None):
        """
        To pixel position from RA and DEC in degrees
        --------------------------------------------
        :param ra: Right ascension (degrees)
        :param dec: Declination (degrees)
        :return: Pixel of position
        """

        from astropy import wcs
        position = coordinates.SkyCoord(
            ra, dec, frame=wcs.utils.wcs_to_celestial_frame(self.wcs).name, unit=(u.degree, u.degree)
        )
        position = np.array(position.to_pixel(self.wcs)).astype(int)

        return position

    def make_cutout(self, pos: tuple = None, size: tuple = (1000, 1000)):
        """
        Make cutout from your image.
        ------------------------------------------------------------
        :param pos: position in pixel size or degrees (RA, DEC)
        :param size: size of your image in pixels
        """

        if type(pos[0]) != int and type(pos[1]) != int:
            pos = self._to_pixel(pos[0], pos[1])
        cutout = Cutout2D(data=self.image_data,
                          position=pos,
                          size=size,
                          wcs=self.wcs,
                          mode="partial")

        self.image_data = cutout.data
        self.header = cutout.wcs.to_header()
        self.wcs = WCS(self.header, naxis=2)

        return self


    def mask_region(self, region_mask):
        """
        Mask region with a ds9 region file
        -----------------------------------
        :param region_mask: ds9 region mask
        """

        r = pyregion.open(region_mask).as_imagecoord(header=self.header)
        mask = np.logical_not(r.get_mask(hdu=self.hdu[0], shape=self.image_data.shape))
        self.image_data *= mask
        return self


    def _get_polylist(self, buff=0, ignore_ra=None, ignore_dec=None):
        """
        Get list with polygons
        --------------------------------------------------------
        :param buff: buffer to enlarge merged polygons in pixels
        :return: list with polygons (source components)
        """

        if self.region_mask is not None:
            self.mask_region(self.region_mask)

        # keep only surface brightness above noise level
        image = np.where((self.image_data > self.rms_island_threshold), self.image_data, 0)

        cs = plt.contour(image, [self.rms_island_threshold], colors='white', linewidths=1)
        cs_list = cs.collections[0].get_paths()

        if ignore_ra is not None and ignore_dec is not None:
            self.sources_to_ignore = self._to_pixel(ignore_ra, ignore_dec).T
            self.sources_to_ignore = np.array(self.sources_to_ignore)
            self.sources_to_ignore = self.sources_to_ignore[np.all((self.sources_to_ignore > 0)
                                                                   & (self.sources_to_ignore < self.image_data.shape[
                0]), axis=1)]
        else:
            self.sources_to_ignore = []

        i = 1
        while len(self.poly_list) == 0 and i < 10:
            if len(self.poly_list) == 0:
                self.poly_list = [Polygon(cs.vertices) for cs in cs_list if (len(cs.vertices) >= 4
                                                             and Polygon(cs.vertices).is_valid == True
                                                             and Polygon(cs.vertices).area > self.beam_pixels / i)]
            i += 1

        self.poly_list = [polygon for polygon in self.poly_list if not any([polygon.contains(Point(p))
                                                                            for p in self.sources_to_ignore])]

        # filter on peak flux
        remove_idx = []
        for k, poly in enumerate(self.poly_list):
            x, y = np.meshgrid(np.arange(len(self.image_data)), np.arange(len(self.image_data)))
            x, y = x.flatten(), y.flatten()
            points = np.vstack((x, y)).T
            grid = Path(poly.exterior.coords).contains_points(points)
            submask = grid.reshape(self.image_data.shape)
            if self.image_data[submask].max() < self.rms_peak_threshold:
                remove_idx.append(k)
        for idx in list(set(remove_idx))[::-1]:
            del self.poly_list[idx]

        if i > 8:
            self.max_n_polygons = 1

        # select max number of polygons based on distance from peak flux
        if self.max_n_polygons is not None:
            if self.peak == 0:
                self.peak = self.peak_flux

            def calc_distance(p1, p2):
                """
                Calculate distance between two points
                :param p1: point 1
                :param p2: point 2
                :return: euclidean distance
                """
                return np.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)

            self.poly_list.sort(key=lambda x: calc_distance(list(x.centroid.coords)[0], self.peak_flux_pos),
                                reverse=False)
            self.poly_list = self.poly_list[0:self.max_n_polygons]

        self.merged_geometry = cascaded_union(self.poly_list)
        if buff > 0:
            self.merged_geometry = self.merged_geometry.buffer(buff)
        return self

    @property
    def _get_polygon_data(self):
        """
        Get image data corresponding to polygons
        ----------------------------------------
        :param image: image data
        :return: mask data outside polygon
        """


        mask = np.zeros(self.image_data.shape)
        for poly in self.poly_list:
            x, y = np.meshgrid(np.arange(len(self.image_data)), np.arange(len(self.image_data)))
            x, y = x.flatten(), y.flatten()
            points = np.vstack((x, y)).T
            grid = Path(poly.exterior.coords).contains_points(points)
            submask = grid.reshape(self.image_data.shape)
            if self.image_data[submask].max() > self.rms_peak_threshold:
                mask += submask
        mask = np.where(mask >= 1, np.ones(mask.shape), 0)
        image = self.image_data * mask
        return image

    @property
    def largest_dist(self):
        """
        Calculate largest distance in convex hull
        -----------------------------------------
        :return: largest distance
        """

        # Function to extract all boundary points from the polygons
        def extract_boundary_points(polygons):
            points = []
            for poly in polygons:
                # Extracting the x, y coordinates of the boundary points
                xs, ys = poly.boundary.xy
                points.extend(zip(xs, ys))
            return points

        def find_max_distance_between_points(points):
            max_distance = 0
            pos = []
            for point1, point2 in combinations(points, 2):
                distance = np.linalg.norm(np.array(point1) - np.array(point2))
                if distance > max_distance:
                    max_distance = distance
                    pos = [point1, point2]
            return max_distance, pos

        # Extract boundary points from all polygons
        boundary_points = extract_boundary_points(self.poly_list)

        # Calculate the maximum distance between boundary points
        max_distance, self.distance_points = find_max_distance_between_points(boundary_points)

        return max_distance * self.header['CDELT2'] * 3600 * u.arcsec

    @property
    def s_code(self):
        """Return s_code similar to pybdsf"""
        if self.peak_flux/self.total_flux.value>0.7 and len(self.poly_list)==1:
            print('S')
            return 'S'
        else:
            return 'M'

    @property
    def total_flux(self):
        """Return total flux"""
        return np.nansum(self._get_convex_hull_data) / self.beam_pixels * u.Jy

    @property
    def _get_convex_hull_data(self):
        """
        :param image: image data
        :return: mask data outside polygon
        """

        polygon = self.merged_geometry.boundary.convex_hull.exterior.coords
        x, y = np.meshgrid(np.arange(len(self.image_data)), np.arange(len(self.image_data)))
        x, y = x.flatten(), y.flatten()
        points = np.vstack((x, y)).T
        grid = Path(polygon).contains_points(points)
        mask = grid.reshape(self.image_data.shape)
        mask = np.where(mask >= 1, np.ones(mask.shape), 0)
        image = self.image_data * mask
        return image

    @property
    def peak_flux(self):
        """
        Peak flux in polygons
        :return: maximum from polygon data
        """

        # if len(self.poly_list) == 0:
        #     self._get_polylist()
        peakflux = np.nanmax(self._get_polygon_data)
        self.peak_flux_pos = np.where(self.image_data == peakflux)
        return peakflux

    @property
    def image_moments(self):
        """
        Calculate the raw image moments up to the second order
        https://en.wikipedia.org/wiki/Image_moment
        """

        image = self._get_polygon_data
        # Create grids of x and y coordinates
        y_indices, x_indices = np.indices(image.shape)

        # Calculate the raw moments
        m00 = np.nansum(image)
        m10 = np.nansum(x_indices * image)
        m01 = np.nansum(y_indices * image)
        # m11 = np.sum(x_indices * y_indices * image)
        # m20 = np.sum(x_indices ** 2 * image)
        # m02 = np.sum(y_indices ** 2 * image)

        cx = m10 / m00
        cy = m01 / m00
        return cx, cy

    @property
    def ra_dec(self):
        """Get ra/dec based on image moment"""

        px, py = self.image_moments
        sky = self.wcs.pixel_to_world(px, py)
        return sky.ra.degree, sky.dec.degree

    def make_plot(self, savefig: str = None):
        """
        Make plot with contour
        ---------------------
        """

        if self.peak == 0:
            self.peak = self.peak_flux

        plt.close()
        imdat = self._get_polygon_data
        fig = plt.figure(figsize=(7, 10))
        plt.subplot(projection=self.wcs)

        plt.contour(imdat, [self.rms_island_threshold], colors='white', linewidths=1)

        if self.merged_geometry.boundary is not None:
            polygon = self.merged_geometry.boundary.convex_hull
            x, y = polygon.exterior.xy

            plt.plot(x, y, label="Merged Polygon", color='lightgreen')

        plt.imshow(self.image_data,
                   norm=PowerNorm(gamma=0.5, vmin=self.rms, vmax=self.rms*9),
                   cmap='RdBu_r')
        plt.scatter(self.peak_flux_pos[1], self.peak_flux_pos[0], color='red', marker='*', zorder=2, s=100)
        WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=self.wcs)
        plt.xlabel('Right Ascension (J2000)', size=15)
        plt.ylabel('Declination (J2000)', size=15)
        # plt.xticks([])
        # plt.yticks([])
        if len(self.sources_to_ignore) > 0:
            plt.scatter(self.sources_to_ignore[:, 0], self.sources_to_ignore[:, 1], color='red')
        im_moments = self.image_moments
        plt.scatter(im_moments[0], im_moments[1], marker='*', color='purple', s=100, zorder=2)
        print(f'Largest distance {self.largest_dist}')
        if len(self.distance_points) > 0:
            plt.plot(np.array(self.distance_points)[:, 0], np.array(self.distance_points)[:, 1], marker='o',
                     linestyle='--', color='darkgreen')
        plt.tight_layout()
        if savefig is not None:
            plt.savefig(savefig, dpi=200)
        else:
            plt.show()
        plt.close()

        return self


def get_region_mask(image, idx):
    regionmask = None
    if 'facet_0' in image:
        if ('0.3' in image and idx == 26) or ('0.6' in image and idx == 15):
            regionmask = 'regionmasks/mask1.reg'
        elif '0.3' in image and idx==27:
            regionmask = 'regionmasks/mask2.reg'
    if 'facet_7' in image:
        if '0.6' in image and idx==2:
            regionmask = 'regionmasks/mask3.reg'
    if 'facet_8' in image:
        if '0.6' in image and idx==24:
            regionmask = 'regionmasks/mask4.reg'
    if 'facet_10' in image:
        if '0.3' in image and idx==14:
            regionmask = 'regionmasks/mask5.reg'
    if 'facet_23' in image:
        if '0.3' in image and idx==31:
            regionmask = 'regionmasks/mask6.reg'
    if 'facet_27' in image:
        if ('0.3' in image and idx==32) or ('0.6' in image and idx==15):
            regionmask = 'regionmasks/mask7.reg'


    return regionmask


def get_source_information(table, image, makeplot, debug_idx=None):
    """
    Get information of source
    :param table: Astropy table
    :param image: image
    :param makeplot: make plot
    :return:
    """

    T = Table.read(table)
    RA, DEC = T["RA"], T['DEC']

    T = T[(T['Peak_flux'] > 5 * T['Isl_rms']) & (T['S_Code'] == 'M')]

    peak_flux = []
    total_flux = []
    largest_dist = []
    s_code = []
    ra_decs = []
    source_ids = []


    for idx in range(len(T)):

        try:

            if debug_idx is not None:
                if int(idx)!=int(debug_idx):
                    continue
                elif int(idx)>int(debug_idx):
                    break

            ra, dec = T[idx]['RA', 'DEC']
            RA = RA[RA != ra]
            DEC = DEC[DEC != dec]

            imsize = 100
            cont = True
            ignoreradec = True

            peakthresh, islandthresh = 5, 1.5
            trys = 0
            while cont \
                and imsize < 1500 \
                and trys < 100:

                regionmask = get_region_mask(image, idx)
                S = MeasureSource(fitsfile=image, rms=T[idx]['Isl_rms'], rms_peak_threshold=peakthresh,
                                  rms_island_threshold=islandthresh, region_mask=regionmask)
                S.make_cutout((ra, dec), (imsize, imsize))
                peakflux = S.peak_flux
                if ignoreradec:
                    S._get_polylist(buff=min(peakflux/S.rms_island_threshold, 10), ignore_ra=RA, ignore_dec=DEC)
                else:
                    S._get_polylist(buff=min(peakflux/S.rms_island_threshold, 10))
                if S.merged_geometry.boundary is None:
                    if islandthresh < 4:
                        islandthresh += 0.2
                        islandthresh = min(4., max(1.5, islandthresh))
                        imsize *= 1.1
                    else:
                        imsize *= 1.05
                        islandthresh -= 0.2
                        islandthresh = min(4., max(1.5, islandthresh))
                        ignoreradec = False
                else:
                    x, y = S.merged_geometry.boundary.convex_hull.boundary.coords.xy
                    if (np.max(x) * 2.5 < imsize and np.max(y) * 2.5 < imsize and np.min(x) > 0 and np.min(y) > 0):
                        cont = False
                    elif (np.max(x) > imsize/1.5 or np.max(y) > imsize/1.5 or np.min(x) < 0 or np.min(y) < 0)\
                            or S.total_flux.value/T[idx]['Total_flux'] > 10:
                        imsize *= 1.1
                    elif S.total_flux.value/T[idx]['Total_flux'] < 10 and (imsize > 100 or trys<10):
                        imsize /= 1.1
                        islandthresh -= 0.2
                        islandthresh = max(1.5, islandthresh)
                        imsize = max(100, imsize)
                    else:
                        cont = False
                try:
                    if S.total_flux.value * 5 < T[idx]['Total_flux']:
                        imsize *= 1.05
                        islandthresh -= 0.2
                        islandthresh = max(1.5, islandthresh)
                    elif S.total_flux.value > T[idx]['Total_flux'] * 5:
                        islandthresh += 0.2
                        islandthresh = min(4, islandthresh)
                        imsize /= 1.05
                        imsize = max(100, imsize)
                except AttributeError:
                    imsize *= 1.05

                trys += 1


            if makeplot:
                if S.peak_flux/T[idx]['Peak_flux'] > 2 \
                    or T[idx]['Peak_flux']/S.peak_flux > 2 \
                    or S.total_flux.value/T[idx]['Total_flux'] > 2.5 \
                    or T[idx]['Total_flux']/S.total_flux.value > 2.5:
                    if debug_idx is not None:
                        imname = 'test.png'
                    else:
                        imname = f'issues/test_resolved_{idx}_{image.split("/")[-1].replace(".fits", "")}_{T["Source_id"][idx]}_{table.split("/")[-2]}.png'
                    print(imname)
                    S.make_plot(
                        savefig=imname)
                    accept = False
                else:
                    if debug_idx is not None:
                        imname = 'test.png'
                    else:
                        imname = f'ok/test_resolved_{idx}_{image.split("/")[-1].replace(".fits", "")}_{T["Source_id"][idx]}_{table.split("/")[-2]}.png'
                    print(imname)
                    S.make_plot(savefig=imname)
                    accept = True

            print('Peak flux image:', str(S.peak_flux), 'Jy/beam')
            print('Peak flux table:', str(T[idx]['Peak_flux']), 'Jy/beam')

            print('Total flux image:', str(S.total_flux))
            print('Total flux table:', str(T[idx]['Total_flux']), 'Jy')

            if accept:
                peak_flux.append(S.peak_flux)
                total_flux.append(S.total_flux)
                largest_dist.append(S.largest_dist)
                s_code.append(S.s_code)
                ra_decs.append(S.ra_dec)
                source_ids.append(T[idx]['Source_id'])
        except:
            pass

    return source_ids, peak_flux, total_flux, largest_dist, s_code, ra_decs

def get_table_index(t, source_id):
    return int(np.argwhere(t['Source_id'] == source_id).squeeze())


def parse_args():
    """
    Parse input arguments
    """

    parser = argparse.ArgumentParser(description='Find flux density, peak flux, and size for resolved sources based on pybdsf table and image')
    parser.add_argument('--table', type=str, help='astropy table')
    parser.add_argument('--image', type=str, help='fits image')
    parser.add_argument('--debug_idx', type=int)
    return parser.parse_args()


def main():
    """
    Main function
    """

    args = parse_args()
    ids, peakflux, totalflux, largestidst, scode, radecs = get_source_information(args.table, args.image, True, args.debug_idx)
    T = Table.read(args.table)
    source_ids = [get_table_index(T, id) for id in ids]
    for n, id in enumerate(source_ids):
        T[id]['Peak_flux'] = peakflux[n]
        T[id]['Total_flux'] = totalflux[n].value
        T[id]['RA'] = float(radecs[n][0])
        T[id]['DEC'] = float(radecs[n][1])
        T[id]['S_Code'] = scode[n]

    os.system('mkdir -p outcat')
    T['RA'] %= 360
    T.write('outcat/'+args.table.replace('_final.fits', '_final_updated.fits').split('/')[-1], format='fits', overwrite=True)



if __name__ == '__main__':
    main()
