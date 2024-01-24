"""
This script can be used to select the best self-calibration cycle from facetselfcal.py see:
https://github.com/rvweeren/lofar_facet_selfcal/blob/main/facetselfcal.py
It will return a few plots and a csv with the statistics for each self-calibration cycle.

You only need to run this script in the folder with your facetselfcal output as
python selfcal_quality.py
"""

__author__ = "Jurjen de Jong (jurjendejong@strw.leidenuniv.nl)"

import functools
import re
import tables
from glob import glob
from scipy.stats import circstd, linregress
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from astropy.io import fits
import csv
from skimage.filters.rank import entropy
from skimage.morphology import disk
import pandas as pd
import sys
from cv2 import bilateralFilter

from typing import Union

class SelfcalQuality:
    def __init__(self, folder: str, station: str):
        """
        Determine quality of selfcal from facetselfcal.py

        :param folder: path to directory where selfcal ran
        :param dutch_only: consider dutch stations only for amp/phase stability
        :param remote_only: consider remote stations only for amp/phase stability
        :param international_only: consider international stations only for amp/phase stability
        """

        # selfcal folder
        self.folder = folder

        # merged selfcal h5parms
        self.h5s = [h5 for h5 in glob(f"{self.folder}/*.h5") if 'linearfulljones' not in h5]
        if len(self.h5s) == 0:
            self.h5s = glob(f"{self.folder}/*.h5")
        if len(self.h5s) == 0:
            raise FileNotFoundError("WARNING: No h5 files found")
        # assert len(self.h5s) != 0, "No h5 files found"

        # select all sources
        sources = []
        for h5 in self.h5s:
            matches = re.findall(r'selfcalcyle\d+_(.*?)\.', h5.split('/')[-1])
            assert len(matches) == 1
            sources.append(matches[0])

        self.sources = set(sources)

        assert len(self.sources) > 0, "No sources found"

        # select all fits images
        fitsfiles = sorted(glob(self.folder + "/*MFS-I-image.fits"))
        if len(fitsfiles) == 0 or '000' not in fitsfiles[0]:
            fitsfiles = sorted(glob(self.folder + "/*MFS-image.fits"))
        self.fitsfiles = [f for f in fitsfiles if 'arcsectaper' not in f]
        assert len(self.fitsfiles) != 0, "No fits files found"

        # for phase/amp evolution
        self.station = station

        # output csv
        self.textfile = open(f'selfcal_performance.csv', 'w')
        self.writer = csv.writer(self.textfile)
        self.writer.writerow(['solutions', 'dirty'] + [str(i) for i in range(len(self.fitsfiles))])

    def get_max_min_pix(self):
        """
        Get max/min pixels from images (measure of dynamic range)
        """

        maxp, minp = 0, 0
        for f in self.fitsfiles:
            with fits.open(f) as fts:
                d = fts[0].data
            if d.min() < minp:
                minp = d.min()
            if d.max() > maxp:
                maxp = d.max()
        return maxp, minp

    @staticmethod
    def get_cycle_num(fitsfile: str = None) -> int:
        """
        Parse cycle number

        :param fitsfile: fits file name
        """
        cycle_num = int(float(re.findall(r"selfcalcyle(\d+)", fitsfile.split('/')[-1])[0]))
        assert cycle_num >= 0
        return cycle_num

    @staticmethod
    def make_utf8(inp=None):
        """
        Convert input to utf8 instead of bytes

        :param inp: string input
        """

        try:
            inp = inp.decode('utf8')
            return inp
        except (UnicodeDecodeError, AttributeError):
            return inp

    @staticmethod
    def make_figure(vals1=None, vals2=None, label1=None, label2=None, plotname=None):
        """
        Make figure (with optionally two axis)

        :param vals1: values 1
        :param vals2: values 2
        :param label1: label corresponding to values 1
        :param label2: label corresponding to values 2
        :param plotname: plot name
        """

        fig, ax1 = plt.subplots()

        color = 'tab:red'
        ax1.set_xlabel('cycle')
        ax1.set_ylabel(label1, color=color)
        ax1.plot([i + 1 for i in range(len(vals1))], vals1, color=color)
        ax1.tick_params(axis='y', labelcolor=color)
        # ax1.set_ylim(0, np.pi/2)
        ax1.set_xlim(1, 11)
        ax1.grid(False)
        ax1.grid('off')
        ax1.grid(None)

        if vals2 is not None:

            ax2 = ax1.twinx()

            color = 'tab:blue'
            ax2.set_ylabel(label2, color=color)
            if 'Amp' in label2:
                ax2.plot([i for i in range(len(vals2))][3:], vals2[3:], color=color)
            else:
                ax2.plot([i for i in range(len(vals2))], vals2, color=color)
            ax2.tick_params(axis='y', labelcolor=color)
            # ax2.set_ylim(0, 2)
            ax2.set_xlim(1, 11)
            ax2.grid(False)
            ax2.grid('off')
            ax2.grid(None)

        fig.tight_layout()

        plt.savefig(plotname, dpi=300)

    def linreg_slope(self, values=None):
        """
        Fit linear regression and return slope

        :param values: Values

        :return: linear regression slope
        """

        return linregress(list(range(len(values))), values).slope

    def image_entropy(self, fitsfile: str = None):
        """
        Calculate entropy of image

        :param fitsfile:

        :return: image entropy value
        """

        with fits.open(fitsfile) as f:
            image = f[0].data

        while image.ndim > 2:
            image = image[0]
        image = np.sqrt((image - self.minp) / (self.maxp - self.minp)) * 255
        image = image.astype(np.uint8)
        val = entropy(image, disk(6)).sum()
        print(f"Entropy: {val}")
        return val

    def get_solution_scores(self, h5_1: str, h5_2: str = None):
        """
        Get solution scores

        :param h5_1: solution file 1
        :param h5_2: solution file 2

        :return: phase_score --> circular std phase difference score
                 amp_score --> std amp difference score
        """

        def extract_data(tables_path):
            with tables.open_file(tables_path) as f:
                return (
                    [self.make_utf8(station) for station in f.root.sol000.antenna[:]['name']],
                    self.make_utf8(f.root.sol000.phase000.val.attrs['AXES']).split(','),
                    f.root.sol000.phase000.pol[:],
                    f.root.sol000.phase000.val[:],
                    f.root.sol000.phase000.weight[:],
                    f.root.sol000.amplitude000.val[:],
                )


        def filter_stations(station_names):
            """Generate indices of filtered stations"""
            if self.station == 'debug':
                return list(range(len(station_names)))

            station_codes = (
                ('CS',) if self.station == 'dutch' else
                ('RS',) if self.station == 'remote' else
                ('RS', 'CS', 'ST')  # i.e.: if self.station == 'international'
            )
            return [
                i for i, station_name in enumerate(station_names)
                if any(station_code in station_name for station_code in station_codes)
            ]

        def filter_params(station_indices, axes, *parameters):
            return tuple(
                np.take(param, station_indices, axes)
                for param in parameters
            )

        def weighted_vals(vals, weights):
            return np.nan_to_num(vals) * weights

        station_names1, axes1, phase_pols1, *params1 = extract_data(h5_1)

        antenna_selection = functools.partial(filter_params, filter_stations(station_names1), axes1.index('ant'))
        phase_vals1, phase_weights1, amps1 = antenna_selection(*params1)

        prep_phase_score = weighted_vals(phase_vals1, phase_weights1)
        prep_amp_score = weighted_vals(amps1, phase_weights1)

        if h5_2 is not None:
            _, _, phase_pols2, *params2 = extract_data(h5_2)
            phase_vals2, phase_weights2, amps2 = antenna_selection(*params2)

            min_length = min(len(phase_pols1), len(phase_pols2))
            assert 0 < min_length < 2

            indices = [0] if min_length == 1 else [0, -1]

            prep_phase_score, prep_amp_score, phase_vals2, phase_weights2, amps2 = filter_params(
                indices, axes1.index('pol'), prep_phase_score, prep_amp_score, phase_vals2, phase_weights2, amps2
            )
            prep_phase_score = np.subtract(prep_phase_score, weighted_vals(phase_vals2, phase_weights2))
            prep_amp_score = np.nan_to_num(
                np.divide(prep_amp_score, weighted_vals(amps2, phase_weights2)),
                posinf=0, neginf=0
            )

        phase_score = circstd(prep_phase_score[prep_phase_score != 0], nan_policy='omit')
        amp_score = np.std(prep_amp_score[prep_amp_score != 0])

        return phase_score, amp_score

    def solution_stability(self):
        """
        #TODO: Under development

        Get solution stability scores and make figure

        :return:    bestcycle --> best cycle according to solutions
                    accept --> accept this selfcal
        """

        # loop over sources to get scores
        assert len(self.sources) > 0
        for k, source in enumerate(self.sources):
            sub_h5s = sorted([h5 for h5 in self.h5s if source in h5])

            phase_scores = []
            amp_scores = []
            for m, sub_h5 in enumerate(sub_h5s):
                number = self.get_cycle_num(sub_h5)

                phase_score, amp_score = self.get_solution_scores(sub_h5, sub_h5s[m - 1] if number > 0 else None)

                phase_scores.append(phase_score)
                amp_scores.append(amp_score)

            if k == 0:
                total_phase_scores = [phase_scores]
                total_amp_scores = [amp_scores]

            total_phase_scores = np.append(total_phase_scores, [phase_scores], axis=0)
            total_amp_scores = np.append(total_amp_scores, [amp_scores], axis=0)

        # plot
        plotname = f'selfcal_stability_{self.station}.png'

        finalphase, finalamp = (np.mean(score, axis=0) for score in (total_phase_scores, total_amp_scores))

        self.make_figure(finalphase, finalamp, 'Phase stability', 'Amplitude stability', plotname)

        bestcycle = np.array(finalphase).argmin()

        self.writer.writerow(['phase', np.nan] + list(finalphase))
        self.writer.writerow(['amp', np.nan] + list(finalamp))

        if len(finalphase) > 3:
            phase_decrease, phase_quality, amp_quality = self.linreg_slope(finalphase[:4]), self.linreg_slope(
                finalphase[-3:]), self.linreg_slope(finalamp[-3:])
            print(phase_decrease, phase_quality, amp_quality)
            accept = phase_decrease < 0 and abs(phase_quality) < 0.05 and abs(amp_quality) < 0.05

            return bestcycle, accept
        else:
            return None, False

    @staticmethod
    def get_rms(inp: Union[str, np.ndarray], maskSup: float = 1e-7):
        """
        find the rms of an array, from Cycil Tasse/kMS

        :param inp: fits file name or numpy array
        :param maskSup: mask threshold

        :return: rms --> rms of image
        """

        if isinstance(inp, str):
            with fits.open(inp) as hdul:
                data = hdul[0].data
        else:
            data = inp

        mIn = np.ndarray.flatten(data)
        m = mIn[np.abs(mIn) > maskSup]
        rmsold = np.std(m)
        diff = 1e-1
        cut = 3.
        med = np.median(m)
        for i in range(10):
            ind = np.where(np.abs(m - med) < rmsold * cut)[0]
            rms = np.std(m[ind])
            if np.abs((rms - rmsold) / rmsold) < diff:
                break
            rmsold = rms


        print(f'rms: {rms}')

        return rms  # jy/beam

    @staticmethod
    def get_minmax(inp: Union[str, np.ndarray]):
        """
        Get min/max value

        :param inp: fits file name or numpy array

        :return: minmax --> pixel min/max value
        """
        if isinstance(inp, str):
            with fits.open(inp) as hdul:
                data = hdul[0].data
        else:
            data = inp

        minmax = np.abs(data.min() / data.max())

        print(f"min/max: {minmax}")
        return minmax

    @staticmethod
    def bilateral_filter(fitsfile=None, sigma_x=1, sigma_y=1, sigma_z=1, general_sigma=None):
        """
        Bilateral filter
        See: https://www.projectpro.io/recipes/what-is-bilateral-filtering-opencv

        :param image: image data
        :param sigma: standard deviation that controls the influence of distant pixels

        :return: Bilateral filter output
        """

        with fits.open(fitsfile) as hdul:
            data = hdul[0].data

        if general_sigma is not None:
            sigma_x = sigma_y = sigma_z = int(general_sigma)
        return bilateralFilter(data, sigma_x, sigma_y, sigma_z)

    @staticmethod
    def select_cycle(cycles=None):
        """
        Select best cycle

        :param cycles: rms or minmax cycles

        :return: best cycle
        """

        b, best_cycle = 0, 0
        for n, c in enumerate(cycles[1:]):
            if c > cycles[n - 1]:
                b += 1
            else:
                b = 0
                best_cycle = n + 1
            if b == 2:
                break
        return best_cycle

    def image_stability(self, bilateral_filter: bool = None):
        """
        Determine image stability

        :param bilateral_filter: use bilateral_filter or not

        :return: bestcycle --> best solution cycle
                 accept    --> accept this selfcal
        """

        if bilateral_filter:
            rmss = [self.get_rms(self.bilateral_filter(fts, general_sigma=40)) * 1000 for fts in self.fitsfiles]
            minmaxs = [self.get_minmax(self.bilateral_filter(fts, general_sigma=40)) for fts in self.fitsfiles]
        else:
            rmss = [self.get_rms(fts) * 1000 for fts in self.fitsfiles]
            minmaxs = [self.get_minmax(fts) for fts in self.fitsfiles]

        self.make_figure(rmss, minmaxs, '$RMS (mJy)$', '$|min/max|$', f'image_stability.png')

        self.writer.writerow(['min/max'] + minmaxs + [np.nan])
        self.writer.writerow(['rms'] + rmss + [np.nan])

        # scores
        best_rms_cycle = self.select_cycle(rmss)
        best_minmax_cycle = self.select_cycle(minmaxs)

        # best cycle
        bestcycle = (best_minmax_cycle + best_rms_cycle) // 2

        # getting slopes for selection
        if len(rmss) > 4:
            rms_slope_start, minmax_slope_start = linregress(list(range(len(rmss[:4]))), rmss[:4]).slope, linregress(
                list(range(len(rmss[:4]))),
                np.array(
                    minmaxs[:4])).slope
        else:
            rms_slope_start, minmax_slope_start = 1, 1

        rms_slope, minmax_slope = linregress(list(range(len(rmss))), rmss).slope, linregress(list(range(len(rmss))),
                                                                                             np.array(
                                                                                                 minmaxs)).slope

        # acceptance criteria
        if minmax_slope > 0 and rms_slope > 0:
            accept = False
        elif minmax_slope_start > 0 or rms_slope_start > 0:
            accept = False
        elif len(self.fitsfiles) < 5:
            accept = False
        elif bestcycle + 2 < len(rmss):
            accept = True
        else:
            accept = True

        return bestcycle - 1, accept


def parse_args():
    """
    Command line argument parser

    :return: parsed arguments
    """

    parser = ArgumentParser(description='Determine selfcal quality')
    parser.add_argument('--selfcal_folder', required=True)
    parser.add_argument('--bilateral_filter', action='store_true')
    parser.add_argument('--station', type=str, default='dutch', choices=['dutch', 'remote', 'international', 'debug'])
    return parser.parse_args()


def main():
    """
    Main function
    """
    args = parse_args()

    sq = SelfcalQuality(args.selfcal_folder, args.station)
    if len(sq.h5s) > 0 or len(sq.fitsfiles) > 0:
        if len(sq.h5s) > 0:
            bestcycle_solutions, accept_solutions = sq.solution_stability()

            print(f"Best cycle according to solutions {bestcycle_solutions}")
            print(f"Accept according to solutions {accept_solutions}")

        if len(sq.fitsfiles) > 0:
            bestcycle_image, accept_image = sq.image_stability(bilateral_filter=args.bilateral_filter)

            print(f"Best cycle according to image {bestcycle_image}")
            print(f"Accept according to image {accept_image}")

        sq.textfile.close()
        df = pd.read_csv(f'selfcal_performance.csv').set_index('solutions').T
        # print(df)
        df.to_csv(f'selfcal_performance.csv', index=False)

    else:
        sq.textfile.close()
        print(f"Need more than 1 h5 or fits file")


if __name__ == '__main__':
    main()
