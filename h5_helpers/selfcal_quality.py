"""
This script is a supporting script to quantify the quality of the self-calibration output from facetselfcal.py from
https://github.com/rvweeren/lofar_facet_selfcal/blob/main/facetselfcal.py
It will return a few plots and a csv with the quality improvements over self-calibration cycle. This might eventually
be useful as a stopping criteria metric.

You only need to run this script in the folder with your facetselfcal output as
python selfcal_quality.py
"""

__author__ = "Jurjen de Jong (jurjendejong@strw.leidenuniv.nl)"

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


class SelfcalQuality:
    def __init__(self, folder: str = None, remote_only: bool = False, international_only: bool = False,
                 dutch_only: bool = False):
        """
        Determine quality of selfcal from facetselfcal.py

        :param folder: path to where selfcal ran
        :param dutch_only: consider dutch stations only for amp/phase stability
        :param remote_only: consider remote stations only for amp/phase stability
        :param international_only: consider international stations only for amp/phase stability
        """

        # selfcal folder
        self.folder = folder

        # merged selfcal h5parms
        self.h5s = [h5 for h5 in glob(f"{self.folder}/merged_selfcalcyle*.h5") if 'linearfulljones' not in h5]
        if len(self.h5s) == 0:
            self.h5s = glob(f"{self.folder}/merged_selfcalcyle*.h5")
        if len(self.h5s) == 0:
            print("WARNING: No h5 files found")
        # assert len(self.h5s) != 0, "No h5 files found"

        # select all sources
        regex = "merged_selfcalcyle\d{3}\_"
        self.sources = set([re.sub(regex, '', h.split('/')[-1]).replace('.ms.copy.phaseup.h5', '') for h in self.h5s])

        # select all fits images
        fitsfiles = sorted(glob(self.folder + "/*MFS-I-image.fits"))
        if len(fitsfiles) == 0 or '000' not in fitsfiles[0]:
            fitsfiles = sorted(glob(self.folder + "/*MFS-image.fits"))
        self.fitsfiles = [f for f in fitsfiles if 'arcsectaper' not in f]
        assert len(self.fitsfiles) != 0, "No fits files found"

        # for phase/amp evolution
        self.remote_only = remote_only
        self.international_only = international_only
        self.dutch_only = dutch_only

        self.textfile = open(f'selfcal_performance.csv', 'w')
        self.writer = csv.writer(self.textfile)
        self.writer.writerow(['solutions', 'dirty'] + [str(i) for i in range(len(self.fitsfiles))])

    def get_max_min_pix(self):
        """
        Get max/min pixels from images
        """

        maxp, minp = 0, 0
        for f in self.fitsfiles:
            fts = fits.open(f)
            d = fts[0].data
            if d.min() < minp:
                minp = d.min()
            if d.max() > maxp:
                maxp = d.max()
        return maxp, minp

    @staticmethod
    def get_cycle_num(fitsfile: str = None):
        """
        Parse cycle number

        :param fitsfile: fits file name
        """
        return int(float(re.findall("\d{3}", fitsfile.split('/')[-1])[0]))

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

        :param phase_scores:
        :param amp_scores:
        :param plotname:
        :return:
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
        file = fits.open(fitsfile)
        image = file[0].data
        file.close()
        while image.ndim > 2:
            image = image[0]
        image = np.sqrt((image - self.minp) / (self.maxp - self.minp)) * 255
        image = image.astype(np.uint8)
        val = entropy(image, disk(6)).sum()
        print(f"Entropy: {val}")
        return val

    @staticmethod
    def euclidean_distance(l1=None, l2=None):
        """
        Take euclidean distance

        :return: euclidean distance
        """
        return np.sqrt(np.sum(np.power(np.subtract(l1, l2), 2)))

    def get_solution_scores(self, h5_1: str = None, h5_2: str = None):
        """
        Get solution scores

        :param h5_1: solution file 1
        :param h5_2: solution file 2

        :return: phasescore --> circular std phase difference score
                 ampscore --> std amp difference score
        """

        # PHASE VALUES
        H = tables.open_file(h5_1)
        axes = self.make_utf8(H.root.sol000.phase000.val.attrs['AXES']).split(',')
        vals1 = H.root.sol000.phase000.val[:]
        weights1 = H.root.sol000.phase000.weight[:]
        pols1 = H.root.sol000.phase000.pol[:]

        if h5_2 is not None:
            F = tables.open_file(h5_2)
            vals2 = F.root.sol000.phase000.val[:]
            weights2 = F.root.sol000.phase000.weight[:]
            pols2 = F.root.sol000.phase000.pol[:]

        if self.dutch_only:
            stations = [i for i, station in enumerate(H.root.sol000.antenna[:]['name']) if
                        ('CS' in self.make_utf8(station))]
            vals1 = np.take(vals1, stations, axis=axes.index('ant'))
            weights1 = np.take(weights1, stations, axis=axes.index('ant'))
            if h5_2 is not None:
                vals2 = np.take(vals2, stations, axis=axes.index('ant'))
                weights2 = np.take(weights2, stations, axis=axes.index('ant'))
        elif self.remote_only:
            stations = [i for i, station in enumerate(H.root.sol000.antenna[:]['name']) if
                        ('RS' in self.make_utf8(station))]
            vals1 = np.take(vals1, stations, axis=axes.index('ant'))
            weights1 = np.take(weights1, stations, axis=axes.index('ant'))
            if h5_2 is not None:
                vals2 = np.take(vals2, stations, axis=axes.index('ant'))
                weights2 = np.take(weights2, stations, axis=axes.index('ant'))
        elif self.international_only:
            stations = [i for i, station in enumerate(H.root.sol000.antenna[:]['name']) if
                        not ('RS' in self.make_utf8(station)
                             or 'CS' in self.make_utf8(station)
                             or 'ST' in self.make_utf8(station))]
            vals1 = np.take(vals1, stations, axis=axes.index('ant'))
            weights1 = np.take(weights1, stations, axis=axes.index('ant'))
            if h5_2 is not None:
                vals2 = np.take(vals2, stations, axis=axes.index('ant'))
                weights2 = np.take(weights2, stations, axis=axes.index('ant'))

        # take circular std from difference of previous and current selfcal cycle
        if h5_2 is not None:
            if len(pols1) != len(pols2):
                if min(len(pols1), len(pols2)) == 1:
                    vals1 = np.take(vals1, [0], axis=axes.index('pol'))
                    vals2 = np.take(vals2, [0], axis=axes.index('pol'))
                    weights1 = np.take(weights1, [0], axis=axes.index('pol'))
                    weights2 = np.take(weights2, [0], axis=axes.index('pol'))
                elif min(len(pols1), len(pols2)) == 2:
                    vals1 = np.take(vals1, [0, -1], axis=axes.index('pol'))
                    vals2 = np.take(vals2, [0, -1], axis=axes.index('pol'))
                    weights1 = np.take(weights1, [0, -1], axis=axes.index('pol'))
                    weights2 = np.take(weights2, [0, -1], axis=axes.index('pol'))
                else:
                    sys.exit("ERROR: SHOULD NOT END UP HERE")
            prepphasescore = np.subtract(np.nan_to_num(vals1) * weights1, np.nan_to_num(vals2) * weights2)
        else:
            prepphasescore = np.nan_to_num(vals1) * weights1
        phasescore = circstd(prepphasescore[prepphasescore != 0], nan_policy='omit')

        # AMP VALUES
        vals1 = H.root.sol000.amplitude000.val[:]
        if h5_2 is not None:
            vals2 = F.root.sol000.amplitude000.val[:]

        if self.remote_only or self.dutch_only or self.international_only:
            vals1 = np.take(vals1, stations, axis=axes.index('ant'))
            if h5_2 is not None:
                vals2 = np.take(vals2, stations, axis=axes.index('ant'))

        # take std from ratio of previous and current selfcal cycle
        if h5_2 is not None:
            if len(pols1) != len(pols2):
                if min(len(pols1), len(pols2)) == 1:
                    vals1 = np.take(vals1, [0], axis=axes.index('pol'))
                    vals2 = np.take(vals2, [0], axis=axes.index('pol'))
                    weights1 = np.take(weights1, [0], axis=axes.index('pol'))
                    weights2 = np.take(weights2, [0], axis=axes.index('pol'))
                elif min(len(pols1), len(pols2)) == 2:
                    vals1 = np.take(vals1, [0, -1], axis=axes.index('pol'))
                    vals2 = np.take(vals2, [0, -1], axis=axes.index('pol'))
                    weights1 = np.take(weights1, [0, -1], axis=axes.index('pol'))
                    weights2 = np.take(weights2, [0, -1], axis=axes.index('pol'))
                else:
                    sys.exit("ERROR: SHOULD NOT END UP HERE")
            prepampscore = np.nan_to_num(
                np.divide(np.nan_to_num(vals1) * weights1, np.nan_to_num(vals2) * weights2),
                posinf=0, neginf=0)
        else:
            prepampscore = np.nan_to_num(vals1) * weights1
        ampscore = np.std(prepampscore[prepampscore != 0])

        H.close()
        if h5_2 is not None:
            F.close()

        return phasescore, ampscore

    def solution_stability(self):
        """
        Get solution stability scores and make figure

        :return:    bestcycle --> best cycle according to solutions
                    accept --> accept this selfcal
        """

        # loop over sources to get scores
        for k, source in enumerate(self.sources):
            sub_h5s = sorted([h5 for h5 in self.h5s if source in h5])
            phase_scores = []
            amp_scores = []
            for m, sub_h5 in enumerate(sub_h5s):
                number = self.get_cycle_num(sub_h5)
                # print(sub_h5, sub_h5s[m - 1])
                if number > 0:
                    phasescore, ampscore = self.get_solution_scores(sub_h5, sub_h5s[m - 1])
                elif number == 0:
                    phasescore, ampscore = self.get_solution_scores(sub_h5, None)
                phase_scores.append(phasescore)
                amp_scores.append(ampscore)
            if k == 0:
                total_phase_scores = [phase_scores]
                total_amp_scores = [amp_scores]

            total_phase_scores = np.append(total_phase_scores, [phase_scores], axis=0)
            total_amp_scores = np.append(total_amp_scores, [amp_scores], axis=0)

        # plot
        if self.dutch_only:
            plotname = f'selfcal_stability_dutch.png'
        elif self.remote_only:
            plotname = f'selfcal_stability_remote.png'
        elif self.international_only:
            plotname = f'selfcal_stability_international.png'
        else:
            plotname = f'selfcal_stability.png'
        finalphase = np.mean(total_phase_scores, axis=0)
        finalamp = np.mean(total_amp_scores, axis=0)

        self.make_figure(finalphase, finalamp, 'Phase stability', 'Amplitude stability', plotname)

        bestcycle = np.array(finalphase).argmin()

        self.writer.writerow(['phase', np.nan] + list(finalphase))
        self.writer.writerow(['amp', np.nan] + list(finalamp))

        if len(finalphase) > 3:
            phase_decrease, phase_quality, amp_quality = self.linreg_slope(finalphase[:4]), self.linreg_slope(
                finalphase[-3:]), self.linreg_slope(finalamp[-3:])
            print(phase_decrease, phase_quality, amp_quality)
            if phase_decrease < 0 and abs(phase_quality) < 0.05 and abs(amp_quality) < 0.05:
                accept = True
            else:
                accept = False
            return bestcycle, accept
        else:
            return None, False

    @staticmethod
    def get_rms(fitsfile: str = None, maskSup: float = 1e-7):
        """
        find the rms of an array, from Cycil Tasse/kMS

        :param fitsfile: fits file name
        :param maskSup: mask threshold

        :return: rms --> rms of image
        """
        hdul = fits.open(fitsfile)
        mIn = np.ndarray.flatten(hdul[0].data)
        m = mIn[np.abs(mIn) > maskSup]
        rmsold = np.std(m)
        diff = 1e-1
        cut = 3.
        med = np.median(m)
        for i in range(10):
            ind = np.where(np.abs(m - med) < rmsold * cut)[0]
            rms = np.std(m[ind])
            if np.abs((rms - rmsold) / rmsold) < diff: break
            rmsold = rms
        hdul.close()

        print(f'rms: {rms}')

        return rms  # jy/beam

    @staticmethod
    def get_minmax(fitsfile: str = None):
        """
        Get min/max value

        :param fitsfile: fits file name

        :return: minmax --> pixel min/max value
        """
        hdul = fits.open(fitsfile)
        data = hdul[0].data
        minmax = np.abs(data.min() / data.max())
        hdul.close()
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

        hdul = fits.open(fitsfile)
        data = hdul[0].data
        hdul.close()

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
            if c>cycles[n-1]:
                b+=1
            else:
                b = 0
                best_cycle = n+1
            if b==2:
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

        # getting slopes
        if len(rmss)>4:
            rms_slope_start, minmax_slope_start = linregress(list(range(len(rmss[:4]))), rmss[:4]).slope, linregress(list(range(len(rmss[:4]))),
                                                                                                 np.array(
                                                                                                     minmaxs[:4])).slope
        else:
            rms_slope_start, minmax_slope_start = 1, 1

        rms_slope, minmax_slope = linregress(list(range(len(rmss))), rmss).slope, linregress(list(range(len(rmss))),
                                                                                             np.array(
                                                                                                 minmaxs)).slope

        # accept direction or not
        if minmax_slope > 0 and rms_slope > 0:
            accept = False
        elif minmax_slope_start > 0 or rms_slope_start > 0:
            accept = False
        elif len(self.fitsfiles) < 5:
            accept = False
        elif bestcycle+2 < len(rmss):
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
    parser.add_argument('--selfcal_folder', default='.')
    parser.add_argument('--bilateral_filter', action='store_true')
    parser.add_argument('--dutch_only', action='store_true', help='Only Dutch stations are considered', default=None)
    parser.add_argument('--remote_only', action='store_true', help='Only remote stations are considered', default=None)
    parser.add_argument('--international_only', action='store_true', help='Only international stations are considered',
                        default=None)
    return parser.parse_args()


def main():
    """
    Main function
    """
    args = parse_args()

    sq = SelfcalQuality(args.selfcal_folder, args.remote_only, args.international_only)
    if len(sq.h5s) > 1 and len(sq.fitsfiles) > 1:
        bestcycle_solutions, accept_solutions = sq.solution_stability()
        bestcycle_image, accept_image = sq.image_stability(bilateral_filter=args.bilateral_filter)
        sq.textfile.close()
        df = pd.read_csv(f'selfcal_performance.csv').set_index('solutions').T
        print(df)
        df.to_csv(f'selfcal_performance.csv', index=False)

        print(f"Best cycle according to solutions {bestcycle_solutions}")
        print(f"Accept according to solutions {accept_solutions}")
        print(f"Best cycle according to image {bestcycle_image}")
        print(f"Accept according to image {accept_image}")
    else:
        sq.textfile.close()
        print(f"Need more than 1 h5 or fits file")


if __name__ == '__main__':
    main()