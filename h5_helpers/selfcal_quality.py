"""
This script is a supporting script to quantify the quality of the self-calibration output from facetselfcal.py from
https://github.com/rvweeren/lofar_facet_selfcal/blob/main/facetselfcal.py

It will return plots quantifying the change of solutions and images over self-calibration cycles.
"""

__author__ = "Jurjen de Jong (jurjendejong@strw.leidenuniv.nl)"

import re
import tables
from glob import glob
from scipy.stats import circstd, linregress, circmean
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from astropy.io import fits
import csv
from skimage.filters.rank import entropy
from skimage.morphology import disk
import pandas as pd


class SelfcalQuality:
    def __init__(self, folder: str = None, remote_only: bool = False, international_only: bool = False):
        """
        Determine quality of selfcal from facetselfcal.py

        :param folder: path to where selfcal ran
        :param remote_only: consider remote stations only for amp/phase stability
        :param international_only: consider international stations only for amp/phase stability
        """
        # selfcal folder
        self.folder = folder
        # merged selfcal h5parms
        self.h5s = glob(f"{self.folder}/merged_selfcalcyle*.h5")
        assert len(self.h5s) != 0, "No h5 files found"

        # select all sources
        regex = "merged_selfcalcyle\d{3}\_"
        self.sources = set([re.sub(regex, '', h.split('/')[-1]).replace('.ms.copy.phaseup.h5', '') for h in self.h5s])
        self.sourcename = list(self.sources)[0].split('_')[-1]

        # select all fits images
        fitsfiles = sorted(glob(self.folder + "/*MFS-I-image.fits"))
        if len(fitsfiles) == 0:
            fitsfiles = sorted(glob(self.folder + "/*MFS-image.fits"))
        self.fitsfiles = [f for f in fitsfiles if 'arcsectaper' not in f]

        # select all fits model images
        modelfiles = sorted(glob(self.folder + "/*MFS-I-model.fits"))
        if len(modelfiles) == 0:
            modelfiles = sorted(glob(self.folder + "/*MFS-model.fits"))
        self.modelfiles = [f for f in modelfiles if 'arcsectaper' not in f]

        # for phase/amp evolution
        self.remote_only = remote_only
        self.international_only = international_only

        self.textfile = open('selfcal_performance.csv', 'w')
        self.writer = csv.writer(self.textfile)
        self.writer.writerow(['solutions', 'dirty'] + [str(i) for i in range(len(self.fitsfiles))])

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

    @staticmethod
    def image_entropy(fitsfile: str = None):
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
        val = entropy(image, disk(10)).sum() / image.max()
        print(val)
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

        if h5_2 is not None:
            F = tables.open_file(h5_2)
            vals2 = F.root.sol000.phase000.val[:]
            weights2 = F.root.sol000.phase000.weight[:]

        if self.remote_only and not self.international_only:
            stations = [i for i, station in enumerate(H.root.sol000.antenna[:]['name']) if
                        ('RS' in self.make_utf8(station))]
            vals1 = np.take(vals1, stations, axis=axes.index('ant'))
            weights1 = np.take(weights1, stations, axis=axes.index('ant'))
            if h5_2 is not None:
                vals2 = np.take(vals2, stations, axis=axes.index('ant'))
                weights2 = np.take(weights2, stations, axis=axes.index('ant'))

        elif self.international_only and not self.remote_only:
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
            prepphasescore = np.subtract(np.nan_to_num(vals1) * weights1, np.nan_to_num(vals2) * weights2)
        else:
            prepphasescore = np.nan_to_num(vals1) * weights1
        phasescore = circstd(prepphasescore[prepphasescore != 0], nan_policy='omit')

        # AMP VALUES
        vals1 = H.root.sol000.amplitude000.val[:]
        if h5_2 is not None:
            vals2 = F.root.sol000.amplitude000.val[:]

        if self.remote_only and not self.international_only:
            vals1 = np.take(vals1, stations, axis=axes.index('ant'))
            if h5_2 is not None:
                vals2 = np.take(vals2, stations, axis=axes.index('ant'))

        elif self.international_only and not self.remote_only:
            vals1 = np.take(vals1, stations, axis=axes.index('ant'))
            if h5_2 is not None:
                vals2 = np.take(vals2, stations, axis=axes.index('ant'))

        # take std from ratio of previous and current selfcal cycle
        if h5_2 is not None:
            prepampscore = np.nan_to_num(np.divide(np.nan_to_num(vals1) * weights1, np.nan_to_num(vals2) * weights2),
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
            print(source)
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
        if self.remote_only:
            plotname = f'selfcal_stability_remote_{self.sourcename}.png'
        elif self.international_only:
            plotname = f'selfcal_stability_international_{self.sourcename}.png'
        else:
            plotname = f'selfcal_stability_{self.sourcename}.png'
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
            return None, None

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

    def image_stability(self):
        """
        Determine image stability

        :return: bestcycle --> best solution cycle
                 accept    --> accept this selfcal
        """

        # cycles = [self.get_cycle_num(fts) for fts in self.fitsfiles]
        rmss = [self.get_rms(fts) * 1000 for fts in self.fitsfiles]
        minmaxs = [self.get_minmax(fts) for fts in self.fitsfiles]
        entropy_images = [self.image_entropy(fts) for fts in self.fitsfiles]
        entropy_models = [self.image_entropy(fts) for fts in self.modelfiles]

        self.make_figure(rmss, minmaxs, '$RMS (mJy)$', '$|min/max|$', f'image_stability_{self.sourcename}.png')
        self.make_figure(vals1=entropy_images, vals2=entropy_models, label1='Entropy image', label2='Entropy model',
                         plotname=f'entropy_{self.sourcename}.png')

        self.writer.writerow(['min/max'] + minmaxs + [np.nan])
        self.writer.writerow(['rms'] + rmss + [np.nan])
        self.writer.writerow(['entropy_image'] + entropy_images + [np.nan])
        self.writer.writerow(['entropy_model'] + entropy_models + [np.nan])

        # SCORING
        best_rms_cycle, best_minmax_cycle, best_entropy_cycle = (np.array(rmss[1:]).argmin() + 1,
                                                                 np.array(minmaxs[1:]).argmin() + 1,
                                                                 min(np.array(entropy_images).argmin(),
                                                                     np.array(entropy_models).argmin()) + 1)
        # using maxmin instead of minmax due to easier slope value to work with
        rms_slope, maxmin_slope = linregress(list(range(len(rmss))), rmss).slope, linregress(list(range(len(rmss))),
                                                                                             1 / np.array(
                                                                                                 minmaxs)).slope

        # accept direction or not
        if maxmin_slope > 10:
            accept = True
        elif maxmin_slope < 0 and rms_slope > 0:
            accept = False
        elif maxmin_slope < 1.5:
            accept = False
        elif maxmin_slope < 2 and rms_slope >= 0:
            accept = False
        else:
            accept = True

        # choose best cycle
        if best_rms_cycle == best_minmax_cycle:
            bestcycle = best_rms_cycle
        elif rmss[best_minmax_cycle] <= 1.1:
            bestcycle = best_minmax_cycle
        elif minmaxs[best_rms_cycle] / minmaxs[0] <= 1:
            bestcycle = best_rms_cycle
        else:
            bestcycle = int(round(np.mean([best_minmax_cycle, best_rms_cycle, best_entropy_cycle]), 0))

        return bestcycle - 1, accept


def parse_args():
    """
    Command line argument parser

    :return: parsed arguments
    """
    parser = ArgumentParser(description='Determine selfcal quality')
    parser.add_argument('--selfcal_folder', default='.')
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
    bestcycle_solutions, accept_solutions = sq.solution_stability()
    bestcycle_image, accept_image = sq.image_stability()
    sq.textfile.close()
    df = pd.read_csv('selfcal_performance.csv').set_index('solutions').T
    print(df)
    df.to_csv(f'selfcal_performance_{sq.sourcename}.csv', index=False)

    print(f"Best cycle according to solutions {bestcycle_solutions} (SCORES IN TESTING PHASE)")
    print(f"Accept according to solutions {accept_solutions} (SCORES IN TESTING PHASE)")
    print(f"Best cycle according to image {bestcycle_image} (SCORES IN TESTING PHASE)")
    print(f"Accept according to image {accept_image} (SCORES IN TESTING PHASE)")


if __name__ == '__main__':
    main()
