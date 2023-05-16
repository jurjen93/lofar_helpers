import tables
import sys
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import normal
from scipy.stats import circstd

# pretty plots
# try:
#     import scienceplots
#     plt.style.use(['science','ieee'])
# except:
#     pass

limit = np.pi

def make_utf8(inp):
    """
    Convert input to utf8 instead of bytes

    :param inp: string input
    :return: input in utf-8 format
    """

    try:
        inp = inp.decode('utf8')
        return inp
    except (UnicodeDecodeError, AttributeError):
        return inp

class GetSolint:
    def __init__(self, h5, optimal_score=0.5, ref_solint=10, station=None):
        """
        Get a score based on the phase difference between XX and YY. This reflects the noise in the observation.
        From this score we can determine an optimal solution interval, by fitting a wrapped normal distribution.

        See:
        - https://en.wikipedia.org/wiki/Wrapped_normal_distribution
        - https://en.wikipedia.org/wiki/Yamartino_method
        - https://en.wikipedia.org/wiki/Directional_statistics

        :param h5: h5parm
        :param optimal_score: score to fit solution interval
        :param ref_solint: reference solution interval
        """
        self.h5 = h5
        self.optimal_score = optimal_score
        self.ref_solint = ref_solint
        self.cstd = 0
        self.C = None
        self.station = station


    def plot_C(self, title=None, saveas=None, extrapoints=None):
        """
        Plot circstd score in function of solint for given C

        :param C: constant defining the noise level
        :param title: title for plot
        """
        # normal_sigmas = [n / 1000 for n in range(1, 10000)]
        # values = [circstd(normal(0, n, 300)) for n in normal_sigmas]
        # x = (self.C*limit**2) / (np.array(normal_sigmas) ** 2) / 2
        bestsolint = self.best_solint
        # plt.plot(x, values, alpha=0.5)
        solints = np.array(range(1, int(max(bestsolint * 200, self.ref_solint * 150))))/100
        plt.plot(solints, [self.theoretical_curve(float(t)) for t in solints], color='green')
        plt.scatter([self.ref_solint], [self.cstd], c='blue', label='measurement', s=80, marker='x')
        plt.scatter([bestsolint], [self.optimal_score], color='red', label='best solint', s=80, marker='x')
        if extrapoints is not None:
            plt.scatter(extrapoints[0], extrapoints[1], color='orange', label='other measurements', s=80, marker='x')
        plt.xlim(0, max(bestsolint * 1.5, self.ref_solint * 1.5))
        # plt.xlim(0, 0.2)
        plt.xlabel("solint (min)")
        plt.ylabel("circstd score")
        plt.legend(frameon=True, loc='upper right', fontsize=10)
        if title is not None:
            plt.title(title)
        if saveas is not None:
            plt.savefig(saveas)
        else:
            plt.show()

        return self


    def _circvar_to_normvar(self, circ_var):
        """
        Convert circular variance to normal variance

        return: circular variance
        """
        if circ_var > limit**2:
            sys.exit('ERROR: optimal score cannot be larger than pi')
        else:
            normvar = -2 * np.log(1 - circ_var / (limit**2))
            return normvar if normvar==normvar else sys.exit('ERROR: variance gives NaN')


    @property
    def _get_C(self):
        """
        Get constant defining the normal circular distribution

        :param cstd: circular standard deviation

        :return: C
        """
        if self.cstd==0:
            self.get_phasediff_score(station=self.station)
        normvar = self._circvar_to_normvar(self.cstd ** 2)
        return normvar * self.ref_solint


    def get_phasediff_score(self, station=False):
        """
        Calculate score for phasediff

        :return: circular standard deviation score
        """
        H = tables.open_file(self.h5)

        stations = [make_utf8(s) for s in list(H.root.sol000.antenna[:]['name'])]

        if not station:
            stations_idx = [stations.index(stion) for stion in stations if
                            ('RS' not in stion) &
                            ('ST' not in stion) &
                            ('CS' not in stion) &
                            ('DE' not in stion)]
        else:
            stations_idx = [stations.index(station)]

        axes = str(H.root.sol000.phase000.val.attrs["AXES"]).replace("b'", '').replace("'", '').split(',')
        axes_idx = sorted({ax: axes.index(ax) for ax in axes}.items(), key=lambda x: x[1], reverse=True)

        phase = H.root.sol000.phase000.val[:] * H.root.sol000.phase000.weight[:]
        H.close()

        phasemod = phase % (2 * np.pi)

        for ax in axes_idx:
            if ax[0] == 'pol':  # YX should be zero
                phasemod = phasemod.take(indices=0, axis=ax[1])
            elif ax[0] == 'dir':  # there should just be one direction
                if phasemod.shape[ax[1]] == 1:
                    phasemod = phasemod.take(indices=0, axis=ax[1])
                else:
                    sys.exit('ERROR: This solution file should only contain one direction, but it has ' +
                             str(phasemod.shape[ax[1]]) + ' directions')
            elif ax[0] == 'freq':  # faraday corrected
                phasemod = np.diff(phasemod, axis=ax[1])
            elif ax[0] == 'ant':  # take only international stations
                phasemod = phasemod.take(indices=stations_idx, axis=ax[1])

        phasemod[phasemod == 0] = np.nan

        self.cstd = circstd(phasemod, nan_policy='omit')

        return circstd(phasemod, nan_policy='omit')


    @property
    def best_solint(self):
        """
        Get optimal solution interval from phasediff, given C

        :return: value corresponding with increase solution interval
        """
        if self.cstd==0:
            self.get_phasediff_score(station=self.station)
        self.C = self._get_C
        optimal_cirvar = self.optimal_score ** 2
        return self.C / (self._circvar_to_normvar(optimal_cirvar))


    def theoretical_curve(self, t):
        """
        Theoretical curve based on circ statistics
        :param t: solution interval
        :return: circular std
        """
        if self.C is None:
            self.C = self._get_C
        return limit * np.sqrt(1 - np.exp(-(self.C / (2 * t))))


if __name__ == "__main__":

    # set std score, for which you want to find the solint
    optimal_score = 2

    # reference solution interval
    ref_solint = 10

    # solution file
    h5 = '../P23872.h5'

    # get solution interval
    S = GetSolint(h5, optimal_score, ref_solint)
    solint = S.best_solint

    # OPTIONAL: plot fit
    S.plot_C("T=" + str(round(solint, 2)) + " min")

