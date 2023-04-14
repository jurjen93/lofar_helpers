import tables
from glob import glob
import numpy as np
import csv
from scipy.stats import circstd
import sys
import pandas as pd

def gonio_score(phasemod):
    """
    Calculate a goniometric score (currently not used)

    :param phasemod: phase value
    :return: score
    """
    phasemod %= (2 * np.pi)
    p = np.zeros(phasemod.shape)
    p += np.where((phasemod < np.pi / 2) | ((phasemod < 2 * np.pi) & (phasemod > 3 * np.pi / 2)),
                  np.abs(np.sin(phasemod)), 0)
    p += np.where((phasemod < 3 * np.pi / 2) & (phasemod > np.pi / 2), 1 + np.abs(np.cos(phasemod)), 0)
    return p

def make_utf8(inp):
    """
    Convert input to utf8 instead of bytes

    :param inp: string input
    """

    try:
        inp = inp.decode('utf8')
        return inp
    except (UnicodeDecodeError, AttributeError):
        return inp


def get_scalarphasediff_score(h5):
    """
    Calculate score for calarphasediff

    :param h5: input h5 file
    :return: circular standard deviation score
    """
    H = tables.open_file(h5)

    stations = [make_utf8(s) for s in list(H.root.sol000.antenna[:]['name'])]
    distant_stations_idx = [stations.index(station) for station in stations if
                            ('RS' not in station) &
                            ('ST' not in station) &
                            ('CS' not in station) &
                            ('DE' not in station)]

    axes = str(H.root.sol000.phase000.val.attrs["AXES"]).replace("b'", '').replace("'", '').split(',')
    axes_idx = sorted({ax: axes.index(ax) for ax in axes}.items(), key=lambda x: x[1], reverse=True)

    phase = H.root.sol000.phase000.val[:]
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
            phasemod = phasemod.take(indices=distant_stations_idx, axis=ax[1])

    return circstd(phasemod, nan_policy='omit')

def rad_to_degree(inp):
    """
    Check if radians and convert to degree

    :param inp: two coordinates (RA, DEC)
    :return: output in degrees
    """
    try:
        if abs(inp[0])<np.pi and abs(inp[1])<np.pi:
            return inp*360/2/np.pi % 360
        else:
            return inp
    except ValueError:
        if abs(inp[0][0])<np.pi and abs(inp[0][1])<np.pi:
            return inp[0]*360/2/np.pi % 360
        else:
            return inp[0]

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--h5', nargs='+', help='selfcal scalarphasediff solutions', default=None)
    args = parser.parse_args()

    h5s = args.h5
    if h5s is None:
        h5s = glob("P*_scalarphasediff/scalarphasediff0*.h5")

    f = open('scalarphasediff_output.csv', 'w')
    writer = csv.writer(f)
    writer.writerow(["Source_id", "spd_score", 'RA', 'DEC'])
    for h5 in h5s:
        print(h5.split("_")[0])
        std = get_scalarphasediff_score(h5)
        print(std)
        H = tables.open_file(h5)
        dir = rad_to_degree(H.root.sol000.source[:]['dir'])
        H.close()
        writer.writerow([h5.split("_")[0], std, dir[0], dir[1]])

    f.close()

    # sort output
    df = pd.read_csv('scalarphasediff_output.csv').sort_values(by='spd_score').to_csv('scalarphasediff_output.csv', index=False)
