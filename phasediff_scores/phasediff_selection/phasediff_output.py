import tables
import numpy as np
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import normal
from scipy.stats import circstd
from glob import glob
import csv

try:
    import scienceplots
    plt.style.use(['science', 'ieee'])
except:
    pass

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
    from ..find_solint import GetSolint

    parser = argparse.ArgumentParser()
    parser.add_argument('--h5', nargs='+', help='selfcal phasediff solutions', default=None)
    parser.add_argument('--station', help='for specific station', default=None)
    parser.add_argument('--all_stations', action='store_true', help='for all stations specifically')
    args = parser.parse_args()

    # set std score, for which you want to find the solint
    optimal_score = 0.5

    # reference solution interval
    ref_solint = 10

    h5s = args.h5
    if h5s is None:
        h5s = glob("P*_phasediff/phasediff0*.h5")

    if args.station is not None:
        station = args.station
    else:
        station = ''

    f = open('phasediff_output.csv', 'w')
    writer = csv.writer(f)
    writer.writerow(["source", "spd_score", "best_solint", 'RA', 'DEC'])
    for h5 in h5s:
        try:
            S = GetSolint(h5, optimal_score, ref_solint)
            if args.all_stations:
                H = tables.open_file(h5)
                stations = [make_utf8(s) for s in list(H.root.sol000.antenna[:]['name'])]
                H.close()

            else:
                stations = [station]
            for station in stations:
                std = S.get_phasediff_score(station=station)
                solint = S.best_solint
                H = tables.open_file(h5)
                dir = rad_to_degree(H.root.sol000.source[:]['dir'])
                writer.writerow([h5 + station, std, solint, dir[0], dir[1]])
                S.plot_C("T=" + str(round(solint, 2)) + " min", saveas=h5+station+'.png')
                H.close()
        except:
            pass


    f.close()

    # sort output
    pd.read_csv('phasediff_output.csv').sort_values(by='spd_score').to_csv('phasediff_output.csv', index=False)
