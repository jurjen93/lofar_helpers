"""
Concat measurement sets with dummies

Example:
    python concat_with_dummies.py --ms *.ms --concat_name concat.ms --data_column CORRECTED_DATA
"""

__author__ = "Jurjen de Jong"

import casacore.tables as ct
import numpy as np
import sys
from glob import glob
import argparse
import os


def get_channels(input_ms):
    """
    Get frequency channels

    :param input: list or string with measurement sets
    :return: sorted concatenated channels
    """

    # if string, make list and sort on name (so give clever names if you want to sort on freq)
    if type(input_ms) == str:
        input_ms = sorted(glob(input_ms))
    elif type(input_ms) == list:
        input_ms = sorted(input_ms)
    else:
        sys.exit("ERROR: wrong input type for " + str(input_ms) + " (only accept strings or lists).")

    if len(input_ms) == 0:
        sys.exit("ERROR: no MS given.")

    # collect all frequencies in 1 numpy array
    for n, ms in enumerate(input_ms):
        t = ct.table(ms + '/::SPECTRAL_WINDOW')

        if n == 0:
            chans = t.getcol("CHAN_FREQ")
        else:
            chans = np.append(chans, t.getcol("CHAN_FREQ"))

        t.close()

    return np.sort(chans), input_ms


def fill_freq_gaps(input, make_dummies, output_name):
    """
    Fill the frequency gaps between sub-blocks with dummies (if requested)
    and return txt file with MS in order

    :param input: list or string with measurement sets
    :return: True if nog gaps, False if gaps
    """

    # get frequency channels
    chans, mslist = get_channels(input)
    # get differences between frequency steps
    chan_diff = np.abs(np.diff(chans, n=2))

    # file to write sorted list of MS
    file = open(output_name, 'w')

    # check gaps in freqs
    if np.sum(chan_diff) != 0:
        dummy_idx = set(
            (np.ndarray.flatten(np.argwhere(chan_diff > 0)) / len(chan_diff) * len(mslist)).round(0).astype(int))
        for n, idx in enumerate(dummy_idx):
            if make_dummies:
                print('dummy_' + str(n) + ' between ' + str(mslist[idx - 1]) + ' and ' + str(mslist[idx]))
                mslist.insert(idx, 'dummy_' + str(n))
            else:
                print('Gap between ' + str(mslist[idx - 1]) + ' and ' + str(mslist[idx]))
        for ms in mslist:
            file.write(ms + '\n')
        file.close()
        return False
    else:
        for ms in mslist:
            file.write(ms + '\n')
        file.close()
        return True


def parse_args():
    """
    Parse input arguments
    """

    parser = argparse.ArgumentParser(
        description='Concattenate measurement sets, while taking into account frequency gaps')
    parser.add_argument('--ms', nargs='+', help='MS', required=True)
    parser.add_argument('--concat_name', help='Concat name', type=str, default='concat.ms')
    parser.add_argument('--parset_name', help='Parset_name', type=str, default='concat.parset')
    parser.add_argument('--data_column', help='Data column', type=str, default='DATA')

    return parser.parse_args()


def make_parset(parset_name, ms, concat_name, data_column):
    """
    Make parset for DP3

    :param parset_name: Name of the parset
    :param ms: input measurement sets
    :param concat_name: name of concattenated measurement sets

    :return: parset
    """

    txtname = parset_name.replace('.parset', '.txt')

    if fill_freq_gaps(input=ms, make_dummies=True, output_name=txtname):
        print('--- SUCCESS: no frequency gaps found ---')

    # write parset
    with open(txtname) as f:
        lines = f.readlines()
    parset = 'msin=' + '[' + ', '.join(lines).replace('\n', '') + ']\n'
    parset += 'msout=' + concat_name
    parset += '\nmsin.datacolumn=' + data_column + \
              '\nmsin.missingdata=True' \
              '\nmsin.orderms=False' \
              '\nmsout.storagemanager=dysco' \
              '\nmsout.writefullresflag=False' \
              '\nsteps=[]'
    with open(parset_name, 'w') as f:
        f.write(parset)

    return parset


def main():
    """
    Main script
    """
    args = parse_args()
    make_parset(args.parset_name, args.ms, args.concat_name, args.data_column)
    os.system('DP3 ' + args.parset_name)


if __name__ == '__main__':
    main()
