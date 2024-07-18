"""
Concat measurement sets with dummies

Example:
    python concat_with_dummies.py --concat_name concat.ms --freq_avg 4 --time_avg 8 *.ms
"""

__author__ = "Jurjen de Jong"

import casacore.tables as ct
import numpy as np
import sys
from glob import glob
import argparse
import os


def get_largest_divider(inp, max=1000):
    """
    Get largest divider

    :param inp: input number
    :param max: max divider

    :return: largest divider from inp bound by max
    """
    for r in range(max)[::-1]:
        if inp % r == 0:
            return r
    sys.exit("ERROR: code should not arrive here.")


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

    input_ms.sort(key=lambda x: x.split('/')[-1])

    # collect all frequencies in 1 numpy array
    for n, ms in enumerate(input_ms):
        t = ct.table(ms + '/::SPECTRAL_WINDOW', ack=False)

        if n == 0:
            chans = t.getcol("CHAN_FREQ")
        else:
            chans = np.append(chans, t.getcol("CHAN_FREQ"))

        t.close()

    return np.sort(chans), input_ms


def fill_freq_gaps(input, make_dummies, output_name, only_basename):
    """
    Fill the frequency gaps between sub-blocks with dummies (if requested)
    and return txt file with MS in order

    :param input: list or string with measurement sets
    :return: True if nog gaps, False if gaps
    """

    # get frequency channels
    chans, mslist = get_channels(input)
    # get differences between frequency steps
    chan_diff_single = np.abs(np.diff(chans, n=1))
    chans_per_ms = len(chan_diff_single)//len(mslist)+1
    chan_diff = np.abs(np.diff(chans, n=2))
    num_missing = [int(i) for i in chan_diff_single[:-1]/chan_diff_single[1:]/chans_per_ms][chans_per_ms-1::chans_per_ms]

    # file to write sorted list of MS
    file = open(output_name, 'w')

    # check gaps in freqs
    if np.sum(chan_diff) != 0:
        dummy_idx = set(
            (np.ndarray.flatten(np.argwhere(chan_diff > 0)) / len(chan_diff) * len(mslist)).round(0).astype(int))
        n = 0
        for idx in dummy_idx:
            if make_dummies:
                for k in range(int(num_missing[idx-1])):
                    print('dummy_' + str(n) + ' between ' + str(mslist[idx - 1]) + ' and ' + str(mslist[idx+n]))
                    mslist.insert(idx+n, 'dummy_' + str(n))
                    n += 1
            else:
                print('Gap between ' + str(mslist[idx - 1]) + ' and ' + str(mslist[idx]))
        for ms in mslist:
            if only_basename:
                ms = ms.split('/')[-1]
            file.write(ms + '\n')
        file.close()
        return False
    else:
        for ms in mslist:
            if only_basename:
                ms = ms.split('/')[-1]
            file.write(ms + '\n')
        file.close()
        return True


def split_ms_phasedir(mslist):
    """
    Separate MS with different phase centers
    """

    d = {}
    for ms in mslist:
        t = ct.table(f"{ms}::FIELD", ack=False)
        pd = ''.join([str(round(i, 6)) for i in t.getcol("PHASE_DIR").squeeze()])
        if pd not in d.keys():
            d.update({pd: [ms]})
        else:
            d[pd] += [ms]
    return d


def parse_args():
    """
    Parse input arguments
    """

    parser = argparse.ArgumentParser(
        description='Concattenate measurement sets, while taking into account frequency gaps')
    parser.add_argument('msin', nargs='+', help='MS')
    parser.add_argument('--msout', help='Concat name', type=str, default=None)
    parser.add_argument('--data_column', help='Data column', type=str, default='DATA')
    parser.add_argument('--phase_center', help='Phase shift to new center', type=str)
    parser.add_argument('--time_avg', help='Time averaging factor', type=int)
    parser.add_argument('--freq_avg', help='Frequency averaging factor', type=int)
    parser.add_argument('--time_res', help='Time resolution (in seconds)', type=int)
    parser.add_argument('--freq_res', help='Frequency resolution', type=str)
    parser.add_argument('--make_only_parset', action='store_true', help='Make only parset')
    parser.add_argument('--only_basename', action='store_true', help='Return only basename of msin')

    return parser.parse_args()


def make_parset(mss, concat_name, data_column, time_avg, freq_avg, time_res, freq_res, phase_center, only_basename):
    """
    Make parset for DP3

    :param ms: input measurement sets
    :param concat_name: name of concattenated measurement sets
    :param data_column: data column
    :param time_avg: time averaging
    :param freq_avg: frequency averaging
    :param only_basename: return only basename

    :return: parset
    """

    ms_dict = split_ms_phasedir(mss)
    parsets = []

    print(ms_dict)

    for dir, ms in ms_dict.items():


        if concat_name is None:
            concatname = ('_'.join([i for i in ms[0].split('_') if 'mhz' not in i.lower()]).
                           replace('mstargetphase','')+'.concat.ms').replace('..', '.').split('/')[-1]
        else:
            concatname = concat_name
        parsetname = concatname.replace('.concat.ms', '.parset')
        txtname = parsetname.replace('.parset', '.txt')

        if fill_freq_gaps(input=ms, make_dummies=True, output_name=txtname, only_basename=only_basename):
            print('--- SUCCESS: no frequency gaps found ---')

        # Write parset
        with open(txtname) as f:
            lines = f.readlines()
        parset = 'msin=' + '[' + ', '.join(lines).replace('\n', '') + ']\n'

        parset += 'msout=' + concatname
        parset += '\nmsin.datacolumn=' + data_column + \
                  '\nmsin.missingdata=True' \
                  '\nmsin.orderms=False' \
                  '\nmsout.storagemanager=dysco' \
                  '\nmsout.writefullresflag=False'
        steps = []

        if phase_center is not None:
            steps.append('ps')
            parset += '\nps.type=phaseshifter'
            parset += f'\nps.phasecenter={phase_center}'

            # Apply beam
            steps.append('beam')
            parset += '\nbeam.type=applybeam'
            parset += '\nbeam.updateweights=True'
            parset += '\nbeam.direction=[]'

        if time_avg is not None or freq_avg is not None or freq_res is not None or time_res is not None:
            steps.append('avg')
            parset += '\navg.type=averager'
            if time_avg is not None and time_res is not None:
                sys.exit("ERROR: choose time averaging or time resolution")
            if freq_avg is not None and freq_res is not None:
                sys.exit("ERROR: choose frequency averaging or frequency resolution")

            if time_avg is not None:
                parset += f'\navg.timestep={time_avg}'
            if time_res is not None:
                parset += f'\navg.timeresolution={time_res}'
            if freq_avg is not None:
                t = ct.table(ms[0] + "::SPECTRAL_WINDOW", ack=False)
                channum = len(t.getcol("CHAN_FREQ")[0])
                t.close()
                freqavg = get_largest_divider(channum, freq_avg + 1)
                if freqavg!=freq_avg:
                    print(f"WARNING: {channum} Channels can only be divided by {freqavg} and not by {freq_avg}")
                parset += f'\navg.freqstep={freqavg}'
            if freq_res is not None:
                parset += f'\navg.freqresolution={freq_res}'

        parset += '\nsteps='+str(steps).replace(" ", "").replace("'", "")
        with open(parsetname, 'w') as f:
            f.write(parset)
        parsets.append(parsetname)

    return parsets


def main():
    """
    Main script
    """
    args = parse_args()
    parsets = make_parset(args.msin, args.msout, args.data_column,
                args.time_avg, args.freq_avg, args.time_res, args.freq_res, args.phase_center, args.only_basename)
    if not args.make_only_parset:
        for parset in parsets:
            os.system('DP3 ' + parset)

if __name__ == '__main__':
    main()
