import tables
import numpy as np
from argparse import ArgumentParser
from sys import exit
import re


def parse_args():
    """
    Command line argument parser

    :return: parsed arguments
    """

    parser = ArgumentParser()
    parser.add_argument('h5', nargs='+', help='Input h5parm')
    parser.add_argument('--freqrange', help='Flag frequency range in MHz (example: 140-158)')
    parser.add_argument('--timerange', type=str, help='Flag time range where each number is followed by s (seconds), m (minutes) or h (hours). (Example: 600s-1258s or 4h30m-5h15m)')
    parser.add_argument('--ant', help='Antenna to flag (example: RS409HBA), if not given it flags for all antennas.')
    return parser.parse_args()


def convert_input_to_time_seconds(time_str):
    """
    Convert time value (5s, 30m, 1h) to a time in seconds
    """
    time_split = re.split(r'(?<=\D)(?=\d)|(?<=\d)(?=\D)', time_str)
    time = 0
    for n, item in enumerate(time_split):
        if item.isdigit():
            if time_split[n + 1] == 's':
                time += float(item)
            elif time_split[n + 1] == 'm':
                time += float(item) * 60
            elif time_split[n + 1] == 'h':
                time += float(item) * 3600
            else:
                exit("Give your --timerange by for instance 600s-1258s or 4h30m-5h15m, where each number ends on time unit 's' (seconds), 'm' (minutes), 'h' (hours).")
    return time


def main():
    """
    Main script
    """

    args = parse_args()
    ant = args.ant

    if args.freqrange is not None:
        min_freq, max_freq = [float(f) * 1e6 for f in args.freqrange.split('-')]
        if ant is not None:
            print(f'Flagging from {min_freq}Hz to {max_freq}Hz for antenna {ant}')
        else:
            print(f'Flagging from {min_freq}Hz to {max_freq}Hz')

    if args.timerange is not None:
        min_time_str, max_time_str = [str(f) for f in args.timerange.split('-')]
        min_time = convert_input_to_time_seconds(min_time_str)
        max_time = convert_input_to_time_seconds(max_time_str)

    if args.timerange is None and args.freqrange is None:
        exit('ERROR: need --timerange or --freqrange')

    for h5 in args.h5:
        print('-----\n' + h5 + '\n-----')
        H = tables.open_file(h5, 'r+')
        for solset in H.root._v_groups.keys():
            print(solset)
            ss = H.root._f_get_child(solset)
            for soltab in ss._v_groups.keys():
                st = ss._f_get_child(soltab)
                new_weights = st.weight[:]
                ants = [a.decode('utf8') for a in st.ant[:]]
                if ant is not None:
                    ant_idxs = [ants.index(ant)]
                else:
                    ant_idxs = list(range(len(ants)))
                AXES = st.weight.attrs["AXES"].decode('utf8').split(',')
                if AXES[0]!='time' or AXES[1]!='freq' or AXES[2]!='ant':
                    exit('ERROR: h5 invalid structure. Fix h5 by running it through h5_merger.py (without merging it)')

                if args.freqrange is not None:
                    freqs = st.freq[:]
                    freq_indices = np.where((freqs >= min_freq) & (freqs <= max_freq))[0]
                else:
                    freq_indices = list(range(len(st.freq[:])+1))

                if args.timerange is not None:
                    time = st.time[:]
                    time_indices = np.where((time-time.min() >= min_time) & (time-time.min() <= max_time))[0]
                else:
                    time_indices = list(range(len(st.time[:])+1))

                for ant_idx in ant_idxs:
                    print(f"Flag for: {ants[ant_idx]}")
                    H.root._f_get_child(solset)._f_get_child(soltab). \
                        weight[time_indices[0]:time_indices[-1], freq_indices[0]:freq_indices[-1], ant_idx, ...] = 0

        H.close()


if __name__ == '__main__':
    main()
