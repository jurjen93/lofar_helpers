#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Jurjen de Jong"

import argparse
import os
import re
import sys
from glob import glob
from pprint import pprint

from casacore.tables import table
import numpy as np


def case_insensitive_replace(text, old, new):
    """
    Case insensivitve_replace function
    """
    return re.sub(re.escape(old), new, text, flags=re.IGNORECASE)


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

    :param input_ms: list or string with MeasurementSets

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
        with table(ms + '/::SPECTRAL_WINDOW', ack=False) as t:

            if n == 0:
                chans = t.getcol("CHAN_FREQ")
            else:
                chans = np.append(chans, t.getcol("CHAN_FREQ"))


    return np.sort(chans), input_ms


def fill_freq_gaps(input: list = None, make_dummies: bool = None, output_name: str= None, only_basename: bool = None):
    """
    Fill the frequency gaps between sub-blocks with dummies (if requested)
    and return txt file with MS in order

    :param input: list or string with MeasurementSets
    :param make_dummies: boolean for making dummies or not
    :param output_name: string with output file name
    :param only_basename: boolean to only use base name

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
    with open(output_name, 'w') as file:

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
            return False
        else:
            for ms in mslist:
                if only_basename:
                    ms = ms.split('/')[-1]
                file.write(ms + '\n')
            return True


def split_ms_phasedir(mslist: list = None):
    """
    Separate MS with different phase centers

    :param mslist: list with MeasurementSets

    :return: phase dirs
    """

    d = {}
    for ms in mslist:
        t = table(f"{ms}::FIELD", ack=False)
        pd = ''.join([str(round(i, 6)) for i in t.getcol("PHASE_DIR").squeeze()])
        if pd not in d.keys():
            d.update({pd: [ms]})
        else:
            d[pd] += [ms]
    return d


def remove_flagged_antennas(msin: str = None):
    """
    Remove antennas that are full flagged (to save storage)

    input:
        - msfile: measurement set name
    """

    # Read antenna names from Measurement Set
    with table(f"{msin}::ANTENNA", ack=False) as ants:
        ants_names = ants.getcol("NAME")

    # Read main tables Measurement Set
    with table(msin, readonly=True, ack=False) as ms:
        # Read the antenna ID columns
        antenna1 = ms.getcol('ANTENNA1')
        antenna2 = ms.getcol('ANTENNA2')

        # Read the FLAG column
        flags = ms.getcol('FLAG')

    # Get the unique antenna indices
    unique_antennas = np.unique(np.concatenate((antenna1, antenna2)))

    # Identify fully flagged antennas
    fully_flagged_antennas = []
    for ant in unique_antennas:
        # Find rows with this antenna
        ant_rows = np.where((antenna1 == ant) | (antenna2 == ant))
        # Check if all data for this antenna is flagged
        if np.all(flags[ant_rows]):
            fully_flagged_antennas.append(ant)

    if len(fully_flagged_antennas) == 0:
        print(f'No flagged antennas for {msin}, move on.')
        return None

    else:
        # Get names of ants to filter
        ants_to_filter = ','.join([ants_names[idx] for idx in fully_flagged_antennas])
        print(f"Filtering fully flagged antennas: {ants_to_filter}")

        # Run DP3
        return f'\nfilter.type=filter\nfilter.remove=true\nfilter.baseline=!{ants_to_filter}'



def parse_args():
    """
    Parse input arguments
    """

    parser = argparse.ArgumentParser(description='Concatenate MeasurementSets or generate corresponding parset files. '
                                                 'This script takes into account frequency gaps.')
    parser.add_argument('--msin', nargs='+', help='Measurement set')
    parser.add_argument('--msout', help='Concat name', type=str, default=None)
    parser.add_argument('--data_column', help='Data column', type=str, default='DATA')
    parser.add_argument('--phase_center', help='Phase shift to new center', type=str)
    parser.add_argument('--time_avg', help='Time averaging factor', type=int)
    parser.add_argument('--freq_avg', help='Frequency averaging factor', type=int)
    parser.add_argument('--time_res', help='Time resolution (in seconds)', type=int)
    parser.add_argument('--freq_res', help='Frequency resolution', type=str)
    parser.add_argument('--remove_flagged_station', action='store_true', help='Remove flagged station (save output)')
    parser.add_argument('--make_only_parset', action='store_true', help='Make only parset')
    parser.add_argument('--only_basename', action='store_true', help='Return only basename of msin')

    return parser.parse_args()


def make_parset(mss: list = None, concat_name: str = None, data_column: str = None,
                time_avg: int= None, freq_avg: int = None, time_res=None, freq_res=None, phase_center: str = None,
                only_basename: bool = None, remove_flagged_station: bool = None):
    """
    Make parset for DP3

    :param mss: input MeasurementSets
    :param concat_name: name of concattenated MeasurementSets
    :param data_column: data column
    :param time_avg: time averaging factor
    :param freq_avg: frequency averaging factor
    :param time_res: time resolution
    :param freq_res: frequency resolution
    :param phase_center: phase center
    :param only_basename: return only basename
    :param remove_flagged_station: remove station when fully flagged

    :return: parset
    """

    ms_dict = split_ms_phasedir(mss)
    parsets = []

    pprint(ms_dict)

    for dir, ms in ms_dict.items():

        if concat_name is None:

            # Parse facet + L-number (if available). Specifically for VLBI pipeline
            matchf = re.search(r'facet_\d{2}-', ms[0].split('/')[-1])
            matchL = re.search(r'L\d{6}', ms[0].split('/')[-1])
            if matchL is None:
                with table(ms[0]+"::OBSERVATION", ack=False) as t:
                    try:
                        matchL = t.getcol("LOFAR_FILENAME")[0].split("_")[0]
                    except RuntimeError:
                        matchL = None

            if matchf is not None and matchL is not None:
                concatname = matchf.group()+matchL.group()+'.concat.ms'

            else:
                concatname = ('_'.join([i for i in ms[0].split('_') if 'mhz' not in i.lower()]).
                               replace('mstargetphase','')+'.concat.ms').replace('..', '.').split('/')[-1]
            parsetname = case_insensitive_replace(concatname, '.concat.ms', '.parset')

        else:
            concatname = concat_name
            parsetname = case_insensitive_replace(concatname, '.ms', '.parset')

        txtname = case_insensitive_replace(parsetname, '.parset', '.txt')

        if fill_freq_gaps(input=ms, make_dummies=True, output_name=txtname, only_basename=only_basename):
            print('--- SUCCESS: no frequency gaps found ---')

        # Write parset
        with open(txtname) as f:
            lines = f.readlines()

        msin = ', '.join(lines).replace('\n', '')

        parset = (
            f"msin=[{msin}]\n"
            f"msout={concatname}\n"
            f"msin.datacolumn={data_column}\n"
            f"msin.missingdata=True\n"
            f"msin.orderms=False\n"
            f"msout.storagemanager=dysco\n"
            f"msout.writefullresflag=False"
        )
        steps = []

        if phase_center is not None:
            steps.append('ps')
            parset += (f'\nps.type=phaseshifter'
                       f'\nps.phasecenter={phase_center}')

            # Apply beam
            steps.append('beam')
            parset += ('\nbeam.type=applybeam'
                       '\nbeam.updateweights=True'
                       '\nbeam.direction=[]')

        if time_avg is not None or freq_avg is not None or freq_res is not None or time_res is not None:
            steps.append('avg')
            parset += '\navg.type=averager'
            if time_avg is not None and time_res is not None:
                sys.exit("ERROR: select only time averaging or time resolution")
            if freq_avg is not None and freq_res is not None:
                sys.exit("ERROR: select only frequency averaging or frequency resolution")

            if time_avg is not None:
                parset += f'\navg.timestep={time_avg}'
            if time_res is not None:
                parset += f'\navg.timeresolution={time_res}'
            if freq_avg is not None:
                with table(ms[0] + "::SPECTRAL_WINDOW", ack=False) as t:
                    channum = len(t.getcol("CHAN_FREQ")[0])
                freqavg = get_largest_divider(channum, freq_avg + 1)
                if freqavg!=freq_avg:
                    print(f"WARNING: {channum} Channels can only be divided by {freqavg} and not by {freq_avg}")
                parset += f'\navg.freqstep={freqavg}'
            if freq_res is not None:
                parset += f'\navg.freqresolution={freq_res}'

        # Remove station when fully flagged
        if remove_flagged_station:
            rmv = remove_flagged_antennas(ms[0])
            if rmv is not None:
                parset += rmv
                steps.append('filter')

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
                args.time_avg, args.freq_avg, args.time_res,
                args.freq_res, args.phase_center, args.only_basename, args.remove_flagged_station)
    if not args.make_only_parset:
        for parset in parsets:
            os.system('DP3 ' + parset)


if __name__ == '__main__':
    main()