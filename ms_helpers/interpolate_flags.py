"""
With this script you can flag data with a specific freq/time resolution and interpolate the flags to a dataset (or datasets) with another freq/time resolution.
Make sure that all datasets have the same antennas with same antenna indices and originate from the same observation.

Strategy:
    1) aoflagger with default_StokesV.lua strategy (Stokes V)
    2) interpolate the new flags to the output measurement set (higher freq/time resolution dataset)

Usage:
    python interpolate_flags.py --msin dataset_lowres.ms dataset_original.ms

    dataset_lowres.ms --> a dataset with a lower time/freq resolution which originates from dataset_original.ms
    dataset_original.ms --> the original dataset

    NOTE: when interpolating to multiple datasets, that the first dataset always the dataset is to run aoflagger on.
"""

from casacore.tables import table
import numpy as np
from argparse import ArgumentParser
from subprocess import call
from scipy.interpolate import griddata
from sys import stdout


def run(command):
    """
    Execute a shell command through subprocess

    Args:
        command (str): the command to execute.
    Returns:
        None
    """

    retval = call(command, shell=True)
    if retval != 0:
        print('FAILED to run ' + command + ': return value is ' + str(retval))
        raise Exception(command)
    return retval


def print_progress_bar(index, total, bar_length=50):
    """
    Prints a progress bar to the console.

    :param::param:
        - index: the current index (0-based) in the iteration.
        - total: the total number of indices.
        - bar_length: the character length of the progress bar (default 50).
    """

    percent_complete = (index + 1) / total
    filled_length = int(bar_length * percent_complete)
    bar = "█" * filled_length + '-' * (bar_length - filled_length)
    stdout.write(f'\rProgress: |{bar}| {percent_complete * 100:.1f}% Complete')
    stdout.flush()  # Important to ensure the progress bar is updated in place

    # Print a new line on completion
    if index == total - 1:
        print()


def runaoflagger(ms, strategy='default_StokesV.lua'):
    """
    Run aoglagger on a Measurement Set.

    Args:
        mslist (list): list of Measurement Sets to iterate over.
    Returns:
        None
    """

    if strategy is not None:
       cmd = 'aoflagger -strategy ' + strategy + ' ' + ms
    else:
       cmd = 'aoflagger ' + ms
    print(cmd)
    run(cmd)
    return


def make_ant_pairs(n_ant, n_time):
    """
    Generate ANTENNA1 and ANTENNA2 arrays for an array with M antennas over N time slots.

    :param:
        - n_ant: Number of antennas in the array.
        - n_int: Number of time slots.

    :return:
        - ANTENNA1
        - ANTENNA2
    """

    # Generate all unique pairs of antennas for one time slot
    antenna_pairs = [(i, j) for i in range(n_ant) for j in range(i + 1, n_ant)]

    # Expand the pairs across n_time time slots
    antenna1 = np.array([pair[0] for pair in antenna_pairs] * n_time)
    antenna2 = np.array([pair[1] for pair in antenna_pairs] * n_time)

    return antenna1, antenna2


def interpolate_flags(flagged_ms, ms, make_backup_flags):
    """
    Args:
        flagged_ms: measurement set from where to interpolate
        ms: the pre-averaged measurement set
        make_backup_flags: make backup of flags
    Returns:
        interpolated flags
    """

    ants = table(flagged_ms + "::ANTENNA", ack=False)
    baselines = np.c_[make_ant_pairs(ants.nrows(), 1)]
    ants.close()

    t1 = table(flagged_ms, ack=False)
    t2 = table(ms, ack=False, readonly=False)

    if make_backup_flags:
        np.save(ms+'.flags.npy', t2.getcol("FLAG"))

    # Get freq axis first table
    t = table(flagged_ms+'::SPECTRAL_WINDOW', ack=False)
    freq_flagged_axis = t.getcol("CHAN_FREQ")[0]
    t.close()

    # Get freq axis second table
    t = table(ms+'::SPECTRAL_WINDOW', ack=False)
    freq_axis = t.getcol("CHAN_FREQ")[0]
    t.close()

    # Loop over baselines
    for n, baseline in enumerate(baselines):
        print_progress_bar(n, len(baselines))
        sub1 = t1.query(f"ANTENNA1={baseline[0]} AND ANTENNA2={baseline[1]}")
        sub2 = t2.query(f"ANTENNA1={baseline[0]} AND ANTENNA2={baseline[1]}")

        time_flagged_axis = sub1.getcol("TIME")
        data_flagged = np.take(sub1.getcol('FLAG'), indices=0, axis=-1).astype(int)

        time_axis = sub2.getcol("TIME")
        data = np.take(sub2.getcol('FLAG'), indices=0, axis=-1).astype(int)

        # Create the grid for interpolation
        grid_x, grid_y = np.meshgrid(freq_flagged_axis, time_flagged_axis)

        # Flatten the grid and data for griddata function
        points = np.column_stack((grid_x.ravel(), grid_y.ravel()))
        values = data_flagged.ravel()

        # Create the interpolation points
        interp_grid_x, interp_grid_y = np.meshgrid(freq_axis, time_axis)

        # Perform the nearest-neighbor interpolation
        new_flags = griddata(points, values, (interp_grid_x, interp_grid_y), method='nearest')

        # Reshape the data to match the original data shape
        new_flags = new_flags.reshape(time_axis.size, freq_axis.size)

        # Apply the new flags to the data
        data += new_flags

        # Store the updated data for the current baseline
        sub2.putcol('FLAG', np.tile(np.expand_dims(np.clip(data, a_min=0, a_max=1), axis=-1), 4).astype(bool))

        sub1.close()
        sub2.close()

    t1.close()
    t2.close()


def parse_args():
    """
    Parse input arguments
    """

    parser = ArgumentParser(description='Flag data from a lower freq/time resolution to a higher one')
    parser.add_argument('--msin', help='MS input from where to interpolate')
    parser.add_argument('--backup_flags', action='store_true', default=None, help='Make backup of flags')
    parser.add_argument('--skip_flagging', action='store_true', default=None, help='Skip flagging')
    parser.add_argument('msout', nargs='+', help='MS output from where to apply new interpolated flags')

    return parser.parse_args()


def main():
    """
    Main script
    """

    args = parse_args()

    # run aoflagger on the input MS
    if not args.skip_flagging:
        runaoflagger(args.msin)

    # interpolate flags
    for ms in args.msout:
        print(f'Interpolate to {ms}')
        interpolate_flags(args.msin, ms, args.backup_flags)


if __name__ == '__main__':
    main()
