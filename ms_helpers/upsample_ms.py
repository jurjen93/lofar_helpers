from argparse import ArgumentParser
import os
from sys import stdout

from casacore.tables import table
from joblib import Parallel, delayed
from numba import njit, prange, set_num_threads
import numpy as np

__author__ = "Jurjen de Jong"

set_num_threads(max(os.cpu_count()-1, 1))


def print_progress_bar(index, total, bar_length=50):
    """
    Prints a progress bar to the console.

    :param:
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


def make_ant_pairs(n_ant, n_time):
    """
    Make ANTENNA1 and ANTENNA2 pair arrays for an array with M antennas over N time slots.

    :param:
        - n_ant: Number of antennas in the array.
        - n_time: Number of time slots.

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


def add_axis(arr, ax_size):
    """
    Add ax dimension with a specific size

    :param:
        - arr: numpy array
        - ax_size: axis size

    :return:
        - output with new axis dimension with a particular size
    """

    or_shape = arr.shape
    new_shape = list(or_shape) + [ax_size]
    return np.repeat(arr, ax_size).reshape(new_shape)


@njit
def find_nearest_index(ref, value):
    """
    Given a sorted 1D array 'ref', find the index of the element nearest to 'value' using a binary search.
    """
    low = 0
    high = len(ref)
    while low < high:
        mid = (low + high) // 2
        if ref[mid] < value:
            low = mid + 1
        else:
            high = mid
    idx = low
    if idx >= len(ref):
        idx = len(ref) - 1
    if idx > 0:
        if abs(ref[idx - 1] - value) < abs(ref[idx] - value):
            idx = idx - 1
    return idx


@njit(parallel=True)
def interpolate_data_single_baseline(time_in, freqs_in, data_in, time_out, freqs_out):
    """
    Standard nearest-neighbor interpolation.

    Parameters:
      time_in, freqs_in : 1D arrays (input time and frequency)
      data_in : 3D array of shape (n_time_in, n_freq_in, 4)
      time_out, freqs_out : 1D arrays (output time and frequency)

    Returns:
      new_data : 3D array of shape (len(time_out), len(freqs_out), 4)
    """

    n_time_out = time_out.shape[0]
    n_freq_out = freqs_out.shape[0]
    n_pol = data_in.shape[2]
    new_data = np.empty((n_time_out, n_freq_out, n_pol), dtype=data_in.dtype)

    for i in prange(n_time_out):
        t_val = time_out[i]
        t_idx = find_nearest_index(time_in, t_val)
        for j in range(n_freq_out):
            f_val = freqs_out[j]
            f_idx = find_nearest_index(freqs_in, f_val)
            for p in range(n_pol):
                new_data[i, j, p] = data_in[t_idx, f_idx, p]
    return new_data


@njit(parallel=True)
def interpolate_data_single_baseline_stokes_i(time_in, freqs_in, data_in, time_out, freqs_out):
    """
    Nearest-neighbor interpolation that computes Stokes I from the input data.
    Stokes I is computed as:

         Stokes I = (data_in[..., 0] + data_in[..., 3]) / 2

    Parameters:
      time_in, freqs_in : 1D arrays (input time and frequency)
      data_in         : 3D array of shape (n_time_in, n_freq_in, 4)
      time_out, freqs_out : 1D arrays (output time and frequency)

    Returns:
      new_data : 2D array of shape (len(time_out), len(freqs_out))
    """

    n_time_out = time_out.shape[0]
    n_freq_out = freqs_out.shape[0]
    new_data = np.empty((n_time_out, n_freq_out), dtype=data_in.dtype)

    for i in prange(n_time_out):
        t_val = time_out[i]
        t_idx = find_nearest_index(time_in, t_val)
        for j in range(n_freq_out):
            f_val = freqs_out[j]
            f_idx = find_nearest_index(freqs_in, f_val)
            new_data[i, j] = data_in[t_idx, f_idx]
    return new_data


@njit(parallel=True)
def compute_stokes_I(data_in):
    """
    Compute Stokes I = (XX + YY) / 2

    Parameters:
      data_in : input data

    Returns:
      A 2D numpy array of shape (n_time, n_freq) with the computed Stokes I.
    """

    n_time, n_freq, n_pol = data_in.shape
    # Create an output array with the same time and frequency dimensions.
    result = np.empty((n_time, n_freq), dtype=data_in.dtype)

    for i in prange(n_time):
        for j in range(n_freq):
            result[i, j] = (data_in[i, j, 0] + data_in[i, j, n_pol - 1]) / 2
    return result


def process_baseline(baseline, from_ms, out_times, out_ant1, out_ant2,
                     freqs_in, freqs_out, column, stokes_i, data_out):
    """
    Process a single baseline:
      - Identify output rows corresponding to the baseline.
      - Open the input measurement set for that baseline.
      - Interpolate the data.

    If no rows are found with the given antenna order, it tries the swapped order (normally not an issue..).
    Returns a tuple (mask, new_data) where 'mask' selects the output rows.
    """

    ant1, ant2 = baseline
    mask = (out_ant1 == ant1) & (out_ant2 == ant2)
    if not np.any(mask):
        return mask, None

    time_out = out_times[mask]

    # Helper function to perform the query.
    def query_ms(query_str):
        with table(from_ms, ack=False) as t_in:
            with t_in.query(query_str) as sub_in:
                time_in = sub_in.getcol("TIME")
                data_in = sub_in.getcol(column)
                if stokes_i:
                    data_in = compute_stokes_I(data_in)
        return time_in, data_in

    # First try the query with the output order.
    query_str = f"ANTENNA1={ant1} AND ANTENNA2={ant2}"
    try:
        time_in, data_in = query_ms(query_str)
    except Exception:
        time_in = None

    # If no rows are returned (or time_in is empty), try the reversed order.
    if time_in is None or time_in.size == 0:
        query_str = f"ANTENNA1={ant2} AND ANTENNA2={ant1}"
        try:
            time_in, data_in = query_ms(query_str)
        except Exception:
            return mask, None

    # Now perform the interpolation.
    if stokes_i:
        new_data = interpolate_data_single_baseline_stokes_i(time_in, freqs_in, data_in, time_out, freqs_out)
    else:
        new_data = interpolate_data_single_baseline(time_in, freqs_in, data_in, time_out, freqs_out)

    data_out[mask, ...] = new_data

    return


@njit(parallel=True)
def subtract_arrays(a, b):
    result = np.empty_like(a)
    # Use parallel range for possible speedup on large arrays.
    for i in prange(a.size):
        result.flat[i] = a.flat[i] - b.flat[i]
    return result


def interpolate(from_ms, to_ms, column: str = "MODEL_DATA", stokes_i: bool = False, tmpdir: str = './', subtract: bool = False):
    """
    Interpolate the given column from an input measurement set to an output measurement set.

    Parameters:
        from_ms (str): Path (or identifier) of the input measurement set.
        to_ms (str): Path (or identifier) of the output measurement set.
        column (str): The name of the column to interpolate.
        stokes_i: Stokes I only interpolation (XX+YY)/2 and XY=YX=0
        tmpdir: temporary directory for memmap storage
        subtract: Subtraction DATA-column in output
    """

    with table(to_ms, ack=False) as t_out:
        # Retrieve frequency axes from the spectral window subtables.
        with table(from_ms + "::SPECTRAL_WINDOW", ack=False) as spw_in:
            freqs_in = spw_in.getcol("CHAN_FREQ")[0]
        with table(to_ms + "::SPECTRAL_WINDOW", ack=False) as spw_out:
            freqs_out = spw_out.getcol("CHAN_FREQ")[0]

        # Get time and antenna info from the output MS
        print("Get TIME, ANTENNA1, ANTENNA2 columns ...")
        out_times = t_out.getcol("TIME")
        out_ant1 = t_out.getcol("ANTENNA1")
        out_ant2 = t_out.getcol("ANTENNA2")

        # Unique baselines from the output.
        baselines = np.unique(np.vstack([out_ant1, out_ant2]).T, axis=0)

        n_rows = out_times.shape[0]
        n_freq_out = freqs_out.shape[0]

        # Create a temporary memmap for the global output data.
        tmp_filename = os.path.join(tmpdir, f"{column}.dat")

        if stokes_i:
            shape = (n_rows, n_freq_out)
        else:
            shape = (n_rows, n_freq_out, 4)

        print(f"Make memmap {tmp_filename} ...")
        global_data = np.memmap(tmp_filename, dtype=np.complex64, mode='w+', shape=shape)
        global_data[:] = 0

        # Process each baseline in parallel.
        print("Interpolate data for all baselines ...")
        Parallel(n_jobs=-1, verbose=10)(
            delayed(process_baseline)(
                baseline, from_ms, out_times, out_ant1, out_ant2, freqs_in, freqs_out, column, stokes_i, global_data
            ) for baseline in baselines
        )

        print("Flush data ...")
        global_data.flush()


    # Write the updated column from the memmap back to the output measurement set.
    print(f"Write results to {column} in chunks...")
    with table(to_ms, ack=False, readonly=False) as t_out:

        # If the column doesn't exist in the output MS, add it.
        if not subtract:
            if column not in t_out.colnames():
                with table(from_ms, ack=False) as t_in:
                    desc = t_in.getcoldesc(column)
                desc['name'] = column
                t_out.addcols(desc)
            else:
                print(column, 'already exists (will overwrite)')

        n_rows = global_data.shape[0]
        chunk_size = min(200000, n_rows)
        if subtract:
            print(f"SUBTRACT ==> DATA = DATA - {column}")
        for start in range(0, n_rows, chunk_size):
            print_progress_bar(start, n_rows)
            end = min(start + chunk_size, n_rows)
            chunk = np.array(global_data[start:end, ...])
            if stokes_i:
                chunk = add_axis(chunk, 4)
                chunk[..., 1] = 0
                chunk[..., 2] = 0
            if subtract:
                chunk = subtract_arrays(t_out.getcol("DATA", startrow=start, nrow=chunk_size), chunk)
                t_out.putcol('DATA', chunk, startrow=start, nrow=chunk_size)
            else:
                t_out.putcol(column, chunk, startrow=start, nrow=chunk_size)
        print_progress_bar(n_rows, n_rows)
    # Clean up the temporary memmap file.
    try:
        os.remove(tmp_filename)
    except Exception:
        pass
    print(f"\nInterpolation Finished ==> {to_ms}::{column}\n")


def parse_args():
    """
    Parse input arguments
    """

    parser = ArgumentParser(description='Interpolate DATA or MODEL_DATA from a specific time/freq resolution to another.')
    parser.add_argument('--from_ms', help='MS input from where to interpolate')
    parser.add_argument('--to_ms', help='MS to interpolate to (your output set)')
    parser.add_argument('--column', help='Column name from in and output', default='MODEL_DATA')
    parser.add_argument('--stokes_i', action='store_true', help='Perform interpolation in Stokes I only (faster)')
    parser.add_argument('--tmpdir', help='Temporary output directory for memmap', default='./')
    parser.add_argument('--subtract', action='store_true', help='Do subtract for DATA = DATA - INPUT_COLUMN.')

    return parser.parse_args()


def main():
    """
    Main script
    """

    args = parse_args()

    print(f'Interpolate from {args.from_ms} to {args.to_ms}')
    interpolate(args.from_ms, args.to_ms, args.column, args.stokes_i, args.tmpdir, args.subtract)

if __name__ == '__main__':
    main()
