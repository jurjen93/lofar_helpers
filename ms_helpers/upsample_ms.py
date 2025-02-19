from argparse import ArgumentParser
import numpy as np
from casacore.tables import table
from numba import njit, set_num_threads
from joblib import Parallel, delayed
import os
from sys import stdout

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


@njit
def interpolate_data_single_baseline(time_in, freqs_in, data_in, time_out, freqs_out):
    """
    For each output time and frequency, pick the nearest input value.
    Parameters:
        time_in   : input times.
        freqs_in  : input frequencies.
        data_in   : data with shape (n_time_in, n_freq_in, 4).
        time_out  : output times (for this baseline).
        freqs_out : output frequencies.

    Returns:
        new_data  : output data with shape (len(time_out), len(freqs_out), 4).
    """

    n_time_out = time_out.shape[0]
    n_freq_out = freqs_out.shape[0]
    n_pol = data_in.shape[2]
    new_data = np.empty((n_time_out, n_freq_out, n_pol), dtype=data_in.dtype)

    for i in range(n_time_out):
        t_val = time_out[i]
        t_idx = find_nearest_index(time_in, t_val)
        for j in range(n_freq_out):
            f_val = freqs_out[j]
            f_idx = find_nearest_index(freqs_in, f_val)
            for p in range(n_pol):
                new_data[i, j, p] = data_in[t_idx, f_idx, p]
    return new_data


def process_baseline(baseline, from_ms, out_times, out_ant1, out_ant2,
                     freqs_in, freqs_out, column):
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
                data_in = sub_in.getcol(column)  # expected shape: (n_time_in, n_freq_in, 4)
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
    new_data = interpolate_data_single_baseline(time_in, freqs_in, data_in, time_out, freqs_out)
    return mask, new_data


def interpolate_flags(from_ms, to_ms, column):
    """
    Interpolate the given column from an input measurement set to an output measurement set.

    Parameters:
        from_ms (str): Path (or identifier) of the input measurement set.
        to_ms (str): Path (or identifier) of the output measurement set.
        column (str): The name of the column to interpolate.
    """

    with table(to_ms, ack=False) as t_out:
        # Retrieve frequency axes from the spectral window subtables.
        with table(from_ms + "::SPECTRAL_WINDOW", ack=False) as spw_in:
            freqs_in = spw_in.getcol("CHAN_FREQ")[0]
        with table(to_ms + "::SPECTRAL_WINDOW", ack=False) as spw_out:
            freqs_out = spw_out.getcol("CHAN_FREQ")[0]

        # Get time and antenna info from the output MS
        out_times = t_out.getcol("TIME")
        out_ant1 = t_out.getcol("ANTENNA1")
        out_ant2 = t_out.getcol("ANTENNA2")

        # Unique baselines from the output.
        baselines = np.unique(np.vstack([out_ant1, out_ant2]).T, axis=0)

        # Assuming global shape: (n_rows, len(freqs_out), 4)
        n_rows = out_times.shape[0]
        n_freq_out = freqs_out.shape[0]
        n_pol = 4  # adjust if needed

        # Create a temporary memmap for the global output data.
        tmp_dir = './'
        tmp_filename = os.path.join(tmp_dir, f"{column}.dat")
        global_data = np.memmap(tmp_filename, dtype=np.complex128, mode='w+',
                                shape=(n_rows, n_freq_out, n_pol))
        global_data[:] = 0

        # Process each baseline in parallel.
        results = Parallel(n_jobs=-1, verbose=10)(
            delayed(process_baseline)(
                baseline, from_ms, out_times, out_ant1, out_ant2, freqs_in, freqs_out, column
            ) for baseline in baselines
        )

        # Update the global memmap with each baseline's new data.
        for mask, new_data in results:
            if new_data is not None:
                global_data[mask, :, :] = new_data
        global_data.flush()

        # If the column doesn't exist in the output MS, add it.
        if column not in t_out.colnames():
            # Get a column description from the input MS's 'DATA' column.
            with table(from_ms, ack=False) as t_in:
                desc = t_in.getcoldesc('DATA')
            print('Creating column', column)
            desc['name'] = column
            # IMPORTANT: add the column to the output MS, not the input.
            t_out.addcols(desc)
        else:
            print(column, 'already exists')

    # Write the updated column from the memmap back to the output measurement set.
    n_rows = global_data.shape[0]
    chunk_size = min(100000, n_rows)
    # Loop over the memmap in chunks
    with table(to_ms, ack=False, readonly=False) as t_out:
        for start in range(0, n_rows, chunk_size):
            print_progress_bar(start, n_rows)
            end = min(start + chunk_size, n_rows)
            chunk = np.array(global_data[start:end, :, :])
            t_out.putcol(column, chunk, startrow=start)

    # Clean up the temporary memmap file.
    try:
        os.remove(tmp_filename)
    except Exception:
        pass


def parse_args():
    """
    Parse input arguments
    """

    parser = ArgumentParser(description='Interpolate data from from a specific time/freq resolution to another.')
    parser.add_argument('--from_ms', help='MS input from where to interpolate')
    parser.add_argument('--to_ms', help='MS to interpolate to (your output set)')
    parser.add_argument('--column', help='Column name from in and output', default='MODEL_DATA')

    return parser.parse_args()


def main():
    """
    Main script
    """

    args = parse_args()

    # interpolate flags
    print(f'Interpolate from {args.from_ms} to {args.to_ms}')
    interpolate_flags(args.from_ms, args.to_ms, args.column)

if __name__ == '__main__':
    main()
