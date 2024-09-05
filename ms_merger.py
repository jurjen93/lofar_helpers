"""
LOFAR SIDEREAL VISIBILITY AVERAGER
This script can be used to average visibilities over sidereal time, when using multiple observations of the same FoV.

Example: python ms_merger.py --msout <MS_NAME> *.ms
The wildcard is in this example combining a collection of measurement sets

Strategy:
    1) Make a template using the 'default_ms' option from casacore.tables (Template class).
       The template inclues all baselines, frequency, and smallest time spacing from all input MS.
       Time is converted to Local Sidereal Time (LST).

    2) Map baselines from input MS to template MS.
        This step makes *baseline_mapping folders with the baseline mappings in json files.

    3) Interpolate new UV data with nearest neighbours.

    4) Make exact mapping between input MS and template MS, using only UV data points.

    5) Average measurement sets in the template (Stack class).
        The averaging is done with a weighted average, using the FLAG and WEIGHT_SPECTRUM columns.
"""

from casacore.tables import table, default_ms, taql
import numpy as np
from os import system as run_command
import os
from shutil import rmtree, move
import sys
from astropy.time import Time
from astropy.coordinates import EarthLocation
import astropy.units as u
from pprint import pprint
from argparse import ArgumentParser
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
import psutil
from math import ceil
from glob import glob
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
from joblib import Parallel, delayed
import tempfile


one_lst_day_sec = 86164.1


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
    sys.stdout.write(f'\rProgress: |{bar}| {percent_complete * 100:.1f}% Complete')
    sys.stdout.flush()  # Important to ensure the progress bar is updated in place

    # Print a new line on completion
    if index == total - 1:
        print()


def is_dysco_compressed(ms):
    """
    Check if MS is dysco compressed

    :param:
        - ms: measurement set
    """

    t = table(ms, readonly=True, ack=False)
    dysco = t.getdesc()["DATA"]['dataManagerGroup'] == 'DyscoData'
    t.close()
    return dysco


def decompress(ms):
    """
    running DP3 to remove dysco compression

    :param:
        - ms: measurement set
    """

    if is_dysco_compressed(ms):

        print('\n----------\nREMOVE DYSCO COMPRESSION\n----------\n')

        if os.path.exists(f'{ms}.tmp'):
            rmtree(f'{ms}.tmp')
        run_command(f"DP3 msin={ms} msout={ms}.tmp steps=[]")
        print('----------')
        return ms + '.tmp'

    else:
        return ms


def make_odd(i):
    """
    Make odd number

    :param:
        - i: digit

    :return: odd digit
    """
    if int(i) % 2 == 0:
        i += 1

    return int(i)


def compress(ms):
    """
    running DP3 to apply dysco compression

    :param:
        - ms: measurement set
    """

    if not is_dysco_compressed(ms):

        print('DYSCO COMPRESSION')

        cmd = f"DP3 msin={ms} msout={ms}.tmp msout.overwrite=true msout.storagemanager=dysco"

        steps = []

        steps = str(steps).replace("'", "").replace(' ','')
        cmd += f' steps={steps}'

        run_command(cmd)

        try:
            t = table(f"{ms}.tmp", ack=False) # test if exists
            t.close()
        except RuntimeError:
            sys.exit(f"ERROR: dysco compression failed (please check {ms})")

        rmtree(ms)
        move(f"{ms}.tmp", ms)

        print('----------')
        return ms

    else:
        return ms


def time_resolution(resolution, fov_diam, time_smearing=0.95):
    """
    Calculate the best time resolution, given a time_smearing allowance

    Using formulas from Bridle & Schwab (1999)

    :params:
        - resolution: resolution in arcseconds
        - fov_diam: longest FoV diameter in degrees
        - time_smearing: allowable time smearing
    :return: integration time in seconds
    """

    # Convert distance from degrees to radians
    distance_from_phase_center_rad = np.deg2rad(fov_diam/2)

    # Calculate angular resolution (radians)
    angular_resolution_rad = resolution*4.8481*1e-6

    int_time = 2.9*10**4*(angular_resolution_rad*np.sqrt(1-time_smearing)/
                          distance_from_phase_center_rad)

    return int_time


def get_ms_content(ms):
    """
    Get MS content

    :param:
        - ms: measurement set

    :return:
        - station names
        - frequency channels
        - total time in seconds
        - delta time
    """

    T = table(ms, ack=False)
    F = table(ms+"::SPECTRAL_WINDOW", ack=False)
    A = table(ms+"::ANTENNA", ack=False)
    L = table(ms+"::LOFAR_ANTENNA_FIELD", ack=False)
    S = table(ms+"::LOFAR_STATION", ack=False)

    # Get all lofar antenna info
    lofar_stations = list(zip(
                        S.getcol("NAME"),
                        S.getcol("CLOCK_ID")
                    ))

    # Get all station information
    stations = list(zip(
                   A.getcol("NAME"),
                   A.getcol("POSITION"),
                   A.getcol("DISH_DIAMETER"),
                   A.getcol("LOFAR_STATION_ID"),
                   A.getcol("LOFAR_PHASE_REFERENCE"),
                   L.getcol("NAME"),
                   L.getcol("COORDINATE_AXES"),
                   L.getcol("TILE_ELEMENT_OFFSET"),
                ))

    chan_num = F.getcol("NUM_CHAN")[0]
    channels = F.getcol("CHAN_FREQ")[0]
    dfreq = np.diff(sorted(set(channels)))[0]
    time = sorted(np.unique(T.getcol("TIME")))
    time_lst = mjd_seconds_to_lst_seconds(T.getcol("TIME"))
    time_min_lst, time_max_lst = time_lst.min(), time_lst.max()
    total_time_seconds = max(time)-min(time)
    dt = np.diff(sorted(set(time)))[0]

    print(f'\nCONTENT from {ms}\n'
          '----------\n'
          f'Stations: {", ".join([s[0] for s in lofar_stations])}\n'
          f'Number of channels: {chan_num}\n'
          f'Channel width: {dfreq} Hz\n'
          f'Total time: {round(total_time_seconds/3600, 2)} hrs\n'
          f'Delta time: {dt} seconds\n'
          f'----------')

    S.close()
    L.close()
    T.close()
    F.close()
    A.close()

    return {'stations': stations,
            'lofar_stations': lofar_stations,
            'channels': channels,
            'dfreq': dfreq,
            'total_time_seconds': total_time_seconds,
            'dt': dt,
            'time_min_lst': time_min_lst,
            'time_max_lst': time_max_lst}


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


def get_station_id(ms):
    """
    Get station with corresponding id number

    :param:
        - ms: measurement set

    :return:
        - antenna names, IDs
    """

    t = table(ms+'::ANTENNA', ack=False)
    ants = t.getcol("NAME")
    t.close()

    t = table(ms+'::FEED', ack=False)
    ids = t.getcol("ANTENNA_ID")
    t.close()

    return ants, ids


def mjd_seconds_to_lst_seconds(mjd_seconds, longitude_deg=6.869837):
    """
    Convert time in modified Julian Date time to LST

    :param:
        - mjd_seconds: modified Julian date time in seconds
        - longitde_deg: longitude telescope in degrees (6.869837 for LOFAR )

    :return:
        time in LST
    """

    # Convert seconds to days for MJD
    mjd_days = mjd_seconds / 86400.0

    # Create an astropy Time object
    time_utc = Time(mjd_days, format='mjd', scale='utc')

    # Define the observer's location using longitude (latitude doesn't affect LST)
    location = EarthLocation(lon=longitude_deg * u.deg, lat=52.915122 * u.deg)

    # Calculate LST in hours
    lst_hours = time_utc.sidereal_time('apparent', longitude=location.lon).hour

    # Convert LST from hours to seconds
    lst_seconds = lst_hours * 3600.0

    return lst_seconds


def same_phasedir(mslist: list = None):
    """
    Have MS same phase center?

    :param:
        - mslist: measurement set list
    """

    for n, ms in enumerate(mslist):
        t = table(ms+'::FIELD', ack=False)
        if n==0:
            phasedir = t.getcol("PHASE_DIR")
        else:
            if not np.all(phasedir == t.getcol("PHASE_DIR")):
                sys.exit("MS do not have the same phase center, check "+ms)


def sort_station_list(zipped_list):
    """
    Sorts a list of lists (or tuples) based on the first element of each inner list or tuple,
    which is necessary for a zipped list with station names and positions

    :param:
        - list_of_lists (list of lists or tuples): The list to be sorted.

    :return:
        sorted list
    """
    return sorted(zipped_list, key=lambda item: item[0])


def unique_station_list(station_list):
    """
    Filters a list of stations only based on first element

    :param:
        - station_list: Stations to be filtered.

    :return:
        - filtered list of stations
    """
    unique_dict = {}
    for item in station_list:
        if item[0] not in unique_dict:
            unique_dict[item[0]] = item
    return list(unique_dict.values())


def n_baselines(n_antennas: int = None):
    """
    Return number of baselines

    :param:
        - n_antennas: number of antennas

    :return: number of baselines
    """

    return n_antennas * (n_antennas - 1) // 2


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


def repeat_elements(original_list, repeat_count):
    """
    Repeat each element in the original list a specified number of times.

    :param:
        - original_list: The original list to be transformed.
        - repeat_count: The number of times each element should be repeated.

    :return:
        - A new list where each element from the original list is repeated.
    """
    return np.array([element for element in original_list for _ in range(repeat_count)])


def find_closest_index(arr, value):
    """
    Find the index of the closest value in the array to the given value.

    :param:
        - arr: A NumPy array of values.
        - value: A float value to find the closest value to.

    :return:
        - The index of the closest value in the array.
    """
    # Calculate the absolute difference with the given value
    differences = np.abs(arr - value)
    print(f"Minimal difference: {min(differences)}")

    # Find the index of the minimum difference
    closest_index = np.argmin(differences)

    return closest_index


def find_closest_index_list(a1, a2):
    """
    Find the indices of the closest values between two arrays.

    :param:
        - a1: first array
        - a2: second array

    :return:
        - The indices of the closest value in the array.
    """
    return abs(np.array(a1)[:, None] - np.array(a2)).argmin(axis=1)


def find_closest_index_multi_array(a1, a2):
    """
    Find the indices of the closest values between two multi-D arrays, preserving nearest neighbor relationships
    even when the arrays are negated.

    :param:
        - a1: first array (shape NxM)
        - a2: second array (shape PxM)

    :return:
        - A list of indices corresponding to the nearest neighbors in a2 for each point in a1.
    """
    # Ensure inputs are numpy arrays
    a1 = np.asarray(a1)
    a2 = np.asarray(a2)

    # Concatenate the original and negated versions of a2
    combined_a2 = np.vstack((a2, -a2))

    # Build a KDTree for the combined array (original + negated)
    tree = cKDTree(combined_a2)

    # Query the tree for the closest points in a1
    distances, indices = tree.query(a1)

    # Handle the indices since we doubled the array by adding -a2
    final_indices = [i if i < len(a2) else i - len(a2) for i in indices]

    return list(final_indices)


def map_array_dict(arr, dct):
    """
    Maps elements of the input_array to new values using a mapping_dict

    :param:
        - arr: numpy array of integers that need to be mapped.
        - dct: dictionary where each key-value pair represents an original value and its corresponding new value.

    :return:
        - An array where each element has been mapped according to the mapping_dict.
        (If an element in the input array does not have a corresponding mapping in the mapping_dict, it will be left unchanged)
    """

    # Vectorize a simple lookup function
    lookup = np.vectorize(lambda x: dct.get(x, x))

    # Apply this vectorized function to the input array
    output_array = lookup(arr)

    return output_array


def get_avg_factor(mslist, less_avg=1):
    """
    Calculate optimized averaging factor

    :param:
        - mslist: measurement set list
        - less_avg: factor to reduce averaging
    :return:
        - averaging factor
    """

    uniq_obs = []
    for ms in mslist:
        obs = table(ms + "::OBSERVATION", ack=False)
        uniq_obs.append(obs.getcol("TIME_RANGE")[0][0])
        obs.close()
    obs_count = len(np.unique(uniq_obs))
    avgfactor = ceil(np.sqrt(obs_count / less_avg))
    if avgfactor < 1:
        return avgfactor
    else:
        return int(avgfactor)


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



def resample_uwv(uvw_arrays, row_idxs, time, time_ref):
    """
    Resample a uvw array to have N rows.

    :param:
        - uvw_arrays: UVW array with shape (num_points, 3)
        - row_idxs: Indices of rows to resample
        - time: Original time array
        - time_ref: Reference time array to resample to

    :return:
        - Resampled UVW array
    """

    # Get the original shape
    num_points, num_coords = uvw_arrays.shape

    if num_coords != 3:
        raise ValueError("Input array must have shape (num_points, 3)")

    # Sort the time array and corresponding UVW arrays
    sorted_indices = np.argsort(time)
    time_sorted = time[sorted_indices]
    uvw_sorted = uvw_arrays[sorted_indices, :]

    # Create a single interpolation function for the entire UVW array
    interp_func = interp1d(time_sorted, uvw_sorted, axis=0, kind='nearest', fill_value='extrapolate')

    # Apply the interpolation function to the reference times
    resampled_array = interp_func(time_ref[row_idxs])

    return resampled_array



def resample_array(data, factor):
    """
    Resamples the input data array such that the number of points increases by a factor.
    The lowest and highest values remain the same, and the spacing between points remains equal.

    :param:
        - data: The input data array to resample.
        - factor: The factor by which to increase the number of data points.

    :return:
        - The resampled data array.
    """
    # Original number of points
    n_points = len(data)

    # Number of points in the resampled array
    new_n_points = factor * (n_points - 1) + 1

    # Generate the new set of equally spaced indices
    original_indices = np.arange(n_points)
    new_indices = np.linspace(0, n_points - 1, new_n_points + 1)

    # Perform linear interpolation
    resampled_data = np.interp(new_indices, original_indices, data)

    return resampled_data


def get_data_arrays(column: str = 'DATA', nrows: int = None, freq_len: int = None):
    """
    Get data arrays (new data and weights)

    :param:
        - column: column name (DATA, WEIGHT_SPECTRUM, WEIGHT, OR UVW)
        - nrows: number of rows
        - freq_len: frequency axis length

    :return:
        - new_data: new data array (empty array with correct shape)
        - weights: weights corresponding to new data array (empty array with correct shape)
    """

    tmpfilename = column.lower()+'.tmp.dat'
    tmpfilename_weights = column.lower()+'_weights.tmp.dat'

    if column in ['UVW']:
        weights = np.memmap(tmpfilename_weights, dtype=np.float16, mode='w+', shape=(nrows, 3))
        weights[:] = 0
    else:
        weights = None

    if column in ['DATA', 'WEIGHT_SPECTRUM']:
        if column == 'DATA':
            dtp = np.complex128
        elif column == 'WEIGHT_SPECTRUM':
            dtp = np.float32
        else:
            dtp = np.float32
        shape = (nrows, freq_len, 4)

    elif column == 'WEIGHT':
        shape, dtp = (nrows, freq_len), np.float32

    elif column == 'UVW':
        shape, dtp = (nrows, 3), np.float32

    else:
        sys.exit("ERROR: Use only DATA, WEIGHT_SPECTRUM, WEIGHT, or UVW")

    new_data = np.memmap(tmpfilename, dtype=dtp, mode='w+', shape=shape)

    return new_data, weights


def load_json(file_path):
    """Load json file"""

    with open(file_path, 'r') as file:
        return json.load(file)


def squeeze_to_intlist(arr):
    """Squeeze array and make list with integers"""

    squeezed = np.squeeze(arr).astype(int)
    if squeezed.ndim == 0:
        return [squeezed.item()]
    elif squeezed.ndim == 1:
        return squeezed.tolist()
    else:
        return squeezed.tolist()


def remove_flagged_entries(input_table):
    """
    Remove flagged entries.
    Note that this corrupts the time axis.
    """
    # Define the output table temporary name
    output_table = input_table + '.copy.tmp'

    # Open the input table
    with table(input_table, ack=False) as tb:
        # Select rows that do not match the deletion criteria
        selected_rows = tb.query('NOT all(WEIGHT_SPECTRUM == 0)')

        # Create a new table with the selected rows
        selected_rows.copy(output_table, deep=True)

    # Overwrite the input table with the new table
    rmtree(input_table)
    move(output_table, input_table)

    # Make new time axis TODO: issues with new time axis for BDA datasets
    # t = table(input_table, ack=False, readonly=False)
    # time_old = np.unique(t.getcol("TIME"))
    # tm = np.linspace(time_old.min(), time_old.max(), t.nrows())
    #
    # ants = table(input_table + "::ANTENNA", ack=False)
    # baselines = np.c_[make_ant_pairs(ants.nrows(), 1)]
    # ants.close()
    #
    # t = repeat_elements(time_range, baseline_count)


def sum_arrays_chunkwise(array1, array2, chunk_size=1000, n_jobs=-1, un_memmap=True):
    """
    Sums two arrays in chunks using joblib for parallel processing.

    :param:
        - array1: np.ndarray or np.memmap
        - array2: np.ndarray or np.memmap
        - chunk_size: int, size of each chunk
        - n_jobs: int, number of jobs for parallel processing (-1 means using all processors)
        - un_memmap: bool, whether to convert memmap arrays to regular arrays if they fit in memory

    :return:
        - np.ndarray or np.memmap: result array which is the sum of array1 and array2
    """

    # Ensure the arrays have the same length
    assert len(array1) == len(array2), "Arrays must have the same length"

    # Check if un-memmap is needed and feasible
    if un_memmap and isinstance(array1, np.memmap):
        try:
            array1 = np.array(array1)
        except MemoryError:
            pass  # If memory error, fall back to using memmap

    if un_memmap and isinstance(array2, np.memmap):
        try:
            array2 = np.array(array2)
        except MemoryError:
            pass  # If memory error, fall back to using memmap

    n = len(array1)
    # Determine the output storage type based on input type
    if isinstance(array1, np.memmap) or isinstance(array2, np.memmap):
        # Create a temporary file to store the result as a memmap
        temp_file = tempfile.NamedTemporaryFile(delete=False)
        result_array = np.memmap(temp_file.name, dtype=array1.dtype, mode='w+', shape=array1.shape)
    else:
        result_array = np.empty_like(array1)

    def sum_chunk_to_result(start, end):
        result_array[start:end] = array1[start:end] + array2[start:end]

    # Create a generator for chunk indices
    chunks = ((i, min(i + chunk_size, n)) for i in range(0, n, chunk_size))

    # Parallel processing with threading preferred for better I/O handling
    Parallel(n_jobs=n_jobs, prefer="threads")(delayed(sum_chunk_to_result)(start, end) for start, end in chunks)

    return result_array


class Template:
    """Make template measurement set based on input measurement sets"""
    def __init__(self, msin: list = None, outname: str = 'empty.ms'):
        self.mslist = msin
        self.outname = outname

    def add_spectral_window(self):
        """
        Add SPECTRAL_WINDOW as sub table
        """

        print("Add table ==> " + self.outname + "::SPECTRAL_WINDOW")

        tnew_spw_tmp = table(self.ref_table.getkeyword('SPECTRAL_WINDOW'), ack=False)
        newdesc = tnew_spw_tmp.getdesc()
        for col in ['CHAN_WIDTH', 'CHAN_FREQ', 'RESOLUTION', 'EFFECTIVE_BW']:
            newdesc[col]['shape'] = np.array([self.channels.shape[-1]])

        tnew_spw = table(self.outname + '::SPECTRAL_WINDOW', newdesc, readonly=False, ack=False)
        tnew_spw.addrows(1)
        chanwidth = np.expand_dims([np.squeeze(np.diff(self.channels))[0]]*self.chan_num, 0)
        tnew_spw.putcol("NUM_CHAN", np.array([self.chan_num]))
        tnew_spw.putcol("CHAN_FREQ", self.channels)
        tnew_spw.putcol("CHAN_WIDTH", chanwidth)
        tnew_spw.putcol("RESOLUTION", chanwidth)
        tnew_spw.putcol("EFFECTIVE_BW", chanwidth)
        tnew_spw.putcol("REF_FREQUENCY", np.nanmean(self.channels))
        tnew_spw.putcol("MEAS_FREQ_REF", np.array([5]))  # Why always 5?
        tnew_spw.putcol("TOTAL_BANDWIDTH", [np.max(self.channels)-np.min(self.channels)-chanwidth[0][0]])
        tnew_spw.putcol("NAME", 'Stacked_MS_'+str(int(np.nanmean(self.channels)//1000000))+"MHz")
        tnew_spw.flush(True)
        tnew_spw.close()
        tnew_spw_tmp.close()

    def add_stations(self):
        """
        Add ANTENNA and FEED tables
        """

        print("Add table ==> " + self.outname + "::ANTENNA")

        stations = [sp[0] for sp in self.station_info]
        st_id = dict(zip(set(
            [stat[0:8] for stat in stations]),
            range(len(set([stat[0:8] for stat in stations])))
        ))
        ids = [st_id[s[0:8]] for s in stations]
        positions = np.array([sp[1] for sp in self.station_info])
        diameters = np.array([sp[2] for sp in self.station_info])
        phase_ref = np.array([sp[4] for sp in self.station_info])
        names = np.array([sp[5] for sp in self.station_info])
        coor_axes = np.array([sp[6] for sp in self.station_info])
        tile_element = np.array([sp[7] for sp in self.station_info])
        lofar_names = np.array([sp[0] for sp in self.lofar_stations_info])
        clock = np.array([sp[1] for sp in self.lofar_stations_info])

        tnew_ant_tmp = table(self.ref_table.getkeyword('ANTENNA'), ack=False)
        newdesc = tnew_ant_tmp.getdesc()
        tnew_ant_tmp.close()

        tnew_ant = table(self.outname + '::ANTENNA', newdesc, readonly=False, ack=False)
        tnew_ant.addrows(len(stations))
        print('Total number of output stations: ' + str(tnew_ant.nrows()))
        tnew_ant.putcol("NAME", stations)
        tnew_ant.putcol("TYPE", ['GROUND-BASED']*len(stations))
        tnew_ant.putcol("POSITION", positions)
        tnew_ant.putcol("DISH_DIAMETER", diameters)
        tnew_ant.putcol("OFFSET", np.array([[0., 0., 0.]] * len(stations)))
        tnew_ant.putcol("FLAG_ROW", np.array([False] * len(stations)))
        tnew_ant.putcol("MOUNT", ['X-Y'] * len(stations))
        tnew_ant.putcol("STATION", ['LOFAR'] * len(stations))
        tnew_ant.putcol("LOFAR_STATION_ID", ids)
        tnew_ant.putcol("LOFAR_PHASE_REFERENCE", phase_ref)
        tnew_ant.flush(True)
        tnew_ant.close()

        print("Add table ==> " + self.outname + "::FEED")

        tnew_ant_tmp = table(self.ref_table.getkeyword('FEED'), ack=False)
        newdesc = tnew_ant_tmp.getdesc()
        tnew_ant_tmp.close()

        tnew_feed = table(self.outname + '::FEED', newdesc, readonly=False, ack=False)
        tnew_feed.addrows(len(stations))
        tnew_feed.putcol("POSITION", np.array([[0., 0., 0.]] * len(stations)))
        tnew_feed.putcol("BEAM_OFFSET", np.array([[[0, 0], [0, 0]]] * len(stations)))
        tnew_feed.putcol("POL_RESPONSE", np.array([[[1. + 0.j, 0. + 0.j], [0. + 0.j, 1. + 0.j]]] * len(stations)).astype(np.complex64))
        tnew_feed.putcol("POLARIZATION_TYPE", {'shape': [len(stations), 2], 'array': ['X', 'Y'] * len(stations)})
        tnew_feed.putcol("RECEPTOR_ANGLE", np.array([[-0.78539816, -0.78539816]] * len(stations)))
        tnew_feed.putcol("ANTENNA_ID", np.array(range(len(stations))))
        tnew_feed.putcol("BEAM_ID", np.array([-1] * len(stations)))
        tnew_feed.putcol("INTERVAL", np.array([28799.9787008] * len(stations)))
        tnew_feed.putcol("NUM_RECEPTORS", np.array([2] * len(stations)))
        tnew_feed.putcol("SPECTRAL_WINDOW_ID", np.array([-1] * len(stations)))
        tnew_feed.putcol("TIME", np.array([5.e9] * len(stations)))
        tnew_feed.flush(True)
        tnew_feed.close()

        print("Add table ==> " + self.outname + "::LOFAR_ANTENNA_FIELD")

        tnew_ant_tmp = table(self.ref_table.getkeyword('LOFAR_ANTENNA_FIELD'), ack=False)
        newdesc = tnew_ant_tmp.getdesc()

        tnew_ant_tmp.close()

        tnew_field = table(self.outname + '::LOFAR_ANTENNA_FIELD', newdesc, readonly=False, ack=False)
        tnew_field.addrows(len(stations))
        tnew_field.putcol("ANTENNA_ID", np.array(range(len(stations))))
        tnew_field.putcol("NAME", names)
        tnew_field.putcol("COORDINATE_AXES", np.array(coor_axes))
        tnew_field.putcol("TILE_ELEMENT_OFFSET", np.array(tile_element))
        tnew_field.putcol("TILE_ROTATION", np.array([0]*len(stations)))
        # tnew_field.putcol("ELEMENT_OFFSET", ???) TODO: fix for primary beam construction
        # tnew_field.putcol("ELEMENT_RCU", ???) TODO: fix for primary beam construction
        # tnew_field.putcol("ELEMENT_FLAG", ???) TODO: fix for primary beam construction
        tnew_field.flush(True)
        tnew_field.close()

        print("Add table ==> " + self.outname + "::LOFAR_STATION")

        tnew_ant_tmp = table(self.ref_table.getkeyword('LOFAR_STATION'), ack=False)
        newdesc = tnew_ant_tmp.getdesc()
        tnew_ant_tmp.close()

        tnew_station = table(self.outname + '::LOFAR_STATION', newdesc, readonly=False, ack=False)
        tnew_station.addrows(len(lofar_names))
        tnew_station.putcol("NAME", lofar_names)
        tnew_station.putcol("FLAG_ROW", np.array([False] * len(lofar_names)))
        tnew_station.putcol("CLOCK_ID", np.array(clock))
        tnew_station.flush(True)
        tnew_station.close()

    def make_template(self, overwrite: bool = True, time_res: int = None, avg_factor: float = 1):
        """
        Make template MS based on existing MS

        :param:
            - overwrite: overwrite output file
            - time_res: time resolution in seconds
            - avg_factor: averaging factor
        """

        if overwrite:
            if os.path.exists(self.outname):
                rmtree(self.outname)

        same_phasedir(self.mslist)

        # Get data columns
        # Initialize variables outside the loop
        unique_stations, unique_channels, unique_lofar_stations = [], [], []
        min_t_lst, min_dt, dfreq_min, max_t_lst = None, None, None, None

        def process_ms(ms):
            """Parallel"""
            mscontent = get_ms_content(ms)
            stations, lofar_stations, channels, dfreq, total_time_seconds, dt, min_t, max_t = mscontent.values()
            return stations, lofar_stations, channels, dfreq, dt, min_t, max_t

        with ThreadPoolExecutor() as executor:
            future_to_ms = {executor.submit(process_ms, ms): ms for ms in self.mslist}
            for future in as_completed(future_to_ms):
                stations, lofar_stations, channels, dfreq, dt, min_t, max_t = future.result()

                if min_t_lst is None:
                    min_t_lst, min_dt, dfreq_min, max_t_lst = min_t, dt, dfreq, max_t
                else:
                    min_t_lst = min(min_t_lst, min_t)
                    min_dt = min(min_dt, dt)
                    dfreq_min = min(dfreq_min, dfreq)
                    max_t_lst = max(max_t_lst, max_t)

                unique_stations.extend(stations)
                unique_channels.extend(channels)
                unique_lofar_stations.extend(lofar_stations)

        self.station_info = unique_station_list(unique_stations)
        self.lofar_stations_info = unique_station_list(unique_lofar_stations)

        chan_range = np.arange(min(unique_channels), max(unique_channels) + dfreq_min, dfreq_min)
        self.channels = np.sort(np.expand_dims(np.unique(chan_range), 0))
        self.chan_num = self.channels.shape[-1]

        if time_res is not None:
            time_range = np.arange(min_t_lst, max_t_lst + min_dt, time_res)

        else:
            time_range = np.arange(min_t_lst, max_t_lst + min_dt, min_dt/avg_factor)

        baseline_count = n_baselines(len(self.station_info))
        nrows = baseline_count*len(time_range)

        # Take one ms for temp usage
        tmp_ms = self.mslist[0]

        # Remove dysco compression
        self.tmpfile = decompress(tmp_ms)
        self.ref_table = table(self.tmpfile, ack=False)

        # Data description
        newdesc_data = self.ref_table.getdesc()

        # Reshape
        for col in ['DATA', 'FLAG', 'WEIGHT_SPECTRUM']:
            newdesc_data[col]['shape'] = np.array([self.chan_num, 4])

        newdesc_data.pop('_keywords_')

        pprint(newdesc_data)
        print()

        # Make main table
        default_ms(self.outname, newdesc_data)
        tnew = table(self.outname, readonly=False, ack=False)
        tnew.addrows(nrows)
        ant1, ant2 = make_ant_pairs(len(self.station_info), len(time_range))
        t = repeat_elements(time_range, baseline_count)
        tnew.putcol("TIME", t)
        tnew.putcol("TIME_CENTROID", t)
        tnew.putcol("ANTENNA1", ant1)
        tnew.putcol("ANTENNA2", ant2)
        tnew.putcol("EXPOSURE", np.array([np.diff(time_range)[0]] * nrows))
        tnew.putcol("FLAG_ROW", np.array([False] * nrows))
        tnew.putcol("INTERVAL", np.array([np.diff(time_range)[0]] * nrows))
        tnew.flush(True)
        tnew.close()

        # Set SPECTRAL_WINDOW info
        self.add_spectral_window()

        # Set ANTENNA/STATION info
        self.add_stations()

        # Set other tables (annoying table locks prevent parallel processing)
        for subtbl in ['FIELD', 'HISTORY', 'FLAG_CMD', 'DATA_DESCRIPTION',
                       'LOFAR_ELEMENT_FAILURE', 'OBSERVATION', 'POINTING',
                       'POLARIZATION', 'PROCESSOR', 'STATE']:
            try:
                print("Add table ==> " + self.outname + "::" + subtbl)

                tsub = table(self.tmpfile+"::"+subtbl, ack=False, readonly=False)
                tsub.copy(self.outname + '/' + subtbl, deep=True)
                tsub.flush(True)
                tsub.close()
            except:
                print(subtbl+" unknown")

        self.ref_table.close()

        # Cleanup
        if 'tmp' in self.tmpfile:
            rmtree(self.tmpfile)

    def make_mapping_lst(self):
        """
        Make mapping json files essential for efficient stacking
        These map LST times from input MS to template MS.
        Note that these are not accurate mappings but good first estimate, which are later corrected with final mappings.
        """

        T = taql(f"SELECT TIME,ANTENNA1,ANTENNA2 FROM {os.path.abspath(self.outname)} ORDER BY TIME")

        ref_time = T.getcol("TIME")
        time_len = ref_time.__len__()
        ref_uniq_time = np.unique(ref_time)
        ref_antennas = np.sort(np.c_[T.getcol("ANTENNA1"), T.getcol("ANTENNA2")])

        T.close()

        def process_antpair(antpair):
            """
            Making json files with antenna pair mappings.
            Mapping INPUT MS idx --> OUTPUT MS idx

            :param:
                - antpair: antenna pair

            """

            # Get idx
            pair_idx = squeeze_to_intlist(np.argwhere(np.all(antennas == antpair, axis=1)))
            ref_pair_idx = squeeze_to_intlist(np.argwhere(np.all(ref_antennas == antpair, axis=1))[time_idxs])

            # Make mapping dict
            mapping = {int(pair_idx[i]): int(ref_pair_idx[i]) for i in range(min(pair_idx.__len__(), ref_pair_idx.__len__()))}

            # Ensure the mapping folder exists
            os.makedirs(mapping_folder, exist_ok=True)
            # Define file path
            file_path = os.path.join(mapping_folder, '-'.join(antpair.astype(str)) + '.json')

            # Write to file
            with open(file_path, 'w') as f:
                json.dump(mapping, f)

        def run_parallel_mapping(uniq_ant_pairs):
            """
            Parallel processing of mapping with unique antenna pairs

            :param:
                - uniq_ant_pairs: unique antenna pairs to loop over in parallel
            """

            # Number of threads in the pool (adjust based on available resources)
            with ThreadPoolExecutor(max_workers=max(os.cpu_count()-3, 1)) as executor:
                # Submit tasks to the executor
                futures = [executor.submit(process_antpair, antpair) for antpair in uniq_ant_pairs]

                # Optionally, gather results or handle exceptions
                for future in futures:
                    future.result()  # This will raise exceptions if any occurred during the execution

        ref_stats, ref_ids = get_station_id(self.outname)

        # Make mapping
        for ms in self.mslist:

            print(f'\nMapping: {ms}')

            # Open MS table
            t = taql(f"SELECT TIME,ANTENNA1,ANTENNA2 FROM {os.path.abspath(ms)} ORDER BY TIME")

            # Make antenna mapping in parallel
            mapping_folder = ms + '_baseline_mapping'

            # Verify if folder exists
            if not check_folder_exists(mapping_folder):
                os.makedirs(mapping_folder, exist_ok=False)

                # Get MS info
                new_stats, new_ids = get_station_id(ms)
                id_map = dict(zip(new_ids, [ref_stats.index(a) for a in new_stats]))

                # Time in LST
                time = mjd_seconds_to_lst_seconds(t.getcol("TIME"))
                uniq_time = np.unique(time)
                time_idxs = find_closest_index_list(uniq_time, ref_uniq_time)

                # Map antenna pairs to same as ref (template)
                antennas = np.sort(np.c_[map_array_dict(t.getcol("ANTENNA1"), id_map), map_array_dict(t.getcol("ANTENNA2"), id_map)])

                # Unique antenna pairs
                uniq_ant_pairs = np.unique(np.sort(antennas), axis=0)

                run_parallel_mapping(uniq_ant_pairs)
            else:
                print(f'{mapping_folder} already exists')

    def make_uvw(self):
        """
        Fill UVW data points
        """

        # Make baseline/time mapping
        self.make_mapping_lst()

        def process_baselines(baseline_indices, baselines, mslist):
            """Process baselines parallel executor"""
            results = []
            for b_idx in baseline_indices:
                baseline = baselines[b_idx]
                c = 0
                # uvw = np.memmap('total_uvw.dat', dtype=np.float32, shape=(1, 3), mode='w+')
                uvw = np.zeros((0, 3))
                # time = np.memmap('total_time.dat', dtype=np.float64, shape=1, mode='w+')
                time = np.array([])
                row_idxs = []
                for ms_idx, ms in enumerate(sorted(mslist)):
                    mappingfolder = ms + '_baseline_mapping'
                    try:
                        mapjson = json.load(open(mappingfolder + '/' + '-'.join([str(a) for a in baseline]) + '.json'))
                    except FileNotFoundError:
                        c += 1
                        continue

                    row_idxs += list(mapjson.values())
                    uvw = np.append(np.memmap(f'{ms}_uvw.tmp.dat', dtype=np.float32).reshape((-1, 3))[
                        [int(i) for i in list(mapjson.keys())]], uvw, axis=0)

                    time = np.append(np.memmap(f'{ms}_time.tmp.dat', dtype=np.float64)[[int(i) for i in list(mapjson.keys())]], time)

                results.append((list(np.unique(row_idxs)), uvw, b_idx, time))
            return results

        # Get baselines
        ants = table(self.outname + "::ANTENNA", ack=False)
        baselines = np.c_[make_ant_pairs(ants.nrows(), 1)]
        ants.close()

        T = table(self.outname, readonly=False, ack=False)
        UVW = np.memmap('UVW.tmp.dat', dtype=np.float32, mode='w+', shape=(T.nrows(), 3))
        TIME = np.memmap('TIME.tmp.dat', dtype=np.float64, mode='w+', shape=(T.nrows()))
        TIME[:] = T.getcol("TIME")

        for ms_idx, ms in enumerate(sorted(self.mslist)):
            with table(ms, ack=False) as f:
                uvw = np.memmap(f'{ms}_uvw.tmp.dat', dtype=np.float32, mode='w+', shape=(f.nrows(), 3))
                time = np.memmap(f'{ms}_time.tmp.dat', dtype=np.float64, mode='w+', shape=(f.nrows()))

                uvw[:] = f.getcol("UVW")
                time[:] = mjd_seconds_to_lst_seconds(f.getcol("TIME"))

        # Determine number of workers
        num_workers = max(os.cpu_count()-5, 1)  # I/O-bound heuristic

        print(f"Using {num_workers} workers for making UVW column and accurate baseline mapping."
              f"\nThis is an expensive operation. So, be patient..")

        batch_size = max(1, len(baselines) // num_workers)  # Ensure at least one baseline per batch

        print("Multithreading...")

        with ThreadPoolExecutor(max_workers=num_workers) as executor:
            future_to_baseline = {
                executor.submit(process_baselines, range(i, min(i + batch_size, len(baselines))), baselines,
                                self.mslist): i
                for i in range(0, len(baselines), batch_size)
            }

            for future in as_completed(future_to_baseline):
                batch_start_idx = future_to_baseline[future]
                try:
                    results = future.result()
                    for row_idxs, uvws, b_idx, time in results:
                        UVW[row_idxs] = resample_uwv(uvws, row_idxs, time, TIME)
                except Exception as exc:
                    print(f'Batch starting at index {batch_start_idx} generated an exception: {exc}')

        UVW.flush()
        T.putcol("UVW", UVW)
        T.close()


        # Make final mapping
        self.make_mapping_uvw()

    def make_mapping_uvw(self):
        """
        Make mapping json files essential for efficient stacking based on UVW points
        """

        def process_baseline(baseline, mslist, UVW):
            """Parallel processing baseline"""
            try:
                folder = '/'.join(mslist[0].split('/')[0:-1])
                if not folder:
                    folder = '.'
                mapping_folder_baseline = sorted(
                    glob(folder + '/*_mapping/' + '-'.join([str(a) for a in baseline]) + '.json'))
                idxs_ref = np.unique(
                    [idx for mapp in mapping_folder_baseline for idx in json.load(open(mapp)).values()])
                uvw_ref = UVW[list(idxs_ref)]
                for mapp in mapping_folder_baseline:
                    idxs = [int(i) for i in json.load(open(mapp)).keys()]
                    ms = glob('/'.join(mapp.split('/')[0:-1]).replace("_baseline_mapping", ""))[0]
                    uvw_in = np.memmap(f'{ms}_uvw.tmp.dat', dtype=np.float32).reshape(-1, 3)[idxs]
                    idxs_new = [int(i) for i in np.array(idxs_ref)[
                        list(find_closest_index_multi_array(uvw_in[:, 0:2], uvw_ref[:, 0:2]))]]
                    with open(mapp, 'w+') as f:
                        json.dump(dict(zip(idxs, idxs_new)), f)
            except Exception as exc:
                print(f'Baseline {baseline} generated an exception: {exc}')

        # Get baselines
        ants = table(self.outname + "::ANTENNA", ack=False)
        baselines = np.c_[make_ant_pairs(ants.nrows(), 1)]
        ants.close()

        UVW = np.memmap('UVW.tmp.dat', dtype=np.float32).reshape(-1, 3)

        num_workers = min(os.cpu_count()-5, len(baselines))

        print('\nMake new mapping based on UVW points')
        with ThreadPoolExecutor(max_workers=num_workers) as executor:
            future_to_baseline = {executor.submit(process_baseline, baseline, self.mslist, UVW): baseline for baseline in
                                  baselines}
            for n, future in enumerate(as_completed(future_to_baseline)):
                baseline = future_to_baseline[future]
                try:
                    future.result()
                except Exception as exc:
                    print(f'Baseline {baseline} generated an exception: {exc}')
                print_progress_bar(n + 1, len(baselines))


class Stack:
    """
    Stack measurement sets in template empty.ms
    """
    def __init__(self, msin: list = None, outname: str = 'empty.ms', chunkmem: float = 1.):
        if not os.path.exists(outname):
            sys.exit(f"ERROR: Template {outname} has not been created or is deleted")
        print("\n\n==== Start stacking ====\n")
        self.template = table(outname, readonly=False, ack=False)
        self.mslist = msin
        self.outname = outname
        self.flag = False

        # Freq
        F = table(self.outname+'::SPECTRAL_WINDOW', ack=False)
        self.ref_freqs = F.getcol("CHAN_FREQ")[0]
        self.freq_len = self.ref_freqs.__len__()
        F.close()

        self.num_cpus = psutil.cpu_count(logical=True)
        total_memory = psutil.virtual_memory().total / (1024 ** 3)  # in GB
        target_chunk_size = total_memory / chunkmem
        self.chunk_size = min(int(target_chunk_size * (1024 ** 3) / np.dtype(np.float128).itemsize/2/self.freq_len), 500_000)
        print(f"\n---------------\nChunk size ==> {self.chunk_size}")


    def smooth_uvw(self):
        """
        Smooth UVW values (EXPERIMENTAL, CURRENTLY NOT USED)
        """

        uvw, _ = get_data_arrays('UVW', self.T.nrows())
        uvw[:] = self.T.getcol("UVW")
        time = self.T.getcol("TIME")

        ants = table(self.outname + "::ANTENNA", ack=False)
        baselines = np.c_[make_ant_pairs(ants.nrows(), 1)]
        ants.close()

        print('Smooth UVW')
        for idx_b, baseline in enumerate(baselines):
            print_progress_bar(idx_b, len(baselines))
            idxs = []
            for baseline_json in glob(f"*baseline_mapping/{baseline[0]}-{baseline[1]}.json"):
                idxs += list(load_json(baseline_json).values())
            sorted_indices = np.argsort(time[idxs])
            for i in range(3):
                uvw[np.array(idxs)[sorted_indices], i] = gaussian_filter1d(uvw[np.array(idxs)[sorted_indices], i], sigma=2)

        self.T.putcol('UVW', uvw)


    def stack_all(self, column: str = 'DATA'):
        """
        Stack all MS

        :param:
            - column: column name (currently only DATA)
        """

        def read_mapping(mapping_folder):
            """
            Read mapping with multi-threads
            """
            # Get the list of JSON files
            json_files = glob(os.path.join(mapping_folder, "*.json"))

            # Load JSON files in parallel
            total_map = {}
            with ThreadPoolExecutor() as executor:
                for result in executor.map(load_json, json_files):
                    total_map.update(result)

            # Convert keys and values to integers and sort
            total_map = {int(k): int(v) for k, v in total_map.items()}
            sorted_total_map = dict(sorted(total_map.items()))

            indices = list(sorted_total_map.keys())
            ref_indices = list(sorted_total_map.values())

            return indices, ref_indices

        if column == 'DATA':
            columns = ['UVW', column, 'WEIGHT_SPECTRUM']
        else:
            sys.exit("ERROR: Only column 'DATA' allowed (for now)")

        # Get template data
        self.T = table(os.path.abspath(self.outname), readonly=False, ack=False)

        # Loop over columns
        for col in columns:

            if col == 'UVW':
                new_data, uvw_weights = get_data_arrays(col, self.T.nrows(), self.freq_len)
            else:
                new_data, _ = get_data_arrays(col, self.T.nrows(), self.freq_len)

            # Loop over measurement sets
            for ms in self.mslist:

                print(f'\nStacking {col}: {ms}')

                # Open MS table
                if col == 'DATA':
                    t = taql(f"SELECT {col} * WEIGHT_SPECTRUM AS DATA_WEIGHTED FROM {os.path.abspath(ms)} ORDER BY TIME")
                elif col == 'UVW':
                    t = taql(f"SELECT {col},WEIGHT_SPECTRUM FROM {os.path.abspath(ms)} ORDER BY TIME")
                else:
                    t = taql(f"SELECT {col} FROM {os.path.abspath(ms)} ORDER BY TIME")

                # Get freqs offset
                if col != 'UVW':
                    f = table(ms+'::SPECTRAL_WINDOW', ack=False)
                    freqs = f.getcol("CHAN_FREQ")[0]
                    freq_idxs = find_closest_index_list(freqs, self.ref_freqs)
                    f.close()

                # Make antenna mapping in parallel
                mapping_folder = ms + '_baseline_mapping'

                print('Read mapping')
                indices, ref_indices = read_mapping(mapping_folder)

                # Chunked stacking!
                chunks = len(indices)//self.chunk_size + 1
                print(f'Stacking in {chunks} chunks')
                for chunk_idx in range(chunks):
                    print_progress_bar(chunk_idx, chunks+1)
                    data = t.getcol(col+"_WEIGHTED" if col == 'DATA' else col,
                                            startrow=chunk_idx * self.chunk_size, nrow=self.chunk_size)

                    row_idxs_new = ref_indices[chunk_idx * self.chunk_size:self.chunk_size * (chunk_idx+1)]
                    row_idxs = [int(i - chunk_idx * self.chunk_size) for i in
                                indices[chunk_idx * self.chunk_size:self.chunk_size * (chunk_idx+1)]]


                    if col == 'UVW':
                        new_data[row_idxs_new, :] = sum_arrays_chunkwise(new_data[row_idxs_new, :], data[row_idxs, :],
                                                                         chunk_size=self.chunk_size//self.num_cpus)

                        uvw_weights[row_idxs_new, :] = sum_arrays_chunkwise(uvw_weights[row_idxs_new, :], np.ones(uvw_weights[row_idxs_new, :].shape),
                                                                         chunk_size=self.chunk_size//self.num_cpus)
                        uvw_weights.flush()
                    else:
                        new_data[np.ix_(row_idxs_new, freq_idxs)] = sum_arrays_chunkwise(new_data[np.ix_(row_idxs_new, freq_idxs)], data[row_idxs, :],
                                                                         chunk_size=self.chunk_size//self.num_cpus)

                    new_data.flush()

                print_progress_bar(chunk_idx, chunks)
                t.close()

            print(f'Put column {col}')
            if col == 'UVW':
                uvw_weights[uvw_weights == 0] = 1
                new_data /= uvw_weights
                new_data[new_data != new_data] = 0.

            for chunk_idx in range(self.T.nrows()//self.chunk_size+1):
                self.T.putcol(col, new_data[chunk_idx * self.chunk_size:self.chunk_size * (chunk_idx+1)],
                              startrow=chunk_idx * self.chunk_size, nrow=self.chunk_size)

        self.T.close()

        # if self.flag:
        #     # ADD FLAG
        print(f'Put column FLAG')
        taql(f'UPDATE {self.outname} SET FLAG = (WEIGHT_SPECTRUM == 0)')
        # else: TODO: FIX FLAGGING
        #     # REMOVE FLAGS
        #     remove_flagged_entries(self.outname)

        # NORM DATA
        print(f'Normalise column DATA')
        taql(f'UPDATE {self.outname} SET DATA = (DATA / WEIGHT_SPECTRUM) WHERE ANY(WEIGHT_SPECTRUM > 0)')

        print("----------\n")


def clean_mapping_files(msin):
    """
    Clean-up mapping files
    """

    for ms in msin:
        rmtree(ms + '_baseline_mapping')


def clean_binary_files():
    """
    Clean-up binary files
    """

    for b in glob('*.tmp.dat'):
        run_command(f'rm {b}')


def check_folder_exists(folder_path):
    """
    Check if folder exists
    """
    return os.path.isdir(folder_path)


def plot_baseline_track(t_final_name: str = None, t_input_names: list = None, baseline='0-1', UV=True, saveas=None):
    """
    Plot baseline track

    :param:
        - t_final_name: table with final name
        - t_input_names: tables to compare with
        - mappingfiles: baseline mapping files
    """


    if len(t_input_names) > 4:
        sys.exit("ERROR: Can just plot 4 inputs")

    colors = ['red', 'green', 'yellow', 'black']

    if not UV:
        print("MAKE UW PLOT")

    ant1, ant2 = baseline.split('-')
    plt.close()

    for n, t_input_name in enumerate(t_input_names):
        print(t_input_name)

        ref_stats, ref_ids = get_station_id(t_final_name)
        new_stats, new_ids = get_station_id(t_input_name)

        id_map = dict(zip([ref_stats.index(a) for a in new_stats], new_ids))

        print(ref_stats[int(float(ant1))], ref_stats[int(float(ant2))])

        with table(t_final_name, ack=False) as f:
            fsub = f.query(f'ANTENNA1={ant1} AND ANTENNA2,={ant2} AND NOT ALL(WEIGHT_SPECTRUM == 0)', columns='UVW')
            uvw1 = fsub.getcol("UVW")

        with table(t_input_name, ack=False) as f:
            fsub = f.query(f'ANTENNA1={id_map[int(ant1)]} AND ANTENNA2={id_map[int(ant2)]} AND NOT ALL(WEIGHT_SPECTRUM == 0)', columns='UVW')
            uvw2 = fsub.getcol("UVW")

        # Scatter plot for uvw1
        if n == 0:
            lbl = 'Final dataset'
        else:
            lbl = None

        if uvw2.ndim>1:

            plt.scatter(uvw1[:, 0], uvw1[:, 2] if UV else uvw1[:, 3], label=lbl, color='blue', edgecolor='black', alpha=0.2, s=130, marker='o')

            # Scatter plot for uvw2
            plt.scatter(uvw2[:, 0], uvw2[:, 2] if UV else uvw2[:, 3], label=f'Dataset {n}', color=colors[n], edgecolor='black', alpha=0.7, s=70, marker='*')


    # Adding labels and title
    plt.xlabel("U (m)", fontsize=14)
    plt.ylabel("V (m)" if UV else "W (m)", fontsize=14)

    # Adding grid
    plt.grid(True, linestyle='--', alpha=0.6)

    # Adding legend
    plt.legend(fontsize=12)

    plt.tight_layout()

    if saveas is None:
        plt.show()
    else:
        plt.savefig(saveas, dpi=150)
        plt.close()


def make_baseline_uvw_plots(tabl, mslist):
    """
    Make baseline plots
    """

    run_command('mkdir -p baseline_plots')

    ants = table(tabl + "::ANTENNA", ack=False)
    baselines = np.c_[make_ant_pairs(ants.nrows(), 1)]
    ants.close()

    for baseline in baselines:
        bl = '-'.join([str(a) for a in baseline])
        plot_baseline_track(tabl, sorted(mslist), bl, saveas=f'baseline_plots/{bl}.png')


def parse_args():
    """
    Parse input arguments
    """

    parser = ArgumentParser(description='Sidereal visibility averaging')
    parser.add_argument('msin', nargs='+', help='Measurement sets to combine')
    parser.add_argument('--msout', type=str, default='empty.ms', help='Measurement set output name')
    parser.add_argument('--time_res', type=float, help='Desired time resolution in seconds')
    parser.add_argument('--resolution', type=float, help='Desired spatial resolution (if given, you also have to give --fov_diam)')
    parser.add_argument('--fov_diam', type=float, help='Desired field of view diameter in degrees. This is used to calculate the optimal time resolution.')
    parser.add_argument('--record_time', action='store_true', help='Record wall-time of stacking')
    parser.add_argument('--chunk_mem', type=float, default=1., help='Additional memory chunk parameter (larger for smaller chunks)')
    parser.add_argument('--no_dysco', action='store_true', help='No Dysco compression of data')
    parser.add_argument('--make_only_template', action='store_true', help='Stop after making empty template')
    parser.add_argument('--keep_mapping', action='store_true', help='Do not remove mapping files')
    parser.add_argument('--plot_uv_baseline_coverage', action='store_true', help='make plots with baseline versus UV')

    return parser.parse_args()


def ms_merger():
    """
    Main function
    """

    # Make template
    args = parse_args()

    # Verify if output exists
    if check_folder_exists(args.msout):
        sys.exit(f"ERROR: {args.msout} already exists! Delete file first if you want to overwrite.")

    avg = 1
    if args.time_res is not None:
        avg = 1
        time_res = args.time_res
        print(f"Use time resolution {time_res} seconds")
    elif args.resolution is not None and args.fov_diam is not None:
        time_res = time_resolution(args.resolution, args.fov_diam)
        print(f"Use time resolution {time_res} seconds")
    elif args.resolution is not None or args.fov_diam is not None:
        sys.exit("ERROR: if --resolution given, you also have to give --fov_diam, and vice versa.")
    else:
        avg = 2
        time_res = None
        print(f"Additional time sampling factor {avg}\n")

    t = Template(args.msin, args.msout)
    t.make_template(overwrite=True, time_res=time_res, avg_factor=avg)
    t.make_uvw()
    print("\n############\nTemplate creation completed\n############")

    # Stack MS
    if not args.make_only_template:
        if args.record_time:
            start_time = time.time()
        s = Stack(args.msin, args.msout, chunkmem=args.chunk_mem)
        s.stack_all()
        if args.record_time:
            end_time = time.time()
            elapsed_time = end_time - start_time
            print(f"Elapsed time for stacking: {elapsed_time//60} minutes")

    if args.plot_uv_baseline_coverage:
        make_baseline_uvw_plots(args.msout, args.msin)

    # Clean up mapping files
    if not args.keep_mapping:
        clean_mapping_files(args.msin)
    clean_binary_files()

    # Apply dysco compression
    if not args.no_dysco:
        compress(args.msout)


if __name__ == '__main__':
    ms_merger()
