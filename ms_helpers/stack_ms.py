"""
LOFAR UV STACKER

Strategy:
    1) Make a template using the 'default_ms' option from casacore.tables (Template class).
       The template inclues all baselines, frequency, and smallest time spacing from all input MS.
       Time is converted to Local Sidereal Time (LST).

    2) Stack measurement sets on the template (Stack class).
        The stacking function maps baseline numbers from the input MS to the template MS and fills
        the template according to this mapping and in LST time.

#TODO: Work in progress!
1) Test on bigger MS
2) Validate issues at sub-arcsecond
3) Improve speed

"""

from casacore.tables import table, tablecopy, default_ms
import numpy as np
import os
import shutil
import sys
from astropy.time import Time
from astropy.coordinates import EarthLocation
import astropy.units as u
from pprint import pprint
from argparse import ArgumentParser


def print_progress_bar(index, total, bar_length=50):
    """
    Prints a progress bar to the console.

    :input:
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

    :input:
        - ms: measurement set
    """

    t = table(ms, readonly=True)
    dysco = t.getdesc()["DATA"]['dataManagerGroup'] == 'DyscoData'
    t.close()

    if dysco:
        return True
    else:
        return False


def decompress(ms):
    """
    running DP3 to remove dysco compression

    :input:
        - ms: measurement set
    """

    if is_dysco_compressed(ms):

        print('\n----------\nREMOVE DYSCO COMPRESSION\n----------\n')

        if os.path.exists(f'{ms}.tmp'):
            shutil.rmtree(f'{ms}.tmp')
        os.system(f"DP3 msin={ms} msout={ms}.tmp steps=[]")
        print('----------')
        return ms + '.tmp'

    else:
        return ms


def compress(ms):
    """
    running DP3 to apply dysco compression

    :input:
        - ms: measurement set
    """

    if not is_dysco_compressed(ms):

        print('\n----------\nDYSCO COMPRESSION\n----------\n')
        os.system(f"DP3 msin={ms} msout={ms} steps=[] msout.overwrite=true")
        print('----------')
        return ms

    else:
        return ms


def get_ms_content(ms):
    """
    Get MS content

    :input:
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
    # time_lst = T.getcol("TIME")
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

    return stations, lofar_stations, channels, dfreq, total_time_seconds, dt, time_min_lst, time_max_lst


def get_station_id(ms):
    """
    Get station with corresponding id number

    :input:
        - ms: measurement set

    :return:
        - antenna names, IDs
    """

    t = table(ms+'/ANTENNA', ack=False)
    ants = t.getcol("NAME")
    t.close()

    t = table(ms+'/FEED', ack=False)
    ids = t.getcol("ANTENNA_ID")
    t.close()

    return ants, ids


def mjd_seconds_to_lst_seconds(mjd_seconds, longitude_deg=52.909):
    """
    Convert time in modified Julian Date time to LST

    :input:
        - mjd_seconds: modified Julian date time in seconds
        - longitde_deg: longitude telescope in degrees (52.909 for LOFAR)

    :return:
        time in LST
    """

    # Convert seconds to days for MJD
    mjd_days = mjd_seconds / 86400.0

    # Create an astropy Time object
    time_utc = Time(mjd_days, format='mjd', scale='utc')

    # Define the observer's location using longitude (latitude doesn't affect LST)
    location = EarthLocation(lon=longitude_deg * u.deg, lat=0 * u.deg)

    # Calculate LST in hours
    lst_hours = time_utc.sidereal_time('mean', longitude=location.lon).hour

    # Convert LST from hours to seconds
    lst_seconds = lst_hours * 3600.0

    return lst_seconds


def same_phasedir(mslist: list = None):
    """
    Have MS same phase center?

    :input:
        - mslist: measurement set list
    """

    for n, ms in enumerate(mslist):
        t = table(ms+'/FIELD', ack=False)
        if n==0:
            phasedir = t.getcol("PHASE_DIR")
        else:
            if not np.all(phasedir == t.getcol("PHASE_DIR")):
                sys.exit("MS do not have the same phase center, check "+ms)


def sort_station_list(zipped_list):
    """
    Sorts a list of lists (or tuples) based on the first element of each inner list or tuple,
    which is necessary for a zipped list with station names and positions

    :input:
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

    :input:
        - n_antennas: number of antennas

    :return: number of baselines
    """

    return n_antennas * (n_antennas - 1) // 2


def make_ant_pairs(n_ant, n_time):
    """
    Generate ANTENNA1 and ANTENNA2 arrays for an array with M antennas over N time slots.

    :input:
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

    :input:
        - original_list: The original list to be transformed.
        - repeat_count: The number of times each element should be repeated.

    :return:
        - A new list where each element from the original list is repeated.
    """
    return np.array([element for element in original_list for _ in range(repeat_count)])


def find_closest_index(arr, value):
    """
    Find the index of the closest value in the array to the given value.

    :input:
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


def map_array_dict(arr, dct):
    """
    Maps elements of the input_array to new values using a mapping_dict

    :input:
        - input_array: numpy array of integers that need to be mapped.
        - mapping_dict: dictionary where each key-value pair represents an original value and its corresponding new value.

    :return:
        - An array where each element has been mapped according to the mapping_dict.
        (If an element in the input array does not have a corresponding mapping in the mapping_dict, it will be left unchanged)
    """

    # Vectorize a simple lookup function
    lookup = np.vectorize(lambda x: dct.get(x, x))

    # Apply this vectorized function to the input array
    output_array = lookup(arr)

    return output_array


class Template:
    """Make template measurement set based on input measurement sets"""
    def __init__(self, msin: list = None, outname: str = 'empty.ms'):
        self.mslist = msin
        self.outname = outname

    def add_spectral_window(self):
        """
        Add SPECTRAL_WINDOW as sub table
        """

        print("\n----------\nADD " + self.outname + "/SPECTRAL_WINDOW\n----------\n")

        tnew_spw_tmp = table(self.ref_table.getkeyword('SPECTRAL_WINDOW'), ack=False)
        newdesc = tnew_spw_tmp.getdesc()
        for col in ['CHAN_WIDTH', 'CHAN_FREQ', 'RESOLUTION', 'EFFECTIVE_BW']:
            newdesc[col]['shape'] = np.array([self.channels.shape[-1]])
        pprint(newdesc)

        tnew_spw = table(self.outname + '/SPECTRAL_WINDOW', newdesc, readonly=False, ack=False)
        tnew_spw.addrows(1)
        chanwidth = np.expand_dims([np.squeeze(np.diff(self.channels))[0]]*self.chan_num, 0)
        tnew_spw.putcol("NUM_CHAN", np.array([self.chan_num]))
        tnew_spw.putcol("CHAN_FREQ", self.channels)
        tnew_spw.putcol("CHAN_WIDTH", chanwidth)
        tnew_spw.putcol("RESOLUTION", chanwidth)
        tnew_spw.putcol("EFFECTIVE_BW", chanwidth)
        tnew_spw.putcol("REF_FREQUENCY", np.mean(self.channels))
        tnew_spw.putcol("MEAS_FREQ_REF", np.array([5])) #TODO: Why always 5?
        tnew_spw.putcol("TOTAL_BANDWIDTH", [np.max(self.channels)-np.min(self.channels)-chanwidth[0][0]])
        tnew_spw.putcol("NAME", 'Stacked_MS_'+str(int(np.mean(self.channels)//1000000))+"MHz")
        tnew_spw.flush(True)
        tnew_spw.close()
        tnew_spw_tmp.close()

    def add_stations(self):
        """
        Add ANTENNA and FEED tables
        """

        print("\n----------\nADD " + self.outname + "/ANTENNA\n----------\n")

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
        pprint(newdesc)
        tnew_ant_tmp.close()

        tnew_ant = table(self.outname + '/ANTENNA', newdesc, readonly=False, ack=False)
        tnew_ant.addrows(len(stations))
        print(len(stations))
        print('Number of stations: ' + str(tnew_ant.nrows()))
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

        print("\n----------\nADD " + self.outname + "/FEED\n----------\n")

        tnew_ant_tmp = table(self.ref_table.getkeyword('FEED'), ack=False)
        newdesc = tnew_ant_tmp.getdesc()
        pprint(newdesc)
        tnew_ant_tmp.close()

        tnew_feed = table(self.outname + '/FEED', newdesc, readonly=False, ack=False)
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

        print("\n----------\nADD " + self.outname + "/LOFAR_ANTENNA_FIELD\n----------\n")

        tnew_ant_tmp = table(self.ref_table.getkeyword('LOFAR_ANTENNA_FIELD'), ack=False)
        newdesc = tnew_ant_tmp.getdesc()
        pprint(newdesc)
        tnew_ant_tmp.close()

        tnew_field = table(self.outname + '/LOFAR_ANTENNA_FIELD', newdesc, readonly=False, ack=False)
        tnew_field.addrows(len(stations))
        tnew_field.putcol("ANTENNA_ID", np.array(range(len(stations))))
        tnew_field.putcol("NAME", names)
        tnew_field.putcol("COORDINATE_AXES", np.array(coor_axes))
        tnew_field.putcol("TILE_ELEMENT_OFFSET", np.array(tile_element))
        tnew_field.putcol("TILE_ROTATION", np.array([0]*len(stations)))
        tnew_field.flush(True)
        tnew_field.close()

        print("\n----------\nADD " + self.outname + "/LOFAR_STATION\n----------\n")

        tnew_ant_tmp = table(self.ref_table.getkeyword('LOFAR_STATION'), ack=False)
        newdesc = tnew_ant_tmp.getdesc()
        pprint(newdesc)
        tnew_ant_tmp.close()

        tnew_station = table(self.outname + '/LOFAR_STATION', newdesc, readonly=False, ack=False)
        tnew_station.addrows(len(lofar_names))
        tnew_station.putcol("NAME", lofar_names)
        tnew_station.putcol("FLAG_ROW", np.array([False] * len(lofar_names)))
        tnew_station.putcol("CLOCK_ID", np.array(clock))
        tnew_station.flush(True)
        tnew_station.close()

    def make_template(self, overwrite: bool =True):
        """
        Make template MS based on existing MS
        """

        if overwrite:
            if os.path.exists(self.outname):
                shutil.rmtree(self.outname)

        same_phasedir(self.mslist)

        # Get data columns
        unique_stations = []
        unique_lofar_stations = []
        unique_channels = []
        for k, ms in enumerate(self.mslist):
            stations, lofar_stations, channels, dfreq, total_time_seconds, dt, min_t, max_t = get_ms_content(ms)
            if k == 0:
                min_t_lst = min_t
                min_dt = dt
                dfreq_min = dfreq
                max_t_lst = max_t
            else:
                min_t_lst = min(min_t_lst, min_t)
                min_dt = max(min_dt, dt)
                dfreq_min = min(dfreq_min, dfreq)
                max_t_lst = max(max_t_lst, max_t)

            unique_stations += list(stations)
            unique_channels += list(channels)
            unique_lofar_stations += list(lofar_stations)
        self.station_info = unique_station_list(unique_stations)
        self.lofar_stations_info = unique_station_list(unique_lofar_stations)

        chan_range = np.arange(min(unique_channels), max(unique_channels)+dfreq_min, dfreq_min)
        self.channels = np.sort(np.expand_dims(np.unique(chan_range), 0))
        self.chan_num = self.channels.shape[-1]
        time_range = np.arange(min_t_lst, max_t_lst+min_dt, min_dt)
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
        newdesc_data.pop("_keywords_")

        pprint(newdesc_data)

        # Make main table
        default_ms(self.outname, newdesc_data)
        tnew = table(self.outname, readonly=False, ack=False)
        tnew.addrows(nrows)
        ant1, ant2 = make_ant_pairs(len(self.station_info), len(time_range))
        t = repeat_elements(time_range, baseline_count)
        tnew.putcol("TIME", t)
        tnew.putcol("TIME_CENTROID", t)
        tnew.putcol("WEIGHT", np.ones((nrows, 4)).astype(np.float32))
        tnew.putcol("SIGMA", np.ones((nrows, 4)).astype(np.float32))
        tnew.putcol("ANTENNA1", ant1)
        tnew.putcol("ANTENNA2", ant2)
        tnew.putcol("EXPOSURE", np.array([np.diff(time_range)[0]] * nrows))
        tnew.putcol("FLAG_ROW", np.array([False] * nrows))
        tnew.putcol("INTERVAL", np.array([np.diff(time_range)[0]] * nrows))
        tnew.putcol("FLAG", np.zeros((nrows, self.chan_num, 4)).astype(bool))
        tnew.putcol("WEIGHT_SPECTRUM", np.ones((nrows, self.chan_num, 4)).astype(np.float32))

        tnew.flush(True)
        tnew.close()

        # Set SPECTRAL_WINDOW info
        self.add_spectral_window()

        # Set ANTENNA/STATION info
        self.add_stations()

        # Set other tables
        for subtbl in ['FIELD', 'HISTORY', 'FLAG_CMD', 'DATA_DESCRIPTION',
                       'LOFAR_ELEMENT_FAILURE', 'OBSERVATION', 'POINTING',
                       'POLARIZATION', 'PROCESSOR', 'STATE']:
            print("----------\nADD " + self.outname + "/" + subtbl + "\n----------")

            tsub = table(self.tmpfile+"/"+subtbl, ack=False, readonly=False)
            tsub.copy(self.outname + '/' + subtbl, deep=True)
            tsub.flush(True)
            tsub.close()

        self.ref_table.close()

        # Cleanup
        if 'tmp' in self.tmpfile:
            shutil.rmtree(self.tmpfile)


class Stack:
    """
    Stack measurement sets in template empty.ms
    """
    def __init__(self, msin: list = None, outname: str = 'empty.ms'):
        if not os.path.exists(outname):
            sys.exit(f"ERROR: Template {outname} has not been created or is deleted")
        self.template = table(outname, readonly=False, ack=False)
        self.mslist = msin
        self.outname = outname
        self.ant_map = None

    def stack_all(self, column: str = 'DATA'):
        """
        Stack all MS

        :input:
            - type: DATA, WEIGHT, WEIGHT_SPECTRUM, WEIGHT
        """

        print(column)

        ref_stats, ref_ids = get_station_id(self.outname)
        T = table(self.outname, ack=False)
        F = table(self.outname+'/SPECTRAL_WINDOW', ack=False)
        ref_freqs = F.getcol("CHAN_FREQ")[0]
        F.close()

        ref_time = T.getcol("TIME")
        ref_uniq_time = np.unique(ref_time)
        ref_antennas = np.c_[T.getcol("ANTENNA1"), T.getcol("ANTENNA2")]

        # Initialize data
        if column in ['DATA', 'WEIGHT_SPECTRUM', 'FLAG', 'WEIGHT']:
            weights = np.zeros((len(ref_time), len(ref_freqs)))
        elif column in ['UVW']:
            weights = np.zeros((len(ref_time)))

        if column in ['DATA', 'WEIGHT_SPECTRUM']:
            if column == 'DATA':
                dtp = np.complex128
            elif column == 'WEIGHT_SPECTRUM':
                dtp = np.float32
            else:
                dtp = np.float32
            new_data = np.zeros((len(ref_time), len(ref_freqs), 4), dtype=dtp)
        elif column == 'FLAG':
            new_data = np.ones((len(ref_time), len(ref_freqs), 4), dtype=bool)

        elif column == 'WEIGHT':
            new_data = np.zeros((len(ref_time), len(ref_freqs)), dtype=np.float32)

        elif column == 'UVW':
            new_data = np.zeros((len(ref_time), 3), dtype=np.float32)

        # Remove ref_time
        ref_time = None

        for ms in self.mslist:
            new_stats, new_ids = get_station_id(ms)
            id_map = dict(zip(new_ids, [ref_stats.index(a) for a in new_stats]))
            t = table(ms, ack=False)
            f = table(ms+'/SPECTRAL_WINDOW', ack=False)
            freqs = f.getcol("CHAN_FREQ")[0]
            f.close()

            # Mapped antenna pairs to same as ref (template) table
            antennas = np.c_[map_array_dict(t.getcol("ANTENNA1"), id_map), map_array_dict(t.getcol("ANTENNA2"), id_map)]

            # Unique antenna pairs
            uniq_ant_pairs = np.array(make_ant_pairs(len(np.unique(t.getcol("ANTENNA1")))+1, 1)).T

            # Time in LST
            time = mjd_seconds_to_lst_seconds(t.getcol("TIME"))
            uniq_time = np.unique(time)
            time_offset = find_closest_index(ref_uniq_time, uniq_time[0])
            # time_offset_last = find_closest_index(ref_uniq_time, uniq_time[-1])

            freq_offset = find_closest_index(ref_freqs, freqs[0])

            # Loop over frequency
            print('\nStacking: '+ms)

            # Loop over antenna pairs
            for m, antpair in enumerate(uniq_ant_pairs):
                print_progress_bar(m, len(uniq_ant_pairs))
                pair_idx = np.squeeze(np.argwhere(np.all(antennas == antpair, axis=1)))
                ref_pair_idx = np.squeeze(np.argwhere(np.all(ref_antennas == antpair, axis=1))[time_offset:len(pair_idx)])
                idx_len = min(len(pair_idx), len(ref_pair_idx))
                if column in ['DATA', 'WEIGHT_SPECTRUM']:
                    new_data[ref_pair_idx[0:idx_len], freq_offset:freq_offset+len(freqs), :] += t.getcol(column)[pair_idx[0:idx_len], :, :]
                    weights[ref_pair_idx[0:idx_len], freq_offset:freq_offset + len(freqs)] += 1
                elif column == 'FLAG':
                    new_data[ref_pair_idx[0:idx_len], freq_offset:freq_offset+len(freqs), :] *= t.getcol(column)[pair_idx[0:idx_len], :, :]
                elif column == 'WEIGHT':
                    new_data[ref_pair_idx[0:idx_len], :] *= t.getcol(column)[pair_idx[0:idx_len], :]
                elif column == 'UVW':
                    new_data[ref_pair_idx[0:idx_len], :] += t.getcol(column)[pair_idx[0:idx_len], :]
                    weights[ref_pair_idx[0:idx_len]] += 1

            t.close()

        # Average
        weights = weights.astype(np.float16)
        weights[weights == 0.] = np.inf
        if column in ['DATA', 'WEIGHT_SPECTRUM']:
            for p in range(4):
                new_data[:, :, p] /= weights
        if column == 'UVW':
            new_data /= np.repeat(weights, 3).reshape(-1, 3)

        T.putcol(column, new_data)

        print("----------\n")

        T.close()


def parse_args():
    """
    Parse input arguments
    """
    parser = ArgumentParser(description='MS stacking')
    parser.add_argument('--msout', type=str, default='empty.ms', help='Measurement set output name')
    parser.add_argument('msin', nargs='+', help='Measurement sets to stack')
    return parser.parse_args()


def main():
    """
    Main function
    """

    # Make template
    args = parse_args()
    t = Template(args.msin, args.msout)
    t.make_template()
    print("############\nTemplate creation completed\n############")

    # Stack MS
    s = Stack(args.msin, args.msout)
    s.stack_all('UVW')
    s.stack_all('DATA')
    s.stack_all('FLAG')
    s.stack_all('WEIGHT_SPECTRUM')

    # Apply dysco compression
    compress(args.msout)


if __name__ == '__main__':
    main()
