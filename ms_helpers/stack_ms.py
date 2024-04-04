"""
CASA for UVW:

fixvis(vis='template.ms', outputvis='fixed.ms', reuse=False)

"""


from casacore.tables import table, tablecopy
import numpy as np
import os
import shutil
import sys
from astropy.time import Time
from astropy.coordinates import EarthLocation
import astropy.units as u
from pprint import pprint


def decompress(ms):
    """
    running DP3 to remove dysco compression

    :input:
        - ms: measurement set
    """

    print('\n----------\nREMOVE DYSCO COMPRESSION (if dysco compressed)\n----------\n')

    if os.path.exists(f'{ms}.tmp'):
        shutil.rmtree(f'{ms}.tmp')
    os.system(f"DP3 msin={ms} msout={ms}.tmp steps=[]")
    print('----------')


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

    T = table(ms)
    F = table(ms+"::SPECTRAL_WINDOW")
    A = table(ms+"::ANTENNA")
    L = table(ms+"::LOFAR_ANTENNA_FIELD")
    S = table(ms+"::LOFAR_STATION")

    # ALL LOFAR STATIONS
    lofar_stations = list(zip(
                        S.getcol("NAME"),
                        S.getcol("CLOCK_ID")
                    ))

    # ALL STATIONS AND SUB STATIONS
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
    time = sorted(np.unique(T.getcol("TIME")))
    time_lst = mjd_seconds_to_lst_seconds(T.getcol("TIME"))
    time_min_lst, time_max_lst = time_lst.min(), time_lst.max()
    total_time_seconds = int(round(max(time)-min(time), 0))
    dt = int(round(np.diff(sorted(set(time)))[0], 0))

    print(f'\nCONTENT from {ms}\n'
          '----------\n'
          f'Stations: {", ".join([s[0] for s in lofar_stations])}\n'
          f'Number of channels: {chan_num}\n'
          f'Total time: {round(total_time_seconds/3600, 2)} hrs\n'
          f'Delta time: {dt} seconds\n'
          f'----------')

    S.close()
    L.close()
    T.close()
    F.close()
    A.close()

    return stations, lofar_stations, channels, total_time_seconds, dt, time_min_lst, time_max_lst


def mjd_seconds_to_lst_seconds(mjd_seconds, longitude_deg=52.909):
    """
    Convert time in modified Julian Date time to LST

    :param mjd_seconds: modified Julian date time in seconds
    :param longitde_deg: longitude telescope in degrees (52.909 for LOFAR)

    :return: LST
    """

    # Convert seconds to days for MJD
    mjd_days = mjd_seconds / 86400.0

    # Create an astropy Time object
    time_utc = Time(mjd_days, format='mjd', scale='utc')

    # Define the observer's location using longitude (latitude doesn't affect LST)
    location = EarthLocation(lon=longitude_deg * u.deg, lat=0*u.deg)

    # Calculate LST in hours
    lst_hours = time_utc.sidereal_time('mean', longitude=location.lon).hour

    # Convert LST from hours to seconds
    lst_seconds = lst_hours * 3600.0

    return lst_seconds


def same_phasedir(mslist: list = None):
    """
    Have MS same phase center?
    """

    for n, ms in enumerate(mslist):
        t = table(ms+'/FIELD')
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

    Returns:
        sorted list
    """
    return sorted(zipped_list, key=lambda item: item[0])


def unique_station_list(station_list):
    """
    Filters a list of stations only based on first element

    Parameters:
        - station_list (list of stations with positions): Stations to be filtered.

    Returns:
        - filtered list of stations
    """
    unique_dict = {}
    for item in station_list:
        if item[0] not in unique_dict:
            unique_dict[item[0]] = item
    return list(unique_dict.values())


def n_baselines(n_antennas: int = None):
    """
    Calculate number of baselines

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


def add_spectral_window(ref_table, template_name, channels):
    """
    Add SPECTRAL_WINDOW as sub table

    :input:
        - ref_table: empty table to pick information from
        - template_name: name of main table
        - channels: numpy array with channels
    """

    print("\n----------\nADD " + template_name + "/SPECTRAL_WINDOW\n----------\n")

    tnew_spw_tmp = table(ref_table.getkeyword('SPECTRAL_WINDOW'))
    newdesc = tnew_spw_tmp.getdesc()
    for col in ['CHAN_WIDTH', 'CHAN_FREQ', 'RESOLUTION']:
        newdesc[col]['shape'] = np.array([channels.shape[-1]])
    tnew_spw_tmp.close()

    tnew_spw = table(template_name + '/SPECTRAL_WINDOW', newdesc, readonly=False)
    tnew_spw.addrows(1)
    channum = channels.shape[-1]
    chanwidth = np.expand_dims([np.squeeze(np.diff(channels))[0]]*channum, 0)
    tnew_spw.putcol("NUM_CHAN", np.array([channum]))
    tnew_spw.putcol("CHAN_FREQ", channels)
    tnew_spw.putcol("CHAN_WIDTH", chanwidth)
    tnew_spw.putcol("RESOLUTION", chanwidth)
    tnew_spw.putcol("EFFECTIVE_BW", chanwidth)
    tnew_spw.putcol("REF_FREQUENCY", np.mean(channels))
    tnew_spw.putcol("MEAS_FREQ_REF", np.array([channum]))
    tnew_spw.putcol("TOTAL_BANDWIDTH", [np.max(channels)-np.min(channels)])
    tnew_spw.flush(True)
    tnew_spw.close()


def add_stations(ref_table, template_name, station_info, lofar_stations_info):
    """
    Add ANTENNA and FEED tables

    :input:
        - ref_table: empty table to pick information from
        - template_name: name of main table
        - station_info: station information
    """

    print("\n----------\nADD " + template_name + "/ANTENNA\n----------\n")

    stations = [sp[0] for sp in station_info]
    positions = np.array([sp[1] for sp in station_info])
    diameters = np.array([sp[2] for sp in station_info])
    ids = np.array([sp[3] for sp in station_info]).astype(np.float32)
    phase_ref = np.array([sp[4] for sp in station_info])
    names = np.array([sp[5] for sp in station_info])
    coor_axes = np.array([sp[6] for sp in station_info])
    tile_element = np.array([sp[7] for sp in station_info])
    lofar_names = np.array([sp[0] for sp in lofar_stations_info])
    clock = np.array([sp[1] for sp in lofar_stations_info])

    tnew_ant_tmp = table(ref_table.getkeyword('ANTENNA'))
    newdesc = tnew_ant_tmp.getdesc()
    tnew_ant_tmp.close()

    tnew_ant = table(template_name + '/ANTENNA', newdesc, readonly=False)
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

    print("\n----------\nADD " + template_name + "/FEED\n----------\n")

    tnew_ant_tmp = table(ref_table.getkeyword('FEED'))
    newdesc = tnew_ant_tmp.getdesc()
    tnew_ant_tmp.close()

    tnew_feed = table(template_name + '/FEED', newdesc, readonly=False)
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

    print("\n----------\nADD " + template_name + "/LOFAR_ANTENNA_FIELD\n----------\n")

    tnew_ant_tmp = table(ref_table.getkeyword('LOFAR_ANTENNA_FIELD'))
    newdesc = tnew_ant_tmp.getdesc()
    tnew_ant_tmp.close()

    tnew_field = table(template_name + '/LOFAR_ANTENNA_FIELD', newdesc, readonly=False)
    tnew_field.addrows(len(stations))
    tnew_field.putcol("ANTENNA_ID", np.array(range(len(stations))))
    tnew_field.putcol("NAME", names)
    tnew_field.putcol("COORDINATE_AXES", np.array(coor_axes))
    tnew_field.putcol("TILE_ELEMENT_OFFSET", np.array(tile_element))
    tnew_field.putcol("TILE_ROTATION", np.array([0]*len(stations)))
    tnew_field.flush(True)
    tnew_field.close()

    print("\n----------\nADD " + template_name + "/LOFAR_STATION\n----------\n")

    tnew_ant_tmp = table(ref_table.getkeyword('LOFAR_STATION'))
    newdesc = tnew_ant_tmp.getdesc()
    tnew_ant_tmp.close()

    tnew_station = table(template_name + '/LOFAR_STATION', newdesc, readonly=False)
    tnew_station.addrows(len(lofar_names))
    tnew_station.putcol("NAME", lofar_names)
    tnew_station.putcol("FLAG_ROW", np.array([False] * len(lofar_names)))
    tnew_station.putcol("CLOCK_ID", np.array(clock))
    tnew_station.flush(True)
    tnew_station.close()


def make_template(mslist: list = None):
    """
    Make template MS based on existing MS

    :input:
        - mslist: list of measurement sets

    """

    # Verify if
    same_phasedir(mslist)

    # Get data columns
    unique_stations = []
    unique_lofar_stations = []
    unique_channels = []
    max_dt = 0
    max_t_lst = 0
    for k, ms in enumerate(mslist):
        stations, lofar_stations, channels, total_time_seconds, dt, min_t, max_t = get_ms_content(ms)
        if k == 0:
            min_t_lst = min_t
        else:
            min_t_lst = min(min_t_lst, min_t)
        max_t_lst = max(max_t_lst, max_t)
        max_dt = max(max_dt, dt)
        unique_stations += list(stations)
        unique_channels += list(channels)
        unique_lofar_stations += list(lofar_stations)
    station_info = unique_station_list(unique_stations)
    lofar_stations_info = unique_station_list(unique_lofar_stations)

    channels = np.expand_dims(np.array(list(set(unique_channels))), 0)
    chan_num = channels.shape[-1]
    time_range = np.arange(min_t_lst, max_t_lst, max_dt)
    baseline_count = n_baselines(len(station_info))
    nrows = baseline_count*len(time_range)

    # Take one ms for temp usage
    tmp_ms = mslist[0]

    # Remove dysco compression
    decompress(tmp_ms)

    # Make empty copy (which later becomes the empty template output)
    outname = 'empty.ms'
    tablecopy(tmp_ms+".tmp", outname, valuecopy=True, copynorows=True)
    ref_table = table(outname)

    # Data description
    newdesc_data = ref_table.getdesc()

    # Reshape
    for col in ['DATA', 'FLAG', 'WEIGHT_SPECTRUM']:
        newdesc_data[col]['shape'] = np.array([chan_num, 4])

    # for key in newdesc_data['_keywords_'].keys():
    #     if key != 'MS_VERSION':
    #         newdesc_data['_keywords_'][key] = newdesc_data['_keywords_'][key].replace('empty.ms', template_name)

    print(pprint(newdesc_data))

    # MAKE MAIN TABLE
    tmp_name = 'tmp.ms'
    tnew = table(tmp_name, newdesc_data, nrow=nrows, _columnnames=ref_table.colnames())
    ant1, ant2 = make_ant_pairs(len(station_info), len(time_range))
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
    tnew.putcol("FLAG", np.zeros((nrows, chan_num, 4)).astype(bool))
    tnew.putcol("WEIGHT_SPECTRUM", np.ones((nrows, chan_num, 4)).astype(np.float32))

    tnew.flush(True)
    tnew.close()

    # SET SPECTRAL WINDOW INFO
    add_spectral_window(ref_table, tmp_name, channels)

    # SET ANTENNA INFO
    add_stations(ref_table, tmp_name, station_info, lofar_stations_info)

    # SET QUALITY_XXX_STATISTIC
    for subt in ['QUALITY_BASELINE_STATISTIC', 'QUALITY_FREQUENCY_STATISTIC',
                 'QUALITY_TIME_STATISTIC', 'QUALITY_KIND_NAME']:
        tablecopy(outname+'/'+subt, tmp_name+'/'+subt)

    # SET OTHER TABLES
    for subtbl in ['FIELD', 'HISTORY', 'FLAG_CMD', 'DATA_DESCRIPTION',
                   'LOFAR_ELEMENT_FAILURE', 'OBSERVATION', 'POINTING',
                   'POLARIZATION', 'PROCESSOR', 'STATE']:
        print("----------\nADD " + tmp_name + "/" + subtbl + "\n----------")

        tsub = table(tmp_ms+".tmp/"+subtbl, ack=False)
        tsub.copy(tmp_name + '/' + subtbl, deep=True)
        tsub.flush(True)
        tsub.close()

    ref_table.close()

    # CLEANUP
    shutil.rmtree(outname)
    shutil.rmtree(tmp_ms+".tmp")
    shutil.move(tmp_name, outname)


#TODO: 1) TEST MULTIPLE STACK 2) ADD MISSING BASELINES 3) ADDING UVW DATA