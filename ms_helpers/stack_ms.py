from casacore.tables import table, tablecopy
import numpy as np
import os
import shutil
import sys


def decompress(ms):
    """
    running DP3 to remove dysco compression

    :input:
        - ms: measurement set
    """

    if os.path.exists(f'{ms}_tmp'):
        shutil.rmtree(f'{ms}_tmp')
    os.system(f"DP3 msin={ms} msout={ms}_tmp steps=[]")


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

    stations = A.getcol("NAME")
    chan_num = F.getcol("NUM_CHAN")[0]
    channels = F.getcol("CHAN_FREQ")[0]
    time = sorted(np.unique(T.getcol("TIME")))
    total_time_seconds = int(round(max(time)-min(time), 0))
    dt = int(round(np.diff(sorted(set(time)))[0], 0))

    print(f'\nCONTENT from {ms}\n'
          '----------\n'
          f'Stations: {", ".join(stations)}\n'
          f'Number of channels: {chan_num}\n'
          f'Total time: {round(total_time_seconds/3600, 2)} hrs\n'
          f'Delta time: {dt} seconds\n'
          f'----------')
    T.close()
    F.close()
    A.close()

    return stations, channels, total_time_seconds, dt


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


def calculate_number_of_baselines(n_antennas: int = None):
    """
    Calculate number of baselines

    :input:
        - n_antennas: number of antennas

    :return: number of baselines
    """

    return n_antennas * (n_antennas - 1) // 2


def add_spectral_window(template_name, channels, desc):
    """
    Add SPECTRAL_WINDOW as sub table

    :input:
        - template_name: name of main table
        - channels: numpy array with channels
        - desc: sub table description
    """
    tnew_spw = table(template_name + '/SPECTRAL_WINDOW', desc)
    tnew_spw.addrows(1)
    channum = channels.shape[-1]
    chanwidth = np.expand_dims([np.squeeze(np.diff(channels))[0]]*channum, 0)
    tnew_spw.putcol("NUM_CHAN", [channum])
    tnew_spw.putcol("CHAN_FREQ", channels)
    tnew_spw.putcol("CHAN_WIDTH", chanwidth)
    tnew_spw.putcol("RESOLUTION", chanwidth)
    tnew_spw.putcol("EFFECTIVE_BW", chanwidth)
    tnew_spw.putcol("REF_FREQUENCY", np.mean(channels))
    tnew_spw.putcol("TOTAL_BANDWIDTH", [np.max(channels)-np.min(channels)])
    tnew_spw.flush(True)


def make_template(mslist: list = None, template_name: str = 'template.ms'):
    """
    Make template MS based on existing MS

    :input:
        - mslist: list of measurement sets

    """

    # verify if
    same_phasedir(mslist)

    # Get data columns
    unique_stations = []
    unique_channels = []
    max_dt = 0
    max_total_time_seconds = 0
    for ms in mslist:
        stations, channels, total_time_seconds, dt = get_ms_content(ms)
        max_total_time_seconds = max(max_total_time_seconds, total_time_seconds)
        max_dt = max(max_dt, dt)
        unique_stations += list(stations)
        unique_channels += list(channels)
    stations = np.array(list(set(unique_stations)))
    channels = np.expand_dims(np.array(list(set(unique_channels))), 0)
    chan_num = channels.shape[-1]

    time_range = np.arange(0, max_total_time_seconds, max_dt)
    baselines = calculate_number_of_baselines(len(channels))

    DATA = np.ones((baselines*len(time_range), len(channels), 4))
    FLAG = np.zeros((baselines*len(time_range), len(channels), 4)).astype(bool)
    UVW = np.zeros((baselines*len(time_range), 3))
    TIME = time_range
    nrows = baselines*len(time_range)

    # Take one ms for temp usage
    tmp_ms = mslist[0]

    # Remove dysco compression
    decompress(tmp_ms)

    # Make empty copy
    tablecopy(tmp_ms+"_tmp", 'empty.ms', valuecopy=True, copynorows=True)

    tnew_data_tmp = table('empty.ms', readonly=False)
    tnew_spw_tmp = table(tnew_data_tmp.getkeyword('SPECTRAL_WINDOW'), readonly=False)

    newdesc_data = tnew_data_tmp.getdesc()
    newdesc_spec = tnew_spw_tmp.getdesc()

    for col in ['DATA', 'FLAG', 'WEIGHT_SPECTRUM']:
        newdesc_data[col]['shape'] = np.array([chan_num, 4])

    for col in ['CHAN_WIDTH', 'CHAN_FREQ', 'RESOLUTION']:
        newdesc_spec[col]['shape'] = np.array([chan_num])

    # MAIN TABLE
    tnew = table(template_name, newdesc_data, nrow=nrows)
    tnew.flush(True)

    # SET SPECTRAL WINDOW
    add_spectral_window(template_name, channels, newdesc_spec)


    #TODO: LOFAR_ANTENNA_FIELD, ANTENNA, FEED, LOFAR_STATION
    for subtbl in ['FIELD', 'HISTORY', 'FLAG_CMD', 'DATA_DESCRIPTION',
                   'LOFAR_ELEMENT_FAILURE', 'OBSERVATION', 'POINTING',
                   'POLARIZATION', 'PROCESSOR', 'STATE']:
        tsub = table(tmp_ms+"_tmp/"+subtbl, ack=False)
        tsub.copy(template_name + '/' + subtbl, deep=True)
        tsub.close()

    tnew_data_tmp.close()
    tnew_spw_tmp.close()
