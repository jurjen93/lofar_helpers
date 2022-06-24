import matplotlib.pyplot as plt
import tables
import numpy as np
from scipy.signal import medfilt
import os
from argparse import ArgumentParser

tables.file._open_files.close_all()

parser = ArgumentParser()
parser.add_argument('--antennas', nargs='+', help='Antenna indices', required=True)
parser.add_argument('--no_plot', action='store_true', help='Do not make plot at the end', default=None)
parser.add_argument('--h5', type=str, help='h5 table', required=True )
args = parser.parse_args()

change_antennas = [int(a) for a in args.antennas]
# change_antennas = [56, 57, 58, 73] # indices of antennas to fix

def overwrite_val(T, table, values, title=None):
    ss = T.root._f_get_child('calibrator')
    ss._f_get_child(table).val._f_remove()
    T.create_array(ss._f_get_child(table), 'val', values)

def flatlist(l):
    return [x for xs in l for x in xs]

newfile = args.h5.replace(".h5", "_smooth.h5")

os.system(f'cp {args.h5} {newfile}')
H = tables.open_file(newfile, 'r+')

freq = H.root.calibrator.bandpass.freq[:]

new_bandpass = np.zeros(H.root.calibrator.bandpass.val.shape)

for bb_pol in range(len(H.root.calibrator.bandpass.pol)):
    for time in range(len(H.root.calibrator.bandpass.time)):
        for antenna in range(len(H.root.calibrator.bandpass.ant)):
            bp = H.root.calibrator.bandpass.val[time, :, antenna, bb_pol]
            if antenna in change_antennas:
                bp_new = bp.copy()
                outliers = np.abs(np.subtract(bp,medfilt(bp, 55)))/np.nanmean(medfilt(bp, 55))>30/np.nanmean(medfilt(bp, 55))
                outliers = sorted(set(flatlist([[j for j in list(range(m - 7, m + 7)) if j<len(outliers)] for m in [i[0] for i in np.argwhere(outliers)]])))
                bp_new[outliers] = medfilt(bp, 55)[outliers]
                new_bandpass[time, :, antenna, bb_pol] = bp_new
            else:
                new_bandpass[time, :, antenna, bb_pol] = bp

overwrite_val(H, 'bandpass', new_bandpass)

if not args.no_plot:
    T = tables.open_file(args.h5)

    fig = plt.figure(figsize=(10, 10))
    rows = len(change_antennas)//2
    columns = 2
    grid = plt.GridSpec(rows, columns, wspace=.25, hspace=.25)
    for i in range(rows*columns):
        exec (f"plt.subplot(grid{[i]})")
        plt.scatter(freq/1000000, T.root.calibrator.bandpass.val[0, :, change_antennas[i], 0], s=1)
        plt.scatter(freq/1000000, H.root.calibrator.bandpass.val[0, :, change_antennas[i], 0], s=1)
        plt.xlabel('freq [MHz]', fontsize=16)
        plt.ylabel('amplitude', fontsize=16)
        plt.xticks(size=14)
        plt.yticks(size=14)
        plt.legend(['old', 'new'], fontsize=14)
        plt.title('Antenna index '+str(change_antennas[i]), size=14)
    plt.show()

    T.close()


H.close()