import tables
import numpy as np
from argparse import ArgumentParser
import re

if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument('--h5', nargs='+', help='Input h5parm', required=True)
    parser.add_argument('--freqrange', help='Flag frequency range in MHz (example: 140-158')
    parser.add_argument('--ant', help='Antenna to flag (example: RS409HBA)')
    parser.add_argument('--ampflag', help='Flag based on amplitude values', action='store_true')
    args = parser.parse_args()

    if args.freqrange is not None:
        min_freq, max_freq = [float(f)*1e6 for f in args.freqrange.split('-')]
        ant = args.ant
        print(f'Flagging from {min_freq}Hz to {max_freq}Hz for antenna {ant}')

    if args.ampflag is not None:
        ampflag = args.ampflag

    for h5 in args.h5:
        print('-----\n'+h5+'\n-----')
        H = tables.open_file(h5, 'r+')
        for solset in H.root._v_groups.keys():
            print(solset)
            ss = H.root._f_get_child(solset)
            for soltab in ss._v_groups.keys():
                st = ss._f_get_child(soltab)
                weights = st.weight[:]

                if args.freqrange is not None:
                    freqs = st.freq[:]
                    ants = [a.decode('utf8') for a in st.ant[:]]
                    ant_idx = ants.index(ant)
                    AXES = st.weight.attrs["AXES"].decode('utf8').split(',')

                    freq_indices = np.where((freqs > min_freq) & (freqs < max_freq))
                    shape = list(weights.shape)
                    shape[AXES.index('freq')] = len(freq_indices)
                    del shape[AXES.index('ant')]
                    newweights = np.zeros(shape)

                    if AXES.index('ant')==2 and AXES.index('freq')==1:
                        print(weights[:, freq_indices[0]:freq_indices[-1], ant_idx, ...].shape, newweights.shape)
                        H.root._f_get_child(solset)._f_get_child(soltab).\
                        weight[:, freq_indices[0]:freq_indices[-1], ant_idx, ...] = newweights

                if args.ampflag is not None and 'amplitude' in soltab:
                    values = st.val[:]
                    flag_vals = (values < 3) & (values > 0.33)
                    flag_percentage = 1-(flag_vals*weights).sum()/weights.sum()
                    print(f'Flagged {round(flag_percentage*100, 4)}%')
                    del values
                    ss._f_get_child(soltab).weight[:] *= flag_vals

                    phasenum = re.findall(r'\d+', soltab)[0]
                    if 'phase'+phasenum in list(ss._v_groups.keys()):
                        ss._f_get_child('phase'+phasenum).weight[:] *= flag_vals


        H.close()


