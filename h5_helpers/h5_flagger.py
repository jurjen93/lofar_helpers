import tables
import numpy as np
from argparse import ArgumentParser

if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument('--h5', help='Input h5parm', required=True)
    parser.add_argument('--freqrange', help='Frequency range in MHz (example: 140-158', required=True)
    parser.add_argument('--ant', help='Antenna (example: RS409HBA', required=True)

    args = parser.parse_args()

    inputh5 = ''
    min_freq, max_freq = [float(f)*1e6 for f in args.freqrange.split('-')]
    ant = args.ant

    print(f'Flagging from {min_freq}Hz to {max_freq}Hz for antenna {ant}')

    H = tables.open_file(args.h5, 'r+')
    for solset in H.root._v_groups.keys():
        print(solset)
        ss = H.root._f_get_child(solset)
        for soltab in ss._v_groups.keys():
            print(soltab)
            st = ss._f_get_child(soltab)
            weights = st.weight[:]
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
                print(weights[:,freq_indices[0]:freq_indices[-1], ant_idx, ...].shape, newweights.shape)
                H.root._f_get_child(solset)._f_get_child(soltab).\
                weight[:,freq_indices[0]:freq_indices[-1], ant_idx, ...] = newweights

    H.close()


