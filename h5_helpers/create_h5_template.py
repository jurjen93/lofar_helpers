import tables
import os
import numpy as np
from argparse import ArgumentParser


class Template:
    """
    Make template based on given h5parm file
    Currently it only sets the phases to 0 and the amplitudes to 1
    """
    def __init__(self, name_in, name_out):
        os.system(' '.join(['cp', name_in, name_out]))
        self.h5 = tables.open_file(name_out, 'r+')

    def make_phase_amplitude_dummies(self):
        for solset in self.h5.root._v_groups.keys():
            ss = self.h5.root._f_get_child(solset)
            for soltab in ss._v_groups.keys():
                st = ss._f_get_child(soltab)
                valtype = str(st._f_get_child('val').dtype)
                attrsaxes = st.val.attrs['AXES']
                if '16' in valtype:
                    atomtype = tables.Float16Atom()
                elif '32' in valtype:
                    atomtype = tables.Float32Atom()
                elif '64' in valtype:
                    atomtype = tables.Float64Atom()
                else:
                    atomtype = tables.Float64Atom()

                if 'phase' in soltab:
                    new_val = np.zeros(st.val[:].shape)
                elif 'amplitude' in soltab:
                    new_val = np.ones(st.val[:].shape)
                else:
                    continue

                st._f_get_child('val')._f_remove()
                self.h5.create_array(st, 'val', new_val.astype(valtype), atom=atomtype)
                st._f_get_child('val').attrs['AXES'] = attrsaxes

        return self


if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument('--input', type=str, help='input name', required=True)
    parser.add_argument('--output', type=str, help='output name', required=True)
    args = parser.parse_args()

    test = Template(args.input, args.output)
    test.make_phase_amplitude_dummies()
    test.h5.close()
