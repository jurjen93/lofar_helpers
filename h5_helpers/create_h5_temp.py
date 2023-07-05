import tables
import os
import numpy as np
from argparse import ArgumentParser


class Template:
    """
    Make template based on given h5parm file
    """
    def __init__(self, name_in, name_out):
        os.system(' '.join(['cp', name_in, name_out]))
        self.h5 = tables.open_file(name_out, 'r+')
        self.axes = ['time', 'freq', 'ant', 'dir', 'pol']

    def update_array(self, st, new_val, arrayname):
        """
        Update array

        :param st: soltab
        :param new_val: new values
        :param array: array name (val, weight, pol, dir, or freq)
        """
        valtype = str(st._f_get_child(arrayname).dtype)
        st._f_get_child(arrayname)._f_remove()
        if 'float' in str(valtype):
            if '16' in valtype:
                atomtype = tables.Float16Atom()
            elif '32' in valtype:
                atomtype = tables.Float32Atom()
            elif '64' in valtype:
                atomtype = tables.Float64Atom()
            else:
                atomtype = tables.Float64Atom()
            self.h5.create_array(st, arrayname, new_val.astype(valtype), atom=atomtype)
        else:
            self.h5.create_array(st, arrayname, new_val.astype(valtype))
        if arrayname=='val' or arrayname=='weight':
            st._f_get_child(arrayname).attrs['AXES'] = bytes(','.join(self.axes), 'utf-8')

        return self

    def make_template(self, shape: tuple = None, polrot: bool = None, freqs = None):
        """
        Make template h5

        :param shape: shape of values and weights solution table
        :param polrot: make rotation matrix to align polarization
        """

        print("Make Template h5parm")

        for solset in self.h5.root._v_groups.keys():
            ss = self.h5.root._f_get_child(solset)
            for soltab in ss._v_groups.keys():
                st = ss._f_get_child(soltab)

                if shape is None and polrot is None:
                    shape = st.val[:].shape
                if polrot is not None:
                    shape = list(st.val[:].shape)
                    shape[0]=1
                    shape[-1]=4
                    if freqs is not None:
                        shape[1] = len(freqs)


                if 'phase' in soltab:
                    new_val = np.zeros(shape)
                elif 'amplitude' in soltab:
                    new_val = np.ones(shape)
                    new_val[..., 1] = 0
                    new_val[..., 2] = 0
                else:
                    continue

                self.update_array(st, new_val, 'val')
                self.update_array(st, np.ones(shape), 'weight')
                self.update_array(st, np.array([st.time[:][0]]), 'time')
                self.update_array(st, np.array(['XX', 'XY', 'YX', 'YY']), 'pol')
                if freqs is not None:
                    self.update_array(st, freqs, 'freq')

        return self

    def rotate(self, rotation_angle):
        """
        Rotate angle by the following matrix:
         /e^(i*rho)  0 \
        |              |
         \ 0         1/

        :param rotation_angle: rotation_angle in radian
        """

        print('Rotate with rotation angle: '+str(rotation_angle) + ' radian')
        for solset in self.h5.root._v_groups.keys():
            ss = self.h5.root._f_get_child(solset)
            for soltab in ss._v_groups.keys():
                if 'phase' in soltab:
                    phaseval = ss._f_get_child(soltab).val[:]
                    phaseval[..., 0] += rotation_angle
                    self.update_array(ss._f_get_child(soltab), phaseval, 'val')

        return self


if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument('--h5_in', type=str, help='input name', required=True)
    parser.add_argument('--h5_out', type=str, help='output name', required=True)
    parser.add_argument('--pol_rotang', type=float, help='polarization rotation angle')
    args = parser.parse_args()

    temp = Template(args.h5_in, args.h5_out)
    if args.pol_rotang is None:
        temp.make_template()
    else:
        temp.make_template(polrot=True)
        temp.rotate(rotation_angle=args.pol_rotang)

    temp.h5.close()
