import tables
import os
import numpy as np
from argparse import ArgumentParser
from scipy.constants import speed_of_light

circ2lin_math = """
Convert circular polarization to linear polarization
-----------------------------
XX = RR + RL + LR + LL
XY = iRR - iRL + iLR - iLL
YX = -iRR - iRL + iLR + iLL
YY = RR - RL - LR + LL
-----------------------------
"""

class Template:
    """
    Make template based on given h5parm file
    """
    def __init__(self, name_in, name_out, freqs = None):
        os.system(' '.join(['cp', name_in, name_out]))
        self.h5 = tables.open_file(name_out, 'r+')
        self.axes = ['time', 'freq', 'ant', 'dir', 'pol']
        if freqs is not None:
            self.freqs = freqs
        else:
            self.freqs = self.h5.root.sol000.phase000.freq[:]

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

    def make_template(self, shape = None, polrot = None):
        """
        Make template h5

        :param shape: shape of values and weights solution table
        :param polrot: make rotation matrix to align polarization
        """

        print("1) MAKE TEMPLATE H5PARM")

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
                    shape[1] = len(self.freqs)


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
                self.update_array(st, self.freqs, 'freq')

        return self

    def rotate(self, intercept, rotation_measure):
        """
        Rotate angle by the following matrix:
         /e^(i*rho)  0 \
        |              |
         \ 0         1/

        :param intercept: intercept
        :param rotation_measure: rotation measure in rad/m^2
        """
        print('2) ADD PHASE ROTATION')

        phaserot = intercept+rotation_measure*(speed_of_light/self.freqs)**2

        mapping = list(zip(list(self.freqs), list(phaserot)))
        print('########################\nFrequency to rotation in radian (circular base):\n------------------------')
        for element in mapping:
            print(str(int(element[0]))+'Hz --> '+str(round(element[1], 3))+'rad')

        # print('Rotate with rotation angle: '+str(intercept) + ' radian')
        for solset in self.h5.root._v_groups.keys():
            ss = self.h5.root._f_get_child(solset)
            for soltab in ss._v_groups.keys():
                if 'phase' in soltab:
                    phaseval = ss._f_get_child(soltab).val[:]
                    # same phaserot for all antennas
                    for ant_idx in range(phaseval.shape[2]):
                        phaseval[0, :, ant_idx, 0, 0] += phaserot
                    self.update_array(ss._f_get_child(soltab), phaseval, 'val')
        print("########################")

        self.circ2lin()

        return self

    def circ2lin(self):
        """
        Convert circular polarization to linear polarization

        XX = RR + RL + LR + LL
        XY = iRR - iRL + iLR - iLL
        YX = -iRR - iRL + iLR + iLL
        YY = RR - RL - LR + LL

        :return: linear polarized solutions
        """

        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'+circ2lin_math+'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

        G = self.h5.root.sol000.amplitude000.val[:] * np.exp(self.h5.root.sol000.phase000.val[:] * 1j)


        XX = (G[..., 0] + G[..., -1])
        XY = 1j * (G[..., 0] - G[..., -1])
        YX = 1j * (G[..., -1] - G[..., 0])
        YY = (G[..., 0] + G[..., -1])

        # XX += (G[..., 2] + G[..., 1])
        # XY += 1j * (G[..., 2] - G[..., 1])
        # YX += 1j * (G[..., 2] - G[..., 1])
        # YY -= (G[..., 1] + G[..., 2])

        XX /= 2
        XY /= 2
        YX /= 2
        YY /= 2

        G_new = np.zeros(G.shape[0:-1] + (4,)).astype(np.complex128)

        G_new[..., 0] += XX
        G_new[..., 1] += XY
        G_new[..., 2] += YX
        G_new[..., 3] += YY

        G_new = np.where(abs(G_new) < 10 * np.finfo(float).eps, 0, G_new)

        phase = np.angle(G_new)
        amplitude = abs(G_new)

        self.update_array(self.h5.root.sol000.phase000, phase, 'val')
        self.update_array(self.h5.root.sol000.amplitude000, amplitude, 'val')

        return G_new


if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument('--h5_in', type=str, help='input name', required=True)
    parser.add_argument('--h5_out', type=str, help='output name', required=True)
    parser.add_argument('--intercept', type=float, help='intercept for rotation angle')
    parser.add_argument('--RM', type=float, help='rotation measure')
    args = parser.parse_args()

    temp = Template(args.h5_in, args.h5_out)
    if args.intercept is None:
        temp.make_template()
    else:
        temp.make_template(polrot=True)
        temp.rotate(intercept=args.intercept, rotation_measure=args.RM)

    temp.h5.close()
