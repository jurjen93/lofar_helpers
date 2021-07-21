"""
LAST UPDATE: 20-7-2021

This script is used to convert from circular to linear polarization and vice versa.
"""

import numpy as np
from argparse import ArgumentParser, ArgumentTypeError
from losoto.h5parm import h5parm
from losoto.lib_operations import reorderAxes
import sys

__author__ = "Jurjen de Jong (jurjendejong@strw.leidenuniv.nl)"

def str2bool(v):
    v = str(v)
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise ArgumentTypeError('Boolean value expected.')

class PolChange:
    """
    This Object is to convert polarization from linear to circular or vice versa.
    """
    def __init__(self, h5_in, h5_out):
        self.h5_in = h5parm(h5_in, readonly=True)
        self.h5_out = h5parm(h5_out, readonly=False)
        self.axes_names = ['time', 'freq', 'ant', 'dir', 'pol']

    @staticmethod
    def lin2circ(G):
        """Convert linear polarization to circular polarization"""
        RR = (G[..., 0] + G[..., -1]).astype(np.complex128)
        LL = (G[..., 0] + G[..., -1]).astype(np.complex128)
        RL = (G[..., 0] - G[..., -1]).astype(np.complex128)
        LR = (G[..., 0] - G[..., -1]).astype(np.complex128)

        if G.shape[-1] == 4:
            RR += 1j * (G[..., 2] - G[..., 1])
            LL += 1j * (G[..., 1] - G[..., 2])
            RL += 1j * (G[..., 2] + G[..., 1])
            LR -= 1j * (G[..., 2] + G[..., 1])

        RR /= 2
        LL /= 2
        RL /= 2
        LR /= 2

        G_new = np.zeros(G.shape[0:-1] + (4,)).astype(np.complex128)
        G_new[..., 0] += RR
        G_new[..., 1] += RL
        G_new[..., 2] += LR
        G_new[..., 3] += LL
        return G_new

    @staticmethod
    def circ2lin(G):
        """Convert circular polarization to linear polarization"""
        XX = (G[..., 0] + G[..., -1]).astype(np.complex128)
        YY = (G[..., 0] + G[..., -1]).astype(np.complex128)
        XY = 1j * (G[..., 0] - G[..., -1]).astype(np.complex128)
        YX = 1j * (G[..., -1] - G[..., 0]).astype(np.complex128)

        if G.shape[-1] == 4:
            XX += (G[..., 2] + G[..., 1]).astype(np.complex128)
            YY -= (G[..., 1] + G[..., 2]).astype(np.complex128)
            XY += 1j * (G[..., 2] - G[..., 1])
            YX += 1j * (G[..., 2] - G[..., 1])

        XX /= 2
        YY /= 2
        XY /= 2
        YX /= 2

        G_new = np.zeros(G.shape[0:-1] + (4,)).astype(np.complex128)
        G_new[..., 0] += XX
        G_new[..., 1] += XY
        G_new[..., 2] += YX
        G_new[..., 3] += YY
        return G_new

    @staticmethod
    def add_polarization(values, dim_pol):
        """
        Add extra polarization if there is no polarization
        :param values: values which need to get a polarization
        :param dim_pol: number of dimensions
        """
        values_temp = np.ones(values.shape+(dim_pol,))
        for i in range(dim_pol):
            values_temp[..., i] = values

        return values_temp


    def make_template(self, soltab):
        """
        Make template of the Gain matrix with only ones
        :param h5_in: h5 file input
        :param soltab: solution table (phase, amplitude)
        """
        self.G, self.axes_vals = np.array([]), {}
        for ss in self.h5_in.getSolsetNames():
            for st in self.h5_in.getSolset(ss).getSoltabNames():
                solutiontable = self.h5_in.getSolset(ss).getSoltab(st)
                if soltab in st:
                    try:
                        if 'pol' in solutiontable.getAxesNames():
                            values = reorderAxes(solutiontable.getValues()[0], solutiontable.getAxesNames(), self.axes_names)
                            self.G = np.ones(values.shape).astype(np.complex128)
                        else:
                            values = reorderAxes(solutiontable.getValues()[0], solutiontable.getAxesNames(), self.axes_names[0:-1])
                            self.G = np.ones(values.shape+(2,)).astype(np.complex128)
                    except:
                        sys.exit('ERROR:\nReceived '+str(solutiontable.getAxesNames())+', but expect at least [time, freq, ant, dir] or [time, freq, ant, dir, pol]')

                    self.axes_vals = {'time': solutiontable.getAxisValues('time'),
                                 'freq': solutiontable.getAxisValues('freq'),
                                 'ant': solutiontable.getAxisValues('ant'),
                                 'dir': solutiontable.getAxisValues('dir'),
                                 'pol': ['XX', 'XY', 'YX', 'YY']}
                    break

        print('Shape of input {shape}'.format(shape=self.G.shape))
        return self

    def add_tec(self, solutiontable):
        tec_axes_names = [ax for ax in self.axes_names if solutiontable.getAxesNames()]
        tec = reorderAxes(solutiontable.getValues()[0], solutiontable.getAxesNames(), tec_axes_names)
        if 'freq' in solutiontable.getAxesNames():
            axes_vals_tec = {'time': solutiontable.getAxisValues('time'),
                             'freq': solutiontable.getAxisValues('freq'),
                             'ant': solutiontable.getAxisValues('ant'),
                             'dir': solutiontable.getAxisValues('dir')}
        else:
            axes_vals_tec = {'dir': solutiontable.getAxisValues('dir'),
                             'ant': solutiontable.getAxisValues('ant'),
                             'time': solutiontable.getAxisValues('time')}
        if 'pol' in solutiontable.getAxesNames():
            if tec.shape[-1] == 2:
                axes_vals_tec.update({'pol': ['XX', 'YY']})
            elif tec.shape[-1] == 4:
                axes_vals_tec.update({'pol': ['XX', 'XY', 'YX', 'YY']})
        axes_vals_tec = [v[1] for v in
                         sorted(axes_vals_tec.items(), key=lambda pair: self.axes_names.index(pair[0]))]
        self.solsetout.makeSoltab('tec', axesNames=tec_axes_names, axesVals=axes_vals_tec, vals=tec,
                             weights=np.ones(tec.shape))

    def make_new_gains(self, lin2circ):
        for ss in self.h5_in.getSolsetNames():

            self.solsetout = self.h5_out.makeSolset(ss)
            solsetin = self.h5_in.getSolset(ss)

            self.solsetout.obj.source.append(solsetin.obj.source[:])

            for st in self.h5_in.getSolset(ss).getSoltabNames():
                solutiontable = self.h5_in.getSolset(ss).getSoltab(st)
                print('Reading {st} from {ss}'.format(ss=ss, st=st))
                if 'phase' in st:
                    if 'pol' in solutiontable.getAxesNames():
                        values = reorderAxes(solutiontable.getValues()[0], solutiontable.getAxesNames(), self.axes_names)
                        self.G *= np.exp(values * 1j)
                    else:
                        values = reorderAxes(solutiontable.getValues()[0], solutiontable.getAxesNames(),
                                             self.axes_names[0:-1])
                        self.G *= np.exp(self.add_polarization(values, 2) * 1j)

                elif 'amplitude' in st:
                    if 'pol' in solutiontable.getAxesNames():
                        values = reorderAxes(solutiontable.getValues()[0], solutiontable.getAxesNames(), self.axes_names)
                        self.G *= values * 1j
                    else:
                        values = reorderAxes(solutiontable.getValues()[0], solutiontable.getAxesNames(),
                                             self.axes_names[0:-1])
                        self.G *= self.add_polarization(values, 2)

                elif 'tec' in st:
                    self.add_tec(solutiontable)
                else:
                    print("Didn't include {st} in this version yet".format(st=st))
                    print("Let me (Jurjen) know if you need to include this.")

            if lin2circ:
                print('Convert linear polarization to circular polarization')
                G_new = self.lin2circ(self.G)
            else:
                print('Convert circular polarization to linear polarization')
                G_new = self.circ2lin(self.G)
            print('Shape of output for amplitude and phase: {shape}'.format(shape=G_new.shape))

            phase = np.angle(G_new)
            amplitude = np.abs(G_new)

            self.axes_vals = [v[1] for v in sorted(self.axes_vals.items(), key=lambda pair: self.axes_names.index(pair[0]))]

            self.solsetout.makeSoltab('phase', axesNames=self.axes_names, axesVals=self.axes_vals, vals=phase,
                                 weights=np.ones(phase.shape))
            print('Created new phase solutions')

            self.solsetout.makeSoltab('amplitude', axesNames=self.axes_names, axesVals=self.axes_vals, vals=amplitude,
                                 weights=np.ones(amplitude.shape))
            print('Created new amplitude solutions')

        self.h5_in.close()
        self.h5_out.close()

        return self

if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument('-h5out', '--h5_file_out', type=str, help='h5 output name')
    parser.add_argument('-h5in', '--h5_file_in', type=str, help='h5 input name')
    parser.add_argument('--lin2circ', type=str2bool, nargs='?', const=True, default=True, help='convert linear to circular')
    args = parser.parse_args()

    """
    C = np.matrix([[1, 1j], [1, -1j]]) / np.sqrt(2)  ----> jones matrix
    C.H is the Hermitian Transpose
    """

    Pol = PolChange(h5_in=args.h5_file_in, h5_out=args.h5_file_out)

    Pol.make_template('phase')
    if len(Pol.G.shape)>1:
        Pol.make_template('amplitude')
        print('Using amplitude as template')
    else:
        print('Using phase as template')

    Pol.make_new_gains(args.lin2circ)