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

    G_new = np.zeros(G.shape[0:-1]) + (4,)
    G_new[..., 0] += RR
    G_new[..., 1] += RL
    G_new[..., 2] += LR
    G_new[..., 3] += LL
    return G_new


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

    G_new = np.zeros(G.shape[0:-1]) + (4,)
    G_new[..., 0] += XX
    G_new[..., 1] += XY
    G_new[..., 2] += YX
    G_new[..., 3] += YY
    return G_new

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

def make_template(h5_in, soltab):
    """
    Make template of the Gain matrix with only ones
    :param h5_in: h5 file input
    :param soltab: solution table (phase, amplitude)
    """
    G, axes_vals = 0, {}
    for ss in h5_in.getSolsetNames():
        for st in h5_in.getSolset(ss).getSoltabNames():
            solutiontable = h5_in.getSolset(ss).getSoltab(st)
            if soltab in st:
                try:
                    if 'pol' in solutiontable.getAxesNames():
                        values = reorderAxes(solutiontable.getValues()[0], solutiontable.getAxesNames(), axes_names)
                        G = np.ones(values.shape)
                    else:
                        values = reorderAxes(solutiontable.getValues()[0], solutiontable.getAxesNames(), axes_names[0:-1])
                        G = np.ones(values.shape+(2,))
                except:
                    sys.exit('ERROR:\nReceived '+str(solutiontable.getAxesNames())+', but expect at least [time, freq, ant, dir] or [time, freq, ant, dir, pol]')
            axes_vals = {'time': solutiontable.getAxisValues('time'),
                         'freq': solutiontable.getAxisValues('freq'),
                         'ant': solutiontable.getAxisValues('ant'),
                         'dir': solutiontable.getAxisValues('dir')}
            if G.shape[-1]==2:
                axes_vals.update({'pol': ['XX', 'YY']})
            elif G.shape[-1]==4:
                axes_vals.update({'pol': ['XX', 'XY', 'YX', 'YY']})
    return G, axes_vals


parser = ArgumentParser()
parser.add_argument('-h5out', '--h5_file_out', type=str, help='h5 output name')
parser.add_argument('-h5in', '--h5_file_in', type=str, help='h5 input name')
parser.add_argument('--lin2circ', type=str2bool, nargs='?', const=True, default=True, help='convert linear to circular')
args = parser.parse_args()

"""
C = np.matrix([[1, 1j], [1, -1j]]) / np.sqrt(2)  ----> jones matrix
C.H is the Hermitian Transpose
"""

h5_in = h5parm(args.h5_file_in, readonly=True)
h5_out = h5parm(args.h5_file_out, readonly=False)

axes_names = ['time', 'freq', 'ant', 'dir', 'pol']

# MAKE TEMPLATE
G, axes_vals = make_template(h5_in, 'phase')
if G==0:
    make_template(h5_in, 'amplitude')

for ss in h5_in.getSolsetNames():

    solsetout = h5_out.makeSolset(ss)
    solsetin = h5_in.getSolset(ss)

    solsetout.obj.source.append(solsetout.obj.source[:])

    for st in h5_in.getSolset(ss).getSoltabNames():
        solutiontable = h5_in.getSolset(ss).getSoltab(st)
        if 'phase' in st:
            if 'pol' in solutiontable.getAxesNames():
                values = reorderAxes(solutiontable.getValues()[0], solutiontable.getAxesNames(), axes_names)
                G *= np.exp(values*1j)
            else:
                values = reorderAxes(solutiontable.getValues()[0], solutiontable.getAxesNames(), axes_names[0:-1])
                G *= np.exp(add_polarization(values, 2)*1j)

        elif 'amplitude' in st:
            if 'pol' in solutiontable.getAxesNames():
                values = reorderAxes(solutiontable.getValues()[0], solutiontable.getAxesNames(), axes_names)
                G *= values*1j
            else:
                values = reorderAxes(solutiontable.getValues()[0], solutiontable.getAxesNames(), axes_names[0:-1])
                G *= add_polarization(values, 2)

        elif 'tec' in st:
            tec_axes_names = [ax for ax in axes_names if solutiontable.getAxesNames()]
            tec = reorderAxes(solutiontable.getValues()[0], solutiontable.getAxesNames(), tec_axes_names)
            if 'freq' in solutiontable.getAxesNames():
                axes_vals = {'time': solutiontable.getAxisValues('time'),
                             'freq': solutiontable.getAxisValues('freq'),
                             'ant': solutiontable.getAxisValues('ant'),
                             'dir': solutiontable.getAxisValues('dir')}
            else:
                axes_vals = {'dir': solutiontable.getAxisValues('dir'),
                             'ant': solutiontable.getAxisValues('ant'),
                             'time': solutiontable.getAxisValues('time')}
            if 'pol' in solutiontable.getAxesNames():
                if tec.shape[-1]==2:
                    axes_vals.update({'pol': ['XX', 'YY']})
                elif tec.shape[-1]==4:
                    axes_vals.update({'pol': ['XX', 'XY', 'YX', 'YY']})
            solsetout.makeSoltab('tec', axesNames=tec_axes_names, axesVals=axes_vals, vals=tec,
                                 weights=np.ones(tec.shape))

    if args.lin2circ:
        G_new = lin2circ(G)
    else:
        G_new = circ2lin(G)

    phase = np.angle(G_new)
    amplitude = np.abs(G_new)

    solsetout.makeSoltab('phase', axesNames=axes_names, axesVals=axes_vals, vals=phase,
                         weights=np.ones(phase.shape))

    solsetout.makeSoltab('amplitude', axesNames=axes_names, axesVals=axes_vals, vals=amplitude,
                         weights=np.ones(amplitude.shape))

h5_in.close()
h5_out.close()