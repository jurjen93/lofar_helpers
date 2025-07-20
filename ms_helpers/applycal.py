#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import tables
from subprocess import call
from argparse import ArgumentParser
from numpy import pi
from sys import exit
from os.path import abspath, basename

__author__ = "Jurjen de Jong"


class ApplyCal:
    def __init__(self, msin: str = None, h5s: str = None, msincol: str = "DATA", msoutcol: str = "CORRECTED_DATA",
                 msout: str = '.', bitrate: int = 10):
        """
        Apply calibration solutions

        :param msin: input measurement set
        :param h5s: solution file(s) to apply
        :param msincol: input column
        :param msoutcol: output column
        :param msout: output measurement set
        :param bitrate: data bitrate
        """

        self.cmd = ['DP3', 'msin=' + abspath(msin)]
        self.cmd += ['msout=' + msout]
        self.cmd += ['msin.datacolumn=' + msincol]
        if msout == '.':
            self.cmd += ['msout.datacolumn=' + msoutcol]
        if bitrate>0:
            self.cmd += [f'msout.storagemanager=dysco msout.storagemanager.databitrate={bitrate}']

        steps = []

        for n, h5 in enumerate(h5s):

            h5 = abspath(h5)

            poldim_num = self.poldim_num(h5)

            # non-scalar
            if poldim_num>1:
                steps.append(f'beam_dir_{n}')
                with tables.open_file(h5) as T:
                    dir = (T.root.sol000.source[:]['dir'][0] * 360 / 2 / pi) % 360  # convert to degree
                    self.cmd += [f'beam_dir_{n}.type=applybeam', f'beam_dir_{n}.direction=[{round(dir[0], 5)}deg,{round(dir[1], 5)}deg]',
                                 f'beam_dir_{n}.updateweights=True']

            # fulljones
            if poldim_num==4:
                steps.append(f'ac_{n}')
                self.cmd += [f'ac_{n}.type=applycal',
                             f'ac_{n}.parmdb=' + h5,
                             f'ac_{n}.correction=fulljones',
                             f'ac_{n}.soltab=[amplitude000,phase000]',
                             f'ac_{n}.updateweights=True',
                             f'ac_{n}.missingantennabehavior=unit']

            # add non-fulljones solutions apply
            else:
                ac_count = 0
                with tables.open_file(h5) as T:
                    for corr in T.root.sol000._v_groups.keys():
                        self.cmd += [f'ac{ac_count}_{n}.type=applycal',
                                     f'ac{ac_count}_{n}.parmdb={h5}',
                                     f'ac{ac_count}_{n}.correction={corr}',
                                     f'ac{ac_count}_{n}.missingantennabehavior=unit',
                                     f'ac{ac_count}_{n}.updateweights=True']
                        steps.append(f'ac{ac_count}_{n}')
                        ac_count += 1

            # non-scalar
            if poldim_num>1:
                # this step inverts the beam at the infield and corrects beam at phase center
                steps.append(f'beam_center_{n}')
                self.cmd += [f'beam_center_{n}.type=applybeam', f'beam_center_{n}.direction=[]',
                             f'beam_center_{n}.updateweights=True']

        self.cmd += ['steps=' + str(steps).replace(" ", "").replace("\'", "")]

    @staticmethod
    def poldim_num(h5: str = None):
        """
        Verify if file is fulljones

        :param h5: h5 file
        """
        with tables.open_file(h5) as T:
            soltab = list(T.root.sol000._v_groups.keys())[0]
            if 'pol' in T.root.sol000._f_get_child(soltab).val.attrs["AXES"].decode('utf8'):
                return T.root.sol000._f_get_child(soltab).pol[:].shape[0]
            else:
                return 0

    def print_cmd(self):
        """Print DP3 command"""
        print('\n'.join(self.cmd))
        return self

    def run(self):
        """Run DP3 command"""
        retval = call(' '.join(self.cmd), shell=True)
        if retval != 0:
            print('FAILED to run ' + ' '.join(self.cmd) + ': return value is ' + str(retval))
            raise Exception(' '.join(self.cmd))
        return retval


def parse_args():
    """Argument parser"""

    parser = ArgumentParser(description='Apply calibration solutions by taking into account beam order corrections.')
    parser.add_argument('msin', nargs='+', type=str, help='Input MeasurementSet(s)')
    parser.add_argument('--msout', type=str, default='.', help='Output MeasurementSet')
    parser.add_argument('--h5', nargs='+', type=str, help='h5parm calibration solution files', required=True)
    parser.add_argument('--colin', type=str, default='DATA', help='Input column name')
    parser.add_argument('--colout', type=str, default="CORRECTED_DATA", help='Output column name')
    parser.add_argument('--bitrate', type=int, help='Number of bits per float used for columns containing visibilities. '
                                                    'Can be set to zero to compress weights only.', default=10)
    return parser.parse_args()


def main():
    """Main function"""

    args = parse_args()

    if len(args.msin) == 1:
        Ac = ApplyCal(msin=args.msin[0], h5s=args.h5, msincol=args.colin, msoutcol=args.colout, msout=args.msout, bitrate=args.bitrate)
    elif len(args.h5)==1:
        for ms in args.msin:
            Ac = ApplyCal(msin=ms, h5s=args.h5, msincol=args.colin, msoutcol=args.colout, msout='applycal_' + basename(ms), bitrate=args.bitrate)
    else:
        exit("ERROR: cannot give multiple MeasurementSets and multiple h5parms.")
    Ac.print_cmd()
    Ac.run()


if __name__ == '__main__':
    main()
