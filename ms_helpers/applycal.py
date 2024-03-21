"""
Use this script to apply solutions in the correct order with beam corrections, based on input from an h5 solution file.
"""

import tables
from subprocess import call
from argparse import ArgumentParser
from numpy import pi
import sys


class ApplyCal:
    def __init__(self, msin: str = None, h5: str = None, msincol: str = "DATA", msoutcol: str = "CORRECTED_DATA",
                 msout: str = '.', dysco: bool = True):
        """
        Apply calibration solutions

        :param msin: input measurement set
        :param h5: solution file to apply
        :param msincol: input column
        :param msoutcol: output column
        :param msout: output measurement set
        :param dysco: compress with dysco
        """

        if len(msin) > 1 and msout != '.':
            sys.exit("ERROR: --msout necessary and unequal to '.' if more than one --msin")

        self.cmd = ['DP3', 'msin=' + msin]
        self.cmd += ['msout=' + msout]
        self.cmd += ['msin.datacolumn=' + msincol]
        if msout == '.':
            self.cmd += ['msout.datacolumn=' + msoutcol]
        if dysco:
            self.cmd += ['msout.storagemanager=dysco']

        steps = []

        steps.append('beam_dir')
        T = tables.open_file(h5)
        dir = (T.root.sol000.source[:]['dir'][0] * 360 / 2 / pi) % 360  # convert to degree
        self.cmd += ['beam_dir.type=applybeam', f'beam_dir.direction=[{round(dir[0], 5)}deg,{round(dir[1], 5)}deg]',
                     'beam_dir.updateweights=True']
        T.close()

        if self.isfulljones(h5):
            steps.append('ac')
            self.cmd += ['ac.type=applycal',
                         'ac.parmdb=' + h5,
                         'ac.correction=fulljones',
                         'ac.soltab=[amplitude000,phase000]']

        # add non-fulljones solutions apply
        else:
            ac_count = 0
            T = tables.open_file(h5)
            for corr in T.root.sol000._v_groups.keys():
                self.cmd += [f'ac{ac_count}.type=applycal',
                             f'ac{ac_count}.parmdb={h5}',
                             f'ac{ac_count}.correction={corr}']
                steps.append(f'ac{ac_count}')
                ac_count += 1

        # this step inverts the beam at the infield and corrects beam at phase center
        steps.append('beam_center')
        self.cmd += ['beam_center.type=applybeam', 'beam_center.direction=[]',
                     'beam_center.updateweights=True']
        self.cmd += ['steps=' + str(steps).replace(" ", "").replace("\'", "")]

    @staticmethod
    def isfulljones(h5: str = None):
        """
        Verify if file is fulljones

        :param h5: h5 file
        """
        T = tables.open_file(h5)
        soltab = list(T.root.sol000._v_groups.keys())[0]
        if 'pol' in T.root.sol000._f_get_child(soltab).val.attrs["AXES"].decode('utf8'):
            if T.root.sol000._f_get_child(soltab).pol[:].shape[0] == 4:
                T.close()
                return True
        T.close()
        return False

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

    parser = ArgumentParser(description='Applycal on MS with H5')
    parser.add_argument('--msin', nargs='+', type=str, help='input measurement set', required=True)
    parser.add_argument('--msout', type=str, default='.', help='output measurement set')
    parser.add_argument('--h5', type=str, help='h5 calibration', required=True)
    parser.add_argument('--colin', type=str, default='DATA', help='input column name')
    parser.add_argument('--colout', type=str, default=None, help='output column name')

    return parser.parse_args()


def main():
    """Main function"""

    args = parse_args()

    if len(args.msin) == 1:
        Ac = ApplyCal(msin=args.msin[0], h5=args.h5, msincol=args.colin, msoutcol=args.colout, msout=args.msout)
    else:
        for ms in args.msin:
            Ac = ApplyCal(msin=ms, h5=args.h5, msincol=args.colin, msoutcol=args.colout, msout='applycal_' + ms)
    Ac.print_cmd()
    Ac.run()


if __name__ == '__main__':
    main()
