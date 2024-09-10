import tables
import os
from astropy.coordinates import SkyCoord
import numpy as np
import sys
from argparse import ArgumentParser
from casacore.tables import table


class FindClosestDir:
    """
    Find closest h5 direction in merged h5
    """
    def __init__(self, h5_in, template_name):
        os.system(' '.join(['cp', h5_in, template_name]))
        print(f'Created {template_name}')

        self.name_out = template_name
        self.h5_in = h5_in
        T = tables.open_file(h5_in)
        # get input directions
        self.dirs = T.root.sol000.source[:]['dir']
        T.close()

    def make_template(self):
        """
        Make template h5 with 1 direction
        """
        # open output h5, which will be modified
        self.h5 = tables.open_file(self.name_out, 'r+')

        # loop over solsets (example: sol000, sol001, ...)
        for solset in self.h5.root._v_groups.keys():
            ss = self.h5.root._f_get_child(solset)
            ss._f_get_child('source')._f_remove()
            values = np.array([(b'Dir00', [0., 0.])], dtype=[('name', 'S128'), ('dir', '<f4', (2,))])
            title = 'Source names and directions'
            self.h5.create_table(ss, 'source', values, title=title)

            # loop over soltabs (example: phase000, amplitude000, ...)
            for soltab in ss._v_groups.keys():
                st = ss._f_get_child(soltab)
                for axes in ['val', 'weight']:
                    AXES = st._f_get_child(axes).attrs['AXES']
                    dir_idx = AXES.decode('utf8').split(',').index('dir')
                    shape = list(st._f_get_child(axes)[:].shape)
                    shape[dir_idx] = 1

                    # only phases are zeros, others are ones
                    if 'phase' in soltab and axes != 'weight':
                        newvals = np.zeros(shape)
                    elif 'amplitude' in soltab or axes == 'weight':
                        newvals = np.ones(shape)
                    else:
                        newvals = np.zeros(shape)

                    # get correct value type
                    valtype = str(st._f_get_child(axes).dtype)
                    if '16' in valtype:
                        atomtype = tables.Float16Atom()
                    elif '32' in valtype:
                        atomtype = tables.Float32Atom()
                    elif '64' in valtype:
                        atomtype = tables.Float64Atom()
                    else:
                        atomtype = tables.Float64Atom()

                    # create new value/weight table
                    st._f_get_child(axes)._f_remove()
                    self.h5.create_array(st, axes, newvals.astype(valtype), atom=atomtype)
                    st._f_get_child(axes).attrs['AXES'] = AXES

                # modify direction axes
                st._f_get_child('dir')._f_remove()
                self.h5.create_array(st, 'dir', np.array([b'Dir00']).astype('|S5'))
        self.h5.close()

        return self

    def closest_dir(self, coor):
        """
        Find closest directions, using the astropy SkyCoord class

        :param coor: coordinate in radian
        """
        min_sep, min_sep_idx = 999, 999
        for n, dir in enumerate(self.dirs):
            c1 = SkyCoord(dir[0], dir[1], unit='radian', frame='icrs')  # your coords
            c2 = SkyCoord(coor[0], coor[1], unit='radian', frame='icrs')
            sep = c1.separation(c2).value
            if sep < min_sep:
                min_sep_idx = n
                min_sep = sep

        # Make sure this function selects a closest direction
        assert min_sep_idx != 999, "ERROR: coordinate comparison went wrong?! (most likely a bug in h5 or code)"
        print(f"Closest direction for {coor} is {self.dirs[min_sep_idx]}")

        return min_sep_idx

    def add_closest_values(self, coor, outcoor):
        """
        Add closest phase, amplitude, and weights to new output h5
        """
        # get closest direction index
        idx = self.closest_dir(coor)

        # open both input and output tables
        self.h5_in = tables.open_file(self.h5_in)
        self.h5_out = tables.open_file(self.name_out, 'r+')

        # loop over solsets (example: sol000, sol001, ...)
        for solset in self.h5_out.root._v_groups.keys():
            ss = self.h5_out.root._f_get_child(solset)

            # modify source table
            ss._f_get_child('source')._f_remove()
            if outcoor in ['msin', 'directions']:
                values = np.array([(b'Dir00', coor)],
                                  dtype=[('name', 'S128'), ('dir', '<f4', (2,))])
            else:
                values = np.array([(b'Dir00', self.h5_in.root.sol000.source[:][idx][1])],
                                  dtype=[('name', 'S128'), ('dir', '<f4', (2,))])
            title = 'Source names and directions'
            self.h5_out.create_table(ss, 'source', values, title=title)

            # loop over soltabs (example: phase000, amplitude000, ...)
            for soltab in ss._v_groups.keys():
                st = ss._f_get_child(soltab)

                st.dir._f_remove()
                self.h5_out.create_array(st, 'dir', np.array([b'Dir00'], dtype='|S5'))

                for axes in ['val', 'weight']:
                    AXES = st._f_get_child(axes).attrs['AXES']
                    dir_idx = AXES.decode('utf8').split(',').index('dir')

                    # select values corresponding to closest direction index (idx variable)
                    allvals = self.h5_in.root._f_get_child(solset)._f_get_child(soltab)._f_get_child(axes)[:]
                    newvals = np.take(allvals, indices=[idx], axis=dir_idx)
                    del allvals  # save memory

                    # get correct value type
                    valtype = str(st._f_get_child(axes).dtype)
                    if '16' in valtype:
                        atomtype = tables.Float16Atom()
                    elif '32' in valtype:
                        atomtype = tables.Float32Atom()
                    elif '64' in valtype:
                        atomtype = tables.Float64Atom()
                    else:
                        atomtype = tables.Float64Atom()

                    # create new value/weight table
                    st._f_get_child(axes)._f_remove()
                    self.h5_out.create_array(st, axes, newvals.astype(valtype), atom=atomtype)
                    st._f_get_child(axes).attrs['AXES'] = AXES
        self.h5_in.close()
        self.h5_out.close()
        return self


def make_list(arglist):
    try:
        return [[float(d) for d in
                 dir.replace('[', '').replace(']', '').replace('(', '').replace(')', '').replace(' ', '').split(',')]
                for dir in arglist]
    except ValueError:
        sys.exit("ERROR: --directions input invalid\nDo not use any spaces, example input: [0.1,0.2] [1.2,4.1]")


def get_dir_from_ms(msin):
    """
    Get phase direction from measurement set
    """

    output = []

    for ms in msin:
        with table(ms+"::FIELD") as t:
            output.append(list(t.getcol("PHASE_DIR").squeeze()))
    return output


def parse_args():
    """
    Command line argument parser

    :return: parsed arguments
    """

    parser = ArgumentParser()
    parser.add_argument('--h5_in', help='Input h5parm', required=True)
    parser.add_argument('--msin', help='MS to get phase center from', nargs="+")
    parser.add_argument('--directions', nargs='+', help='directions to find the closest h5_in direction to (Example: (0.1, 0.2) (1.2, 3.1) (3.5, 1.2)', default=None)
    parser.add_argument('--outcoor', help='Output coordinates from msin, directions, h5_in (only touch this if you know what you are doing)', type=str, default='h5_in')

    return parser.parse_args()


def main():
    """Main function"""

    args = parse_args()
    inputh5 = args.h5_in
    if args.directions is not None:
        dirs = make_list(args.directions)
    elif args.msin is not None:
        dirs = get_dir_from_ms(args.msin)
    else:
        sys.exit('ERROR: give --msin or --directions')

    os.system('mkdir output_h5s')
    print('Make folder --> output_h5s')

    for n, dir in enumerate(dirs):
        T = FindClosestDir(inputh5, f'output_h5s/source_{n}.h5')
        T.make_template()
        T.add_closest_values(dir, args.outcoor)


if __name__ == '__main__':
    main()
