import tables
from argparse import ArgumentParser
import os
import numpy as np


def parse_args():
    """
    Command line argument parser

    :return: parsed arguments
    """

    parser = ArgumentParser()
    parser.add_argument('--h5_in', type=str, help='input h5 (from which to extract the frequencies and antennas)',
                        required=True)
    parser.add_argument('--h5_dirs', type=str, help='h5 with directions to copy from', required=True)
    parser.add_argument('--h5_out', type=str, help='output h5')
    return parser.parse_args()


def main():
    """Main function"""

    args = parse_args()

    os.system(f'cp {args.h5_in} {args.h5_out}')

    print("Make " + args.h5_out)

    H = tables.open_file(args.h5_out, "r+")
    T = tables.open_file(args.h5_dirs)

    H.root.sol000.source._f_remove()
    H.create_table(H.root.sol000, 'source', T.root.sol000.source[:], title='Source names and directions')

    for table in ['phase000', 'amplitude000']:
        print(table)

        for axes in ['val', 'weight']:
            print(axes)

            vals = H.root.sol000._f_get_child(table)._f_get_child(axes)
            AXES = vals.attrs["AXES"]
            dirindex = AXES.decode('utf8').split(',').index('dir')
            newshape = list(vals.shape)
            newshape[dirindex] = len(T.root.sol000.source[:])
            newvals = np.zeros(newshape)

            for d in range(newshape[dirindex]):
                newvals[:, :, :, d, :] = vals[:, :, :, 0, :]
            H.root.sol000._f_get_child(table)._f_get_child(axes)._f_remove()
            H.create_array(H.root.sol000._f_get_child(table), axes, newvals)
            H.root.sol000._f_get_child(table)._f_get_child(axes).attrs['AXES'] = AXES

        H.root.sol000._f_get_child(table).dir._f_remove()
        H.create_array(H.root.sol000._f_get_child(table), 'dir', T.root.sol000._f_get_child(table).dir[:])


if __name__ == '__main__':
    main()
