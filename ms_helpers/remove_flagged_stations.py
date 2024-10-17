from casacore.tables import table
import numpy as np
import os
from shutil import rmtree, move
from argparse import ArgumentParser
from sys import exit


def remove_flagged_antennas(msin: str = None, msout: str = None, overwrite: bool = False):
    """
    Remove antennas that are full flagged (to save storage)

    input:
        - msfile: measurement set name
    """

    # Cannot both overwrite and give an output name
    if msout is not None and overwrite:
        exit('ERROR: You specified an --msout and ask to --overwrite. Please give only one of both.')

    # Set name for output if not given
    if msout is None:
        msout = f"flagged_{msin.split('/')[-1]}"

    # Read antenna names from Measurement Set
    with table(f"{msin}::ANTENNA", ack=False) as ants:
        ants_names = ants.getcol("NAME")

    # Read main tables Measurement Set
    with table(msin, readonly=True, ack=False) as ms:
        # Read the antenna ID columns
        antenna1 = ms.getcol('ANTENNA1')
        antenna2 = ms.getcol('ANTENNA2')

        # Read the FLAG column
        flags = ms.getcol('FLAG')

    # Get the unique antenna indices
    unique_antennas = np.unique(np.concatenate((antenna1, antenna2)))

    # Identify fully flagged antennas
    fully_flagged_antennas = []
    for ant in unique_antennas:
        # Find rows with this antenna
        ant_rows = np.where((antenna1 == ant) | (antenna2 == ant))
        # Check if all data for this antenna is flagged
        if np.all(flags[ant_rows]):
            fully_flagged_antennas.append(ant)

    if len(fully_flagged_antennas) == 0:
        print(f'No flagged antennas for {msin}, move on.')
        return

    # Get names of ants to filter
    ants_to_filter = ','.join([ants_names[idx] for idx in fully_flagged_antennas])
    print(f"Filtering fully flagged antennas: {ants_to_filter}")

    # Run DP3
    dp3_cmd = f'DP3 msin={msin} msout={msout} msout.storagemanager=dysco steps=[filter] \
    filter.type=filter filter.remove=true filter.baseline=!{ants_to_filter}'

    os.system(dp3_cmd)

    # Overwrite input
    if overwrite:
        rmtree(msin)
        move(msout, msin)


def parse_args():
    """
    Parse input arguments
    """

    parser = ArgumentParser(description='Remove stations from the MS that are fully flagged')
    parser.add_argument('msin', type=str, help='Input Measurement Set')
    parser.add_argument('--msout', type=str, default=None, help='Output Measurement Set')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite input Measurement Set')

    return parser.parse_args()


def main():
    """
    Main function
    """

    args = parse_args()
    remove_flagged_antennas(args.msin, args.msout, args.overwrite)


if __name__ == '__main__':
    main()
