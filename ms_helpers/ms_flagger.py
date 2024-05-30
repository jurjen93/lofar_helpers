from argparse import ArgumentParser
import os


def parse_args():
    """
    Command line argument parser

    :return: parsed arguments
    """
    parser = ArgumentParser()
    parser.add_argument('--ms', nargs='+', help='Input ms', required=True)
    parser.add_argument('--freqrange', help='Frequency range in MHz (example: 140-158)')
    parser.add_argument('--timerange', help='Time range in hours from the start (example: 1:30:0-2:30:0)')
    parser.add_argument('--ant', help='Antenna (example: RS409HBA)', required=True)
    return parser.parse_args()


def main():
    """
    Main function
    """

    args = parse_args()

    if args.freqrange is not None:
        min_freq, max_freq = [float(f) for f in args.freqrange.split('-')]
    if args.timerange is not None:
        min_time, max_time = [f for f in args.timerange.split('-')]

    for ms in args.ms:
        command = ['DP3',
                   'msout.storagemanager=dysco',
                   f'msin={ms}',
                   f'msout=flagged_{ms.split("/")[-1]}',
                   'steps=[flag]',
                   'flag.type=preflagger',
                   f'flag.baseline={args.ant}']
        # frequency flagging
        if args.freqrange is not None:
            command += [f"flag.freqrange='[{min_freq} .. {max_freq} MHz]'"]
        # time flagging
        if args.timerange is not None:
            command += [f"flag.reltime='[{min_time} .. {max_time}]'"]

        print(' '.join(command))
        os.system(' '.join(command))

if __name__ == '__main__':
    main()
