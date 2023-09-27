from argparse import ArgumentParser
import os

if __name__ == '__main__':

    parser = ArgumentParser()
    parser.add_argument('--ms', nargs='+', help='Input ms', required=True)
    parser.add_argument('--freqrange', help='Frequency range in MHz (example: 140-158)', required=True)
    parser.add_argument('--ant', help='Antenna (example: RS409HBA)', required=True)
    args = parser.parse_args()

    min_freq, max_freq = [float(f) for f in args.freqrange.split('-')]

    for ms in args.ms:
        command = ['DP3',
                   'msout.storagemanager=dysco',
                   f'msin={ms}',
                   f'msout=flagged_{ms.split("/")[-1]}',
                   'steps=[flag]',
                   'flag.type=preflagger',
                   f'flag.baseline={args.ant}',
                   f"flag.freqrange='[{min_freq} .. {max_freq} MHz]'"]

        print(' '.join(command))
        os.system(' '.join(command))