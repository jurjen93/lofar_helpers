"""
Script to flag freq
"""

__author__ = "Jurjen de Jong (jurjendejong@strw.leidenuniv.nl)"

from argparse import ArgumentParser
import os

parser = ArgumentParser()
parser.add_argument('-ff', '--frequency_flag', nargs='+', help='freq flag: channels', required=True, type=str)
parser.add_argument('-ms', '--measurement_set', help='measurement set', required=True)
args = parser.parse_args()

os.system("DPPP msin={MSIN} msout={MSIN}.goodfreq msout.storagemanager=dysco steps=[pf] pf.type=preflagger pf.chan={FREQCHAN}".format(MSIN=args.measurement_set, FREQCHAN=args.frequency_flag))