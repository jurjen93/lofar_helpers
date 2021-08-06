"""
Script to flag freq
"""

__author__ = "Jurjen de Jong (jurjendejong@strw.leidenuniv.nl)"

from argparse import ArgumentParser
import os

parser = ArgumentParser()
parser.add_argument('-ff', '--frequency_flag', nargs='+', help='freq flag: channels', required=True, type=str)
parser.add_argument('-msin', '--ms_in', help='measurement set', required=True)
parser.add_argument('-msout', '--ms_out', help='measurement set', required=False, default='')
args = parser.parse_args()

if not args.ms_out:
    msout = args.ms_in+'.goodfreq'
else:
    msout = args.ms_out

os.system("DPPP msin={MSIN} msout={MSOUT} msout.storagemanager=dysco steps=[pf] pf.type=preflagger pf.chan={FREQCHAN}".format(MSIN=args.ms_in, MSOUT=msout, FREQCHAN=args.frequency_flag))