"""
Script to flag time
"""

__author__ = "Jurjen de Jong (jurjendejong@strw.leidenuniv.nl)"

from argparse import ArgumentParser
import pyrap.tables as pt

parser = ArgumentParser()
parser.add_argument('-tf', '--time_flag', nargs='+', help='flag time: start_time end_time', required=True)
parser.add_argument('-msin', '--ms_in', help='measurement set', required=True)
parser.add_argument('-msout', '--ms_out', help='measurement set', required=False, default='')
args = parser.parse_args()

if not args.msout:
    msout = args.ms_in+'.goodtimes'
else:
    msout = args.ms_in

pt.taql('SELECT FROM {MSIN} WHERE TIME IN (SELECT DISTINCT TIME FROM {MSIN} OFFSET {time[0]} LIMIT {time[1]}) GIVING {MSOUT} AS PLAIN'.format(MSIN=args.ms_in, MSOUT=msout, time=args.time_flag))