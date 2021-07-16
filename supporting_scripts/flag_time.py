from argparse import ArgumentParser
import pyrap.tables as pt

parser = ArgumentParser()
parser.add_argument('-tf', '--time_flag', nargs='+', help='flag time: start_time end_time', required=True)
parser.add_argument('-ms', '--measurement_set', help='measurement set', required=True)
args = parser.parse_args()

pt.taql('SELECT FROM {MS} WHERE TIME IN (SELECT DISTINCT TIME FROM {MS} OFFSET {time[0]} LIMIT {time[1]}) GIVING {MS}.goodtimes AS PLAIN'.format(MS=args.measurement_set, time=args.time_flag))
# os.system(SINGULARITY+' DPPP msin={MS} msout.storagemanager=dysco msout={MS}.goodtimes msin.ntimes=1500 steps=[]'.format(LOCATION=LOCATION, MS=MS))