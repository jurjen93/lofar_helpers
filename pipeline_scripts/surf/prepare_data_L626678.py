from argparse import ArgumentParser
from glob import glob
import os
import casacore.tables as ct
import time
import sys

parser = ArgumentParser()
parser.add_argument('--box', type=str, help='Measurement set input')
args = parser.parse_args()

TO = "/project/lofarvwf/Share/jdejong/output/A399"

BOX = 'box_' + args.box

ms_archives = sorted([b.split('/')[-1] for b in glob(TO + '/extract/' + BOX+'/*' + BOX + '.dysco.sub.shift.avg.weights.ms.archive0')])

# starting times for measurement sets that have to be cutted for freq
CUTFREQS = [5021107868.011121]

os.system("mkdir -p " + TO + "/selfcal/" + BOX)
for MS in ms_archives:
    if ct.table(TO + "/extract/" + BOX + '/' + MS).getcol('TIME')[0] in CUTFREQS:
        print('Cutting freq for ' + MS)
        os.system("python /home/lofarvwf-jdejong/scripts/lofar_helpers/supporting_scripts/flag_freq.py -ff='[15..19]' -msin " + TO + "/extract/" + BOX + '/' + MS+" -msout " + TO + "/selfcal/" + BOX + '/' + MS + '.goodfreq')
        os.system("python /home/lofarvwf-jdejong/scripts/lofar_helpers/supporting_scripts/flag_time.py -tf 0 1500 -msin " + TO + "/selfcal/" + BOX + '/' + MS + '.goodfreq' + " -msout " + TO + "/selfcal/" + BOX + '/' + MS + '.goodfreq.goodtimes')
        os.system("rm -rf " + TO + "/selfcal/" + BOX + '/' + MS + '.goodfreq')
    else:
        print('Cutting time for '+MS)
        os.system("python /home/lofarvwf-jdejong/scripts/lofar_helpers/supporting_scripts/flag_time.py -tf 0 1500 -msin " + TO + "/extract/" + BOX + '/' + MS + " -msout " + TO + "/selfcal/" + BOX + '/' + MS + '.goodtimes')

start = time.time()
MS = [ms.split('/')[-1] for ms in glob(TO + '/selfcal/' + BOX + '/*' + BOX + '.dysco.sub.shift.avg.weights.ms.archive*.goodtimes')]
while len(MS) != 1: # important to wait until everything is ready before moving on to the next script --> selfcal
    MS = [ms.split('/')[-1] for ms in glob(TO + '/selfcal/' + BOX + '/*' + BOX + '.dysco.sub.shift.avg.weights.ms.archive*.goodtimes')]
    if (time.time()-start)/3600>3:
        sys.exit('Time cutting did not finish on time (took longer than 3h). Please check why.')
print('DATA CUT AND PREPARED TOOK '+str(round((time.time()-start)/3600, 2))+' HOURS.')