from argparse import ArgumentParser
from glob import glob
import time
import os
import casacore.tables as ct

parser = ArgumentParser()
parser.add_argument('--box', type=str, help='Measurement set input')
args = parser.parse_args()

BOX = 'box_' + args.box

ms_archives = sorted([b.split('/')[-1] for b in glob('*.dysco.sub.shift.avg.weights.ms.archive*')])

while len(ms_archives) != 6:
    time.sleep(5)

#starting times fom measurement sets that have to be cutted for time
CUTTIMES = [5019387068.011121, 5017577408.011121, 5020506668.011121]

#starting times for measurement sets that have to be cutted for freq
CUTFREQS = [5021107868.011121]

os.system("mkdir -p selfcal")
for MS in ms_archives:
    if ct.table(MS).getcol('TIME')[0] in CUTTIMES:
        print('Cutting time for '+MS)
        os.system("python /home/jurjendejong/scripts/lofar_helpers/supporting_scripts/flag_time.py -tf 0 1500 -msin " + MS + " -msout selfcal/" + MS + '.goodtimes')
    elif ct.table(MS).getcol('TIME')[0] in CUTFREQS:
        print('Cutting freq for ' + MS)
        os.system("python /home/jurjendejong/scripts/lofar_helpers/supporting_scripts/flag_freq.py -ff='[15..19]' -msin " + MS+" -msout selfcal/" + MS + '.goodfreq')
    else:
        os.system("cp -r " + MS + " selfcal/")

MS = [ms.split('/')[-1] for ms in glob('selfcal/*' + BOX + '.dysco.sub.shift.avg.weights.ms.archive*')]
while len(MS) != 6:#important to wait until everything is ready before moving on to the next script --> selfcal
    MS = [ms.split('/')[-1] for ms in glob('selfcal/*' + BOX + '.dysco.sub.shift.avg.weights.ms.archive*')]
print('DATA CUT AND PREPARED')