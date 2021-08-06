from argparse import ArgumentParser
from glob import glob
import time
import os
import casacore.tables as ct

parser = ArgumentParser()
parser.add_argument('--box', type=str, help='Measurement set input')
args = parser.parse_args()

TO = "/project/lofarvwf/Share/jdejong/output/A399"

BOX = 'box_' + args.box

SING_IMAGE = "/home/lofarvwf-jdejong/singularities/lofar_sksp_fedora31_ddf.sif"
SING_BIND = "/project/lofarvwf/Share/jdejong,/home/lofarvwf-jdejong/scripts"

ms_archives = sorted([b.split('/')[-1] for b in glob(TO + '/extract/*' + BOX + '.dysco.sub.shift.avg.weights.ms.archive*')])

while len(ms_archives) != 6:
    time.sleep(5)

#starting times fom measurement sets that have to be cutted for time
CUTTIMES = [5019387068.011121, 5017577408.011121, 5020506668.011121]

#starting times for measurement sets that have to be cutted for freq
CUTFREQS = [5021107868.011121]

os.system("mkdir " + TO+"/selfcal/" + BOX)
for MS in ms_archives:
    if ct.table(TO+"/extract/" + MS).getcol('TIME')[0] in CUTTIMES:
        print('Cutting time for '+MS)
        os.system("python /home/lofarvwf-jdejong/scripts/lofar_helpers/supporting_scripts/flag_time.py -tf 0 1500 -msin "+TO+"/extract/" + MS+" -msout "+TO+"/selfcal/"+BOX+'/'+MS+'.goodtimes')
    elif ct.table(TO+"/extract/" + MS).getcol('TIME')[0] in CUTFREQS:
        print('Cutting freq for '+MS)
        os.system("python /home/lofarvwf-jdejong/scripts/lofar_helpers/supporting_scripts/flag_freq.py -ff='[15..19]' -msin "+TO+"/extract/" + MS+" -msout "+TO+"/selfcal/"+BOX+'/'+MS+'.goodfreq')
    else:
        os.system("cp -r "+TO+"/extract/" + MS + " " + TO+"/selfcal/" + BOX)


os.system("cd "+TO+"/selfcal/"+BOX+" && python /home/lofarvwf-jdejong/scripts/runwscleanLBautoR.py -b "+TO+"/boxes/" + BOX + ".reg --auto --imager=DDFACET --helperscriptspath="+args.script_path+" --autofrequencyaverage-calspeedup='True' "+BOX + '.dysco.sub.shift.avg.weights.ms.archive*')