from argparse import ArgumentParser
from glob import glob
import time
import os
import casacore.tables as ct

parser = ArgumentParser()
parser.add_argument('--box', type=str, help='Measurement set input')
parser.add_argument('--source', type=str, help='Source name')
parser.add_argument('--script_path', type=str, help='Path to scripts')
args = parser.parse_args()

TO = "/project/lofarvwf/Share/jdejong/output/" + args.source

BOX = 'box_' + args.box

SING_IMAGE = "/home/lofarvwf-jdejong/singularities/lofar_sksp_fedora31_ddf.sif"
SING_BIND = "/project/lofarvwf/Share/jdejong,/home/lofarvwf-jdejong/scripts"

box_archives = sorted([b.split('/')[-1] for b in glob(TO + '/extract/*' + BOX + '.dysco.sub.shift.avg.weights.ms.archive*')])

while len(box_archives) != 6:
    time.sleep(5)

#starting times fom measurement sets that have to be cutted for time
CUTTIMES = [5019387068.011121, 5017577408.011121, 5020506668.011121]

#starting times for measurement sets that have to be cutted for freq
CUTFREQS = [5021107868.011121]


if len(box_archives) == 6:
    for N, SUBBOX in enumerate(box_archives):
        N = str(N + 1)
        if ct.table(TO+"/extract/" + SUBBOX).getcol('TIME')[0] in CUTTIMES:
            print('Cutting time for '+SUBBOX)
            cml = [
                "singularity exec -B " + SING_BIND + " " + SING_IMAGE + " python "+args.script_path+"/supporting_scripts/flag_time.py -tf 0 1500 -msin "+TO+"/extract/" + SUBBOX+" -msout "+TO+"/selfcal/"+BOX + '.' + N+'/'+SUBBOX,
                "cd "+TO+"/selfcal/" + BOX + '.' + N,
                "singularity exec -B " + SING_BIND + " " + SING_IMAGE + " python "+args.script_path+"/runwscleanLBautoR.py -b "+TO+"/boxes/" + BOX + ".reg --auto --imager=DDFACET --helperscriptspath="+args.script_path+" --autofrequencyaverage-calspeedup='True' "+SUBBOX+'.goodtimes']
        elif ct.table(TO+"/extract/" + SUBBOX).getcol('TIME')[0] in CUTFREQS:
            print('Cutting freq for '+SUBBOX)
            cml = [
                "singularity exec -B " + SING_BIND + " " + SING_IMAGE + " python "+args.script_path+"/supporting_scripts/flag_freq.py -ff='[15..19]' -msin "+TO+"/extract/" + SUBBOX+" -msout "+TO+"/selfcal/"+BOX + '.' + N+'/'+SUBBOX,
                "cd "+TO+"/selfcal/" + BOX + '.' + N,
                "singularity exec -B " + SING_BIND + " " + SING_IMAGE + " python "+args.script_path+"/runwscleanLBautoR.py -b "+TO+"/boxes/" + BOX + ".reg --auto --imager=DDFACET --helperscriptspath="+args.script_path+" --autofrequencyaverage-calspeedup='True' "+SUBBOX+'.goodtimes']
        else:
            cml = [
                "cp -r "+TO+"/extract/" + SUBBOX + " " + TO+"/selfcal/" + BOX + '.' + N,
                "cd "+TO+"/selfcal/" + BOX + '.' + N,
                "singularity exec -B " + SING_BIND + " " + SING_IMAGE + " python "+args.script_path+"/runwscleanLBautoR.py -b "+TO+"/boxes/" + BOX + ".reg --auto --imager=DDFACET --helperscriptspath="+args.script_path+" --autofrequencyaverage-calspeedup='True' "+SUBBOX]
        os.system("mkdir " + TO+"/selfcal/" + BOX + '.' + N)
        with open(TO+"/selfcal/" + BOX + '.' + N + "/command.sh", "w+") as f:
            f.write("#!/bin/bash\n")
            f.write("\n".join(cml))
        os.system("chmod u+x "+TO+"/selfcal/" + BOX + '.' + N + "/command.sh")
else:
    raise ValueError("SOMETHING WENT WRONG WITH SELFCALLING " + BOX)
