from argparse import ArgumentParser
from glob import glob
import time
import os

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
    print('file not found')

if len(box_archives) == 6:
    for N, SUBBOX in enumerate(box_archives):
        print(N)
        N = str(N + 1)
        cml = [
            "cp -r "+TO+"/extract/" + SUBBOX + " " + TO+"/selfcal/" + BOX + '.' + N,
            "cd "+TO+"/selfcal/" + TO+"/selfcal/" + BOX + '.' + N,
            "singularity exec -B +" + SING_BIND + " " + SING_IMAGE + " " + " python "+args.script_path+"/runwscleanLBautoR.py -b "+TO+"/boxes/" + BOX + ".reg --auto --imager=DDFACET --helperscriptspath=SCRIPT_PATH/ --autofrequencyaverage-calspeedup='True' "+TO+"/selfcal/" + BOX + "/" +SUBBOX,
            "rm -rf "+TO+"/selfcal/"+SUBBOX
        ]
        cml = " && ".join([c.replace('TO', TO) for c in cml])
        print("FOLLOWING COMMAND READY TO BE EXECUTED:\n" + cml)
        os.system("mkdir " + TO+"/selfcal/" + BOX + '.' + N)
        with open(TO+"/selfcal/" + BOX + '.' + N + "/command.txt", "w+") as f:
            f.write(cml)
else:
    raise ValueError("SOMETHING WENT WRONG WITH SELFCALLING " + BOX)
