from argparse import ArgumentParser
from glob import glob
import time

parser = ArgumentParser()
parser.add_argument('--box', type=str, help='Measurement set input')
parser.add_argument('--source', type=str, help='Source name')
parser.add_argument('--script_path', type=str, help='Path to scripts')
args = parser.parse_args()

TO = "/project/lofarvwf/Share/jdejong/output/" + args.field

BOX = 'box_' + args.box

SING_IMAGE = "/home/lofarvwf-jdejong/singularities/lofar_sksp_fedora31_ddf.sif"
SING_BIND = "/project/lofarvwf/Share/jdejong,/home/lofarvwf-jdejong/scripts"

box_archives = sorted(glob(TO + '/extract/*' + BOX + '.dysco.sub.shift.avg.weights.ms.archive*'))

while len(box_archives) != 6:
    time.sleep(5)
    print('FILE NOT FOUND IN WRITE SCRIPT')

if len(box_archives) == 6:
    for N, SUBBOX in enumerate(box_archives):
        N = str(N + 1)
        print(N)
        SELFCAL_FOLDER = "TO/selfcal/" + BOX + '.' + N
        cml = [
            "cd TO/selfcal/",
            "mkdir " + SELFCAL_FOLDER,
            "cp -r " + SUBBOX + " " + SELFCAL_FOLDER,
            "cd TO/selfcal/" + SELFCAL_FOLDER,
            "singularity exec -B +" + SING_BIND + " " + SING_IMAGE + " " + " python SCRIPT_PATH/runwscleanLBautoR.py -b TO/boxes/" + BOX + ".reg --auto --imager=DDFACET --helperscriptspath=SCRIPT_PATH/ --autofrequencyaverage-calspeedup='True' TO/selfcal/" + BOX + "/" +
            SUBBOX.split('/')[-1]
        ]
        cml = " && ".join([c.replace('TO', TO).replace('SCRIPT_PATH', args.script_path) for c in cml])
        print("FOLLOWING COMMAND READY TO BE EXECUTED:\n" + cml)
        with open(SELFCAL_FOLDER + "/command.txt", "w") as f:
            f.write(cml)
else:
    raise ValueError("SOMETHING WENT WRONG WITH SELFCALLING " + BOX)
