"""
Example:
    python ~/scripts/lofar_helpers/make_new_image.py
        -from /disks/paradata/shimwell/LoTSS-DR2/archive_other/L626678
        -to /net/tussenrijn/data2/jurjendejong/L626678/result
        -tf 0 3000
"""

import os
from argparse import ArgumentParser
from glob import glob
from os import path

parser = ArgumentParser()
parser.add_argument('-from', '--from_where', type=str, help='directory where data is originally from', required=True)
parser.add_argument('-to', '--to_where', type=str, help='destination directory', required=True)
parser.add_argument('-tf', '--time_flag', nargs='+', help='flag time: start_time end_time', required=False)
args = parser.parse_args()

LOCATION=args.to_where#/net/tussenrijn/data2/jurjendejong/L626678
FROM=args.from_where #/disks/paradata/shimwell/LoTSS-DR2/archive_other/L626678
SING_IMAGE='/net/rijn/data2/rvweeren/data/pill-latestMay2021.simg'
SING_BIND='/net/tussenrijn'
SINGULARITY=' '.join(['singularity exec -B', SING_BIND, SING_IMAGE])
# SINGULARITY=' '.join(['singularity exec ', SING_IMAGE])

#CREATE DESTINATION DIRECTORY IF NOT EXISTS
# if not path.exists(LOCATION):
#     os.system('mkdir {LOCATION}'.format(LOCATION=LOCATION))
#     print('Created {LOCATION}'.format(LOCATION=LOCATION))

#CLEAN CACHE
os.system('{SINGULARITY} CleanSHM.py'.format(SINGULARITY=SINGULARITY))

#MOVE FILES
print('Moving files to '+LOCATION)
# # os.system("scp -r lofarvwf-jdejong@spider.surfsara.nl:/project/lofarvwf/Share/jdejong/output/L626678/selfcal/all_directions.h5 "+LOCATION)
os.system('cp -r '+FROM+'/image_full_ampphase_di_m.NS.mask01.fits '+LOCATION)
os.system('cp -r '+FROM+'/image_full_ampphase_di_m.NS.DicoModel '+LOCATION)
# os.system('cp -r '+FROM+'/*_uv.pre-cal_*.pre-cal.ms.archive '+LOCATION)
print('Finished moving files')

#FLAG TIME
if args.time_flag:
    print('Making goodtimes')
    for MS in glob('{LOCATION}/*_uv.pre-cal_*.pre-cal.ms.archive'.format(LOCATION=LOCATION)):
        os.system('{SINGULARITY} python /home/jurjendejong/scripts/lofar_helpers/supporting_scripts/flag_time.py -tf {time} -ms {MS}'.format(SINGULARITY=SINGULARITY, MS=MS, time=' '.join(args.time_flag)))
        print('Created ' + MS.split('/')[-1] + '.goodtimes')

#MAKE LIST WITH MEASUREMENT SETS
os.system('ls -1d {LOCATION}/*.goodtimes > {LOCATION}/big-mslist.txt'.format(LOCATION=LOCATION))

with open('/home/jurjendejong/scripts/lofar_helpers/DDF_scripts/ddf.txt') as f:
    lines = [l.replace('\n','') for l in f.readlines()]
    lines+=['--Data-MS={LOCATION}/big-mslist.txt'.format(LOCATION=LOCATION)]
    lines+=['--Predict-InitDicoModel={LOCATION}/image_full_ampphase_di_m.NS.DicoModel'.format(LOCATION=LOCATION)]
    lines+=['--DDESolutions-DDSols={LOCATION}/complete_merged.h5:sol000/amplitude000+phase000'.format(LOCATION=LOCATION)]
    lines+=['--Mask-External={LOCATION}/image_full_ampphase_di_m.NS.mask01.fits'.format(LOCATION=LOCATION)]

#RUN DDF COMMAND
print('Running DDF COMMAND')
os.system(' '.join(['cd', LOCATION, '&&', SINGULARITY] + lines))
print('Finished making new image')