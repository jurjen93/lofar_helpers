import os
from argparse import ArgumentParser
from glob import glob
from os import path


parser = ArgumentParser()
parser.add_argument('-from', '--from_where', type=str, help='directory where data is from', required=True)
parser.add_argument('-to', '--to_where', type=str, help='destination directory', required=True)
parser.add_argument('-tf', '--time_flag', nargs='+', help='flag time: start_time end_time', required=False)
args = parser.parse_args()

LOCATION=args.to_where#/net/tussenrijn/data2/jurjendejong/L626678
FROM=args.from_where #/disks/paradata/shimwell/LoTSS-DR2/archive_other/L626678
SING_IMAGE='/net/rijn/data2/rvweeren/data/pill-latestMay2021.simg'
SING_BIND='/tmp,/dev/shm,/disks/paradata,/data1,/net/lofar1,/net/rijn,/net/nederrijn/,/net/bovenrijn,/net/botlek,/net/para10,/net/lofar2,/net/lofar3,/net/lofar4,/net/lofar5,/net/lofar6,/net/lofar7,/disks/ftphome,/net/krommerijn,/net/voorrijn,/net/achterrijn,/net/tussenrijn,/net/ouderijn,/net/nieuwerijn,/net/lofar8,/net/lofar9,/net/rijn8,/net/rijn7,/net/rijn5,/net/rijn4,/net/rijn3,/net/rijn2'
SINGULARITY=' '.join(['singularity exec -B', SING_BIND, SING_IMAGE])

#CREATE DESTINATION DIRECTORY IF NOT EXISTS
if not path.exists(LOCATION):
    os.system('mkdir {LOCATION}'.format(LOCATION=LOCATION))
    print('Created {LOCATION}'.format(LOCATION=LOCATION))

#CLEAN CACHE
os.system('{SINGULARITY} CleanSHM.py'.format(SINGULARITY=SINGULARITY))

#MOVE FILES
print('Moving files to '+LOCATION)
# os.system('scp -r lofarvwf-jdejong@spider.surfsara.nl:/project/lofarvwf/Share/jdejong/output/L626678/selfcal/all_directions.h5 '+LOCATION)
os.system('cp -r '+FROM+'/image_full_ampphase_di_m.NS.mask01.fits '+LOCATION)
os.system('cp -r '+FROM+'/image_full_ampphase_di_m.NS.DicoModel '+LOCATION)
# os.system('cp -r '+FROM+'/*_uv.pre-cal_*.pre-cal.ms.archive '+LOCATION)
print('Finished moving files')

#FLAG TIME
if args.time_flag:
    print('Making goodtimes')
    for MS in glob('{LOCATION}/*_uv.pre-cal_*.pre-cal.ms.archive'.format(LOCATION=LOCATION)):
        os.system('{SINGULARITY} python supporting_scripts/flag_time.py -tf {time[0]} {time[1]} -ms {MS}'.format(SINGULARITY=SINGULARITY, MS=MS, time=args.flag_time))

#MAKE LIST WITH MEASUREMENT SETS
os.system('ls -1d {LOCATION}/*.goodtimes > {LOCATION}/big-mslist.txt'.format(LOCATION=LOCATION))

with open('/home/jurjendejong/scripts/lofar_helpers/DDF_scripts/ddf.txt') as f:
    lines = [l.replace('\n','') for l in f.readlines()]

#RUN DDF COMMAND
print('Running DDF COMMAND')
os.system('cd '+LOCATION)
os.system(' '.join([SINGULARITY] + lines))
print('Finished making new image')