import os
from argparse import ArgumentParser
from glob import glob

parser = ArgumentParser()
parser.add_argument('-from', '--from_where', type=str, help='data from where')
parser.add_argument('-to', '--to_where', type=str, help='to where')
args = parser.parse_args()

LOCATION=args.to_where+'/result'#/net/tussenrijn/data2/jurjendejong/L626678
LOCATION_BACKUP=args.to_where+'/result_backup'
FROM=args.from_where #/disks/paradata/shimwell/LoTSS-DR2/archive_other/L626678
SING_IMAGE='/net/rijn/data2/rvweeren/data/pill-latestMay2021.simg'
SING_BIND='/tmp,/dev/shm,/disks/paradata,/data1,/net/lofar1,/net/rijn,/net/nederrijn/,/net/bovenrijn,/net/botlek,/net/para10,/net/lofar2,/net/lofar3,/net/lofar4,/net/lofar5,/net/lofar6,/net/lofar7,/disks/ftphome,/net/krommerijn,/net/voorrijn,/net/achterrijn,/net/tussenrijn,/net/ouderijn,/net/nieuwerijn,/net/lofar8,/net/lofar9,/net/rijn8,/net/rijn7,/net/rijn5,/net/rijn4,/net/rijn3,/net/rijn2'
SINGULARITY=' '.join(['singularity exec -B', SING_BIND, SING_IMAGE])

os.system('{SINGULARITY} CleanSHM.py'.format(SINGULARITY=SINGULARITY))

#MAKE RESULT FOLDER
# os.system('mkdir '+LOCATION)
# os.system('mkdir '+LOCATION_BACKUP)
print('Made '+LOCATION)

#MOVE FILES
print('Moving files to '+LOCATION)
# os.system('scp -r lofarvwf-jdejong@spider.surfsara.nl:/project/lofarvwf/Share/jdejong/output/L626678/selfcal/all_directions.h5 '+LOCATION)
os.system('cp -r '+FROM+'/image_full_ampphase_di_m.NS.mask01.fits '+LOCATION)
os.system('cp -r '+FROM+'/image_full_ampphase_di_m.NS.DicoModel '+LOCATION)
os.system('cp -r '+FROM+'/*_uv.pre-cal_*.pre-cal.ms.archive '+LOCATION)
print('Finished moving files')

#CUT TIME FOR MESSY END PART (ONLY FOR THIS CASE APPLICABLE)
print('Making goodtimes')
for MS in glob('{LOCATION}/*_uv.pre-cal_*.pre-cal.ms.archive'.format(LOCATION=LOCATION)):
    os.system(SINGULARITY+' DPPP msin={MS} msout.storagemanager=dysco msout={MS}.goodtimes msin.ntimes=1500 steps=[]'.format(LOCATION=LOCATION, MS=MS))
    os.system('mv {MS} {LOCATION_BACKUP}'.format(LOCATION=LOCATION, MS=MS, LOCATION_BACKUP=LOCATION_BACKUP))
    print('Made '+MS+'.goodtimes')

#MAKE LIST WITH MEASUREMENT SETS
os.system('ls -1d {LOCATION}/*.goodtimes > {LOCATION}/big-mslist.txt'.format(LOCATION=LOCATION))

with open('/home/jurjendejong/scripts/lofar_helpers/pipeline_scripts/strw/ddf.txt') as f:
    lines = [l.replace('\n','') for l in f.readlines()]

#RUN DDF COMMAND
# print('Running DDF COMMAND')
# os.system('cd '+LOCATION)
# os.system(' '.join([SINGULARITY] + lines))
# print('Finished making new image')