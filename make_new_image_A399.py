"""
Example:
    python ~/scripts/lofar_helpers/make_new_image.py
        -from /disks/paradata/shimwell/LoTSS-DR2/archive_other/L626678
        -to /net/tussenrijn/data2/jurjendejong/L626678/result
        -h5 complete_merged.h5
"""

__author__ = "Jurjen de Jong (jurjendejong@strw.leidenuniv.nl)"

import os
from os import path
import tables
import casacore.tables as ct
from glob import glob
from numpy import unique
import pathlib
import operator

def get_DDS3(folder):
    """
    This function returns from a folder the last DDS file for all the measurement sets within that same folder.
    :param folder: string of folder path
    """
    ms = glob(folder+'/*.pre-cal.ms.archive*') # get all measurement sets in folder
    ms_observations = unique([m.split('/')[-1].split('_')[0] for m in ms]) # get unique observation names
    for i in '321': # get correct DDS*.npz files [should be DDS3*]
        DDS = glob(folder+'/DDS'+i+'*')
        if DDS: break

    time_modification = {D: pathlib.Path(D).stat().st_mtime for D in DDS} # get dict with DDS and time of modification
    time_observation = {D: int(D.split('full_')[1].split('.')[0]) for D in DDS if 'slow' not in D} # get dict with DDS and time of observation
    single_m = [[m for m in ms if night in m][0] for night in ms_observations] # get one measurement set per observation

    DDS_output = [] # DDS output
    DDS_dict = {}
    for observation in ms_observations:
        for sm in single_m:
            if observation in sm:
                table = ct.table(sm) # open table
                t = table.getcol('TIME')[0] # get first time element from measurement set
                table.close() # close table
                diff = lambda ob_time : abs(ob_time - t) #f formula to compare difference between first time element of ms and observation times
                closest_value = min(list(time_observation.values()), key=diff) # get closest value with lambda function
                DDS_options = {D: time_modification[D] for D in
                               [f for f, time in time_observation.items() if closest_value == time]} # get optional DDS files
                try: # python 2
                    correct_DDS = max(DDS_options.iteritems(), key=operator.itemgetter(1))[0] # get correct DDS
                except: # python 3
                    correct_DDS = max(DDS_options.items(), key=operator.itemgetter(1))[0] # get correct DDS
                DDS_output.append([D for D in DDS if correct_DDS.split('full_')[1].split('_smoothed')[0] in D]) # append right DDS
                DDS_dict.update({sm : [D for D in DDS if correct_DDS.split('full_')[1].split('_smoothed')[0] in D]})

    return [file for sublist in DDS_output for file in sublist], DDS_dict

TO='/net/tussenrijn/data2/jurjendejong/A399/result'
FROM='/net/tussenrijn/data2/jurjendejong/A399_DEEP'
SING_IMAGE='/net/rijn/data2/rvweeren/data/pill-latestMay2021.simg'
SING_BIND='/net/tussenrijn'
SINGULARITY=' '.join(['singularity exec -B', SING_BIND, SING_IMAGE])
# SINGULARITY=' '.join(['singularity exec ', SING_IMAGE])

#CREATE DESTINATION DIRECTORY IF NOT EXISTS
if not path.exists(TO):
    os.system('mkdir {LOCATION}'.format(LOCATION=TO))
    print('Created {LOCATION}'.format(LOCATION=TO))

#CLEAN CACHE
os.system('CleanSHM.py')

#MOVE FILES
print('Moving files to '+TO)
command = 'cp -r '+FROM+'/image_full_ampphase_di_m.NS.mask01.fits '+TO+ ' && '+\
        'cp -r '+FROM+'/image_full_ampphase_di_m.NS.DicoModel '+TO+' && wait'
        # 'cp -r '+FROM+'/*_uv.pre-cal_*.pre-cal.ms.archive '+TO+' && wait'

# os.system(command)
print('Finished moving files')


#----------------------------------------------------------------------------------------------------------------------

#FLAG TIME AND FREQUENCY

#starting times fom measurement sets that have to be cutted for time
CUTTIMES = [5019387068.011121, 5019387064.005561, 5017577408.011121, 5017577404.005561, 5020506668.011121, 5020506664.005561]

#starting times for measurement sets that have to be cutted for freq
CUTFREQS = [5021107868.011121, 5021107864.005561]

for MS in glob('/net/tussenrijn/data2/jurjendejong/A399_DEEP/*.ms.archive'):
    t = ct.table(MS)
    time = t.getcol('TIME')[0]
    t.close()
    if time in CUTTIMES:
        print('Cutting time for '+MS)
        # os.system("python /home/jurjendejong/scripts/lofar_helpers/supporting_scripts/flag_time.py -tf 0 1500 -msin " + MS + " -msout " + TO + '/' + MS.split('/')[-1] + '.goodtimes')
    elif time in CUTFREQS:
        print('Cutting freq for ' + MS)
        # os.system("python /home/jurjendejong/scripts/lofar_helpers/supporting_scripts/flag_freq.py -ff='[15..19]' -msin " + MS+" -msout " + TO + '/' + MS.split('/')[-1] + '.goodfreq')
        # os.system("cp -r " + MS + " " + TO)
    else:
        print('Copying for ' + MS)
        # os.system("cp -r " + MS + " " + TO)

# important to wait until everything is ready before moving on
while len(glob('/net/tussenrijn/data2/jurjendejong/A399_DEEP/*.ms.archive')) != len(glob(TO+'/*.pre-cal.ms.archive*')):
    print('TIME AND FREQUENCY FLAGGING')
#----------------------------------------------------------------------------------------------------------------------

#MERGE LOTSS OUTER EDGE


DDS3, DDS3_dict = get_DDS3('/net/tussenrijn/data2/jurjendejong/A399_DEEP')

soltable_times = {}
for soltable in glob('/net/tussenrijn/data2/jurjendejong/A399/result/all_directions*.h5'):
    tab = tables.open_file(soltable)
    t = tab.root.sol000.phase000.time[0]
    soltable_times.update({t: soltable})
    tab.close()  # close table
print(soltable_times)

# os.system('mkdir /net/tussenrijn/data2/jurjendejong/A399/result_filtered')
# os.system(' && '.join(['cp '+s+' /net/tussenrijn/data2/jurjendejong/A399/result_filtered' for s in DDS3]))
command = []
for ms in DDS3_dict.items():
    new_h5=[]
    for npz in ms[1]:
        command.append('KillMS2H5parm.py ' + npz.split('/')[-1].replace('npz','h5 ') + npz + ' --nofulljones')
        new_h5.append(npz.split('/')[-1].replace('npz','h5 '))

    table = ct.table(ms[0])  # open table
    t = table.getcol('TIME')[0]  # get first time element from measurement set
    table.close()  # close table
    diff = lambda ob_time: abs(ob_time - t)  # formula to compare difference between first time element of ms and observation times
    closest_value = min(list(soltable_times.keys()), key=diff)  # get closest value with lambda function
    h5 = soltable_times[closest_value]
    command.append('python /home/jurjendejong/scripts/lofar_helpers/h5_merger.py -out final_lotss_'+str(closest_value)+'.h5 -in '+' '.join(new_h5) + '--convert_tec 0')
    command.append('python /home/jurjendejong/scripts/lofar_helpers/supporting_scripts/h5_filter.py -f /net/tussenrijn/data2/jurjendejong/A399_DEEP/image_full_ampphase_di_m.NS.app.restored.fits -ac 2.5 -in false -h5out lotss_full_merged_filtered_'+str(closest_value)+'.h5 -h5in final_lotss_'+str(closest_value)+'.h5')
    command.append('python /home/jurjendejong/scripts/lofar_helpers/h5_merger.py -out complete_merged_'+str(closest_value)+'.h5 -in lotss_full_merged_filtered_'+str(closest_value)+'.h5 ' + soltable_times[closest_value]+' --convert_tec 0')
print('cd /net/tussenrijn/data2/jurjendejong/A399/result_filtered && '+' && '.join(command))
os.system('cd /net/tussenrijn/data2/jurjendejong/A399/result_filtered && '+' && '.join(command))
"""
OUTPUT_FOLDER=${FOLDER}/result_filtered
mkdir ${OUTPUT_FOLDER}

cd ${OUTPUT_FOLDER}

#converging and merging .npz to .h5 files from LoTSS output
singularity exec -B ${SING_BIND} ${SING_IMAGE} killMS2H5parm.py lotss_slow.h5 ${FOLDER}/extract/DDS3_full_slow*merged.npz --nofulljones
singularity exec -B ${SING_BIND} ${SING_IMAGE} killMS2H5parm.py lotss_smoothed.h5 ${FOLDER}/extract/DDS3_full*smoothed.npz --nofulljones
singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/h5_merger.py -out lotss_full_merged.h5 -in lotss_*.h5 -ms '/net/tussenrijn/data2/jurjendejong/L626678/result/*.goodtimes' --convert_tec 0

#h5 filter
singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/supporting_scripts/h5_filter.py -f ${FOLDER}/extract/image_full_ampphase_di_m.NS.app.restored.fits -ac 2.5 -in false -h5out lotss_full_merged_filtered.h5 -h5in lotss_full_merged.h5
#singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/supporting_scripts/h5_filter.py -f ${FOLDER}/extract/image_full_ampphase_di_m.NS.app.restored.fits -ac 2.5 -in true -h5out all_directions_filtered.h5 -h5in ${FOLDER}/result/all_directions.h5

#merging h5 files
singularity exec -B ${SING_BIND} ${SING_IMAGE} python ${SCRIPT_FOLDER}/h5_merger.py -out complete_merged.h5 -in lotss_full_merged_filtered.h5 all_directions_filtered.h5 --convert_tec 0
"""

#----------------------------------------------------------------------------------------------------------------------

#MAKE LIST WITH MEASUREMENT SETS
os.system('ls -1d /net/tussenrijn/data2/jurjendejong/A399/result/*.pre-cal.ms.archive* > /net/tussenrijn/data2/jurjendejong/A399/result/big-mslist.txt'.format(LOCATION=TO))

#----------------------------------------------------------------------------------------------------------------------

#MAKE DDF COMMAND
# with open('/home/jurjendejong/scripts/lofar_helpers/DDF_scripts/ddf.txt') as f:
#     lines = [l.replace('\n','') for l in f.readlines()]
#     lines+=['--Data-MS=/net/tussenrijn/data2/jurjendejong/A399/result/big-mslist.txt']
#     lines+=['--Predict-InitDicoModel=/net/tussenrijn/data2/jurjendejong/A399/result/image_full_ampphase_di_m.NS.DicoModel']
#     lines+=['--DDESolutions-DDSols=/net/tussenrijn/data2/jurjendejong/A399/result/complete_merged*.h5:sol000/amplitude000+phase000']
#     lines+=['--Mask-External=/net/tussenrijn/data2/jurjendejong/A399/result/image_full_ampphase_di_m.NS.mask01.fits']
#
# #RUN DDF COMMAND
# print('Running DDF COMMAND')
# os.system(' '.join(['cd', TO, '&&', SINGULARITY] + lines))
# print('Finished making new image')