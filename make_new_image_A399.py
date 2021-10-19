"""
Example:
    python ~/scripts/lofar_helpers/make_new_image.py
        -from /disks/paradata/shimwell/LoTSS-DR2/archive_other/L626678
        -to /net/tussenrijn/data2/jurjendejong/L626678/result
        -h5 complete_merged.h5
"""

__author__ = "Jurjen de Jong (jurjendejong@strw.leidenuniv.nl)"

import os
from glob import glob
from os import path
import casacore.tables as ct
import tables
import numpy as np

TO='/net/tussenrijn/data2/jurjendejong/A399/result'
FROM='/net/rijn/data2/jdejong/A399_DEEP'

#CREATE DESTINATION DIRECTORY IF NOT EXISTS
# if not path.exists(TO):
#     os.system('mkdir {LOCATION}'.format(LOCATION=TO))
#     print('Created {LOCATION}'.format(LOCATION=TO))

#CLEAN CACHE
os.system('CleanSHM.py')

#MOVE FILES
print('Moving files to '+TO)
# command = 'cp -r '+FROM+'/image_full_ampphase_di_m.NS.mask01.fits '+TO+ ' && '+\
#         'cp -r '+FROM+'/image_full_ampphase_di_m.NS.DicoModel '+TO+' && '+\
#         'scp lofarvwf-jdejong@spider.surfsara.nl:/project/lofarvwf/Share/jdejong/output/A399/selfcal/all_directions*.h5 '+TO+' && wait'
        # 'cp -r '+FROM+'/*_uv.pre-cal_*.pre-cal.ms.archive '+TO+' && wait'

# os.system(command)
print('Finished moving files')

#----------------------------------------------------------------------------------------------------------------------

#FLAG TIME AND FREQUENCY

# starting times fom measurement sets that have to be cutted for time
# CUTTIMES = [5019387068.011121, 5019387064.005561, 5017577408.011121, 5017577404.005561, 5020506668.011121, 5020506664.005561]
#
# #starting times for measurement sets that have to be cutted for freq
# CUTFREQS = [5021107868.011121, 5021107864.005561]
#
# for MS in glob(FROM+'/*.ms.archive'):
#     t = ct.table(MS)
#     time = t.getcol('TIME')[0]
#     t.close()
#     if time in CUTTIMES:
#         print('Cutting time for '+MS)
#         os.system("python /home/jurjendejong/scripts/lofar_helpers/supporting_scripts/flag_time.py -tf 0 3000 -msin " + MS + " -msout " + TO + '/' + MS.split('/')[-1] + '.goodtimes')
#     elif time in CUTFREQS:
#         print('Cutting freq for ' + MS)
#         if '127' not in MS:
#             os.system("cp -r " + MS + " " + TO)
#     else:
#         print('Copying for ' + MS)
#         os.system("cp -r " + MS + " " + TO)
#
# # important to wait until everything is ready before moving on
# while len(glob(FROM+'/*.ms.archive')) != len(glob(TO+'/*.pre-cal.ms.archive*'))+1:
#     print('TIME AND FREQUENCY FLAGGING')
#----------------------------------------------------------------------------------------------------------------------

#MERGE LOTSS OUTER EDGE

from supporting_scripts.get_DDS3 import get_DDS3

DDS3, DDS3_dict = get_DDS3(FROM)

soltable_times = {}
for soltable in glob('/net/tussenrijn/data2/jurjendejong/A399/result/all_directions*.h5'):
    tab = tables.open_file(soltable)
    t = tab.root.sol000.phase000.time[0]
    soltable_times.update({t: soltable})
    tab.close()  # close table

os.system('mkdir /net/tussenrijn/data2/jurjendejong/A399/H5files')
os.system(' && '.join(['cp '+s+' /net/tussenrijn/data2/jurjendejong/A399/H5files' for s in DDS3]))
command = []
for ms in DDS3_dict.items():
    print(ms)
    new_h5=[]
    for npz in ms[1]:
        command.append('killMS2H5parm.py ' + npz.split('/')[-1].replace('npz','h5 ') + npz + ' --nofulljones')
        new_h5.append(npz.split('/')[-1].replace('npz','h5 '))

    table = ct.table(ms[0])  # open table
    t = np.array(sorted(np.unique(table.getcol('TIME'))))[0]  # get first time element from measurement set
    table.close()  # close table
    diff = lambda ob_time: abs(ob_time - t)  # formula to compare difference between first time element of ms and observation times
    closest_value = min(list(soltable_times.keys()), key=diff)  # get closest value with lambda function
    h5 = soltable_times[closest_value]
    print(h5)
    print(closest_value)
    command.append('python /home/jurjendejong/scripts/lofar_helpers/h5_merger.py -out final_lotss_' + str(closest_value) + '.h5 -in ' + ' '.join(new_h5) + '--convert_tec 0')
    command.append('python /home/jurjendejong/scripts/lofar_helpers/supporting_scripts/h5_filter.py -f /net/rijn/data2/jdejong/A399_DEEP/image_full_ampphase_di_m.NS.app.restored.fits -ac 2.5 -in false -h5out lotss_full_merged_filtered_' + str(closest_value) + '.h5 -h5in final_lotss_' + str(closest_value) + '.h5')
    command.append('python /home/jurjendejong/scripts/lofar_helpers/h5_merger.py -out complete_merged_' + str(closest_value)+'.h5 -in ' + soltable_times[closest_value] + ' lotss_full_merged_filtered_' + str(closest_value) + '.h5 --convert_tec 0')
print('cd /net/tussenrijn/data2/jurjendejong/A399/H5files && '+' && '.join(command))
os.system('cd /net/tussenrijn/data2/jurjendejong/A399/H5files && '+' && '.join(command))
# os.system('mv /net/tussenrijn/data2/jurjendejong/A399/H5files/complete_merged*.h5 /net/tussenrijn/data2/jurjendejong/A399/result/ && wait')
# os.system('rm -rf /net/tussenrijn/data2/jurjendejong/A399/H5files')

#----------------------------------------------------------------------------------------------------------------------

#MAKE LIST WITH MEASUREMENT SETS
# os.system('ls -1d /net/tussenrijn/data2/jurjendejong/A399/result/*.pre-cal.ms.archive* > /net/tussenrijn/data2/jurjendejong/A399/result/big-mslist.txt'.format(LOCATION=TO))
#
# #----------------------------------------------------------------------------------------------------------------------
#
# #MAKE DDF COMMAND
# with open('/home/jurjendejong/scripts/lofar_helpers/DDF_scripts/ddf.txt') as f:
#     lines = [l.replace('\n','') for l in f.readlines()]
#     lines+=['--Data-MS=/net/tussenrijn/data2/jurjendejong/A399/result/big-mslist.txt']
#     lines+=['--Predict-InitDicoModel=/net/tussenrijn/data2/jurjendejong/A399/result/image_full_ampphase_di_m.NS.DicoModel']
#     lines+=['--DDESolutions-DDSols=/net/tussenrijn/data2/jurjendejong/A399/result/complete_merged*.h5:sol000/amplitude000+phase000']
#     lines+=['--Mask-External=/net/tussenrijn/data2/jurjendejong/A399/result/image_full_ampphase_di_m.NS.mask01.fits']
#
# #RUN DDF COMMAND
# print('Running DDF COMMAND')
# os.system(' '.join(['cd', TO, '&&'] + lines))
# print('Finished making new image')