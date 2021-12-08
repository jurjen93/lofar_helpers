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

TO='/net/vdesk/data2/jurjendejong/L626678_test'
FROM='/net/tussenrijn/data2/jurjendejong/A399_extracted'

#CREATE DESTINATION DIRECTORY IF NOT EXISTS
if not path.exists(TO):
    os.system('mkdir {LOCATION}'.format(LOCATION=TO))
    print('Created {LOCATION}'.format(LOCATION=TO))

#CLEAN CACHE
os.system('CleanSHM.py')

#MOVE FILES
# print('Moving files to '+TO)
# command = 'cp -r '+FROM+'/dicoMask.fits '+TO+ ' && '+\
#         'cp -r '+FROM+'/extr.DicoModel '+TO+' && wait'#+\
        # 'scp lofarvwf-jdejong@spider.surfsara.nl:/project/lofarvwf/Share/jdejong/output/A399/selfcal/all_directions*.h5 '+TO+' && wait'
        # 'cp -r '+FROM+'/*_uv.pre-cal_*.pre-cal.ms.archive '+TO+' && wait'

# os.system(command)
# print('Finished moving files')

#----------------------------------------------------------------------------------------------------------------------

#FLAG TIME AND FREQUENCY

# starting times fom measurement sets that have to be cutted for time
# CUTTIMES = [5019387068.011121, 5019387064.005561, 5017577408.011121, 5017577404.005561, 5020506668.011121, 5020506664.005561]

#starting times for measurement sets that have to be cutted for freq
CUTFREQS = [5021107868.011121, 5021107864.005561]

for MS in glob(FROM+'/extr_L626678*.ms'):
    t = ct.table(MS)
    time = t.getcol('TIME')[0]
    t.close()
    if not (time in CUTFREQS and '127' in MS):
        print('Making goodtimes for'+MS)
        os.system("python /home/jurjendejong/scripts/lofar_helpers/supporting_scripts/flag_time.py -tf 0 3000 -msin " + MS + " -msout " + TO + '/' + MS.split('/')[-1] + '.goodtimes')

#----------------------------------------------------------------------------------------------------------------------

#MERGE LOTSS OUTER EDGE

# from supporting_scripts.get_DDS3 import get_DDS3
#
# DDS3, DDS3_dict = get_DDS3(FROM)
#
# soltable_times = {}
# for soltable in glob(TO+'/all_directions*.h5'):
#     tab = tables.open_file(soltable)
#     t = tab.root.sol000.phase000.time[0]
#     soltable_times.update({t: soltable})
#     tab.close()  # close table
#
# os.system('mkdir '+TO+'/H5files')
# os.system(' && '.join(['cp '+s+' '+TO+'/H5files' for s in DDS3]))
# command = []
# for ms in DDS3_dict.items():
#     new_h5=[]
#     for npz in ms[1]:
#         command.append('killMS2H5parm.py ' + npz.split('/')[-1].replace('npz','h5 ') + npz + ' --nofulljones')
#         new_h5.append(npz.split('/')[-1].replace('npz','h5 '))
#
#     table = ct.table(ms[0])  # open table
#     t = np.array(sorted(np.unique(table.getcol('TIME'))))[0]  # get first time element from measurement set
#     table.close()  # close table
#     diff = lambda ob_time: abs(ob_time - t)  # formula to compare difference between first time element of ms and observation times
#     closest_value = min(list(soltable_times.keys()), key=diff)  # get closest value with lambda function
#     h5 = soltable_times[closest_value]
#     command.append('python /home/jurjendejong/scripts/lofar_helpers/h5_merger.py -out final_lotss_' + str(closest_value) + '.h5 -in ' + ' '.join(new_h5) + '--convert_tec 0')
#     command.append('python /home/jurjendejong/scripts/lofar_helpers/supporting_scripts/h5_filter.py -f /net/rijn/data2/jdejong/A399_DEEP/image_full_ampphase_di_m.NS.app.restored.fits -ac 2.5 -in false -h5out lotss_full_merged_filtered_' + str(closest_value) + '.h5 -h5in final_lotss_' + str(closest_value) + '.h5')
#     command.append('python /home/jurjendejong/scripts/lofar_helpers/h5_merger.py -out complete_merged_' + str(closest_value)+'.h5 -in ' + soltable_times[closest_value] + ' lotss_full_merged_filtered_' + str(closest_value) + '.h5 --convert_tec 0')
# print('cd ' + TO + '/H5files' + ' && ' + '\n'.join(command))
# os.system('cd ' + TO + '/H5files' + ' && '+' && '.join(command))
# os.system('mv ' + TO + '/H5files/complete_merged*.h5 '+TO+' && wait')
# os.system('rm -rf ' + TO + '/H5files')

#----------------------------------------------------------------------------------------------------------------------

#MAKE LIST WITH MEASUREMENT SETS
os.system('ls -1d '+TO+'/extr_L626678*.pre-cal.ms* > '+TO+'/big-mslist.txt'.format(LOCATION=TO))

#----------------------------------------------------------------------------------------------------------------------

#MAKE DDF COMMAND
with open('/'.join(__file__.split('/')[0:-1])+'/DDF_scripts/ddf.txt') as f:
    lines = [l.replace('\n','') for l in f.readlines()]
    lines+=['--Data-MS='+TO+'/big-mslist.txt']
    # lines+=['--Predict-InitDicoModel='+TO+'/extr.DicoModel']
    lines+=['--DDESolutions-DDSols='+TO+'/all_directions0.h5:sol000/amplitude000+phase000']
    # lines+=['--Mask-External='+TO+'/dicoMask.fits']
    lines+=['--Weight-ColName=WEIGHT_SPECTRUM']

#RUN DDF COMMAND
print('Running DDF COMMAND')
os.system(' '.join(['cd', TO, '&&'] + lines))
print('Finished making new image')