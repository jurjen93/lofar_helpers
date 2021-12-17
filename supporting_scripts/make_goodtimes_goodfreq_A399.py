import os
from glob import glob
import casacore.tables as ct


TO='/net/tussenrijn/data2/jurjendejong/A399_extracted_avg'
FROM='/net/tussenrijn/data2/jurjendejong/A399_extracted'


CUTFREQS = [5021107868.011121, 5021107864.005561]


for MS in glob(FROM+'/Abell399-401_extr*.ms.archive*.avg'):
    if ct.table(MS).getcol('TIME')[0] in CUTFREQS:
        print('Cutting freq for ' + MS)
        os.system("python /home/jurjendejong/scripts/lofar_helpers/supporting_scripts/flag_freq.py -ff='[15..19]' -msin " + MS+" -msout " + TO + '/' + MS.split('/')[-1] + '.goodfreq')
        os.system("python /home/jurjendejong/scripts/lofar_helpers/supporting_scripts/flag_time.py -tf 0 1500 -msin " + MS + " -msout " + TO + '/' + MS.split('/')[-1] + '.goodfreq.goodtimes')
        os.system("rm -rf " + TO + '/' + MS.split('/')[-1] + '.goodfreq')
    else:
        print('Cutting time for '+MS)
        os.system("python /home/jurjendejong/scripts/lofar_helpers/supporting_scripts/flag_time.py -tf 0 1500 -msin " + MS + " -msout " + TO + '/' + MS.split('/')[-1] + '.goodtimes')