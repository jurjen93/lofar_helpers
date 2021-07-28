import os
from argparse import ArgumentParser
import casacore.tables as ct
from numpy import unique
import pathlib
import operator
from glob import glob

def get_DDS3(folder):
    """
    This function returns from a folder the last DDS file for all the measurement sets within that same folder.
    :param folder: string of folder path
    """
    ms = glob(folder+'/*.ms.archive') # get all measurement sets in folder
    ms_observations = unique([m.split('/')[-1].split('_')[0] for m in ms]) # get unique observation names
    for i in '321': # get correct DDS*.npz files [should be DDS3*]
        DDS = glob(folder+'/DDS'+i+'*')
        if DDS: break

    time_modification = {D: pathlib.Path(D).stat().st_mtime for D in DDS} # get dict with DDS and time of modification
    time_observation = {D: int(D.split('full_')[1].split('.')[0]) for D in DDS if 'slow' not in D} # get dict with DDS and time of observation
    single_m = [[m for m in ms if night in m][0] for night in ms_observations] # get one measurement set per observation

    DDS_output = [] # DDS output
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

    return [file for sublist in DDS_output for file in sublist]

parser = ArgumentParser()
parser.add_argument('--frm', type=str, help='from where')
parser.add_argument('--to', type=str, help='to where')
args = parser.parse_args()

# FROM --> /disks/paradata/shimwell/LoTSS-DR2/archive_other/L626678
# TO --> /project/lofarvwf/Share/jdejong/data/L626678

files = [D.split('/')[-1] for D in get_DDS3(args.frm)] + \
        ['image_full_ampphase_di_m.NS.tessel.reg',
         'image_full_ampphase_di_m.NS.mask01.fits',
         'image_full_ampphase_di_m.NS.DicoModel',
         'image_dirin_SSD_m.npy.ClusterCat.npy',
         'SOLSDIR'] + \
        glob('*.ms.archive')

command_1 = 'tar -C {FROM} -cvzf /net/tussenrijn/data2/jurjendejong/data_archive.tar.gz --absolute-names '.format(FROM=args.frm) + ' '.join(files)
os.system(command_1)
command_2 = 'scp -r /net/tussenrijn/data2/jurjendejong/data_archive.tar.gz lofarvwf-jdejong@spider.surfsara.nl:{TO}'.format(TO=args.to)
os.system(command_2)
