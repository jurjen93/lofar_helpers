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
                t = ct.table(sm).getcol('TIME')[0] # get first time element from measurement set
                diff = lambda ob_time : abs(ob_time - t) #f formula to compare difference between first time element of ms and observation times
                closest_value = min(list(time_observation.values()), key=diff) # get closest value with lambda function
                DDS_options = {D: time_modification[D] for D in
                               [f for f, time in time_observation.items() if closest_value == time]} # get optional DDS files
                try: # python 2
                    correct_DDS = max(DDS_options.iteritems(), key=operator.itemgetter(1))[0] # get correct DDS
                except: # python 3
                    correct_DDS = max(DDS_options.items(), key=operator.itemgetter(1))[0] # get correct DDS
                DDS_output.append([D for D in DDS if correct_DDS.split('full_')[1].split('_smoothed')[0] in D]) # append right DDS

    return DDS_output

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-f', '--folder', type=str, help='folder with the files')
    args = parser.parse_args()
    print(get_DDS3(args.folder))