## LOFAR helper scripts

These scripts are helper script for LOFAR data reduction.\
Scripts were originally used for imaging and recalibrating A399-401 (see paper https://arxiv.org/abs/2209.13930) but can also be used for other purposes for reducing LOFAR data as well.

Clone repo with: ```git clone https://github.com/jurjen93/lofar_helpers.git``` \
If you are only interested in ```h5_merger.py```, use:\
```wget "https://raw.githubusercontent.com/jurjen93/lofar_helpers/master/h5_merger.py"```

The main scripts are described below:

-------------------------------
## 1. Merge solutions with ```h5_merger.py```

With ```h5_merger.py``` it is possible to merge H5parm solution files (often ending on .h5), which are used in LOFAR (self-)calibration.
These files contain 'sol000' or 'sol001' or other solution set names. Within these solution sets there has to be a source and antenna table. 
They should also contain 'amplitude000' and/or 'phase000' solution tables. The script will most likely crash or not perform well when one of the above is not satisfied.
If you want to see the structure of your input solution file, you can run ```h5_helpers/test_input.py```.
This script has also other functionalities, such as time and frequency averaging or linear to circular conversion and vice versa. 
These can be independently used from the merging functionality. Find all the options below.

```h5_merger.py``` is compatible with both Python 2 and Python 3.

### Special requirements

```h5_mergers.py``` uses the following libraries:
* ```losoto``` (install with ```pip install lososto```)
* ```casacore.tables``` (install with ```pip install python-casacore```)
* ```numpy```
* ```scipy```
* ```tables```

#### Usage

There are two ways to use this script, by importing directly into your python script or by using the command line. 

#### Python import
Import the main function with: \
```from h5_merger import merge_h5```\
\
You can use these variables:
###### REQUIRED
* ```h5_out``` --> Solution file output name. This name cannot be in the list of input solution files.
* ```h5_tables``` --> Solution input files (can be both given as list with wildcard or string).
###### OPTIONAL
* ```ms_files``` --> Measurement set input files (can be both given as list with wildcard or string).
* ```h5_time_freq``` --> Solution file to use time and frequency arrays from. This is useful if the input solution files do not have the preferred time/freq resolution.
* ```convert_tec``` --> Convert TEC to phase.
* ```merge_all_in_one``` --> Merge all solutions into one single direction.
* ```lin2circ``` --> Transform linear polarization to circular.
* ```circ2lin``` --> Transform circular polarization to linear.
* ```add_directions``` --> Add direction with amplitude 1 and phase 0 (example: --add_direction [0.73,0.12]).
* ```single_pol``` --> Return only a single polarization axis if both polarizations are the same.
* ```no_pol``` --> Remove polarization axis if both polarizations are the same.
* ```use_solset``` --> Choose a solset to merge from your input solution files. Default is sol000.
* ```filtered_dir``` --> Filter out a list of indexed directions from your solution file. Only lists allowed (example: --filter_directions [2,3]).
* ```add_cs``` --> Add core stations to antenna output from measurement set (needs --ms).
* ```use_ants_from_ms``` --> Use only antenna stations from measurement set (needs --ms). Note that this is different to --add_cs, as it does not keep the international stations if these are not in the measurement set.
* ```time_av``` --> Time averaging factor.
* ```freq_av``` --> Frequency averaging factor.
* ```check_flagged_station``` --> Check if input stations are flagged, if so flag same stations in output.
* ```propagate_flags``` --> Interpolate weights and return in output file.
* ```no_antenna_check``` --> Do not compare antennas.
* ```check_output``` --> Check if the output has all the correct output information.
* ```output_summary``` --> Give summary of solution file content.
* ```merge_diff_freq``` --> Merging tables with different frequencies.

You can also import the functionalities to perform the output check and to get the summary of the solution file separately:\
```from h5_merger import h5_check```\
or \
```from h5_merger import output_check``` \
where you can insert a solution file in the function.

#### Command line
Use the script with the command line with the following parameters with ```python h5_merger.py <PARAM>```:
###### REQUIRED
* ```--h5_out``` --> Solution file output name. This name cannot be in the list of input H5 files.
* ```--h5_tables``` --> Solution input files (can be both given as list with wildcard or string).
###### OPTIONAL
* ```--ms``` --> Measurement set input files (can be both given as list with wildcard or string).
* ```--h5_time_freq``` --> Solution file to use time and frequency arrays from. This is useful if the input solution files do not have the preferred time/freq resolution.
* ```--time_av``` --> Time averaging factor.
* ```--freq_av``` --> Frequency averaging factor.
* ```--keep_tec``` --> Do not convert TEC to phase.
* ```--merge_all_in_one``` --> Merge all solutions into one single direction.
* ```--lin2circ``` --> Transform linear polarization to circular.
* ```--circ2lin``` --> Transform circular polarization to linear.
* ```--add_direction``` --> Add direction with amplitude 1 and phase 0 (example: --add_direction [0.73,0.12]).
* ```--single_pol``` --> Return only a single polarization axis if both polarizations are the same.
* ```--no_pol``` --> Remove polarization axis if both polarizations are the same.
* ```--usesolset``` --> Choose a solset to merge from your input H5 files. Default is sol000.
* ```--filter_directions``` --> Filter out a list of indexed directions from your H5 file. Only lists allowed (example: --filter_directions [2,3]).
* ```--add_cs``` --> Add core stations to antenna output from MS (needs --ms).
* ```--use_ants_from_ms``` --> Use only antenna stations from measurement set (needs --ms). Note that this is different to --add_cs, as it does not keep the international stations if these are not in the MS.
* ```--check_output``` --> Check if the output has all the correct output information.
* ```--not_flagstation``` --> Do not flag any station if station is flagged in input h5.
* ```--propagate_flags``` --> Interpolate weights and return in output file.
* ```--no_antenna_check``` --> Do not compare antennas.
* ```--output_summary``` --> Give output summary.
* ```--check_output``` --> Check if the output has all the correct output information.
* ```--merge_diff_freq``` --> Merging tables with different frequencies.

Or use the help function to list all the funcationalities in your version of ```h5_merger.py``` with:\
```python h5_merger.py -h```

This script is being used in:\
https://github.com/rvweeren/lofar_facet_selfcal \
https://github.com/lmorabit/lofar-vlbi

-------------------------------

#### Contact
Let me know if you are using this script or other scripts and have any issues or suggestions for improvements.
You can contact me on: jurjendejong(AT)strw.leidenuniv.nl
