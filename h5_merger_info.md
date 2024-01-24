## Merge solutions with ```h5_merger.py```

-------------------------------

With ```h5_merger.py``` it is possible to merge H5parm solution files (typically with extension .h5).
Within these solution sets there has to be a source, antenna, and solXXX tables. 
They should also contain 'amplitude000' and/or 'phase000' solution tables.
This script has also other functionalities, such as time and frequency averaging or linear to circular conversion and vice versa. 
These can be independently used from the merging functionality. You can find all the options below.

```h5_merger.py``` is compatible with both Python 2 and Python 3 but is likely most reliable with Python 3.

### Requirements

```h5_mergers.py``` uses the following libraries:
* ```losoto``` (install with ```pip install lososto```)
* ```casacore.tables``` (install with ```pip install python-casacore```)
* ```numpy```
* ```scipy```
* ```tables```

#### Usage

Import the main function in Python with \
```from h5_merger import merge_h5```\
Or use the script on the command line with \
```python h5_merger.py <PARAM>```\
List all the options in your version of ```h5_merger.py``` with:\
```python h5_merger.py -h```

You can also import the functionalities to perform the output check and to get the summary of the solution file separately:\
```from h5_merger import h5_check```\
or run the script on the command line\
```from h5_merger import output_check``` \
where you can insert a solution file in the function.

This script is being used in:\
https://github.com/rvweeren/lofar_facet_selfcal \
https://github.com/lmorabit/lofar-vlbi

-------------------------------

# Examples

Examples below demonstrate how to use h5_merger.py for specific cases:

#### 1) Merge series of h5 files:

```python h5_merger.py -in *.h5 -out out.h5```\
By adding ```--propagate_flags``` it also propagate the weights from the input h5 file, which usually corresponds with the flagged data.

#### 2) Merge series of h5 files with frequency and time axis from specific measurement set:

```python h5_merger.py --propagate_flags -in *.h5 -out out.h5 -ms <YOUR_MS>```\
The input MS can also be multiple measurement sets (where the script will generate one time and freq axis convering all these MS).

#### 3) Merge series of h5 files with specific frequency and time axis from one h5:

```python h5_merger.py --propagate_flags -in *.h5 -out out.h5 --h5_time_freq <H5>```\
This will take the time and freq axis from H5 to interpolate to.
You can also give a boolean to ```--h5_time_freq``` which means that all input h5 files will be used to generate the time and freq axis from (this is the default if ```-ms``` is not used).

#### 4) Merge all h5 files in one direction:

```python h5_merger.py --propagate_flags -in *.h5 -out out.h5 --merge_all_in_one```\
The script merges by default the same directions with each other, but if the user wishes to merge all h5 files in one direction, you can add the ```--merge_all_in_one``` option.

#### 5) Merge all h5 files and do not convert tec:

```python h5_merger.py --propagate_flags -in *.h5 -out out.h5 --keep_tec```\
The script converts by default TEC input to phases (see Equation 1 in Sweijen et al +22), but user can also turn this conversion off.

#### 6) Convert circular to linear polarization and vice versa:

```python h5_merger.py --propagate_flags -in *.h5 -out out.h5 --circ2lin```\
or\
```python h5_merger.py --propagate_flags -in *.h5 -out out.h5 --lin2circ```\
```--circ2lin``` and ```-lin2circ``` can be added to any type of merge, as this conversion will be done at the end of the algorithm.

#### 7) Merge h5 solutions with different freq and time axis:

```python h5_merger.py --propagate_flags -in *.h5 -out out.h5 --merge_diff_freq```\
This is for example useful if you wish to merge solution files from different frequency bands together into 1 big solution file.

#### 8) Add Core Stations back to output:

```python h5_merger.py --propagate_flags -in *.h5 -out out.h5 -ms <YOUR_MS> --add_cs```\
This functionality will replace the super station (ST001) in the h5 file with the core stations from a given MS file.

#### 9) Return only stations from specific MS file:

```python h5_merger.py --propagate_flags -in *.h5 -out out.h5 -ms <YOUR_MS> --add_ms_stations```\
This functionality will return in the output H5 file only the stations that are present in the MS file. 
For international stations default values will be returned (amplitudes=1 for diagonals, amplitudes=0 for off-diagonals, phases=0).

Combinations of all the above are possible to merge more exotic cases, such as for example a merge like this:\
```python h5_merger.py -in *.h5 -out test.h5 -ms <YOUR_MS> --h5_time_freq=true --add_ms_stations --propagate_flags --merge_diff_freq --merge_all_in_one --keep_tec --circ2lin```

-------------------------------

#### Contact
Let me know if you are using this script or other scripts and have any issues or suggestions for improvements.
You can contact me on: jurjendejong(AT)strw.leidenuniv.nl