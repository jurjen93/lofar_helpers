## LOFAR helper scripts

These scripts are helper script for LOFAR data reduction.\
Scripts were originally used for imaging and recalibrating A399-401 (see paper https://arxiv.org/abs/2209.13930) but can also be used for other purposes for reducing LOFAR data as well.

Clone repo with: ```git clone https://github.com/jurjen93/lofar_helpers.git``` \
If you are only interested in ```h5_merger.py```, use:\
```wget "https://raw.githubusercontent.com/jurjen93/lofar_helpers/master/h5_merger.py"```

The main scripts are described below:

-------------------------------
## 1. Merge solutions with ```h5_merger.py```

With ```h5_merger.py``` it is possible to merge H5parm solution files (typically with extension .h5).
Within these solution sets there has to be a source, antenna, and solXXX tables. 
They should also contain 'amplitude000' and/or 'phase000' solution tables.
This script has also other functionalities, such as time and frequency averaging or linear to circular conversion and vice versa. 
These can be independently used from the merging functionality. You can find all the options below.

```h5_merger.py``` is compatible with both Python 2 and Python 3 but is likely most reliable with Python 3.

### Special requirements

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

#### Contact
Let me know if you are using this script or other scripts and have any issues or suggestions for improvements.
You can contact me on: jurjendejong(AT)strw.leidenuniv.nl
