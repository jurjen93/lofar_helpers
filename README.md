## LOFAR helper scripts

Note that these scripts are a part of the LOFAR pipeline, which is still in development.

Clone this repo with: ```git clone https://github.com/jurjen93/lofar_helpers.git```

### Make boxes for LOFAR self-calibration

With ```make_boxes.py``` it is possible to make boxes around sources that can be extracted and self-calibrated.\
Note that this script only works with Python 3.

#### Usage

Use ```make_boxes.py``` as a standalone script by running on the command line:\
```python make_boxes.py <FLAGS>``` \
You can use the following flags: \
* ```--file``` --> followed by the fits file name (and path)
* ```--location``` --> followed by the location (path) to store the data \
* ```--no_images``` --> don't save the images locally \
* ```--ds9``` --> interactive mode to validate the box selection in ds9
* ```--max_boxes``` --> max number of boxes

The script returns the following:
* directory with .reg region boxes.
* directory with box images, to check the quality of the boxes.

#### Details
The following steps are taken in the algorithm:
* Determine peak fluxes (>0.07). Close flux peaks are considered from the same source.
* Loop through all peak fluxes (sources) and check if these are already included in other images. If so, we stop here.
* Check if the source is 'interesting' by comparing the outlier fluxes in the images with the outliers in a gaussian filtered version. If not, we stop here.
* Reposition the image by taking into account the borders where we don't want to have high flux peaks or other sources. This is done by moving the center of the image and resizing.
* We also look for islands of flux, which we find by binning 100x100 pixels and summing over these pixels. When the flux exceeds 30, we extract a box from this position.
* If there are multiple peak fluxes (sources) we flag these as done from our list.
* Now we save the image of the box and the .reg.

### Merge solutions

With ```h5_merger.py``` it is possible to merge h5 solutions, which are a result from self-calibrating the different directions that are extracted with the boxes.\
This script is compatible with both Python 2 and Python 3.

### Special requirements

```h5_mergers.py``` relies heavily on losoto: https://github.com/revoltek/losoto
which works really well to work with .h5 files.

#### Usage

There are two ways to use this script. To use it within your script, you can
just import the most important function with:
```from h5_merger import merge_h5```

And use the following parameters:
* ```h5_out``` --> output file name
* ```h5_tables``` --> h5 tables to merge in list form or string (with/without wildcard)
* ```ms_files``` --> ms files to use in list form or string (with/without wildcard)
* ```convert_tec``` --> convert tec to phase (boolean)
* ```merge_all_in_one``` --> choose to merge al directions into one direction
* ```lin2circ``` --> linear to circular polarization transformation
* ```circ2lin``` --> circular to linear polarization transformation
* ```add_direction``` --> add default direction (phase 0 and amplitude 1)

Or you can run the script on the command line with
```python h5_merger.py <FLAGS>```
* ```-out``` --> output file name
* ```-in``` --> h5 tables to merge in list form or string (with/without wildcard)
* ```-ms``` --> ms files to use in list form or string (with/without wildcard)
* ```-ct``` --> convert tec to phase
* ```-merge_all_in_one``` --> choose to merge al directions into one direction
* ```--lin2circ``` --> linear to circular polarization transformation
* ```--circ2lin``` --> circular to linear polarization transformation
* ```--add_direction``` --> add default direction (phase 0 and amplitude 1)

#### Contact
Let me know if you are using this script and have any issues or suggestions for improvements.
You can contact me on: jurjendejong(AT)strw.leidenuniv.nl
