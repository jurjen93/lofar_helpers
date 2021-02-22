### Make boxes for LOFAR self-calibration

With this script it is possible to make boxes around sources that can be extracted and self-calibrated.
Note that this script is just a part of the LOFAR pipeline, which is still in development.

#### Usage

Clone this repo.
Use ```make_boxes.py``` as a standalone script by running on the command line:\
```python make_boxes.py -f <fits_file>``` \
The script returns the following:
* directory with .reg region boxes.
* directory with box images, to check the quality of the boxes.

#### Details
The following steps are taken in the algorithm:
* Determine peak fluxes (>0.07). Close flux peaks are considered from the same source.
* Loop through all peak fluxes (sources) and check if these are already included in other images. If so, we stop here.
* Check if the source is 'interesting' by comparing the outlier fluxes in the images with the outliers in a gaussian filtered version. If not, we stop here.
* Reposition the image by taking into account the borders where we don't want to have high flux peaks or other sources. This is done by moving the center of the image and resizing.
* If the flux is >0.07 and <0.1 we accept the box only if there are other sources with flux >0.07.
* If there are multiple peak fluxes (sources) we flag these as done from our list.
* Now we save the image of the box and the .reg.

#### Contact
Let me know if you are using this script and have any issues.
You can contact me on: jurjendejong(AT)strw.leidenuniv.nl
