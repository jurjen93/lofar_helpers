## LOFAR helper scripts

These scripts were developed for reducing LOFAR data in the following peer-reviewed publications:
- Study of the bridge between galaxy clusters A399 and A401: https://arxiv.org/abs/2209.13930
- Producing the deepest sub-arcsecond wide-field image of ELAIS-N1: https://arxiv.org/pdf/2407.13247v1

Some of these scripts have been integrated in the LOFAR-VLBI pipeline: https://github.com/LOFAR-VLBI \
You are free to use the scripts for your own work.

Install with: ```pip install -v git+https://github.com/jurjen93/lofar_helpers.git``` \
This gives you the following command line functions:
- ```split_polygon_facets```
- ```crop_nan_boundaries```
- ```cut_fits_with_region```
- ```close_h5```
- ```concat_with_dummies```
- ```remove_flagged_stations```
- `applycal`
- `subtract_with_wsclean`

Use ```--help``` to get more details about the functionalities.

#### Contact
Let me know if you have any issues or suggestions for improvements.
You can contact me on: jurjendejong(AT)strw.leidenuniv.nl


#### Acknowledgments
This repository is part of the project CORTEX (NWA.1160.18.316) of the research programme NWA-ORC which is (partly) financed by the Dutch Research Council (NWO). 
