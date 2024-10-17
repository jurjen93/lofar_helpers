## LOFAR helper scripts

These scripts were developed for reducing LOFAR data for the following studies:
- Study of the bridge between galaxy clusters A399 and A401: https://arxiv.org/abs/2209.13930
- Producing the deepest sub-arcsecond wide-field image of ELAIS-N1: https://arxiv.org/pdf/2407.13247v1

Although the scripts were developed for specific purposes, 
I have tried to write them such they can be used for more general use-cases as well.
Some of these scripts have been integrated in data reduction pipelines as well.
You are free to use the scripts for your own work.

Clone the repo with: ```git clone https://github.com/jurjen93/lofar_helpers.git```

Or install with: ```pip install -v git+https://github.com/jurjen93/lofar_helpers.git```
This gives you the following command line functionalities:
- split_polygon_facets
- crop_nan_boundaries
- cut_fits_with_region
- close_h5
- concat_with_dummies
- remove_flagged_stations
- applycal
- subtract_with_wsclean

#### Contact
Let me know if you have any issues or suggestions for improvements.
You can contact me on: jurjendejong(AT)strw.leidenuniv.nl


#### Acknowledgments
This repository is part of the project CORTEX (NWA.1160.18.316) of the research programme NWA-ORC which is (partly) financed by the Dutch Research Council (NWO). 
