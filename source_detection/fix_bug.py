from astropy.table import Table
import numpy as np

fixtable = "facet_19_source_catalog_final.fits"

def fix(fixtable):

    T = Table.read(fixtable, format='fits')
    Told = Table.read(fixtable.replace("_final",""), format='fits')

    subT = Told['Source_id', 'Peak_flux']

    def get_table_index(t, source_id):
        return int(np.argwhere(t['Source_id'] == source_id).squeeze())

    for sourceid, peakflux in subT:
        try:
            idx = get_table_index(T, sourceid)
            if T[idx]['Peak_flux'] != peakflux:
                T[idx]['Peak_flux'] = peakflux
        except:
            pass

    T.write(fixtable, format='fits', overwrite=True)
