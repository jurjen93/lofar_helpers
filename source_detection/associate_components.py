from astropy.io import fits
from astropy.table import Table
import numpy as np

"""
ds9 <FITS_IMAGE> -asinh -region components.reg -scale limits -0.000001 0.001 -cmap ch05m151008
"""

def error_prop(errors):
    return np.sqrt(np.sum(np.power(errors, 2)))

def associate(associate_components, table):
    """
    associate components by merging components

    :param associate_components: Example: [{4: [1, 2, 3, 4, 5]}, {10: [9, 10, 11]}]. This will merge component 1,2,3,4,5 with 4 as middle component
    :param fits table
    """

    t = Table.read(table, format='fits')
    columns = ['Source_id', 'RA', 'E_RA', 'DEC', 'E_DEC', 'Total_flux', 'E_Total_flux', 'Peak_flux', 'E_Peak_flux',
               'S_Code', 'Isl_rms']
    # Take Source_id, RA, E_RA, DEC, E_DEC main component
    # Total_flux sum of all components
    # E_Total_flux uncertainty all components
    # Peak flux max of all components
    # S_Code = 'M'

    t = t[columns]

    to_delete = []
    to_not_delete = []

    for p in associate_components:
        main_ID = list(p.keys())[0]
        to_not_delete.append(main_ID)
        t[main_ID]['Total_flux'] = t[p[main_ID]]['Total_flux'].sum()
        t[main_ID]['E_Total_flux'] = error_prop(t[p[main_ID]]['E_Total_flux'])
        t[main_ID]['Peak_flux'] = t[p[main_ID]]['E_Peak_flux'].max()
        t[main_ID]['S_Code'] = 'M'
        t[main_ID]['Isl_rms'] = t[p[main_ID]]['Isl_rms'].mean()
        for i in p[main_ID]:
            to_delete.append(i)

    for i in sorted(to_delete)[::-1]:
        if i not in to_not_delete:
            del t[i]

    t.write(table.replace('.fits', '_final.fits'), format='fits')

