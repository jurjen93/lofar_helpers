from astropy.table import Table
import numpy as np

"""
ds9 <FITS_IMAGE> -asinh -region components.reg -scale limits -0.000001 0.001 -cmap ch05m151008
"""


def error_prop(errors):
    return np.sqrt(np.sum(np.power(errors, 2)))


def get_table_index(t, source_id):
    return int(np.argwhere(t['Source_id'] == source_id).squeeze())


def associate(associate_components, table):
    """
    associate components by merging components

    :param associate_components: Example: [{4: [1, 2, 3, 4, 5]}, {10: [9, 10, 11]}]. This will merge component 1,2,3,4,5 with 4 as middle component
    :param fits table
    """

    t = Table.read(table, format='fits')
    # columns = ['Source_id', 'RA', 'E_RA', 'DEC', 'E_DEC', 'Total_flux', 'E_Total_flux', 'Peak_flux', 'E_Peak_flux',
    #            'S_Code', 'Isl_rms']
    # Take Source_id, RA, E_RA, DEC, E_DEC main component
    # Total_flux sum of all components
    # E_Total_flux uncertainty all components
    # Peak flux max of all components
    # S_Code = 'M'

    # t = t[columns]

    to_delete = []
    to_not_delete = []

    for p in associate_components:
        if type(p) == dict:
            main_ID = list(p.keys())[0]
            main_ID = get_table_index(t, main_ID)
            to_not_delete.append(main_ID)
            ids = [get_table_index(t, x) for x in list(set(p[main_ID]))]
            t[main_ID]['Total_flux'] = t[ids]['Total_flux'].sum()
            t[main_ID]['E_Total_flux'] = error_prop(t[ids]['E_Total_flux'])
            t[main_ID]['Peak_flux'] = t[ids]['Peak_flux'].max()
            t[main_ID]['E_Peak_flux'] = error_prop(t[ids]['E_Peak_flux'])
            t[main_ID]['S_Code'] = 'M'
            for i in ids:
                to_delete.append(i)
        if type(p) == int:
            to_delete.append(get_table_index(t, p))

    to_delete = list(set(to_delete))
    to_not_delete = list(set(to_not_delete))
    for i in sorted(to_delete)[::-1]:
        if i not in to_not_delete:
            del t[i]

    t.write(table.replace('.fits', '_final.fits'), format='fits', overwrite=True)
