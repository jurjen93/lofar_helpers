from astropy.table import Table
import numpy as np

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
            t[main_ID]['Maj'] = t[ids]['Maj'].min()
            t[main_ID]['Min'] = t[ids]['Min'].min()
            t[main_ID]['E_Maj'] = error_prop(t[ids]['E_Maj'])
            t[main_ID]['E_Min'] = error_prop(t[ids]['E_Min'])
            t[main_ID]['PA'] = t[ids]['PA'].min()
            t[main_ID]['E_PA'] = error_prop(t[ids]['E_PA'])
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
