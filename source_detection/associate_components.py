from astropy.io import fits
from astropy.table import Table

def associate(associate_components, table):
    """
    associate components by merging components

    :param associate_components: Example: [{4: [1, 2, 3, 4, 5]}, {10: [9, 10, 11]}]. This will merge component 1,2,3,4,5 with 4 as middle component
    :param fits table
    """

    t = Table.read(table, format='fits')
    columns = ['Source_id', 'RA', 'E_RA', 'DEC', 'E_DEC', 'Total_flux', 'E_Total_flux', 'Peak_flux', 'E_Peak_flux',
               'S_Code', 'Isl_rms']
    t = t[columns]

    for p in associate_components:
        main_ID = list(p.keys())[0]
        subt = t[[t.Source_id==main_ID]]