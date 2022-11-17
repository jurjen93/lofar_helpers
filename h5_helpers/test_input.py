"""
Script to give table structure
Written for busy week 17/11/2022 by Jurjen
"""

import tables

output = []

if __name__ == '__main__':
    h5 = tables.open_file(input('Give H5 file to check:\n'))

    solsets = list(h5.root._v_groups.keys())

    for solset in solsets:

        sols = h5.root._f_get_child(solset)

        soltabs = list(sols._v_children.keys())

        if not 'antenna' in soltabs:
            print('WARNING: There is no antenna table in '+solset)
        if not 'source' in soltabs:
            print('WARNING: There is no source table in '+solset)

        for soltab in soltabs:

            solt = sols._f_get_child(soltab)

            try:
                st_list = list(solt._v_children.keys())
                for st in st_list:
                    output.append('/'.join([solset, soltab, st]))
            except:
                output.append('/'.join([solset, soltab]))

print('Table structure is:\n'+'\n'.join(output))