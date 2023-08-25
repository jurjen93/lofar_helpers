import tables
import casacore.tables as ct
from numpy import array

def add_antenna_from_ms(ms, h5, solset=None):
    """
    Insert antenna table from MS in given h5
    WARNING: this overwrites the antenna table from the h5

    :param ms: Measurement set
    :param h5: hdf5 solution file
    :param solset: specific solset (sol000, sol001,...)

    :return:
    """

    t = ct.table(ms + "::ANTENNA", ack=False)
    new_antlist = t.getcol('NAME')
    new_antpos = t.getcol('POSITION')
    t.close()
    values = array([list(zip(*(new_antlist, new_antpos)))], dtype=[('name', 'S16'), ('position', '<f4', (3,))])

    T = tables.open_file(h5, 'r+')
    if solset:
        T.root._f_get_child(solset).antenna._f_remove()
        T.create_table(T.root._f_get_child(solset), 'antenna', values, title='Antenna names and positions')

    else:
        for solset in T.root._v_groups.keys():
            T.root._f_get_child(solset).antenna._f_remove()
            T.create_table(T.root._f_get_child(solset), 'antenna', values, title='Antenna names and positions')
    T.close()