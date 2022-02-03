import tables
import casacore.tables as ct
from numpy import array

def add_antenna_from_ms(ms, h5):
    """
    Add antenna table to h5 from corresponding MS
    WARNING: this overwrites the antenna table from the h5

    :param ms: Measurement set
    :param h5: hdf5 solution file
    """

    t = ct.table(ms + "::ANTENNA", ack=False)
    new_antlist = t.getcol('NAME')
    new_antpos = t.getcol('POSITION')
    t.close()
    values = array([list(zip(*(new_antlist, new_antpos)))], dtype=[('name', 'S16'), ('position', '<f4', (3,))])

    T = tables.open_file(h5, 'r+')
    T.root.sol000.antenna._f_remove()
    T.create_table(T.root.sol000, 'antenna', values, title='Antenna names and positions')
    T.close()