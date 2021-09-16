import losoto.h5parm as h5parm
from h5_merger import merge_h5

testfile = '*.h5'
merge_h5(h5_tables=testfile, h5_out='test'+testfile, single_pol=True)

h5 = h5parm.h5parm(testfile)
h5.getSolset('sol000').getSoltab('phase000')
h5.getSolset('sol000').getSoltab('amplitude000')