import sys
sys.path.append("..")
from h5_merger import merge_h5
from glob import glob
import tables
import numpy as np

fulljones = glob('fulljonestest/*.h5')

scalartest = glob('scalartest/*.h5')

tectest = glob('tectest/*.h5')


#fulljones
# merge_h5(h5_out='testfulljones.h5', h5_tables=fulljones, ms_files=None,
#              lin2circ=False, circ2lin=False, add_directions=None, single_pol=None, no_pol=None,
#              filtered_dir=None, add_cs=None, check_output=None, freq_av=None, time_av=None,
#              check_flagged_station=True, propagate_flags=None)
# H = tables.open_file('testfulljones.h5')
# values = H.root.sol000.phase000.val[0,0,0,0,...]
# if np.all(np.array([ 0.0953,  0.5581, -2.9505, -0.0253]) == values.round(4)):
#     print('CORRECT TEST')
# else:
#     sys.exit('Values fulljones diff phase:\n'+str(values.round(4))+'\nVS\n'+str(np.array([ 0.0953,  0.5581, -2.9505, -0.0253])))
#
# values = H.root.sol000.amplitude000.val[0,0,0,0,...]
# if np.all(np.array([0.9597, 0.0432, 0.0124, 0.9722]) == values.round(4)):
#     print('CORRECT TEST')
# else:
#     sys.exit('Values fulljones diff amp:\n'+str(values.round(4))+'\nVS\n'+str(np.array([0.9597, 0.0432, 0.0124, 0.9722])))
#
# H.close()
#
# #fulljones circ2lin
# merge_h5(h5_out='testfulljones.h5', h5_tables=fulljones, ms_files=None,
#              lin2circ=False, circ2lin=True, add_directions=None, single_pol=None, no_pol=None,
#              filtered_dir=None, add_cs=None, check_output=None, freq_av=None, time_av=None,
#              check_flagged_station=True, propagate_flags=None)
# H = tables.open_file('testfulljones.h5')
# values = H.root.sol000.phase000.val[0,0,0,0,...]
# if np.all(np.array([ 0.0447, -2.5162, -0.2252,  0.0243]) == values.round(4)):
#     print('CORRECT TEST')
# else:
#     sys.exit('Values scalartest diff circ2lin phase:\n'+str(values.round(4))+'\nVS\n'+str(np.array([ 0.0447, -2.5162, -0.2252,  0.0243])))
#
# values = H.root.sol000.amplitude000.val[0,0,0,0,...]
# if np.all(np.array([0.9768, 0.0559, 0.0724, 0.9517]) == values.round(4)):
#     print('CORRECT TEST')
# else:
#     sys.exit('Values scalartest diff circ2lin amp:\n'+str(values.round(4))+'\nVS\n'+str(np.array([0.9768, 0.0559, 0.0724, 0.9517])))
#
# H.close()
#
# #fulljones lin2circ
# merge_h5(h5_out='testfulljones.h5', h5_tables=fulljones, ms_files=None,
#              lin2circ=True, circ2lin=False, add_directions=None, single_pol=None, no_pol=None,
#              filtered_dir=None, add_cs=None, check_output=None, freq_av=None, time_av=None,
#              check_flagged_station=True, propagate_flags=None)
# H = tables.open_file('testfulljones.h5')
# values = H.root.sol000.phase000.val[0,0,0,0,...]
# if np.all(np.array([0.0091, 1.8292, 1.5274, 0.0607]) == values.round(4)):
#     print('CORRECT TEST')
# else:
#     sys.exit('Values scalartest diff lin2circ phase:\n'+str(values.round(4))+'\nVS\n'+str(np.array([0.0091, 1.8292, 1.5274, 0.0607])))
# H.close()

#fulljones lin2circ
merge_h5(h5_out='doublefulljones.h5',
         h5_tables=['testfulljones.h5']+glob('fulljonestest/*.h5'),
         ms_files=None,
             lin2circ=False, circ2lin=False, add_directions=None, single_pol=None, no_pol=None,
             filtered_dir=None, add_cs=None, check_output=None, freq_av=None, time_av=None,
             check_flagged_station=True, propagate_flags=None)
H = tables.open_file('doublefulljones.h5')
values = H.root.sol000.phase000.val[0,0,0,0,...]
if np.all(np.array([0.0091, 1.8292, 1.5274, 0.0607]) == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values scalartest diff lin2circ phase:\n'+str(values.round(4))+'\nVS\n'+str(np.array([0.0091, 1.8292, 1.5274, 0.0607])))
H.close()

#scalartest
merge_h5(h5_out='scalartest.h5', h5_tables=scalartest, ms_files=None,
             lin2circ=False, circ2lin=False, add_directions=None, single_pol=None, no_pol=None,
             filtered_dir=None, add_cs=None, check_output=None, freq_av=None, time_av=None,
             check_flagged_station=True, propagate_flags=None)
H = tables.open_file('scalartest.h5')
values = H.root.sol000.phase000.val[0,0,0,0,...]
if np.all(-0.0087 == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values scalartest diff:\n'+str(values.round(4))+'\nVS\n'+str(np.array([ 0.203 , -0.0087])))
H.close()

#tectest
merge_h5(h5_out='tectest.h5', h5_tables=tectest, ms_files=None,
             lin2circ=False, circ2lin=False, add_directions=None, single_pol=None, no_pol=None,
             filtered_dir=None, add_cs=None, check_output=None, freq_av=None, time_av=None,
             check_flagged_station=True, propagate_flags=None)
H = tables.open_file('tectest.h5')
values = H.root.sol000.phase000.val[0,-4,0,0,...]
if np.all(29.424 == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values tec diff:\n'+str(values.round(4))+'\nVS\n'+str(29.424))
H.close()

#no convert
merge_h5(h5_out='tecnoconverttest.h5', h5_tables=tectest, ms_files=None, convert_tec=False,
             lin2circ=False, circ2lin=False, add_directions=None, single_pol=None, no_pol=None,
             filtered_dir=None, add_cs=None, check_output=None, freq_av=None, time_av=None,
             check_flagged_station=True, propagate_flags=None)
H = tables.open_file('tecnoconverttest.h5')
values = H.root.sol000.tec000.val[0,-4,0,0,...]
if np.all(-0.5744 == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values tec diff:\n'+str(values.round(4))+'\nVS\n'+str(-0.5744))
H.close()
