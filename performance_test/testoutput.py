import sys
sys.path.append("..")
from h5_merger import merge_h5
from glob import glob
import tables
import numpy as np
from random import shuffle

fulljones = glob('fulljonestest/*.h5')

scalartest = glob('scalartest/*.h5')

tectest = glob('tectest/*.h5')

directiontest = glob('merging_directions/*.h5')


#direction merge
merge_h5(h5_out='directions.h5', h5_tables=directiontest,
             lin2circ=False, circ2lin=False, add_directions=None, single_pol=None, no_pol=None,
             filtered_dir=None, add_cs=None, check_output=True, freq_av=None, time_av=None,
             check_flagged_station=True, propagate_flags=True, output_summary=True)
H = tables.open_file('directions.h5')
values = H.root.sol000.phase000.val[0,0,0,0,...]
if np.all(np.array([0.0389, 0.0389]) == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values fulljones diff phase:\n'+str(values.round(4))+'\nVS\n'+str(np.array([0.0389, 0.0389])))

values = H.root.sol000.amplitude000.val[0,0,0,0,...]
if np.all(np.array([1.0285, 1.0285]) == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values fulljones diff amp:\n'+str(values.round(4))+'\nVS\n'+str(np.array([1.0285, 1.0285])))

H.close()

# MS test
merge_h5(h5_out='mstest.h5', h5_tables='merging_directions/merged_selfcalcyle011_L816272_P43479.ms.copy.phaseup.h5', ms_files='merging_directions/L816272_149MHz_P43479.ms',
             lin2circ=False, circ2lin=False, add_directions=None, single_pol=None, no_pol=None,
             filtered_dir=None, add_cs=None, check_output=True, freq_av=None, time_av=None,
             check_flagged_station=True, propagate_flags=True, output_summary=True)
H = tables.open_file('mstest.h5')
values = H.root.sol000.phase000.val[0,0,0,0,...]
if np.all(np.array([-1.2622, -1.2622]) == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values fulljones diff phase:\n'+str(values.round(4))+'\nVS\n'+str(np.array([-1.2622, -1.2622])))

values = H.root.sol000.amplitude000.val[0,0,0,0,...]
if np.all(np.array([1, 1]) == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values fulljones diff amp:\n'+str(values.round(4))+'\nVS\n'+str(np.array([1, 1])))

H.close()

merge_h5(h5_out='mstest.h5', h5_tables='merging_directions/merged_selfcalcyle011_L816272_P43479.ms.copy.phaseup.h5', ms_files='merging_directions/L816272_149MHz_P43479.ms',
             lin2circ=False, circ2lin=False, add_directions=None, single_pol=None, no_pol=None, h5_time_freq='merging_directions/merged_selfcalcyle011_L816272_P43479.ms.copy.phaseup.h5',
             filtered_dir=None, add_cs=None, check_output=True, freq_av=None, time_av=None,
             check_flagged_station=True, propagate_flags=True, output_summary=True)
H = tables.open_file('mstest.h5')
values = H.root.sol000.phase000.val[0,0,0,0,...]
if np.all(np.array([0.3111, 0.3111]) == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values fulljones diff phase:\n'+str(values.round(4))+'\nVS\n'+str(np.array([0.3111, 0.3111])))

values = H.root.sol000.amplitude000.val[0,0,0,0,...]
if np.all(np.array([1, 1]) == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values fulljones diff amp:\n'+str(values.round(4))+'\nVS\n'+str(np.array([1, 1])))

H.close()

#fulljones
merge_h5(h5_out='testfulljones.h5', h5_tables=fulljones, ms_files=None,
             lin2circ=False, circ2lin=False, add_directions=None, single_pol=None, no_pol=None,
             filtered_dir=None, add_cs=None, check_output=True, freq_av=None, time_av=None,
             check_flagged_station=True, propagate_flags=True, output_summary=True)
H = tables.open_file('testfulljones.h5')
values = H.root.sol000.phase000.val[0,0,0,0,...]
if np.all(np.array([ 0.0953, 0.5753, -2.9758, -0.0253]) == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values fulljones diff phase:\n'+str(values.round(4))+'\nVS\n'+str(np.array([0.0953, 0.5753, -2.9758, -0.0253])))

values = H.root.sol000.amplitude000.val[0,0,0,0,...]
if np.all(np.array([0.9597, 0.0415, 0.0121, 0.9722]) == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values fulljones diff amp:\n'+str(values.round(4))+'\nVS\n'+str(np.array([0.9597, 0.0415, 0.0121, 0.9722])))

H.close()

#fulljones circ2lin
merge_h5(h5_out='testfulljones.h5', h5_tables=fulljones, ms_files=None,
             lin2circ=False, circ2lin=True, add_directions=None, single_pol=None, no_pol=None,
             filtered_dir=None, add_cs=None, check_output=True, freq_av=None, time_av=None,
             check_flagged_station=True, propagate_flags=True, output_summary=True)
H = tables.open_file('testfulljones.h5')
values = H.root.sol000.phase000.val[0,0,0,0,...]
if np.all(np.array([ 0.0447, -2.5354, -0.2116,  0.0242]) == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values scalartest diff circ2lin phase:\n'+str(values.round(4))+'\nVS\n'+str(np.array([ 0.0447, -2.5354, -0.2116,  0.0242])))

values = H.root.sol000.amplitude000.val[0,0,0,0,...]
if np.all(np.array([0.976, 0.0556, 0.0719, 0.9525]) == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values scalartest diff circ2lin amp:\n'+str(values.round(4))+'\nVS\n'+str(np.array([0.976, 0.0556, 0.0719, 0.9525])))

H.close()

#fulljones lin2circ
merge_h5(h5_out='testdoublefulljones.h5', h5_tables=fulljones, ms_files=None,
             lin2circ=True, circ2lin=False, add_directions=None, single_pol=None, no_pol=None,
             filtered_dir=None, add_cs=None, check_output=True, freq_av=None, time_av=None,
             check_flagged_station=True, propagate_flags=True, output_summary=True)
H = tables.open_file('testdoublefulljones.h5')
values = H.root.sol000.phase000.val[0,0,0,0,...]
if np.all(np.array([0.0102, 1.8323, 1.5276, 0.0596]) == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values scalartest diff lin2circ phase:\n'+str(values.round(4))+'\nVS\n'+str(np.array([0.0102, 1.8323, 1.5276, 0.0596])))
H.close()

#double fulljones
merge_h5(h5_out='doublefulljones.h5',
         h5_tables=glob('doublefulljonestest/*dummy*.h5')+glob('fulljonestest/*.h5'),
         ms_files=None,
             lin2circ=False, circ2lin=False, add_directions=None, single_pol=None, no_pol=None,
             filtered_dir=None, add_cs=None, check_output=True, freq_av=None, time_av=None,
             check_flagged_station=True, propagate_flags=None, output_summary=True)
H = tables.open_file('doublefulljones.h5')
values = H.root.sol000.phase000.val[0,0,0,0,...]
if np.all(np.array([1.094, 0.9983, 1.094, 0.9983]) == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values scalartest double fulljones phase:\n'+str(values.round(4))+'\nVS\n'+str(np.array([1.094, 0.9983, 1.094, 0.9983])))
H.close()

#double fulljones
merge_h5(h5_out='doublefulljones.h5',
         h5_tables=glob('doublefulljonestest/*dummy*.h5'),
         ms_files=None,
             lin2circ=False, circ2lin=False, add_directions=None, single_pol=None, no_pol=None,
             filtered_dir=None, add_cs=None, check_output=True, freq_av=None, time_av=None,
             check_flagged_station=True, propagate_flags=None, output_summary=True)
H = tables.open_file('doublefulljones.h5')
values = H.root.sol000.phase000.val[0,0,0,0,...]
if np.all(np.array([1, 1, 1, 1]) == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values scalartest double fulljones phase:\n'+str(values.round(4))+'\nVS\n'+str(np.array([1, 1, 1, 1])))
H.close()

#double fulljones
merge_h5(h5_out='doublefulljones.h5',
         h5_tables=glob('doublefulljonestest/*dummy*.h5'),
         ms_files=None,
             lin2circ=False, circ2lin=False, add_directions=None, single_pol=None, no_pol=None,
             filtered_dir=None, add_cs=None, check_output=True, freq_av=None, time_av=None,
             check_flagged_station=True, propagate_flags=None, output_summary=True)
H = tables.open_file('doublefulljones.h5')
values = H.root.sol000.phase000.val[0,0,0,0,...]
if np.all(np.array([1, 1, 1, 1]) == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values scalartest double fulljones phase:\n'+str(values.round(4))+'\nVS\n'+str(np.array([1, 1, 1, 1])))
H.close()

#triple fulljones
merge_h5(h5_out='doublefulljones.h5',
         h5_tables=glob('doublefulljonestest/*.h5'),
         ms_files=None,
             lin2circ=False, circ2lin=False, add_directions=None, single_pol=None, no_pol=None,
             filtered_dir=None, add_cs=None, check_output=True, freq_av=None, time_av=None,
             check_flagged_station=True, propagate_flags=None, output_summary=True)
H = tables.open_file('doublefulljones.h5')
values = H.root.sol000.phase000.val[0,0,0,0,...]
if np.all(np.array([1.1318, 1.1184, 1.1318, 1.1184]) == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values scalartest double fulljones phase:\n'+str(values.round(4))+'\nVS\n'+str(np.array([1.1318, 1.1184, 1.1318, 1.1184])))
H.close()

#scalartest
merge_h5(h5_out='scalartest.h5', h5_tables=['scalartest/scalarphase1_selfcalcyle011_LBCS_120_168MHz_averaged.ms.avg.phaseup.h5'], ms_files=None,
             lin2circ=False, circ2lin=False, add_directions=None, single_pol=True, no_pol=None,
             filtered_dir=None, add_cs=None, check_output=True, freq_av=None, time_av=None,
             check_flagged_station=True, propagate_flags=None, output_summary=True)

merge_h5(h5_out='scalartest.h5', h5_tables=['scalartest/scalarphase1_selfcalcyle011_LBCS_120_168MHz_averaged.ms.avg.phaseup.h5'], ms_files=None,
             lin2circ=False, circ2lin=False, add_directions=None, single_pol=None, no_pol=True,
             filtered_dir=None, add_cs=None, check_output=True, freq_av=None, time_av=None,
             check_flagged_station=True, propagate_flags=None, output_summary=True)

merge_h5(h5_out='scalartest.h5', h5_tables=scalartest, ms_files=None,
             lin2circ=False, circ2lin=False, add_directions=None, single_pol=None, no_pol=None,
             filtered_dir=None, add_cs=None, check_output=True, freq_av=None, time_av=None,
             check_flagged_station=True, propagate_flags=True, output_summary=True)
H = tables.open_file('scalartest.h5')
values = H.root.sol000.phase000.val[0,0,0,0,...]
if np.all(-0.0087 == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values scalartest diff:\n'+str(values.round(4))+'\nVS\n'+str(-0.0087))
H.close()

shuffle(scalartest)

merge_h5(h5_out='scalartest.h5', h5_tables=scalartest, ms_files=None,
             lin2circ=False, circ2lin=False, add_directions=None, single_pol=None, no_pol=None,
             filtered_dir=None, add_cs=None, check_output=True, freq_av=None, time_av=None,
             check_flagged_station=True, propagate_flags=True, output_summary=True)
H = tables.open_file('scalartest.h5')
values = H.root.sol000.phase000.val[0,0,0,0,...]
if np.all(-0.0087 == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values scalartest diff:\n'+str(values.round(4))+'\nVS\n'+str(-0.0087))
H.close()

shuffle(scalartest)

merge_h5(h5_out='scalartest.h5', h5_tables=scalartest, ms_files=None,
             lin2circ=False, circ2lin=False, add_directions=None, single_pol=None, no_pol=None,
             filtered_dir=None, add_cs=None, check_output=True, freq_av=None, time_av=None,
             check_flagged_station=True, propagate_flags=True, output_summary=True)
H = tables.open_file('scalartest.h5')
values = H.root.sol000.phase000.val[0,0,0,0,...]
if np.all(-0.0087 == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values scalartest diff:\n'+str(values.round(4))+'\nVS\n'+str(-0.0087))
H.close()

#tectest
merge_h5(h5_out='tectest.h5', h5_tables=tectest[0], ms_files='tectest/basenametestP17*',
             lin2circ=False, circ2lin=False, add_directions=None, single_pol=None, no_pol=None,
             filtered_dir=None, add_cs=None, check_output=True, freq_av=8, time_av=None,
             check_flagged_station=True, propagate_flags=True)
H = tables.open_file('tectest.h5')
values = H.root.sol000.phase000.val[0,-4,0,0,...]
if np.all(8.6688 == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values tec diff:\n'+str(values.round(4))+'\nVS\n'+str(8.6688))
H.close()

merge_h5(h5_out='tectest.h5', h5_tables=tectest[1], ms_files=None,
             lin2circ=False, circ2lin=False, add_directions=None, single_pol=None, no_pol=None,
             filtered_dir=None, add_cs=None, check_output=True, freq_av=8, time_av=None,
             check_flagged_station=True, propagate_flags=True)
H = tables.open_file('tectest.h5')
values = H.root.sol000.phase000.val[0,-4,0,0,...]
if np.all(29.424 == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values tec diff:\n'+str(values.round(4))+'\nVS\n'+str(29.424))
H.close()

#no convert
merge_h5(h5_out='tecnoconverttest.h5', h5_tables=tectest[1], ms_files=None, convert_tec=False,
             lin2circ=False, circ2lin=False, add_directions=None, single_pol=None, no_pol=None,
             filtered_dir=None, add_cs=None, check_output=True, freq_av=3, time_av=2,
             check_flagged_station=True, propagate_flags=True)
H = tables.open_file('tecnoconverttest.h5')
values = H.root.sol000.tec000.val[0,-4,0,0,...]
if np.all(-0.5744 == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values tec diff:\n'+str(values.round(4))+'\nVS\n'+str(-0.5744))
H.close()
