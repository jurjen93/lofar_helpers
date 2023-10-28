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

weighttest = ['weighttest/scalarphase3_skyselfcalcyle000_L816272_120_168MHz_averaged.ms.avg.h5',
              'weighttest/scalarcomplexgain4_skyselfcalcyle000_L816272_120_168MHz_averaged.ms.avg.h5',
              'weighttest/fulljones5_skyselfcalcyle000_L816272_120_168MHz_averaged.ms.avg.h5',
              'weighttest/scalarcomplexgain6_skyselfcalcyle000_L816272_120_168MHz_averaged.ms.avg.h5']

freqtimetest = glob('different_freqs/*.h5')

#time/freq coverage
merge_h5(h5_out='diff.h5', h5_tables=freqtimetest, ms_files='different_freqs/L693725_SB303_uv_12D4EA9ADt_132MHz.msdpppconcat.avg',
         h5_time_freq=True, lin2circ=False, circ2lin=False, add_directions=None, single_pol=None, no_pol=None,
         filtered_dir=None, add_cs=None, check_output=True, freq_av=None, time_av=None,
         check_flagged_station=True, propagate_flags=True, output_summary=True, add_ms_stations=True, merge_diff_freq=True)

H = tables.open_file('diff.h5')
shape = H.root.sol000.phase000.weight.shape
if shape[0]!=1618 or len(H.root.sol000.antenna[:])!=75 or np.all(H.root.sol000.amplitude000.val[:, : , -1, :, :]!=1):
    sys.exit('Freq/Time propagation wrong')
else:
    print("CORRECT TEST")

H.close()

#time/freq coverage
merge_h5(h5_out='diff.h5', h5_tables=freqtimetest,
         h5_time_freq=True, lin2circ=False, circ2lin=False, add_directions=None, single_pol=None, no_pol=None,
         filtered_dir=None, add_cs=None, check_output=True, freq_av=None, time_av=None,
         check_flagged_station=True, propagate_flags=True, output_summary=True, merge_diff_freq=True)

H = tables.open_file('diff.h5')
shape = H.root.sol000.phase000.weight.shape
if shape[0]!=1618 or len(H.root.sol000.antenna[:])!=61 or round(H.root.sol000.phase000.val[:][0,0,0,0,0],3)!=0.043:
    sys.exit('Freq/Time propagation wrong')
else:
    print("CORRECT TEST")

#weight merge
merge_h5(h5_out='weight.h5', h5_tables=weighttest,
             lin2circ=False, circ2lin=False, add_directions=None, single_pol=None, no_pol=None,
             filtered_dir=None, add_cs=None, check_output=True, freq_av=None, time_av=None,
             check_flagged_station=True, propagate_flags=True, output_summary=True)
H = tables.open_file('weight.h5')
weights = H.root.sol000.phase000.weight[-1, -1, -1, -1,:]
if np.sum(weights)!=0:
    sys.exit('Propagate weights incorrect')
else:
    print("CORRECT TEST")

H.close()

merge_h5(h5_out='weight.h5', h5_tables='circ2lin_weights/merged.h5',
             lin2circ=False, circ2lin=True, add_directions=None, single_pol=None, no_pol=None,
             filtered_dir=None, add_cs=None, check_output=None, freq_av=None, time_av=None,
             check_flagged_station=True, propagate_flags=True, output_summary=True)
H = tables.open_file('weight.h5')
F = tables.open_file('circ2lin_weights/merged.h5')
if not np.all(H.root.sol000.phase000.weight[:, 0, 0, 0, 0] == F.root.sol000.phase000.weight[:, 0, 0, 0, 0]):
    sys.exit('Propagate weights incorrect')
else:
    print("CORRECT TEST")

F.close()
H.close()

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
if np.all(np.array([0.0953,  0.5581, -2.9505, -0.0253]) == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values fulljones diff phase:\n'+str(values.round(4))+'\nVS\n'+str(np.array([0.0953,  0.5581, -2.9505, -0.0253])))

values = H.root.sol000.amplitude000.val[0,0,0,0,...]
if np.all(np.array([0.9597, 0.0432, 0.0124, 0.9722]) == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values fulljones diff amp:\n'+str(values.round(4))+'\nVS\n'+str(np.array([0.9597, 0.0432, 0.0124, 0.9722])))

H.close()

#fulljones circ2lin
merge_h5(h5_out='testfulljones.h5', h5_tables=fulljones, ms_files=None,
             lin2circ=False, circ2lin=True, add_directions=None, single_pol=None, no_pol=None,
             filtered_dir=None, add_cs=None, check_output=True, freq_av=None, time_av=None,
             check_flagged_station=True, propagate_flags=True, output_summary=True)
H = tables.open_file('testfulljones.h5')
values = H.root.sol000.phase000.val[0,0,0,0,...]
if np.all(np.array([ 0.0447, -2.5162, -0.2252,  0.0243]) == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values scalartest diff circ2lin phase:\n'+str(values.round(4))+'\nVS\n'+str(np.array([ 0.0447, -2.5162, -0.2252,  0.0243])))

values = H.root.sol000.amplitude000.val[0,0,0,0,...]
if np.all(np.array([0.9768, 0.0559, 0.0724, 0.9517]) == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values scalartest diff circ2lin amp:\n'+str(values.round(4))+'\nVS\n'+str(np.array([0.9768, 0.0559, 0.0724, 0.9517])))

H.close()

#fulljones lin2circ
merge_h5(h5_out='testdoublefulljones.h5', h5_tables=fulljones, ms_files=None,
             lin2circ=True, circ2lin=False, add_directions=None, single_pol=None, no_pol=None,
             filtered_dir=None, add_cs=None, check_output=True, freq_av=None, time_av=None,
             check_flagged_station=True, propagate_flags=True, output_summary=True)
H = tables.open_file('testdoublefulljones.h5')
values = H.root.sol000.phase000.val[0,0,0,0,...]
if np.all(np.array([0.0091, 1.8292, 1.5274, 0.0607]) == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values scalartest diff lin2circ phase:\n'+str(values.round(4))+'\nVS\n'+str(np.array([0.0091, 1.8292, 1.5274, 0.0607])))
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
             lin2circ=False, circ2lin=False, add_directions=None, single_pol=None, no_pol=True,
             filtered_dir=None, add_cs=None, check_output=True, freq_av=1, time_av=None,
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

merge_h5(h5_out='tectest.h5', h5_tables=glob('tectest2/*.h5'), ms_files=None,
             lin2circ=False, circ2lin=False, add_directions=None, single_pol=None, no_pol=None,
             filtered_dir=None, add_cs=None, check_output=True, freq_av=8, time_av=None,
             check_flagged_station=True, propagate_flags=True)
H = tables.open_file('tectest.h5')
values = H.root.sol000.phase000.val[0,-4,0,0,0]
if np.all(0.2857 == values.round(4)):
    print('CORRECT TEST')
else:
    sys.exit('Values tec diff:\n'+str(values.round(4))+'\nVS\n'+str(0.2857))
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
