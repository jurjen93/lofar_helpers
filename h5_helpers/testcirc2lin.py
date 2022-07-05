import pyrap.tables as pt
import numpy as np

t = pt.table('LBCS_120_168MHz_averaged.ms/')

datalin = t.getcol('DATA_LIN')
linapplied = t.getcol('LIN_APPLIED')

diff = datalin-linapplied
diff_amp = np.abs(datalin)-np.abs(linapplied)
diffangle = np.angle(datalin)-np.angle(linapplied)

print("Complex diff:\n"+str(diff))
print("Amplitude diff:\n"+str(diff_amp))
print("Phase diff:\n"+str(diffangle))