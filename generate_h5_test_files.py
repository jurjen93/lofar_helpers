import numpy as np
import tables
from losoto.h5parm import h5parm



def make_soltab(solax={'pol': 2, 'dir': 1, 'ant': 10, 'freq': 20, 'time': 100}, soltab='phase000', solset='sol000', name='test.h5'):
  
   h5_out = h5parm(name, readonly=False)
   if solset in h5_out.getSolsetNames():
       solsetout = h5_out.getSolset(solset)
   else:
       solsetout = h5_out.makeSolset(solset)
   if 'amplitude' in soltab:
      values = np.ones(list(solax.values))-0.5
   else:
      values = np.zeros(list(solax.values))+0.5
   print(values)
   ax_vals=[]
   for sa in solax.keys():
      if sa=='pol':
         ax_vals += ['XX', 'YY']
      if sa=='dir':
         ax_vals += ['dir_1']
      if sa=='ant':
         ants=[] 
         for a in range(solax['ant']):
            ants+='antenna_{num}'.format(num=a)
         ax_vals+=ants
      if sa=='freq':
         ax_vals+=list(np.linspace(0, 1, solax['freq']))
      if sa=='time':
         ax_vals+=list(np.linspace(0,100, solax['time']))
   
   solsetout.makeSoltab(soltab, axesNames=list(solax.keys()), axesVals=ax_vals, vals=values,
                                 weights=np.ones(values.shape))
   h5_out.close()
   
   
make_soltab(soltab='phase000', solset='sol000', name='test_test.h5')