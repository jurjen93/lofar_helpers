#!/usr/bin/env python


# do not predict sky second time in pertubation solve?
# only do scalarphasediff solve once?
# to do: log command into the FITS header
# auto solints tune, HBA-international slow
# solnorm fulljones fix?
# fulljones flagging and medamps not working correctly
# for LBA low S/N cases switch to [tec,tec] automatically
# avg to units of herz and seconds? (for input data that hass different averaging)
# BLsmooth not for gain solves opttion
# BLsmooth constant smooth for gain solves
# Look into Wiener/Kalman smoothing
# only trigger HBA upper band selection for sources outside the FWHM?
# add 1 Jy source in phase center option


# example:
# python /net/rijn/data2/rvweeren/LoTSS_ClusterCAL/runwscleanLBautoR.py -b box_18.reg --forwidefield --usewgridder --avgfreqstep=2 --avgtimestep=2 --smoothnessconstraint-list="[0.0,0.0,5.0]" --antennaconstraint-list="['core']" --solint-list=[1,20,120] --soltypecycles-list="[0,1,3]" --soltype-list="['tecandphase','tecandphase','scalarcomplexgain']" test.ms


#import logging
#logging.basicConfig(filename='selfcal.log', format='%(levelname)s:%(asctime)s ---- %(message)s', datefmt='%m/%d/%Y %H:%M:%S', #level=logging.DEBUG)

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(filename='selfcal.log', format='%(levelname)s:%(asctime)s ---- %(message)s', datefmt='%m/%d/%Y %H:%M:%S')
logger.setLevel(logging.DEBUG)


import matplotlib
matplotlib.use('Agg')
import os, sys
import numpy as np
import losoto
import losoto.lib_operations
import glob, time
from astropy.io import fits
import astropy.stats
from astroquery.skyview import SkyView
import pyrap.tables as pt
import os.path
from losoto import h5parm
import bdsf
import pyregion
import argparse
import pickle
import aplpy
import fnmatch
import tables
from astropy.io import ascii
import multiprocessing
import ast
from lofar.stationresponse import stationresponse
from itertools import product

# this function does not work, for some reason cannot modify the source table
#def copy_over_sourcedirection_h5(h5ref, h5):
   #Href = tables.open_file(h5ref, mode='r') 
   #ssdir = np.copy(Href.root.sol000.source[0]['dir']) 
   #Href.close()
   #H = tables.open_file(h5, mode='a') 
   #print(ssdir, H.root.sol000.source[0]['dir'])
   #H.root.sol000.source[0]['dir'] = np.copy(ssdir)
   #H.flush()
   #print(ssdir, H.root.sol000.source[0]['dir'])
   #H.close()
   #return

def create_mergeparmdbname(mslist, selfcalcycle):
   parmdblist = mslist[:]
   for ms_id, ms in enumerate(mslist):
     parmdblist[ms_id] = 'merged_selfcalcyle' + str(selfcalcycle).zfill(3) + '_' + ms + '.avg.h5'
   print('Created parmdblist', parmdblist)
   return parmdblist

def preapply(H5filelist, mslist):
   for ms in mslist:
      parmdb = time_match_mstoH5(H5filelist, ms)
      applycal(ms, parmdb, msincol='DATA',msoutcol='CORRECTED_DATA')
      os.system("taql 'update " + ms + " set DATA=CORRECTED_DATA'")
   return

def time_match_mstoH5(H5filelist, ms):
   t = pt.table(ms)
   timesms = np.unique(t.getcol('TIME'))
   t.close()
   H5filematch = None
  
   for H5file in H5filelist:
      H = tables.open_file(H5file, mode='r')    

      try:
        times = H.root.sol000.amplitude000.time[:]
      except:
        pass
      try:
        times = H.root.sol000.rotation000.time[:]
      except:
        pass 
      try:
        times = H.root.sol000.phase000.time[:]
      except:
        pass      
      try:
        times = H.root.sol000.tec000.time[:]
      except:
        pass
      if np.median(times) >= np.min(timesms) and np.median(times) <= np.max(timesms):
         print(H5file, 'overlaps in time with', ms)
         H5filematch = H5file
 
      H.close()

   if H5filematch == None:
      print('Cannot find matching H5file and ms')
      sys.exit()

   return H5filematch





def logbasicinfo(args, fitsmask, mslist, version, inputsysargs):
   
   logger.info(' '.join(map(str,inputsysargs)))
   
   logger.info('Version:                   ' + str(version))
   logger.info('Imsize:                    ' + str(args['imsize']))
   logger.info('Pixelscale:                ' + str(args['pixelscale']))
   logger.info('Niter:                     ' + str(args['niter']))
   logger.info('Uvmin:                     ' + str(args['uvmin']  ))
   logger.info('Multiscale:                ' + str(args['multiscale']))
   logger.info('No beam correction:        ' + str(args['no_beamcor']))
   logger.info('IDG:                       ' + str(args['idg']))
   logger.info('Widefield:                 ' + str(args['forwidefield']))
   logger.info('Flagslowamprms:            ' + str(args['flagslowamprms']))
   logger.info('flagslowphaserms:          ' + str(args['flagslowphaserms']))
   logger.info('Do linear:                 ' + str(args['dolinear']))
   logger.info('Do circular:               ' + str(args['docircular']))
   if args['boxfile'] != None:
     logger.info('Bobxfile:                  ' + args['boxfile'])
   logger.info('Mslist:                    ' + ' '.join(map(str,mslist)))
   logger.info('User specified clean mask: ' + str(fitsmask))
   logger.info('Threshold for MakeMask:    ' + str(args['maskthreshold']))
   logger.info('Briggs robust:             ' + str(args['robust']))
   logger.info('Imagename prefix:          ' + args['imagename'] + '_')

    
   return    

def max_area_of_island(grid):
    rlen, clen = len(grid), len(grid[0])
    def neighbors(r, c):
        """
        Generate the neighbor coordinates of the given row and column that
        are within the bounds of the grid.
        """
        for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            if (0 <= r + dr < rlen) and (0 <= c + dc < clen):
                yield r + dr, c + dc

    visited = [[False] * clen for _ in range(rlen)]
    def island_size(r, c):
        """
        Find the area of the land connected to the given coordinate.
        Return 0 if the coordinate is water or if it has already been
        explored in a previous call to island_size().
        """
        if grid[r][c] == 0 or visited[r][c]:
            return 0
        area = 1
        stack = [(r, c)]
        visited[r][c] = True
        while stack:
            for r, c in neighbors(*stack.pop()):
                if grid[r][c] and not visited[r][c]:
                    stack.append((r, c))
                    visited[r][c] = True
                    area += 1
        return area

    return max(island_size(r, c) for r, c in product(range(rlen), range(clen)))

def getlargestislandsize(fitsmask):
   hdulist = fits.open(fitsmask)
   data = hdulist[0].data
   max_area = max_area_of_island(data[0,0,:,:])
   hdulist.close()
   return max_area    

#print(getlargestislandsize('image_7-MFS-image.fits.mask.fits'))


def create_phase_slope(inmslist, incol='DATA', outcol='DATA_PHASE_SLOPE', ampnorm=False):
   if not isinstance(inmslist,list):
      inmslist = [inmslist] 
   for ms in inmslist:
     t = pt.table(ms, readonly=False, ack=True)    
     if outcol not in t.colnames():
       print('Adding',outcol,'to',ms)            
       desc = t.getcoldesc(incol)
       newdesc = pt.makecoldesc(outcol, desc)
       newdmi = t.getdminfo(incol)
       newdmi['NAME'] = 'Dysco' + outcol
       t.addcols(newdesc, newdmi) 
     data = t.getcol(incol)
     dataslope = np.copy(data)
     for ff in range(data.shape[1]-1):
       if ampnorm:
         dataslope[:,ff,0] = np.copy(np.exp(1j * (np.angle(data[:,ff,0])-np.angle(data[:,ff+1,0]))))
         dataslope[:,ff,3] = np.copy(np.exp(1j * (np.angle(data[:,ff,3])-np.angle(data[:,ff+1,3]))))
       else:
         dataslope[:,ff,0] = np.copy(np.abs(data[:,ff,0])*np.exp(1j * (np.angle(data[:,ff,0])-np.angle(data[:,ff+1,0]))))
         dataslope[:,ff,3] = np.copy(np.abs(data[:,ff,3])*np.exp(1j * (np.angle(data[:,ff,3])-np.angle(data[:,ff+1,3]))))
       
     #last freq set to second to last freq because difference reduces length of freq axis with one
     dataslope[:,-1,:] = np.copy(dataslope[:,-2,:])
     t.putcol(outcol, dataslope) 
     t.close()
     #print( np.nanmedian(np.abs(data)))
     #print( np.nanmedian(np.abs(dataslope)))
     del data, dataslope
   return


def create_phasediff_column(inmslist, incol='DATA', outcol='DATA_CIRCULAR_PHASEDIFF'):
   if not isinstance(inmslist,list):
      inmslist = [inmslist] 
   for ms in inmslist:
     t = pt.table(ms, readonly=False, ack=True)    
     if outcol not in t.colnames():
       print('Adding',outcol,'to',ms)            
       desc = t.getcoldesc(incol)
       newdesc = pt.makecoldesc(outcol, desc)
       newdmi = t.getdminfo(incol)
       newdmi['NAME'] = 'Dysco' + outcol
       t.addcols(newdesc, newdmi) 
     
    
     data = t.getcol(incol)
     phasediff =  np.copy(np.angle(data[:,:,0]) - np.angle(data[:,:,3])) #RR - LL
     data[:,:,0] = 0.5*np.exp(1j * phasediff) # because I = RR+LL/2 (this is tricky because we work with phase diff)
     data[:,:,3] = data[:,:,0]
     t.putcol(outcol, data) 
     t.close()
     del data
     del phasediff
     
     
     if False:
        #data = t.getcol(incol)
        #t.putcol(outcol, data)
        #t.close()
        
        time.sleep(2)
        cmd = "taql 'update " + ms + " set " 
        cmd += outcol + "[,0]=0.5*EXP(1.0i*(PHASE(" + incol + "[,0])-PHASE(" + incol + "[,3])))'"
        #cmd += outcol + "[,3]=" + outcol + "[,0],"
        #cmd += outcol + "[,1]=0+0i,"
        #cmd += outcol + "[,2]=0+0i'"
        print(cmd)
        os.system(cmd)
        cmd = "taql 'update " + ms + " set " 
        #cmd += outcol + "[,0]=0.5*EXP(1.0i*(PHASE(" + incol + "[,0])-PHASE(" + incol + "[,3]))),"
        cmd += outcol + "[,3]=" + outcol + "[,0]'"
        #cmd += outcol + "[,1]=0+0i,"
        #cmd += outcol + "[,2]=0+0i'"
        print(cmd)
        os.system(cmd)
   return

def create_phase_column(inmslist, incol='DATA', outcol='DATA_PHASEONLY'):
   if not isinstance(inmslist,list):
      inmslist = [inmslist] 
   for ms in inmslist:
     t = pt.table(ms, readonly=False, ack=True)    
     if outcol not in t.colnames():
       print('Adding',outcol,'to',ms)            
       desc = t.getcoldesc(incol)
       newdesc = pt.makecoldesc(outcol, desc)
       newdmi = t.getdminfo(incol)
       newdmi['NAME'] = 'Dysco' + outcol
       t.addcols(newdesc, newdmi) 
     data = t.getcol(incol)
     data[:,:,0] = np.copy(np.exp(1j * np.angle(data[:,:,0]))) # because I = xx+yy/2
     data[:,:,3] = np.copy(np.exp(1j * np.angle(data[:,:,3]))) # because I = xx+yy/2
     t.putcol(outcol, data) 
     t.close()
     del data
   return

def create_MODEL_DATA_PDIFF(inmslist):
   if not isinstance(inmslist,list):
      inmslist = [inmslist] 
   for ms in inmslist:
     os.system('DPPP msin=' + ms + ' msout=. msout.datacolumn=MODEL_DATA_PDIFF steps=[]')
     os.system("taql" + " 'update " + ms + " set MODEL_DATA_PDIFF[,0]=(0.5+0i)'") # because I = RR+LL/2 (this is tricky because we work with phase diff)
     os.system("taql" + " 'update " + ms + " set MODEL_DATA_PDIFF[,3]=(0.5+0i)'") # because I = RR+LL/2 (this is tricky because we work with phase diff)
     os.system("taql" + " 'update " + ms + " set MODEL_DATA_PDIFF[,1]=(0+0i)'")
     os.system("taql" + " 'update " + ms + " set MODEL_DATA_PDIFF[,2]=(0+0i)'")

def fulljonesparmdb(h5):
    H=tables.open_file(h5) 
    try:
        phase = H.root.sol000.phase000.val[:]
        amplitude = H.root.sol000.amplitude000.val[:]
        if phase.shape[-1] == 4 and amplitude.shape[-1] == 4:
            fulljones = True
        else:
            fulljones = False
    except:
        fulljones = False 
    H.close()
    return fulljones

def reset_gains_noncore(h5parm, keepanntennastr='CS'):
   fulljones = fulljonesparmdb(h5parm) # True/False
   hasphase = True
   hasamps  = True
   hasrotatation = True
   hastec = True
   
   H=tables.open_file(h5parm, mode='a')
   # figure of we have phase and/or amplitude solutions
   try:
     antennas = H.root.sol000.amplitude000.ant[:]
     axisn = H.root.sol000.amplitude000.val.attrs['AXES'].decode().split(',')
   except: 
      hasamps = False
   try:
     antennas = H.root.sol000.phase000.ant[:]
     axisn = H.root.sol000.phase000.val.attrs['AXES'].decode().split(',')
   except:
     hasphase = False
   try:
     antennas = H.root.sol000.tec000.ant[:]
     axisn = H.root.sol000.tec000.val.attrs['AXES'].decode().split(',')
   except:
     hastec = False
   try:
     antennas = H.root.sol000.rotation000.ant[:]
     axisn = H.root.sol000.rotation000.val.attrs['AXES'].decode().split(',')
   except:
     hasrotatation = False     

   if hasphase:
     phase = H.root.sol000.phase000.val[:]
   if hasamps:  
     amp = H.root.sol000.amplitude000.val[:]
   if hastec:
     tec = H.root.sol000.tec000.val[:]  
   if hasrotatation:
     rotation = H.root.sol000.rotation000.val[:]       
     
     
   for antennaid,antenna in enumerate(antennas):
     if antenna[0:2] != keepanntennastr:
       if hasphase:
         antennaxis = axisn.index('ant')  
         axisn = H.root.sol000.phase000.val.attrs['AXES'].decode().split(',')
         print('Resetting phase', antenna, 'Axis entry number', axisn.index('ant'))
         #print(phase[:,:,antennaid,...])
         if antennaxis == 0:
           phase[antennaid,...] = 0.0
         if antennaxis == 1:
           phase[:,antennaid,...] = 0.0
         if antennaxis == 2:
           phase[:,:,antennaid,...] = 0.0
         if antennaxis == 3:
           phase[:,:,:,antennaid,...] = 0.0  
         if antennaxis == 4:
           phase[:,:,:,:,antennaid,...] = 0.0
         #print(phase[:,:,antennaid,...])  
       if hasamps:
         antennaxis = axisn.index('ant')  
         axisn = H.root.sol000.amplitude000.val.attrs['AXES'].decode().split(',')
         print('Resetting amplitude', antenna, 'Axis entry number', axisn.index('ant'))
         if antennaxis == 0:
           amp[antennaid,...] = 1.0
         if antennaxis == 1:
           amp[:,antennaid,...] = 1.0
         if antennaxis == 2:
           amp[:,:,antennaid,...] = 1.0
         if antennaxis == 3:
           amp[:,:,:,antennaid,...] = 1.0  
         if antennaxis == 4:
           amp[:,:,:,:,antennaid,...] = 1.0
         if fulljones:  
           amp[...,1] = 0.0 # XY, assumpe pol is last axis
           amp[...,2] = 0.0 # YX, assume pol is last axis
           
       if hastec:
         antennaxis = axisn.index('ant')  
         axisn = H.root.sol000.tec000.val.attrs['AXES'].decode().split(',')
         print('Resetting TEC', antenna, 'Axis entry number', axisn.index('ant'))
         if antennaxis == 0:
           tec[antennaid,...] = 0.0
         if antennaxis == 1:
           tec[:,antennaid,...] = 0.0
         if antennaxis == 2:
           tec[:,:,antennaid,...] = 0.0
         if antennaxis == 3:
           tec[:,:,:,antennaid,...] = 0.0  
         if antennaxis == 4:
           tec[:,:,:,:,antennaid,...] = 0.0                         
       if hasrotatation:
         antennaxis = axisn.index('ant')  
         axisn = H.root.sol000.rotation000.val.attrs['AXES'].decode().split(',')
         print('Resetting rotation', antenna, 'Axis entry number', axisn.index('ant'))
         if antennaxis == 0:
           rotation[antennaid,...] = 0.0
         if antennaxis == 1:
           rotation[:,antennaid,...] = 0.0
         if antennaxis == 2:
           rotation[:,:,antennaid,...] = 0.0
         if antennaxis == 3:
           rotation[:,:,:,antennaid,...] = 0.0  
         if antennaxis == 4:
           rotation[:,:,:,:,antennaid,...] = 0.0     
   # fill values back in
   if hasphase:
     H.root.sol000.phase000.val[:] = np.copy(phase)
   if hasamps:  
     H.root.sol000.amplitude000.val[:] = np.copy(amp)
   if hastec:
     H.root.sol000.tec000.val[:] = np.copy(tec) 
   if hasrotatation:
     H.root.sol000.rotation000.val[:] = np.copy(rotatation)       
     
   H.flush()
   H.close()
   return

#reset_gains_noncore('merged_selfcalcyle11_testquick260.ms.avg.h5')
#sys.exit()

def phaseup(msinlist,datacolumn='DATA',superstation='core', parmdbmergelist=None):
  msoutlist = []
  for ms in msinlist:
    msout=ms + '.phaseup'
    msoutlist.append(msout)
    if os.path.isdir(msout):
      os.system('rm -rf ' + msout)
    cmd = "DPPP msin=" + ms + " msout.storagemanager=dysco steps=[add,filter] msout.writefullresflag=False "
    cmd += "msout=" + msout + " msin.datacolumn=" + datacolumn + " "
    cmd += "filter.type=filter filter.remove=True "
    cmd += "add.type=stationadder "
    if superstation == 'core':
      cmd += "add.stations={ST001:'CS*'} filter.baseline='!CS*&&*' "
    if superstation == 'superterp':
      cmd += "add.stations={ST001:'CS00[2-7]*'} filter.baseline='!CS00[2-7]*&&*' "  

    print(cmd)
    os.system(cmd)
  
  return msoutlist

def findfreqavg(ms, imsize, bwsmearlimit=1.0):
    
  t = pt.table(ms + '/SPECTRAL_WINDOW',ack=False)
  bwsmear = bandwidthsmearing(np.median(t.getcol('CHAN_WIDTH')), \
            np.min(t.getcol('CHAN_FREQ')[0]), np.float(imsize), verbose=False)
  nfreq = len(t.getcol('CHAN_FREQ')[0])
  t.close()
  avgfactor = 0
  
  for count in range(2,21): # try average values between 2 to 20
     if bwsmear  < (bwsmearlimit/np.float(count)): # factor X avg
        if nfreq % count == 0:
           avgfactor = count

  #if bwsmear  < (bwsmearlimit/2.): # factor 2 avg
    #if nfreq % 2 == 0:
      #avgfactor = 2

  #if bwsmear  < (bwsmearlimit/3.): # factor 3 avg
    #if nfreq % 3 == 0:
      #avgfactor = 3

  #if bwsmear  < (bwsmearlimit/4.): # factor 4 avg
    #if nfreq % 4 == 0:
      #avgfactor = 4

  #if bwsmear  < (bwsmearlimit/5.): # factor 5 avg
    #if nfreq % 5 == 0:
      #avgfactor = 5

  #if bwsmear  < (bwsmearlimit/6.): # factor 6 avg
    #if nfreq % 6 == 0:
      #avgfactor = 6      
      
  #if bwsmear  < (bwsmearlimit/7.): # factor 7 avg
    #if nfreq % 7 == 0:
      #avgfactor = 7      

  #if bwsmear  < (bwsmearlimit/8.): # factor 8 avg
    #if nfreq % 8 == 0:
      #avgfactor = 8

  #if bwsmear  < (bwsmearlimit/9.): # factor 9 avg
    #if nfreq % 9 == 0:
      #avgfactor = 9
      
  #if bwsmear  < (bwsmearlimit/10.): # factor 10 avg
    #if nfreq % 10 == 0:
      #avgfactor = 10      

  return avgfactor



def ntimesH5(H5file):
   # function to return number of timeslots in H5 solution
   H=tables.open_file(H5file, mode='r')
   try:
     times= H.root.sol000.amplitude000.time[:]
   except: # apparently no slow amps available
     try:
       times= H.root.sol000.phase000.time[:]
     except:
       print('No amplitude000 or phase000 solution found')  
       sys.exit()
   H.close()
   return len(times)


def create_backup_flag_col(ms, flagcolname='FLAG_BACKUP'):
    cname = 'FLAG'
    flags = []
    t = pt.table(ms, readonly=False, ack=True)
    if flagcolname not in t.colnames():
      flags = t.getcol('FLAG')  
      print('Adding flagging column',flagcolname,'to',ms)            
      desc = t.getcoldesc(cname)
      newdesc = pt.makecoldesc(flagcolname, desc)
      newdmi = t.getdminfo(cname)
      newdmi['NAME'] = flagcolname
      t.addcols(newdesc, newdmi)  
      t.putcol(flagcolname, flags)     

    t.close()
    del flags
    return    
    

def checklongbaseline(ms):
    t   = pt.table(ms + '/ANTENNA',ack=False)
    antennasms = list(t.getcol('NAME'))
    t.close()
    substr = 'DE' # to check if a German station is present, if yes assume this is long baseline data
    haslongbaselines =  any(substr in mystring for mystring in antennasms)
    print('Contains long baselines?', haslongbaselines)
    return haslongbaselines

def average(mslist, freqstep, timestep=None, start=0, msinnchan=None):
    # sanity check
    if len(mslist) != len(freqstep):
      print('Hmm, made a mistake with freqstep?')
      sys.exit()
    
    outmslist = []
    for ms_id, ms in enumerate(mslist):
      if (freqstep[ms_id] > 0) or (timestep != None) or (msinnchan != None): # if this is True then average
        msout = ms + '.avg'  
        cmd = 'DPPP msin=' + ms + ' msout.storagemanager=dysco steps=[av] av.type=averager '
        cmd+= 'msout='+ msout + ' msin.weightcolumn=WEIGHT_SPECTRUM msout.writefullresflag=False '
        if freqstep[ms_id] != None:
          cmd+='av.freqstep=' + str(freqstep[ms_id]) + ' '
        if timestep != None:  
          cmd+='av.timestep=' + str(timestep) + ' '
        if msinnchan != None:
           cmd+='msin.nchan=' + str(msinnchan) + ' ' 
        if start == 0:
          print('Average with default WEIGHT_SPECTRUM:', cmd)
          if os.path.isdir(msout):
            os.system('rm -rf ' + msout)
          os.system(cmd)

        msouttmp = ms + '.avgtmp'  
        cmd = 'DPPP msin=' + ms + ' msout.storagemanager=dysco steps=[av] av.type=averager '
        cmd+= 'msout='+ msouttmp + ' msin.weightcolumn=WEIGHT_SPECTRUM_SOLVE msout.writefullresflag=False '
        if freqstep[ms_id] != None:
          cmd+='av.freqstep=' + str(freqstep[ms_id]) + ' '
        if timestep != None:  
          cmd+='av.timestep=' + str(timestep) + ' '
        if msinnchan != None:
           cmd+='msin.nchan=' + str(msinnchan) + ' '           
        if start == 0:
          t = pt.table(ms)
          if 'WEIGHT_SPECTRUM_SOLVE' in t.colnames(): # check if present otherwise this is not needed
            t.close()   
            print('Average with default WEIGHT_SPECTRUM_SOLVE:', cmd)
            if os.path.isdir(msouttmp):
              os.system('rm -rf ' + msouttmp)
            os.system(cmd)
          
            # Make a WEIGHT_SPECTRUM from WEIGHT_SPECTRUM_SOLVE
            t  = pt.table(msout, readonly=False)
            print('Adding WEIGHT_SPECTRUM_SOLVE')
            desc = t.getcoldesc('WEIGHT_SPECTRUM')
            desc['name']='WEIGHT_SPECTRUM_SOLVE'
            t.addcols(desc)

            t2 = pt.table(msouttmp, readonly=True)
            imweights = t2.getcol('WEIGHT_SPECTRUM')
            t.putcol('WEIGHT_SPECTRUM_SOLVE', imweights)

            # Fill WEIGHT_SPECTRUM with WEIGHT_SPECTRUM from second ms
            t2.close()
            t.close() 

            # clean up
            os.system('rm -rf ' + msouttmp)
          
          
        outmslist.append(msout)
      else:
        outmslist.append(ms)  # so no averaging happened
    
    return outmslist

def makeh5templates(mslist, parmdb, phasesoltype, slowsoltype, solint_phase, solint_ap, nchan_phase, nchan_ap):
  for msnumber, ms in enumerate(mslist):
    if phasesoltype == 'scalarphase':
       runDPPP(ms, np.int(solint_ap[msnumber]), np.int(solint_phase[msnumber]), \
               np.int(nchan_phase[msnumber]), np.int(nchan_ap[msnumber]), \
               ms + parmdb + str(0).zfill(3) + '_polversion.h5' ,args['phase_soltype'], \
               preapplyphase=False, TEC=TEC, puretec=args['pure_tec'], maxiter=1)
    if slowsoltype != 'fulljones':
       runDPPP(ms, np.int(solint_ap[msnumber]), np.int(solint_phase[msnumber]), \
               np.int(nchan_phase[msnumber]), np.int(nchan_ap[msnumber]), \
               ms + parmdb + str(0).zfill(3) + '_slowgainversion.h5' ,'complexgain', \
               preapplyphase=False, TEC=False, puretec=False, maxiter=1)
    else:  
       runDPPP(ms, np.int(solint_ap[msnumber]), np.int(solint_phase[msnumber]), \
               np.int(nchan_phase[msnumber]), np.int(nchan_ap[msnumber]), \
               ms + parmdb + str(0).zfill(3) + '_slowgainversion.h5' ,'fulljones', \
               preapplyphase=False, TEC=False, puretec=False, maxiter=1)
    resetgains(ms + parmdb + str(0).zfill(3) + '_slowgainversion.h5')
  return


def tecandphaseplotter(h5, ms, outplotname='plot.png'):
    if not os.path.isdir('plotlosoto%s'  % ms): # needed because if this is the first plot this directory does not yet exist
      os.system('mkdir plotlosoto%s'  % ms)
    cmd = 'python plot_tecandphase.py  '
    cmd += '--H5file=' + h5 + ' --outfile=plotlosoto%s/%s_nolosoto.png' % (ms,outplotname)
    print(cmd)
    os.system(cmd)
    return

def runaoflagger(mslist):
    for ms in mslist:
       cmd = 'aoflagger ' + ms
       os.system(cmd)
    return


def applycal(ms, inparmdblist, msincol='DATA',msoutcol='CORRECTED_DATA'):

    # to allow both a list or a single file (string)
    if not isinstance(inparmdblist,list):
     inparmdblist = [inparmdblist]    
    
    cmd = 'DPPP numthreads='+ str(multiprocessing.cpu_count()) + ' msin=' + ms +' msout=. '
    cmd +='msin.datacolumn=' + msincol + ' '
    cmd += 'msout.datacolumn=' + msoutcol + ' msout.storagemanager=dysco '
    count = 0
    for parmdb in inparmdblist:
      if fulljonesparmdb(parmdb):
        cmd += 'ac' + str(count) +'.parmdb='+parmdb + ' '
        cmd += 'ac' + str(count) +'.type=applycal '  
        cmd += 'ac' + str(count) +'.correction=fulljones '
        cmd += 'ac' + str(count) +'.soltab=[phase000,amplitude000] '  
        count = count + 1
      else:  
        H=tables.open_file(parmdb) 
        try:
          phase = H.root.sol000.phase000.val[:]
          cmd += 'ac' + str(count) +'.parmdb='+parmdb + ' '
          cmd += 'ac' + str(count) +'.type=applycal '  
          cmd += 'ac' + str(count) +'.correction=phase000 '
          count = count + 1
        except:
          pass  
        
        try:
          phase = H.root.sol000.tec000.val[:]
          cmd += 'ac' + str(count) +'.parmdb='+parmdb + ' '
          cmd += 'ac' + str(count) +'.type=applycal '  
          cmd += 'ac' + str(count) +'.correction=tec000 '
          count = count + 1        
        except:
          pass  
        
        try:
          phase = H.root.sol000.rotation000.val[:]
          cmd += 'ac' + str(count) +'.parmdb='+parmdb + ' '
          cmd += 'ac' + str(count) +'.type=applycal '  
          cmd += 'ac' + str(count) +'.correction=rotation000 '
          count = count + 1        
        except:
          pass  
        
        try:
          phase = H.root.sol000.amplitude000.val[:]
          cmd += 'ac' + str(count) +'.parmdb='+parmdb + ' '
          cmd += 'ac' + str(count) +'.type=applycal '  
          cmd += 'ac' + str(count) +'.correction=amplitude000 '
          count = count + 1        
        except:
          pass  
      
        H.close()
    
    if count < 1:
        print('Something went wrong, cannot build the applycal command. H5 file is valid?')
        sys.exit(1)
    # build the steps command    
    cmd += 'steps=['
    for i in range(count):
      cmd += 'ac'+ str(i)
      if i < count-1: # to avoid last comma in the steps list
        cmd += ','
    cmd += ']'

    print('DPPP applycal:', cmd)
    os.system(cmd) 
    return


def inputchecker(args):

  if not os.path.isfile('lib_multiproc.py'):
    print('Cannot find lib_multiproc.py, file does not exist')
    sys.exit(1)
  if not os.path.isfile('h5_merger.py'):
    print('Cannot find h5_merger.py, file does not exist')
    sys.exit(1)
  if not os.path.isfile('plot_tecandphase.py'):
    print('Cannot find plot_tecandphase.py, file does not exist')
    sys.exit(1)
  if not os.path.isfile('lin2circ.py'):
    print('Cannot find lin2circ.py, file does not exist')
    sys.exit(1)    
  if not os.path.isfile('BLsmooth.py'):
    print('Cannot find BLsmooth.py, file does not exist')
    sys.exit(1)

  if not args['no_beamcor'] and args['idg']:
    print('beamcor=True and IDG=True is not possible')
    sys.exit(1)
  
  for antennaconstraint in args['antennaconstraint_list']:
    if antennaconstraint not in ['superterp', 'coreandfirstremotes','core', 'remote',\
                                 'all', 'international', 'alldutch', 'core-remote','coreandallbutmostdistantremotes'] \
                         and antennaconstraint != None:
      print('Invalid input, antennaconstraint can only be core, superterp, coreandfirstremotes, remote, alldutch, international, or all')
      sys.exit(1)

  for soltype in args['soltype_list']:
    if soltype not in ['complexgain','scalarcomplexgain','scalaramplitude','amplitudeonly', 'phaseonly',\
                       'fulljones', 'rotation', 'rotation+diagonal','tec','tecandphase','scalarphase',\
                       'scalarphasediff', 'phaseonly_phmin', 'rotation_phmin', 'tec_phmin',\
                       'tecandphase_phmin','scalarphase_phmin','scalarphase_slope','phaseonly_slope']:
      print('Invalid soltype input')
      sys.exit(1)    

  if args['boxfile'] != None:
    if not (os.path.isfile(args['boxfile'])):
      print('Cannot find boxfile, file does not exist')
      sys.exit(1)
      
  if args['fitsmask'] != None:
    if not (os.path.isfile(args['fitsmask'])):
      print('Cannot find fitsmask, file does not exist')
      sys.exit(1)      

  if args['skymodel'] != None:
    if not (os.path.isfile(args['skymodel'])) and not (os.path.isdir(args['skymodel'])):
      print('Cannot find skymodel, file does not exist')
      sys.exit(1)

  if args['docircular'] and args['dolinear']:
      print('Conflicting input, docircular and dolinear used')
      sys.exit(1)

  if which('DPPP') == None:
    print('Cannot find DPPP, forgot to source lofarinit.[c]sh?')
    sys.exit(1)

  # Check boxfile and imsize settings
  if args['boxfile'] == None and args['imsize'] == None:
    print('Incomplete input detected, either boxfile or imsize is required')
    sys.exit(1)
  if args['boxfile'] != None and args['imsize'] != None:
    print('Wrong input detected, both boxfile and imsize are set')
    sys.exit(1)

  if args['imager'] not in ['DDFACET', 'WSCLEAN']: 
    print('Wrong input detected for option --imager, should be DDFACET or WSCLEAN')
    sys.exit(1)  

  if args['phaseupstations'] != None:
    if args['phaseupstations'] not in ['core', 'superterp']:    
      print('Wrong input detected for option --phaseupstations, should be core or superterp')
      sys.exit(1)  

  if args['soltypecycles_list'][0] != 0:
     print('Wrong input detected for option --soltypecycles-list should always start with 0') 
     sys.exit(1)

  if len(args['soltypecycles_list']) != len(args['soltype_list']): 
     print('Wrong input detected, length soltypecycles-list does not match that of soltype-list') 
     sys.exit(1)
 
  for soltype_id, soltype in enumerate(args['soltype_list']):
    wronginput = False
    if soltype in ['tecandphase', 'tec', 'tec_phmin', 'tecandphase_phmin']:    
      try: # in smoothnessconstraint_list is not filled by the user
        if args['smoothnessconstraint_list'][soltype_id] > 0.0:
          print('smoothnessconstraint should be 0.0 for a tec-like solve')
          wronginput = True
      except:
        pass    
      if wronginput:
       sys.exit(1)    
  
  for smoothnessconstraint in args['smoothnessconstraint_list']:
    if smoothnessconstraint < 0.0:
      print('Smoothnessconstraint must be equal or larger than 0.0')
      sys.exit(1)
  for smoothnessreffrequency in args['smoothnessreffrequency_list']:
    if smoothnessreffrequency < 0.0:
      print('Smoothnessreffrequency must be equal or larger than 0.0')
      sys.exit(1)
  
  if (args['skymodel'] != None) and (args['skymodelpointsource']) !=None:
    print('Wrong input, you cannot use a separate skymodel file and then also set skymodelpointsource')
    sys.exit(1)
  if (args['skymodelpointsource'] != None):
    if (args['skymodelpointsource'] <= 0.0):
      print('Wrong input, flux density provided for skymodelpointsource is <= 0.0')
      sys.exit(1)
  if (args['msinnchan'] != None):
    if (args['msinnchan'] <= 0):
      print('Wrong input for msinnchan, must be larger than zero')
      sys.exit(1)

  if (args['skymodelpointsource'] != None) and (args['predictskywithbeam']):
    print('Combination of skymodelpointsource and predictskywithbeam not supported')
    print('Provide a skymodel file to predict the sky with the beam')
    sys.exit(1)

  if (args['wscleanskymodel'] != None) and (args['skymodelpointsource']) !=None:
    print('Wrong input, you cannot use a wscleanskymodel and then also set skymodelpointsource')
    sys.exit(1)

  if (args['wscleanskymodel'] != None) and (args['skymodel']) !=None:
    print('Wrong input, you cannot use a wscleanskymodel and then also set skymodel')
    sys.exit(1)

  if (args['wscleanskymodel'] != None) and (args['predictskywithbeam']):
    print('Combination of wscleanskymodel and predictskywithbeam not supported')
    print('Provide a skymodel component file to predict the sky with the beam')
    sys.exit(1)

  if (args['wscleanskymodel'] != None) and (args['imager'] == 'DDFACET'):
    print('Combination of wscleanskymodel and DDFACET as an imager is not supported')
    sys.exit(1)
  if (args['wscleanskymodel'] != None): 
    if len(glob.glob(args['wscleanskymodel'] + '-????-model.fits')) < 2:
      print('Not enough WSClean channel model images found')
      print(glob.glob(args['wscleanskymodel'] + '-????-model.fits'))
      sys.exit(1)
    if (args['wscleanskymodel'].find('/') != -1):
      print('wscleanskymodel contains a slash, not allowed, needs to be in pwd') 
      sys.exit(1)
    if (args['wscleanskymodel'].find('..') != -1):
      print('wscleanskymodel contains .., not allowed, needs to be in pwd')      
      sys.exit(1)  
  return
  

def get_uvwmax(ms):
    t = pt.table(ms)
    uvw = t.getcol('UVW')
    ssq = np.sqrt(np.sum(uvw**2, axis=1))
    print(uvw.shape)
    t.close()
    return np.max(ssq)    

def makeBBSmodelforTGSS(boxfile, fitsimage=None, pixelscale=None, imsize=None):
    r = pyregion.open(boxfile)
    if len(r[:]) > 1:
       print('Composite region file, not allowed') 
       sys.exit()
    tgsspixsize = 6.2
    xs = np.ceil((r[0].coord_list[2])*3600./tgsspixsize)
    ys = np.ceil((r[0].coord_list[3])*3600./tgsspixsize)
        

    phasecenter = getregionboxcenter(boxfile)
    phasecenterc = phasecenter.replace('deg','')
    print(xs, phasecenterc)
 
    if fitsimage == None:
        filename = SkyView.get_image_list(position=phasecenterc,survey='TGSS ADR1', pixels=np.int(xs))
        print(filename)
        if os.path.isfile(filename[0].split('/')[-1]):
          os.system('rm -f ' + filename[0].split('/')[-1])
        time.sleep(10)
        os.system('wget ' + filename[0])
        filename = filename[0].split('/')[-1]
        print(filename)
    else:
        filename = fitsimage  
    
    img = bdsf.process_image(filename,mean_map='zero', rms_map=True, rms_box = (100,10), \
                             frequency=150e6, beam=(25./3600,25./3600,0.0) )
    img.write_catalog(format='bbs', bbs_patches='single', outfile='tgss.skymodel', clobber=True)
    #bbsmodel = 'bla.skymodel'
    del img
    return 'tgss.skymodel'

def getregionboxcenter(regionfile, standardbox=True):
    """
    Extract box center of a DS9 box region. 
    Input is regionfile Return NDPPP compatible string for phasecenter shifting
    """
    r = pyregion.open(regionfile)
    
    if len(r[:]) > 1:
      print('Only one region can be specified, your file contains', len(r[:]))
      sys.exit() 
    
    if r[0].name != 'box':
      print('Only box region supported')
      sys.exit()
    
    ra  = r[0].coord_list[0]
    dec = r[0].coord_list[1]
    boxsizex = r[0].coord_list[2]
    boxsizey = r[0].coord_list[3]
    angle = r[0].coord_list[4]
    
    if standardbox:
      if boxsizex != boxsizey:
        print('Only a square box region supported, you have these sizes:', boxsizex, boxsizey)
        sys.exit()
      if np.abs(angle) > 1:
        print('Only normally oriented sqaure boxes are supported, your region is oriented under angle:', angle)
        sys.exit()   
    
    regioncenter =  ('{:12.8f}'.format(ra) + 'deg,' + '{:12.8f}'.format(dec) + 'deg').replace(' ', '')
    return regioncenter



def bandwidthsmearing(chanw, freq, imsize, verbose=True):

  R =  (chanw/freq)*(imsize/6.) # asume we have used 3 pixels per beam
  if verbose:
    print('R value for bandwidth smearing is:', R)
    if R > 1.:
      print('Warning, try to increase your frequency resolution, or lower imsize, to reduce the R value below 1')
  
  return R

def number_freqchan_h5(h5parm):
    H=tables.open_file(h5parm)
    
    try:
       freq = H.root.sol000.phase000.freq[:]
       #print('You solutions do not contain phase values')
    except:    
       pass
    
    try:
        freq = H.root.sol000.amplitude000.freq[:] # apparently we only have amplitudes
    except:
        pass
    
    try:
        freq = H.root.sol000.rotation000.freq[:] # apparently we only have rotatioon
    except:
        pass

    try:
        freq = H.root.sol000.tec000.freq[:] # apparently we only have rotatioon
    except:
        pass

    H.close()
    print('Number of frequency channels in this solutions file is:', len(freq))
    return len(freq)


def calculate_restoringbeam(mslist, LBA):
    
    if LBA: # so we have LBA
      restoringbeam = 15.
    else : # so we have HBA
      restoringbeam = 6.  
    
    return restoringbeam



def print_title(version):
    print("""
              __        ______    _______    ___      .______      
             |  |      /  __  \  |   ____|  /   \     |   _  \     
             |  |     |  |  |  | |  |__    /  ^  \    |  |_)  |    
             |  |     |  |  |  | |   __|  /  /_\  \   |      /     
             |  `----.|  `--'  | |  |    /  _____  \  |  |\  \----.
             |_______| \______/  |__|   /__/     \__\ | _| `._____|
                                                                   
               _______    ___       ______  _______ .___________.
              |   ____|  /   \     /      ||   ____||           |
              |  |__    /  ^  \   |  ,----'|  |__   `---|  |----`
              |   __|  /  /_\  \  |  |     |   __|      |  |     
              |  |    /  _____  \ |  `----.|  |____     |  |     
              |__|   /__/     \__\ \______||_______|    |__|     
                                                                 
         _______. _______  __       _______   ______      ___       __      
        /       ||   ____||  |     |   ____| /      |    /   \     |  |     
       |   (----`|  |__   |  |     |  |__   |  ,----'   /  ^  \    |  |     
        \   \    |   __|  |  |     |   __|  |  |       /  /_\  \   |  |     
    .----)   |   |  |____ |  `----.|  |     |  `----. /  _____  \  |  `----.
    |_______/    |_______||_______||__|      \______|/__/     \__\ |_______|
                                                                              
    
                      Reinout van Weeren (2021, A&A, in press)

                              Starting.........
          """)

    print('\n\nVERSION: ' + version + '\n\n')
    return

def makemslist(mslist):
    os.system('rm -rf mslist.txt')
    f=open('mslist.txt', 'w')
    for ms in mslist:
       f.write(str(ms)+'\n')
    f.close()
    return

def antennaconstraintstr(ctype, antennasms, HBAorLBA):
    antennasms = list(antennasms)
    #print(antennasms)
    if ctype != 'superterp' and ctype != 'core' and ctype != 'coreandfirstremotes' and \
       ctype != 'remote' and ctype != 'alldutch' and ctype != 'all' and \
       ctype != 'international' and ctype != 'core-remote' and ctype != 'coreandallbutmostdistantremotes' :
        print('Invalid input, ctype can only be "superterp" or "core"')
        sys.exit(1)
    if HBAorLBA == 'LBA':  
      if ctype == 'superterp':  
        antstr=['CS002LBA','CS003LBA','CS004LBA','CS005LBA','CS006LBA','CS007LBA','ST001']
      if ctype == 'core':
        antstr=['CS001LBA','CS002LBA','CS003LBA','CS004LBA','CS005LBA','CS006LBA','CS007LBA', \
                'CS011LBA','CS013LBA','CS017LBA','CS021LBA','CS024LBA','CS026LBA','CS028LBA', \
                'CS030LBA','CS031LBA','CS032LBA','CS101LBA','CS103LBA','CS201LBA','CS301LBA', \
                'CS302LBA','CS401LBA','CS501LBA','ST001']
      if ctype == 'coreandfirstremotes':
        antstr=['CS001LBA','CS002LBA','CS003LBA','CS004LBA','CS005LBA','CS006LBA','CS007LBA', \
                'CS011LBA','CS013LBA','CS017LBA','CS021LBA','CS024LBA','CS026LBA','CS028LBA', \
                'CS030LBA','CS031LBA','CS032LBA','CS101LBA','CS103LBA','CS201LBA','CS301LBA', \
                'CS302LBA','CS401LBA','CS501LBA','RS503LBA','RS305LBA','RS205LBA','RS306LBA', \
                'RS106LBA','ST001']
      if ctype == 'coreandallbutmostdistantremotes':
        antstr=['CS001LBA','CS002LBA','CS003LBA','CS004LBA','CS005LBA','CS006LBA','CS007LBA', \
                'CS011LBA','CS013LBA','CS017LBA','CS021LBA','CS024LBA','CS026LBA','CS028LBA', \
                'CS030LBA','CS031LBA','CS032LBA','CS101LBA','CS103LBA','CS201LBA','CS301LBA', \
                'CS302LBA','CS401LBA','CS501LBA','RS503LBA','RS305LBA','RS205LBA','RS306LBA', \
                'RS106LBA','RS307LBA','RS406LBA','RS407LBA','ST001']
      if ctype == 'remote':
        antstr=['RS503LBA','RS305LBA','RS205LBA','RS306LBA', 'RS310LBA','RS406LBA','RS407LBA',\
                'RS106LBA','RS307LBA','RS208LBA','RS210LBA', 'RS409LBA','RS508LBA','RS509LBA']
      if ctype == 'alldutch':
        antstr=['CS001LBA','CS002LBA','CS003LBA','CS004LBA','CS005LBA','CS006LBA','CS007LBA', \
                'CS011LBA','CS013LBA','CS017LBA','CS021LBA','CS024LBA','CS026LBA','CS028LBA', \
                'CS030LBA','CS031LBA','CS032LBA','CS101LBA','CS103LBA','CS201LBA','CS301LBA', \
                'CS302LBA','CS401LBA','CS501LBA','RS503LBA','RS305LBA','RS205LBA','RS306LBA', \
                'RS310LBA','RS406LBA','RS407LBA','RS106LBA','RS307LBA','RS208LBA','RS210LBA', \
                'RS409LBA','RS508LBA','RS509LBA', 'ST001']
      if ctype == 'all':
        antstr=['CS001LBA','CS002LBA','CS003LBA','CS004LBA','CS005LBA','CS006LBA','CS007LBA', \
                'CS011LBA','CS013LBA','CS017LBA','CS021LBA','CS024LBA','CS026LBA','CS028LBA', \
                'CS030LBA','CS031LBA','CS032LBA','CS101LBA','CS103LBA','CS201LBA','CS301LBA', \
                'CS302LBA','CS401LBA','CS501LBA','RS503LBA','RS305LBA','RS205LBA','RS306LBA', \
                'RS310LBA','RS406LBA','RS407LBA','RS106LBA','RS307LBA','RS208LBA','RS210LBA', \
                'RS409LBA','RS508LBA','RS509LBA', \
                'DE601LBA','DE602LBA','DE603LBA','DE604LBA', 'DE605LBA','DE609LBA','FR606LBA', \
                'SE607LBA','UK608LBA','PL610LBA','PL611LBA', 'PL612LBA','IE613LBA','LV614LBA','ST001']          
      if ctype == 'international':
        antstr=['DE601LBA','DE602LBA','DE603LBA','DE604LBA', 'DE605LBA','DE609LBA','FR606LBA', \
                'SE607LBA','UK608LBA','PL610LBA','PL611LBA', 'PL612LBA','IE613LBA','LV614LBA']    
      if ctype == 'core-remote':
        antstr1=['CS001LBA','CS002LBA','CS003LBA','CS004LBA','CS005LBA','CS006LBA','CS007LBA', \
                'CS011LBA','CS013LBA','CS017LBA','CS021LBA','CS024LBA','CS026LBA','CS028LBA', \
                'CS030LBA','CS031LBA','CS032LBA','CS101LBA','CS103LBA','CS201LBA','CS301LBA', \
                'CS302LBA','CS401LBA','CS501LBA','ST001']
        antstr2=['RS503LBA','RS305LBA','RS205LBA','RS306LBA', 'RS310LBA','RS406LBA','RS407LBA',\
                'RS106LBA','RS307LBA','RS208LBA','RS210LBA', 'RS409LBA','RS508LBA','RS509LBA']
          

    if HBAorLBA == 'HBA':    
      if ctype == 'superterp': 
         antstr=['CS002HBA0','CS003HBA0','CS004HBA0','CS005HBA0','CS006HBA0','CS007HBA0', \
                 'CS002HBA1','CS003HBA1','CS004HBA1','CS005HBA1','CS006HBA1','CS007HBA1','ST001']
      if ctype == 'remote':
        antstr=['RS503HBA','RS305HBA','RS205HBA','RS306HBA', 'RS310HBA','RS406HBA','RS407HBA', \
                'RS106HBA','RS307HBA','RS208HBA','RS210HBA', 'RS409HBA','RS508HBA','RS509HBA']
      if ctype == 'core':
        antstr=['CS001HBA0','CS002HBA0','CS003HBA0','CS004HBA0','CS005HBA0','CS006HBA0','CS007HBA0', \
                'CS011HBA0','CS013HBA0','CS017HBA0','CS021HBA0','CS024HBA0','CS026HBA0','CS028HBA0', \
                'CS030HBA0','CS031HBA0','CS032HBA0','CS101HBA0','CS103HBA0','CS201HBA0','CS301HBA0', \
                'CS302HBA0','CS401HBA0','CS501HBA0',\
                'CS001HBA1','CS002HBA1','CS003HBA1','CS004HBA1','CS005HBA1','CS006HBA1','CS007HBA1', \
                'CS011HBA1','CS013HBA1','CS017HBA1','CS021HBA1','CS024HBA1','CS026HBA1','CS028HBA1', \
                'CS030HBA1','CS031HBA1','CS032HBA1','CS101HBA1','CS103HBA1','CS201HBA1','CS301HBA1', \
                'CS302HBA1','CS401HBA1','CS501HBA1','ST001']
      if ctype == 'coreandfirstremotes':
        antstr=['CS001HBA0','CS002HBA0','CS003HBA0','CS004HBA0','CS005HBA0','CS006HBA0','CS007HBA0', \
                'CS011HBA0','CS013HBA0','CS017HBA0','CS021HBA0','CS024HBA0','CS026HBA0','CS028HBA0', \
                'CS030HBA0','CS031HBA0','CS032HBA0','CS101HBA0','CS103HBA0','CS201HBA0','CS301HBA0', \
                'CS302HBA0','CS401HBA0','CS501HBA0',\
                'CS001HBA1','CS002HBA1','CS003HBA1','CS004HBA1','CS005HBA1','CS006HBA1','CS007HBA1', \
                'CS011HBA1','CS013HBA1','CS017HBA1','CS021HBA1','CS024HBA1','CS026HBA1','CS028HBA1', \
                'CS030HBA1','CS031HBA1','CS032HBA1','CS101HBA1','CS103HBA1','CS201HBA1','CS301HBA1', \
                'CS302HBA1','CS401HBA1','CS501HBA1','RS503HBA' ,'RS305HBA' ,'RS205HBA' ,'RS306HBA',  \
                'RS106HBA','ST001']
      if ctype == 'coreandallbutmostdistantremotes':
        antstr=['CS001HBA0','CS002HBA0','CS003HBA0','CS004HBA0','CS005HBA0','CS006HBA0','CS007HBA0', \
                'CS011HBA0','CS013HBA0','CS017HBA0','CS021HBA0','CS024HBA0','CS026HBA0','CS028HBA0', \
                'CS030HBA0','CS031HBA0','CS032HBA0','CS101HBA0','CS103HBA0','CS201HBA0','CS301HBA0', \
                'CS302HBA0','CS401HBA0','CS501HBA0',\
                'CS001HBA1','CS002HBA1','CS003HBA1','CS004HBA1','CS005HBA1','CS006HBA1','CS007HBA1', \
                'CS011HBA1','CS013HBA1','CS017HBA1','CS021HBA1','CS024HBA1','CS026HBA1','CS028HBA1', \
                'CS030HBA1','CS031HBA1','CS032HBA1','CS101HBA1','CS103HBA1','CS201HBA1','CS301HBA1', \
                'CS302HBA1','CS401HBA1','CS501HBA1','RS503HBA' ,'RS305HBA' ,'RS205HBA' ,'RS306HBA',  \
                'RS106HBA','RS307HBA','RS406HBA','RS407HBA','ST001']
      if ctype == 'alldutch':
        antstr=['CS001HBA0','CS002HBA0','CS003HBA0','CS004HBA0','CS005HBA0','CS006HBA0','CS007HBA0', \
                'CS011HBA0','CS013HBA0','CS017HBA0','CS021HBA0','CS024HBA0','CS026HBA0','CS028HBA0', \
                'CS030HBA0','CS031HBA0','CS032HBA0','CS101HBA0','CS103HBA0','CS201HBA0','CS301HBA0', \
                'CS302HBA0','CS401HBA0','CS501HBA0',\
                'CS001HBA1','CS002HBA1','CS003HBA1','CS004HBA1','CS005HBA1','CS006HBA1','CS007HBA1', \
                'CS011HBA1','CS013HBA1','CS017HBA1','CS021HBA1','CS024HBA1','CS026HBA1','CS028HBA1', \
                'CS030HBA1','CS031HBA1','CS032HBA1','CS101HBA1','CS103HBA1','CS201HBA1','CS301HBA1', \
                'CS302HBA1','CS401HBA1','CS501HBA1', \
                'RS503HBA','RS305HBA','RS205HBA','RS306HBA', 'RS310HBA','RS406HBA','RS407HBA', \
                'RS106HBA','RS307HBA','RS208HBA','RS210HBA', 'RS409HBA','RS508HBA','RS509HBA','ST001']
      if ctype == 'all':
        antstr=['CS001HBA0','CS002HBA0','CS003HBA0','CS004HBA0','CS005HBA0','CS006HBA0','CS007HBA0', \
                'CS011HBA0','CS013HBA0','CS017HBA0','CS021HBA0','CS024HBA0','CS026HBA0','CS028HBA0', \
                'CS030HBA0','CS031HBA0','CS032HBA0','CS101HBA0','CS103HBA0','CS201HBA0','CS301HBA0', \
                'CS302HBA0','CS401HBA0','CS501HBA0',\
                'CS001HBA1','CS002HBA1','CS003HBA1','CS004HBA1','CS005HBA1','CS006HBA1','CS007HBA1', \
                'CS011HBA1','CS013HBA1','CS017HBA1','CS021HBA1','CS024HBA1','CS026HBA1','CS028HBA1', \
                'CS030HBA1','CS031HBA1','CS032HBA1','CS101HBA1','CS103HBA1','CS201HBA1','CS301HBA1', \
                'CS302HBA1','CS401HBA1','CS501HBA1', \
                'RS503HBA','RS305HBA','RS205HBA','RS306HBA', 'RS310HBA','RS406HBA','RS407HBA', \
                'RS106HBA','RS307HBA','RS208HBA','RS210HBA', 'RS409HBA','RS508HBA','RS509HBA', \
                'DE601HBA','DE602HBA','DE603HBA','DE604HBA', 'DE605HBA','DE609HBA','FR606HBA', \
                'SE607HBA','UK608HBA','PL610HBA','PL611HBA', 'PL612HBA','IE613HBA','LV614HBA','ST001']
      if ctype == 'international':
        antstr=['DE601HBA','DE602HBA','DE603HBA','DE604HBA', 'DE605HBA','DE609HBA','FR606HBA', \
                'SE607HBA','UK608HBA','PL610HBA','PL611HBA', 'PL612HBA','IE613HBA','LV614HBA']        
      if ctype == 'core-remote':
        antstr1=['CS001HBA0','CS002HBA0','CS003HBA0','CS004HBA0','CS005HBA0','CS006HBA0','CS007HBA0', \
                'CS011HBA0','CS013HBA0','CS017HBA0','CS021HBA0','CS024HBA0','CS026HBA0','CS028HBA0', \
                'CS030HBA0','CS031HBA0','CS032HBA0','CS101HBA0','CS103HBA0','CS201HBA0','CS301HBA0', \
                'CS302HBA0','CS401HBA0','CS501HBA0',\
                'CS001HBA1','CS002HBA1','CS003HBA1','CS004HBA1','CS005HBA1','CS006HBA1','CS007HBA1', \
                'CS011HBA1','CS013HBA1','CS017HBA1','CS021HBA1','CS024HBA1','CS026HBA1','CS028HBA1', \
                'CS030HBA1','CS031HBA1','CS032HBA1','CS101HBA1','CS103HBA1','CS201HBA1','CS301HBA1', \
                'CS302HBA1','CS401HBA1','CS501HBA1','ST001']
        antstr2=['RS503HBA','RS305HBA','RS205HBA','RS306HBA', 'RS310HBA','RS406HBA','RS407HBA', \
                'RS106HBA','RS307HBA','RS208HBA','RS210HBA', 'RS409HBA','RS508HBA','RS509HBA']


    if ctype != 'core-remote':
        antstrtmp = list(antstr) # important to use list here, otherwise it's not a copy(!!) and antstrtmp refers to antstr
        for ant in antstr:
            if ant not in antennasms:
                antstrtmp.remove(ant)    
        antstr = ','.join(map(str, antstrtmp))
        antstr = '[[' + antstr + ']]'
    else:
        antstrtmp1 = list(antstr1)
        for ant in antstr1:
            if ant not in antennasms:
                antstrtmp1.remove(ant)    
        antstr1 = ','.join(map(str, antstrtmp1))

        antstrtmp2 = list(antstr2)
        for ant in antstr2:
            if ant not in antennasms:
                antstrtmp2.remove(ant)    
        antstr2 = ','.join(map(str, antstrtmp2))        
        
        antstr =  '[[' + antstr1 + '],[' + antstr2 + ']]'

    return antstr    


def makephasediffh5(phaseh5): 
    #note for scalarphase/phaseonly solve, does not work for tecandphase as freq axis is missing there for phase000
    H5pol = tables.open_file(phaseh5,mode='a')

    phase_pol = H5pol.root.sol000.phase000.val[:] # time, freq, ant, dir, pol
    phase_pol_tmp = np.copy(phase_pol)
    #antenna   = H5pol.root.sol000.phase000.ant[:]
    print('Shape to make phase diff array', phase_pol.shape)
    
    #for ant in range(len(antenna)):
    phase_pol[:, :, :, :,0]  = phase_pol_tmp[:, :, :, :,0] # XX
    phase_pol[:, :, :, :,-1] = 0.0*phase_pol_tmp[:, :, :, :,0] # YY


    H5pol.root.sol000.phase000.val[:] = phase_pol
    H5pol.flush()
    H5pol.close()
    return

def makephaseCDFh5(phaseh5): 
    #note for scalarphase/phaseonly solve, does not work for tecandphase as freq axis is missing there for phase000
    H5 = tables.open_file(phaseh5,mode='a')

    phaseCDF = H5.root.sol000.phase000.val[:] # time, freq, ant, dir, pol
    phaseCDF_tmp = np.copy(phaseCDF)
    print('Shape to make phase CDF array', phaseCDF.shape)
    nfreq = len(H5.root.sol000.phase000.freq[:])
    for ff in range(nfreq-1):
      #reverse order so phase increase towards lower frequnecies
      phaseCDF[:,nfreq-ff-2, ...]  = np.copy(phaseCDF[:,nfreq-ff-2, ...] + phaseCDF[:, nfreq-ff-1, ...])
    
    print(phaseCDF.shape)
    H5.root.sol000.phase000.val[:] = phaseCDF
    H5.flush()
    H5.close()
    return


def copyoverscalarphase(scalarh5, phasexxyyh5): 
    #note for scalarphase/phaseonly solve, does not work for tecandphase as freq axis is missing there for phase000
    H5    = tables.open_file(scalarh5, mode='r')
    H5pol = tables.open_file(phasexxyyh5,mode='a')

    phase     = H5.root.sol000.phase000.val[:] # time, freq, ant, dir
    phase_pol = H5pol.root.sol000.phase000.val[:] # time, freq, ant, dir, pol
    antenna   = H5.root.sol000.phase000.ant[:]
    print('Shapes for pol copy', phase.shape, phase_pol.shape)
    
    for ant in range(len(antenna)):
      phase_pol[:, :, ant, :,0] = phase[:, :, ant, :] # XX
      phase_pol[:, :, ant, :,-1] = phase[:, :, ant, :] # YY

    H5pol.root.sol000.phase000.val[:] = phase_pol[:,:,:]
    H5pol.flush()


    H5.close()
    H5pol.close()
    return

def copyovergain(gaininh5,gainouth5, soltype):
    H5in    = tables.open_file(gaininh5, mode='r')
    H5out   = tables.open_file(gainouth5,mode='a')
    antenna   = H5in.root.sol000.amplitude000.ant[:]
   
    h5 = h5parm.h5parm(gaininh5)
    ss = h5.getSolset('sol000')
    st = ss.getSoltab('amplitude000')    
    axesnames = st.getAxesNames()
    h5.close()
   
    if 'pol' in axesnames:

        if soltype != 'scalaramplitude' and soltype != 'amplitudeonly':
            phase     = H5in.root.sol000.phase000.val[:] 
            H5out.root.sol000.phase000.val[:] = phase
        else:    
            H5out.root.sol000.phase000.val[:] = 0.0
        
        amplitude = H5in.root.sol000.amplitude000.val[:] 
        print('Shapes for gain copy with polarizations', amplitude.shape)
        H5out.root.sol000.amplitude000.val[:] = amplitude
        
    else:
        if soltype != 'scalaramplitude' and soltype != 'amplitudeonly':
            phase     = H5in.root.sol000.phase000.val[:]
            phase_pol = H5out.root.sol000.phase000.val[:] # time, freq, ant, dir, pol
           
        amplitude = H5in.root.sol000.amplitude000.val[:] 
        amplitude_pol   = H5out.root.sol000.amplitude000.val[:] # time, freq, ant, dir, pol
        print('Shapes for gain copy 1 pol', amplitude.shape)
        
        for ant in range(len(antenna)):
            if soltype != 'scalaramplitude' and soltype != 'amplitudeonly':
                phase_pol[:, :, ant, :,0] = phase[:, :, ant, :] # XX
                phase_pol[:, :, ant, :,-1] = phase[:, :, ant, :] # YY    
            amplitude_pol[:, :, ant, :,0] = amplitude[:, :, ant, :] # XX
            amplitude_pol[:, :, ant, :,-1] = amplitude[:, :, ant, :] # YY   
        if soltype != 'scalaramplitude' and soltype != 'amplitudeonly':
            H5out.root.sol000.phase000.val[:] = phase_pol[:,:,:]
        else:
            H5out.root.sol000.phase000.val[:] = 0.0
        H5out.root.sol000.amplitude000.val[:] = amplitude_pol[:,:,:]
    
    H5out.flush()
    H5in.close()
    H5out.close()
    return    

def resetgains(parmdb):
   H5 = tables.open_file(parmdb, mode='a')
   H5.root.sol000.phase000.val[:] = 0.0
   H5.root.sol000.amplitude000.val[:] = 1.0
   H5.flush()
   H5.close()
   return

def removenans(parmdb, soltab):
   H5 = h5parm.h5parm(parmdb, readonly=False)
   vals =H5.getSolset('sol000').getSoltab(soltab).getValues()[0]
   weights = H5.getSolset('sol000').getSoltab(soltab).getValues(weight=True)[0]   
   
   idxnan  = np.where((~np.isfinite(vals))) 
   
   #print(idxnan)
   #print('Found some NaNs', vals[idxnan])
   print('Found some NaNs, flagging them....')

   if H5.getSolset('sol000').getSoltab(soltab).getType() == 'phase':
       vals[idxnan] = 0.0
   if H5.getSolset('sol000').getSoltab(soltab).getType() == 'amplitude':
       vals[idxnan] = 1.0            
   if H5.getSolset('sol000').getSoltab(soltab).getType() == 'rotation':
       vals[idxnan] = 0.0
   
   weights[idxnan] = 0.0

   H5.getSolset('sol000').getSoltab(soltab).setValues(weights,weight=True)
   H5.getSolset('sol000').getSoltab(soltab).setValues(vals)
   H5.close()
   return


def losotolofarbeam(parmdb, soltabname, ms, inverse=False, useElementResponse=True, useArrayFactor=True, useChanFreq=True):
    """
    Do the beam correction via this imported losoto operation
    """

    H5 = h5parm.h5parm(parmdb, readonly=False)
    soltab = H5.getSolset('sol000').getSoltab(soltabname)

    #t = pt.table(ms)
    sr = stationresponse(ms, inverse, useElementResponse, useArrayFactor, useChanFreq)

    numants = pt.taql('select gcount(*) as numants from '+ms+'::ANTENNA').getcol('numants')[0]
    times = soltab.getAxisValues('time')

    for vals, coord, selection in soltab.getValuesIter(returnAxes=['ant','time','pol','freq'], weight=False):
        vals = losoto.lib_operations.reorderAxes( vals, soltab.getAxesNames(), ['ant','time','freq','pol'] )

        for stationnum in range(numants):
            logger.debug('Working on station number %i' % stationnum)
            for itime, time in enumerate(times):
                beam = sr.evaluateStation(time=time, station=stationnum)
                # Reshape from [nfreq, 2, 2] to [nfreq, 4]
                beam = beam.reshape(beam.shape[0], 4)

                if soltab.getAxisLen('pol') == 2:
                    beam = beam[:,[0,3]] # get only XX and YY
                   
                if soltab.getType() == 'amplitude':
                    vals[stationnum, itime, :, :] = np.abs(beam)
                elif soltab.getType() == 'phase':
                    vals[stationnum, itime, :, :] = np.angle(beam)
                else:
                    logger.error('Beam prediction works only for amplitude/phase solution tables.')
                    return 1

        vals = losoto.lib_operations.reorderAxes( vals, ['ant','time','freq','pol'], [ax for ax in soltab.getAxesNames() if ax in ['ant','time','freq','pol']] )
        soltab.setValues(vals, selection)
    
    H5.close()
    return
 

#losotolofarbeam('P214+55_PSZ2G098.44+56.59.dysco.sub.shift.avg.weights.ms.archive_templatejones.h5', 'amplitude000', 'P214+55_PSZ2G098.44+56.59.dysco.sub.shift.avg.weights.ms.archive', inverse=False, useElementResponse=False, useArrayFactor=True, useChanFreq=True)


def cleanup(mslist):
    
    for ms in mslist:
        os.system('rm -rf ' + ms)
        
    os.system('rm -f *first-residual.fits')    
    os.system('rm -f *psf.fits') 
    os.system('rm -f *-00*-*.fits')
    os.system('rm -f *dirty.fits')
    os.system('rm -f solintimage*model.fits')
    os.system('rm -f solintimage*residual.fits')
    os.system('rm -f *pybdsf.log')
    return


def flagms_startend(ms, tecsolsfile, tecsolint):

    taql = 'taql'
       
    msout = ms + '.cut'
    
    H5 =h5parm.h5parm(tecsolsfile)
    tec = H5.getSolset('sol000').getSoltab('tec000').getValues() 
    tecvals = tec[0]
    
    axis_names = H5.getSolset('sol000').getSoltab('tec000').getAxesNames()
    time_ind = axis_names.index('time')
    ant_ind = axis_names.index('ant')
    
    #['time', 'ant', 'dir', 'freq']
    reftec = tecvals[:,0,0,0] 
    
    #print np.shape( tecvals[:,:,0,0]), np.shape( reftec[:,None]), 
    tecvals = tecvals[:,:,0,0] - reftec[:,None] # reference to zero
    
    times   = tec[1]['time']
    
    #print tecvals[:,0]
    
    goodtimesvec = []
    
    for timeid, time in enumerate(times):
    
      tecvals[timeid,:]
      
      #print timeid, np.count_nonzero( tecvals[timeid,:])
      goodtimesvec.append(np.count_nonzero( tecvals[timeid,:]))



    goodstartid = np.argmax (np.array(goodtimesvec) > 0)
    goodendid   = len(goodtimesvec) - np.argmax (np.array(goodtimesvec[::-1]) > 0)
    
    print('First good solutionslot,', goodstartid, ' out of', len(goodtimesvec))
    print('Last good solutionslot,', goodendid, ' out of', len(goodtimesvec))    
    H5.close()
    
    if (goodstartid != 0) or (goodendid != len(goodtimesvec)): # only do if needed to save some time
    
        cmd = taql + " ' select from " + ms + " where TIME in (select distinct TIME from " + ms 
        cmd+= " offset " + str(goodstartid*np.int(tecsolint)) 
        cmd+= " limit " + str((goodendid-goodstartid)*np.int(tecsolint)) +") giving " 
        cmd+= msout + " as plain'"
        
        print(cmd)
        os.system(cmd)
        
        os.system('rm -rf ' + ms)
        os.system('mv ' + msout + ' ' + ms)
    return


#flagms_startend('P215+50_PSZ2G089.52+62.34.dysco.sub.shift.avg.weights.ms.archive','phaseonlyP215+50_PSZ2G089.52+62.34.dysco.sub.shift.avg.weights.ms.archivesolsgrid_9.h5', 2)
#sys.exit()




def removestartendms(ms, starttime=None, endtime=None):

    # chdeck if output is already there and remove    
    if os.path.isdir(ms + '.cut'):
          os.system('rm -rf ' + ms + '.cut')  
    if os.path.isdir(ms + '.cuttmp'):
          os.system('rm -rf ' + ms + '.cuttmp')  

        
    cmd = 'DPPP msin=' + ms + ' ' + 'msout.storagemanager=dysco msout=' + ms + '.cut '
    cmd+=  'msin.weightcolumn=WEIGHT_SPECTRUM steps=[] msout.writefullresflag=False ' 
    if starttime is not None:
      cmd+= 'msin.starttime=' + starttime + ' '
    if endtime is not None:  
      cmd+= 'msin.endtime=' + endtime   + ' '   
    print(cmd)  
    os.system(cmd)
    
    cmd = 'DPPP msin=' + ms + ' ' + 'msout.storagemanager=dysco msout=' + ms + '.cuttmp '
    cmd+= 'msin.weightcolumn=WEIGHT_SPECTRUM_SOLVE steps=[] msout.writefullresflag=False '  
    if starttime is not None:
      cmd+= 'msin.starttime=' + starttime + ' '
    if endtime is not None:  
      cmd+= 'msin.endtime=' + endtime   + ' '
    print(cmd)
    os.system(cmd)    


    # Make a WEIGHT_SPECTRUM from WEIGHT_SPECTRUM_SOLVE
    t  = pt.table(ms + '.cut' , readonly=False)

    print('Adding WEIGHT_SPECTRUM_SOLVE')
    desc = t.getcoldesc('WEIGHT_SPECTRUM')
    desc['name']='WEIGHT_SPECTRUM_SOLVE'
    t.addcols(desc)

    t2 = pt.table(ms + '.cuttmp' , readonly=True)
    imweights = t2.getcol('WEIGHT_SPECTRUM')
    t.putcol('WEIGHT_SPECTRUM_SOLVE', imweights)

    # Fill WEIGHT_SPECTRUM with WEIGHT_SPECTRUM from second ms
    t2.close()
    t.close() 

    # clean up
    os.system('rm -rf ' + ms + '.cuttmp')



  
    return

#removestartendms('P219+50_PSZ2G084.10+58.72.dysco.sub.shift.avg.weights.ms.archive',endtime='16-Apr-2015/02:14:47.0')
#removestartendms('P223+50_PSZ2G084.10+58.72.dysco.sub.shift.avg.weights.ms.archive',starttime='24-Feb-2015/22:16:00.0')
#removestartendms('P223+52_PSZ2G088.98+55.07.dysco.sub.shift.avg.weights.ms.archive',starttime='19-Feb-2015/22:40:00.0')
#removestartendms('P223+55_PSZ2G096.14+56.24.dysco.sub.shift.avg.weights.ms.archive',starttime='31-Mar-2015/20:11:00.0')
#removestartendms('P227+53_PSZ2G088.98+55.07.dysco.sub.shift.avg.weights.ms.archive',starttime='19-Feb-2015/22:40:00.0')



def which(file_name):
    for path in os.environ["PATH"].split(os.pathsep):
        full_path = os.path.join(path, file_name)
        if os.path.exists(full_path) and os.access(full_path, os.X_OK):
            return full_path
    return None



def plotimage(fitsimagename, outplotname, mask=None, rmsnoiseimage=None):
  
  #image noise for plotting
  if rmsnoiseimage == None:
    hdulist = fits.open(fitsimagename)
  else:
    hdulist = fits.open(rmsnoiseimage)   
  imagenoise = findrms(np.ndarray.flatten(hdulist[0].data))
  hdulist.close() 
  
  #image noise info
  hdulist = fits.open(fitsimagename) 
  imagenoiseinfo = findrms(np.ndarray.flatten(hdulist[0].data))
  hdulist.close()   
  
  f = aplpy.FITSFigure(fitsimagename, slices=[0, 0])
  f.show_colorscale(vmax=16*imagenoise, vmin=-6*imagenoise, cmap='bone')
  f.set_title(fitsimagename+' (noise = {} mJy/beam)'.format(round(imagenoiseinfo*1e3, 3)))
  try: # to work around an aplpy error
    f.add_beam()
    f.beam.set_frame(True)
    f.beam.set_color('white')
    f.beam.set_edgecolor('black')
    f.beam.set_linewidth(1.)
  except:
    pass

  f.add_grid()
  f.grid.set_color('white')
  f.grid.set_alpha(0.5)
  f.grid.set_linewidth(0.2)
  f.add_colorbar()
  f.colorbar.set_axis_label_text('Flux (Jy beam$^{-1}$)')
  if mask is not None:
    try:  
      f.show_contour(mask, colors='red', levels=[0.1*imagenoise], filled=False, smooth=1, alpha=0.6, linewidths=1)
    except:
      pass  
  if os.path.isfile(outplotname + '.png'):
      os.system('rm -f ' + outplotname + '.png')
  f.save(outplotname, dpi=120, format='png')
  logger.info(fitsimagename + ' RMS noise: ' + str(imagenoiseinfo))
  return



def archive(mslist, outtarname, regionfile, fitsmask, imagename):
  path = '/disks/ftphome/pub/vanweeren'
  for ms in mslist:
    msout = ms + '.calibrated'
    if os.path.isdir(msout):
      os.system('rm -rf ' + msout)
    cmd  ='DPPP numthreads='+ str(multiprocessing.cpu_count()) +' msin=' + ms + ' msout=' + msout + ' '
    cmd +='msin.datacolumn=CORRECTED_DATA msout.storagemanager=dysco msout.writefullresflag=False steps=[]'
    os.system(cmd)
 

  msliststring = ' '.join(map(str, glob.glob('*.calibrated') ))
  cmd = 'tar -zcf ' + outtarname + ' ' + msliststring + ' selfcal.log ' +  imagename + ' '
     
  if fitsmask != None:  # add fitsmask to tar if it exists
    if os.path.isfile(fitsmask):
      cmd +=  fitsmask + ' '

  if regionfile != None:  # add box regionfile to tar if it exists
    if os.path.isfile(regionfile):
      cmd +=  regionfile + ' '
     
 
  if os.path.isfile(outtarname):
      os.system('rm -f ' + outtarname)
  logger.info('Creating archived calibrated tarball: ' + outtarname)    
  os.system(cmd)
  
  for ms in mslist:
    msout = ms + '.calibrated'   
    os.system('rm -rf ' + msout)
  return


#def reweight(mslist, pixsize, imsize, channelsout, niter, robust, multiscale=False, fitsmask=None):
   #"""
   #determine the solution time and frequency intervals based on the amount of compact source flux
   #"""
   
   #rmslist = []

   #logger.info('Adjusting weights')

   #for ms in mslist:
          #imageout =  'rmsimage' + ms.split('.ms')[0] 
          #makeimage([ms], imageout, pixsize, imsize, channelsout, np.int(niter/(len(mslist)**(1./3.))), robust, multiscale=multiscale, predict=False,fitsmask=fitsmask)
          
          #hdulist = fits.open(imageout + '-MFS-image.fits')
          #imagenoise = findrms(np.ndarray.flatten(hdulist[0].data))
          #hdulist.close() 
          #rmslist.append(imagenoise)
          
   #weightslist = []       
   #return 


def setinitial_solint(mslist, soltype_list, longbaseline, LBA,\
                      innchan_list, insolint_list, insmoothnessconstraint_list, \
                      insmoothnessreffrequency_list, inantennaconstraint_list, insoltypecycles_list):
   """
   take user input solutions,nchan,smoothnessconstraint,antennaconstraint and expand them to all ms
   these list can then be updated later with values from auto_determinesolints for example
   """

      
   if os.path.isfile('nchan.p') and os.path.isfile('solint.p'):
    
      f = open('nchan.p', 'rb') 
      nchan_list = pickle.load(f)        
      f.close()   
  
      f = open('solint.p', 'rb') 
      solint_list = pickle.load(f)        
      f.close()   

      f = open('antennaconstraint.p', 'rb') 
      antennaconstraint_list = pickle.load(f)        
      f.close()   

      f = open('smoothnessconstraint.p', 'rb') 
      smoothnessconstraint_list = pickle.load(f)        
      f.close()
      
      f = open('smoothnessreffrequency.p', 'rb') 
      smoothnessreffrequency_list = pickle.load(f)        
      f.close()

      f = open('soltypecycles.p', 'rb') 
      soltypecycles_list = pickle.load(f)        
      f.close()   
      
  
   else:
      nchan_list  = [] # list with len(soltype_list)
      solint_list = [] # list with len(soltype_list)
      smoothnessconstraint_list = [] # nested list with len(soltype_list), inner list is for ms)
      smoothnessreffrequency_list = [] # nested list with len(soltype_list), inner list is for ms)
      antennaconstraint_list = [] # nested list with len(soltype_list), inner list is for ms)
      soltypecycles_list = []  # nested list with len(soltype_list), inner list is for ms)

      for soltype_id, soltype in enumerate(soltype_list):
        nchan_ms   = [] # list with len(mslist)
        solint_ms  = [] # list with len(mslist)
        antennaconstraint_list_ms   = [] # list with len(mslist)
        smoothnessconstraint_list_ms  = [] # list with len(mslist)
        smoothnessreffrequency_list_ms  = [] # list with len(mslist)
        soltypecycles_list_ms = [] # list with len(mslist)

        for ms in mslist:
          # use try statement in case the user did not provide all the info for certain soltypes
          
          #solint 
          try:
            solint = insolint_list[soltype_id]
          except:
            solint = 1  
         
          # nchan 
          try:
            nchan = innchan_list[soltype_id]
          except:
            nchan = 10     

          # smoothnessconstraint 
          try:
            smoothnessconstraint = insmoothnessconstraint_list[soltype_id]
          except:
            smoothnessconstraint = 0.0

          # smoothnessreffrequency 
          try:
            smoothnessreffrequency = insmoothnessreffrequency_list[soltype_id]
          except:
            smoothnessreffrequency = 0.0
            
          # antennaconstraint 
          try:
            antennaconstraint = inantennaconstraint_list[soltype_id]
          except:
            antennaconstraint = None

          # soltypecycles
          soltypecycles = insoltypecycles_list[soltype_id]

          # force nchan 1 for tec(andphase) solve and in case smoothnessconstraint is invoked
          if soltype == 'tec' or  soltype == 'tecandphase' or smoothnessconstraint > 0.0:
            nchan  = 1

 
          nchan_ms.append(nchan)
          solint_ms.append(solint)
          smoothnessconstraint_list_ms.append(smoothnessconstraint)
          smoothnessreffrequency_list_ms.append(smoothnessreffrequency)
          antennaconstraint_list_ms.append(antennaconstraint)
          soltypecycles_list_ms.append(soltypecycles)
          #logging.info('MS, NCHAN_PHASE - SOLINT_PHASE ||| NCHAN_SLOW - SOLINT_SLOW: ' + str(ms) + ':: ' + str(nchan_phase) + \
          #             ' - ' + str (solint_phase) + ' ||| ' + str(nchan_ap) + ' - ' + str(solint_ap))
        
        nchan_list.append(nchan_ms)   # list of lists
        solint_list.append(solint_ms) # list of lists
        antennaconstraint_list.append(antennaconstraint_list_ms)   # list of lists
        smoothnessconstraint_list.append(smoothnessconstraint_list_ms) # list of lists
        smoothnessreffrequency_list.append(smoothnessreffrequency_list_ms) # list of lists
        
        soltypecycles_list.append(soltypecycles_list_ms)

      f = open('nchan.p', 'wb') 
      pickle.dump(nchan_list,f)        
      f.close()   
  
      f = open('solint.p', 'wb') 
      pickle.dump(solint_list,f)        
      f.close()  
      
      f = open('smoothnessconstraint.p', 'wb') 
      pickle.dump(smoothnessconstraint_list,f)        
      f.close()  

      f = open('smoothnessreffrequency.p', 'wb') 
      pickle.dump(smoothnessreffrequency_list,f)        
      f.close()  
      
      f = open('antennaconstraint.p', 'wb') 
      pickle.dump(antennaconstraint_list,f)        
      f.close()        

      f = open('soltypecycles.p', 'wb') 
      pickle.dump(soltypecycles_list,f)        
      f.close()   
      
      
   print('soltype:',soltype_list, mslist)   
   print('nchan:',nchan_list)
   print('solint:',solint_list)
   print('smoothnessconstraint:',smoothnessconstraint_list)
   print('smoothnessreffrequency:',smoothnessreffrequency_list)
   print('antennaconstraint:',antennaconstraint_list)
   print('soltypecycles:',soltypecycles_list)
   return nchan_list, solint_list, smoothnessconstraint_list, smoothnessreffrequency_list, antennaconstraint_list, soltypecycles_list

def getmsmodelinfo(ms, modelcolumn, fastrms=False):
   t = pt.table(ms + '/SPECTRAL_WINDOW')
   chanw = np.median(t.getcol('CHAN_WIDTH'))
   freq = np.median(t.getcol('CHAN_FREQ'))
   nfreq = len(t.getcol('CHAN_FREQ')[0])
   t.close()
   uvdismod = get_uvwmax(ms)*0.333 # take range [0.333uvmax - 1.0uvmax]
   
   HBA_upfreqsel = 0.75 # select only freqcencies above 75% of the available bandwidth
   freqct = 1000e6
   # the idea is that for HBA if you are far out in the beam the noise gets up much more at the higher freqs and the model flux goes down due to the spectral index. In this way we get more conservative solints
          
          
   t = pt.taql('SELECT ' + modelcolumn + ',DATA,UVW,TIME,FLAG FROM ' + ms + ' WHERE SQRT(SUMSQR(UVW[:2])) > '+ str(uvdismod) )
   model = np.abs(t.getcol(modelcolumn))
   flags = t.getcol('FLAG')
   data  = t.getcol('DATA')
   print('Compute visibility noise of the dataset with robust sigma clipping', ms)
   if fastrms:    # take only every fifth element of the array to speed up the computation
     if freq > freqct: # HBA
        noise = astropy.stats.sigma_clipping.sigma_clipped_stats(data[0:data.shape[0]:5,np.int(np.floor(np.float(nfreq)*HBA_upfreqsel)):-1,1:3],\
        mask=flags[0:data.shape[0]:5,np.int(np.floor(np.float(nfreq)*HBA_upfreqsel)):-1,1:3])[2] # use XY and YX
     else:   
        noise = astropy.stats.sigma_clipping.sigma_clipped_stats(data[0:data.shape[0]:5,:,1:3],\
        mask=flags[0:data.shape[0]:5,:,1:3])[2] # use XY and YX
   else:
     if freq > freqct: # HBA
        noise = astropy.stats.sigma_clipping.sigma_clipped_stats(data[:,np.int(np.floor(np.float(nfreq)*HBA_upfreqsel)):-1,1:3],\
        mask=flags[:,np.int(np.floor(np.float(nfreq)*HBA_upfreqsel)):-1,1:3])[2] # use XY and YX
     else:
        noise = astropy.stats.sigma_clipping.sigma_clipped_stats(data[:,:,1:3],\
        mask=flags[:,:,1:3])[2] # use XY and YX         
   
   model = np.ma.masked_array(model, flags)
   if freq > freqct: # HBA:
      flux  = np.ma.mean((model[:,np.int(np.floor(np.float(nfreq)*HBA_upfreqsel)):-1,0] + model[:,np.int(np.floor(np.float(nfreq)*HBA_upfreqsel)):-1,3])*0.5) # average XX and YY (ignore XY and YX, they are zero, or nan, in other words this is Stokes I)
   else:
      flux  = np.ma.mean((model[:,:,0] + model[:,:,3])*0.5) # average XX and YY (ignore XY and YX, they are zero, or nan)
   time  = np.unique(t.getcol('TIME'))
   tint  = np.abs(time[1]-time[0])
   print('Integration time visibilities', tint)
   t.close()

   del data, flags, model
   print('Noise visibilities:', noise, 'Jy')
   print('Flux in model', flux, 'Jy')
   print('UV-selection to compute model flux', str(uvdismod/1e3), 'km')

   
   return noise, flux, tint, chanw

def return_soltype_index(soltype_list, soltype, occurence=1, onetectypeoccurence=False):
   
   if onetectypeoccurence:
     if soltype == 'tecandphase' or soltype == 'tec':
       soltype_list = [sol.replace('tecandphase', 'tec') for sol in soltype_list]   
       soltype = 'tec'
       
   sol_index = None    
   count = 0
   for sol_id, sol in enumerate(soltype_list):
     if sol == soltype:
       sol_index = sol_id
       count = count + 1
       if occurence == count:
         return sol_index
       else:
         sol_index = None  
   return sol_index

def auto_determinesolints(mslist, soltype_list, longbaseline, LBA,\
                          innchan_list=None, insolint_list=None,\
                          uvdismod=None, modelcolumn='MODEL_DATA', redo=False,\
                          insmoothnessconstraint_list=None, insmoothnessreffrequency_list=None, \
                          inantennaconstraint_list=None, \
                          insoltypecycles_list=None, tecfactorsolint=1.0):
   """
   determine the solution time and frequency intervals based on the amount of compact source flux and noise
   """
   # 1 find the first tec/tecandphase in the soltype_list
   # set the solints there based on this code
   # update antennaconstraint if needed
   
   # find the first (scalar)complexgain
   # compute solints
   
   # based on scalarcomplexgain determine what to do
   # proceed normally, as is now
   # do an 'all' constrained solve
   # no complexgain solves at all, how to do this? update soltype_list here, or set 
   
   # soltype_list == 'tec'

   for ms_id, ms in enumerate(mslist):          
      noise, flux, tint, chanw = getmsmodelinfo(ms, modelcolumn)
      for soltype_id, soltype in enumerate(soltype_list):          
          
          ######## TEC or TECANDPHASE ######
          ######## for first occurence of tec(andphase) #######
          if soltype in ['tec', 'tecandphase'] and \
              ((soltype_id == return_soltype_index(soltype_list, 'tec', occurence=1, onetectypeoccurence=True)) or \
              (soltype_id == return_soltype_index(soltype_list, 'tecandphase', occurence=1, onetectypeoccurence=True))) :
          
             if LBA: 
               if longbaseline:
                 solint_sf = 3.0e-3*tecfactorsolint # untested
               else: #for -- LBA dutch --
                 solint_sf = 4.0e-2*tecfactorsolint #0.5e-3 # for tecandphase and coreconstraint
          
             else: # for -- HBA --
               if longbaseline:
                 solint_sf = 3.0e-3*tecfactorsolint # for tecandphase, no coreconstraint          
               else: #for -- HBA dutch --
                 solint_sf = 4.0e-2*tecfactorsolint # for tecandphase, no coreconstraint          
          
             if soltype == 'tec':
               solint_sf = solint_sf/np.sqrt(2.) # tec and coreconstraint
 
 
             # trigger antennaconstraint_phase core if solint > tint 
             if not longbaseline and (tint*solint_sf* ((noise/flux)**2) * (chanw/390.625e3) > tint):
               print(tint*solint_sf* ((noise/flux)**2) * (chanw/390.625e3))
               solint_sf = solint_sf/30. 
               print('Trigger_antennaconstraint core:', soltype, ms)
               inantennaconstraint_list[soltype_id][ms_id] = 'core'
               # do another pertubation, a slow solve of the core stations
               #if (tint*solint_sf* ((noise/flux)**2) * (chanw/390.625e3) < 360.0): # less than 6 min now, also doing constraint remote
               if (tint*solint_sf* ((noise/flux)**2) * (chanw/390.625e3) < 720.0): # less than 12 min now, also doing
                 inantennaconstraint_list[soltype_id+1][ms_id] = 'remote' # or copy over input ??
                 insoltypecycles_list[soltype_id+1][ms_id] = insoltypecycles_list[soltype_id][ms_id] # do + 1 here??
                 insolint_list[soltype_id+1][ms_id] = np.int(np.rint(10.*solint_sf* ((noise/flux)**2) * (chanw/390.625e3) ))
                 if insolint_list[soltype_id+1][ms_id] < 1:
                   insolint_list[soltype_id+1][ms_id] = 1    
               else:
                 insoltypecycles_list[soltype_id+1][ms_id] = 999                 
             
             else:  
               inantennaconstraint_list[soltype_id][ms_id] = None # or copy over input            
             
             # round to nearest integer  
             solint = np.rint(solint_sf* ((noise/flux)**2) * (chanw/390.625e3) )
             # frequency scaling is need because if we avearge in freqeuncy the solint should not change for a tec(andphase) solve
             if solint < 1:
                solint = 1        
             if (np.float(solint)*tint/3600.) > 0.5: # so check if larger than 30 min
               print('Warning, it seems there is not enough flux density on the longer baselines for solving')
               solint = np.rint(0.5*3600./tint) # max is 30 min 

             print(solint_sf*((noise/flux)**2)*(chanw/390.625e3), 'Using tec(andphase) solint:', solint)
             print('Using tec(andphase) solint [s]:', np.float(solint)*tint)
          
             insolint_list[soltype_id][ms_id] = np.int(solint)
             innchan_list[soltype_id][ms_id] = 1
             

          ######## COMPLEXGAIN or SCALARCOMPLEXGAIN or AMPLITUDEONLY or SCALARAMPLITUDE ######
          ######## requires smoothnessconstraint
          ######## for first occurence of (scalar)complexgain 
          if soltype in ['complexgain', 'scalarcomplexgain'] and (insmoothnessconstraint_list[soltype_id][ms_id] > 0.0) and \
              ((soltype_id == return_soltype_index(soltype_list, 'complexgain', occurence=1)) or \
              (soltype_id == return_soltype_index(soltype_list, 'scalarcomplexgain', occurence=1))):

             thr_disable_gain = 64. # 32. #  72.
             thr_SM15Mhz = 4.
             thr_gain_trigger_allantenna =  32. # 16. # 8.
             
             tgain_max = 4. # do not allow ap solves that are more than 4 hrs
             tgain_min = 0.3333  # check if less than 20 min, min solint is 20 min
             
             innchan_list[soltype_id][ms_id] = 1
 
             if LBA: 
               if longbaseline:
                 solint_sf = 0.4 # untested
               else: #for -- LBA dutch --
                 solint_sf = 10.0 
          
             else: # for -- HBA --
               if longbaseline:
                 solint_sf = 0.4 # for tecandphase    
               else: #for -- HBA dutch --
                 solint_sf = 0.8 #10.0 # for tecandphase         

             solint = np.rint(solint_sf*((noise/flux)**2)*(chanw/390.625e3)) 
             print(solint_sf*((noise/flux)**2)*(chanw/390.625e3), 'Computes gain solint:', solint, ' ')
             print('Computes gain solint [hr]:', np.float(solint)*tint/3600.)

             # do not allow very short ap solves
             if ((solint_sf*((noise/flux)**2)*(chanw/390.625e3))*tint/3600.) < tgain_min: #  check if less than tgain_min (20 min)
               solint = np.rint(tgain_min*3600./tint) # minimum tgain_min is 20 min 
               print('Setting gain solint to 20 min (the min value allowed):', np.float(solint)*tint/3600.)

             # do not allow ap solves that are more than tgain_max (4) hrs
             if ((solint_sf*((noise/flux)**2)*(chanw/390.625e3))*tint/3600.) > tgain_max: # so check if larger than 30 min
               print('Warning, it seems there is not enough flux density for gain solving')
               solint = np.rint(tgain_max*3600./tint) # max is tgain_max (4) hrs  

             # trigger 15 MHz smoothnessconstraint 
             #print('TEST:', ((solint_sf*((noise/flux)**2)*(chanw/390.625e3))*tint/3600.))
             if ((solint_sf*((noise/flux)**2)*(chanw/390.625e3))*tint/3600.) < thr_SM15Mhz: # so check if larger than 30 min
               insmoothnessconstraint_list[soltype_id][ms_id] = 5.0
             else:
               print('Increasing smoothnessconstraint to 15 MHz')   
               insmoothnessconstraint_list[soltype_id][ms_id] = 15.0

             # trigger nchan=0 solve because not enough S/N
             if not longbaseline and (((solint_sf*((noise/flux)**2)*(chanw/390.625e3))*tint/3600.) > thr_gain_trigger_allantenna):
               inantennaconstraint_list[soltype_id][ms_id] = 'all'
               solint = np.rint(2.0*3600./tint) # 2 hrs nchan=0 solve (do not do bandpass because slope can diverge)
               innchan_list[soltype_id][ms_id] = 0 # no frequency dependence, smoothnessconstraint will be turned of in runDPPPbase
               print('Triggering antennaconstraint all:', soltype, ms)
             else:  
               inantennaconstraint_list[soltype_id][ms_id] = None
          
             # completely disable slow solve if the solints get too long, target is too faint
             if not longbaseline and (((solint_sf*((noise/flux)**2)*(chanw/390.625e3))*tint/3600.) > thr_disable_gain):
               insoltypecycles_list[soltype_id][ms_id] = 999
               print('Disabling solve:', soltype, ms)
             else:  
               insoltypecycles_list[soltype_id][ms_id] = 3 # set to user input value? problem because not retained now
  
             insolint_list[soltype_id][ms_id] = np.int(solint)
  
             # --------------- NCHAN --------------------- not used for now

             if insmoothnessconstraint_list[soltype_id][ms_id] == 0.0 and innchan_list[soltype_id][ms_id] != 0: # DOES NOT GET HERE BECAUSE smoothnessconstraint > 0 test above

                if LBA: 
                  if longbaseline:
                    print('Not supported')
                    sys.exit()
                  else: #for -- LBA dutch, untested --
                    nchan_sf = 0.75 # for tecandphase and coreconstraint
          
                else: # for -- HBA --
                  if longbaseline:
                    nchan_sf = 0.0075 #   
                  else: #for -- HBA dutch --
                    nchan_sf = 0.75 #          
          
                nchan = np.rint(nchan_sf*(noise/flux)**2)

                # do not allow very low nchan solves
                if (np.float(nchan)*chanw/1e6) < 2.0: #  check if less than 2 MHz
                  nchan = np.rint(2.0*1e6/chanw) # 2 MHz
               
                # do not allow nchan solves that are more than 15 MHz
                if (np.float(nchan)*chanw/1e6) > 15.0: 
                  print('Warning, it seems there is not enough flux density on the longer baselines for solving')
                  nchan = np.rint(15*1e6/chanw) # 15 MHz  
 
                print(nchan_sf*(noise/flux)**2, 'Using gain nchan:', nchan)
                print('Using gain nchan [MHz]:', np.float(nchan)*chanw/1e6) 

                innchan_list[soltype_id][ms_id] = np.int(nchan)


   f = open('nchan.p', 'wb') 
   pickle.dump(innchan_list,f)        
   f.close()   
  
   f = open('solint.p', 'wb') 
   pickle.dump(insolint_list,f)        
   f.close()   

   f = open('smoothnessconstraint.p', 'wb') 
   pickle.dump(insmoothnessconstraint_list,f)        
   f.close()         

   f = open('smoothnessreffrequency.p', 'wb') 
   pickle.dump(insmoothnessreffrequency_list,f)        
   f.close()   
  
   f = open('antennaconstraint.p', 'wb') 
   pickle.dump(inantennaconstraint_list,f)        
   f.close()   

   f = open('soltypecycles.p', 'wb') 
   pickle.dump(insoltypecycles_list,f)        
   f.close()         

   print('soltype:',soltype_list, mslist)   
   print('nchan:',innchan_list)
   print('solint:',insolint_list)
   print('smoothnessconstraint:',insmoothnessconstraint_list)
   print('smoothnessreffrequency:',insmoothnessreffrequency_list)
   print('antennaconstraint:',inantennaconstraint_list)
   print('soltypecycles:',insoltypecycles_list)
      
   return innchan_list, insolint_list, insmoothnessconstraint_list, insmoothnessreffrequency_list, inantennaconstraint_list, insoltypecycles_list



def create_beamcortemplate(ms):
  """
  create a DPPP gain H5 template solutution file that can be filled with losoto
  """
  H5name = ms + '_templatejones.h5'   

  cmd = 'DPPP numthreads='+str(multiprocessing.cpu_count())+ ' msin=' + ms + ' msin.datacolumn=DATA msout=. '
  cmd += 'msin.modelcolumn=DATA '
  cmd += 'steps=[ddecal] ddecal.type=ddecal '
  cmd += 'ddecal.maxiter=1 ddecal.usemodelcolumn=True ddecal.nchan=1 '
  cmd += 'ddecal.mode=complexgain ddecal.h5parm=' + H5name  + ' '
  cmd += 'ddecal.solint=10'

  print(cmd)
  os.system(cmd)

  return H5name

def create_losoto_beamcorparset(ms, refant='CS003HBA0'):
    """
    Create a losoto parset to fill the beam correction values'.
    """
    parset = 'losotobeam.parset'
    os.system('rm -f ' + parset)
    f=open(parset, 'w')

    f.write('pol = [XX,YY]\n')
    f.write('soltab = [sol000/*]\n\n\n')

    f.write('[plotphase]\n')
    f.write('operation = PLOT\n')
    f.write('soltab = [sol000/phase000]\n')
    f.write('axesInPlot = [time,freq]\n')
    f.write('axisInTable = ant\n')
    f.write('minmax = [-0.5,0.5]\n')
    f.write('prefix = plotlosoto%s/phases_beam\n' % ms)
    f.write('refAnt = %s\n\n\n' % refant)

    f.write('[plotamp]\n')
    f.write('operation = PLOT\n')
    f.write('soltab = [sol000/amplitude000]\n')
    f.write('axesInPlot = [time,freq]\n')
    f.write('axisInTable = ant\n')
    f.write('minmax = [0.2,1]\n')
    f.write('prefix = plotlosoto%s/amplitudes_beam\n' %ms)

    f.close()
    return parset

def create_losoto_tecandphaseparset(ms, refant='CS003HBA0', outplotname='fasttecandphase'):
    parset = 'losoto_plotfasttecandphase.parset'
    os.system('rm -f ' + parset)
    f=open(parset, 'w')
  
    f.write('pol = []\n')
    f.write('Ncpu = 0\n\n\n')

    f.write('[plottecandphase]\n')
    f.write('operation = PLOT\n')
    f.write('soltab = [sol000/phase000]\n')
    f.write('axesInPlot = [time]\n')
    f.write('axisInTable = ant\n')
    f.write('minmax = [-3.14,3.14]\n')
    f.write('soltabToAdd = tec000\n')
    f.write('figSize=[120,20]\n')
    f.write('prefix = plotlosoto%s/fasttecandphase\n' % ms)
    f.write('refAnt = %s\n' % refant)
  
    f.close()
    return parset

def create_losoto_tecparset(ms, refant='CS003HBA0', outplotname='fasttec'):
    parset = 'losoto_plotfasttec.parset'
    os.system('rm -f ' + parset)
    f=open(parset, 'w')
  
    f.write('pol = []\n')
    f.write('Ncpu = 0\n\n\n')

    f.write('[plottec]\n')
    f.write('operation = PLOT\n')
    f.write('soltab = [sol000/tec000]\n')
    f.write('axesInPlot = [time]\n')
    f.write('axisInTable = ant\n')
    f.write('minmax = [-0.2,0.2]\n')
    f.write('figSize=[120,20]\n')
    f.write('prefix = plotlosoto%s/%s\n' % (ms,outplotname))
    f.write('refAnt = %s\n' % refant)
  
    f.close()
    return parset



def create_losoto_rotationparset(ms, refant='CS003HBA0', onechannel=False, outplotname='rotatation'):
    parset = 'losoto_plotrotation.parset'
    os.system('rm -f ' + parset)
    f=open(parset, 'w')

    f.write('pol = [XX,YY]\n')
    f.write('soltab = [sol000/*]\n')
    f.write('Ncpu = 0\n\n\n')

    f.write('[plotrotation]\n')
    f.write('operation = PLOT\n')
    f.write('soltab = [sol000/rotation000]\n')
    if onechannel:
      f.write('axesInPlot = [time]\n')      
    else:
      f.write('axesInPlot = [time,freq]\n')
    f.write('axisInTable = ant\n')
    f.write('minmax = [-3.14,3.14]\n')
    f.write('figSize=[120,20]\n')
    f.write('prefix = plotlosoto%s/%s\n' % (ms,outplotname))
    f.write('refAnt = %s\n' % refant)
    f.close()
    return parset


def create_losoto_fastphaseparset(ms, refant='CS003HBA0', onechannel=False, onepol=False, outplotname='fastphase'):
    parset = 'losoto_plotfastphase.parset'
    os.system('rm -f ' + parset)
    f=open(parset, 'w')

    f.write('pol = [XX,YY]\n')
    f.write('soltab = [sol000/*]\n')
    f.write('Ncpu = 0\n\n\n')

    f.write('[plotphase]\n')
    f.write('operation = PLOT\n')
    f.write('soltab = [sol000/phase000]\n')
    if onechannel:
      f.write('axesInPlot = [time]\n')
      if not onepol:
        f.write('axisInCol = pol\n')
      
    else:
      f.write('axesInPlot = [time,freq]\n')
    f.write('axisInTable = ant\n')
    f.write('minmax = [-3.14,3.14]\n')
    f.write('figSize=[120,20]\n')
    f.write('prefix = plotlosoto%s/%s\n' % (ms,outplotname))
    f.write('refAnt = %s\n' % refant)

    if not onepol:
      f.write('[plotphasediff]\n')
      f.write('operation = PLOT\n')
      f.write('soltab = [sol000/phase000]\n')
      if onechannel:
        f.write('axesInPlot = [time]\n')
      else:
        f.write('axesInPlot = [time,freq]\n')
      f.write('axisInTable = ant\n')
      f.write('minmax = [-3.14,3.14]\n')
      f.write('figSize=[120,20]\n')
      f.write('prefix = plotlosoto%s/%spoldiff\n' % (ms,outplotname))
      f.write('refAnt = %s\n' % refant)  
      f.write('axisDiff=pol\n')
 
        
    f.close()
    return parset


def create_losoto_flag_apgridparset(ms, flagging=True, maxrms=7.0, maxrmsphase=7.0, includesphase=True, \
                                    refant='CS003HBA0', onechannel=False, medamp=2.5, flagphases=True, \
                                    onepol=False, outplotname='slowamp'):

    parset= 'losoto_flag_apgrid.parset'
    os.system('rm -f ' + parset)
    f=open(parset, 'w')

    #f.write('pol = []\n')
    f.write('soltab = [sol000/*]\n')
    f.write('Ncpu = 0\n\n\n')
   
    f.write('[plotamp]\n')
    f.write('operation = PLOT\n')
    f.write('soltab = [sol000/amplitude000]\n')
    if onechannel:
      f.write('axesInPlot = [time]\n')
      if not onepol:
        f.write('axisInCol = pol\n')
    else:
      f.write('axesInPlot = [time,freq]\n')   
    f.write('axisInTable = ant\n')
    #if longbaseline:
    #  f.write('minmax = [0,2.5]\n')        
    #else:    
    f.write('minmax = [%s,%s]\n' % (str(medamp/4.0), str(medamp*2.5)))
    #f.write('minmax = [0,2.5]\n')
    f.write('prefix = plotlosoto%s/%samp\n\n\n' % (ms,outplotname))

    if includesphase:   
        f.write('[plotphase]\n')
        f.write('operation = PLOT\n')
        f.write('soltab = [sol000/phase000]\n')
        if onechannel:
          f.write('axesInPlot = [time]\n')
          if not onepol:
            f.write('axisInCol = pol\n')
        else:
           f.write('axesInPlot = [time,freq]\n')    
        f.write('axisInTable = ant\n')
        f.write('minmax = [-3.14,3.14]\n')
        f.write('prefix = plotlosoto%s/%sphase\n' % (ms,outplotname))
        f.write('refAnt = %s\n\n\n' % refant)

    if flagging:
        f.write('[flagamp]\n')
        f.write('soltab = [sol000/amplitude000]\n')
        f.write('operation = FLAG\n')
        if onechannel:
          f.write('axesToFlag = [time]\n')
        else:
          f.write('axesToFlag = [time,freq]\n')
        f.write('mode = smooth\n')
        f.write('maxCycles = 3\n')
        f.write('windowNoise = 7\n')
        f.write('maxRms = %s\n' % str(maxrms))
        if onechannel:
          f.write('order  = [5]\n\n\n')
        else:
          f.write('order  = [5,5]\n\n\n')  
    
        if includesphase and flagphases:
            f.write('[flagphase]\n')
            f.write('soltab = [sol000/phase000]\n')
            f.write('operation = FLAG\n')
            if onechannel:
              f.write('axesToFlag = [time]\n')
            else:
              f.write('axesToFlag = [time,freq]\n')
            f.write('mode = smooth\n')
            f.write('maxCycles = 3\n')
            f.write('windowNoise = 7\n')
            f.write('maxRms = %s\n' % str(maxrmsphase))
            if onechannel:
              f.write('order  = [5]\n\n\n')
            else:
              f.write('order  = [5,5]\n\n\n')  

        f.write('[plotampafter]\n')
        f.write('operation = PLOT\n')
        f.write('soltab = [sol000/amplitude000]\n')
        if onechannel:
          f.write('axesInPlot = [time]\n')
          if not onepol:
            f.write('axisInCol = pol\n')
        else:
          f.write('axesInPlot = [time,freq]\n')   
        f.write('axisInTable = ant\n')
        #f.write('minmax = [0,2.5]\n')
        f.write('minmax = [%s,%s]\n' % (str(medamp/4.0), str(medamp*2.5)))
        f.write('prefix = plotlosoto%s/%sampfl\n\n\n' % (ms,outplotname))

        if includesphase and flagphases:
            f.write('[plotphase_after]\n')
            f.write('operation = PLOT\n')
            f.write('soltab = [sol000/phase000]\n')
            if onechannel:
              f.write('axesInPlot = [time]\n')
              if not onepol:
                f.write('axisInCol = pol\n')
            else:
              f.write('axesInPlot = [time,freq]\n')   
            f.write('axisInTable = ant\n')
            f.write('minmax = [-3.14,3.14]\n')
            f.write('prefix = plotlosoto%s/%sphasefl\n' % (ms,outplotname))
            f.write('refAnt = %s\n' % refant)
  
  
    f.close()
    return parset

def create_losoto_mediumsmoothparset(ms, boxsize, longbaseline, includesphase=True, refant='CS003HBA0',\
                                     onechannel=False, outplotname='runningmedian'):
    parset= 'losoto_mediansmooth.parset'
    os.system('rm -f ' + parset)
    f=open(parset, 'w')

    f.write('pol = []\n')
    f.write('soltab = [sol000/*]\n')
    f.write('Ncpu = 0\n\n\n')

    if includesphase:
        f.write('[smoothphase]\n')
        f.write('soltab = [sol000/phase000]\n')
        f.write('operation= SMOOTH\n')
        if onechannel:
          f.write('axesToSmooth = [time]\n')
          f.write('size = [%s]\n' % (boxsize, boxsize))
        else:
          f.write('axesToSmooth = [freq,time]\n')  
          f.write('size = [%s,%s]\n' % (boxsize, boxsize))
        f.write('mode = runningmedian\n\n\n')

    f.write('[smoothamp]\n')
    f.write('soltab = [sol000/amplitude000]\n')
    f.write('operation= SMOOTH\n')
    if onechannel:
      f.write('axesToSmooth = [time]\n')
      f.write('size = [%s]\n' % (boxsize, boxsize))
    else:
      f.write('axesToSmooth = [freq,time]\n')  
      f.write('size = [%s,%s]\n' % (boxsize, boxsize))
    f.write('mode = runningmedian\n\n\n')

    f.write('[plotamp_after]\n')
    f.write('operation = PLOT\n')
    f.write('soltab = [sol000/amplitude000]\n')
    if onechannel:
      f.write('axesInPlot = [time]\n')
    else:
      f.write('axesInPlot = [time,freq]\n')   
    f.write('axisInTable = ant\n')
    if longbaseline:
      f.write('minmax = [0,2.5]\n')        
    else:    
      f.write('minmax = [0,2.5]\n')
    f.write('prefix = plotlosoto%s/amps_smoothed\n\n\n' % ms)

    if includesphase:
        f.write('[plotphase_after]\n')
        f.write('operation = PLOT\n')
        f.write('soltab = [sol000/phase000]\n')
        if onechannel:
          f.write('axesInPlot = [time]\n')
        else:
          f.write('axesInPlot = [time,freq]\n')
        f.write('axisInTable = ant\n')
        f.write('minmax = [-3.14,3.14]\n')
        f.write('prefix = plotlosoto%s/phases_smoothed\n\n\n' % ms)
        f.write('refAnt = %s\n' % refant)


        f.write('[plotphase_after1rad]\n')
        f.write('operation = PLOT\n')
        f.write('soltab = [sol000/phase000]\n')
        if onechannel:
          f.write('axesInPlot = [time]\n')
        else:
          f.write('axesInPlot = [time,freq]\n')
        f.write('axisInTable = ant\n')
        f.write('minmax = [-1,1]\n')
        f.write('prefix = plotlosoto%s/phases_smoothed1rad\n' % ms)
        f.write('refAnt = %s\n' % refant)

    f.close()
    return parset

def check_phaseup(H5name):
    H5 =  tables.open_file(H5name, mode='r')
    try:
      ants   = H5.root.sol000.phase000.ant[:]
    except:
      pass

    try:
      ants   = H5.root.sol000.amplitude000.ant[:]
    except:
      pass
  
    try:
      ants   = H5.root.sol000.rotation000.ant[:]
    except:
      pass
  
    try:
      ants   = H5.root.sol000.tec000.ant[:]
    except:
      pass
  
    #H5 = h5parm.h5parm(H5name, readonly=False)
    #ants = H5.getSolset('sol000').getAnt().keys()
    H5.close()
    if 'ST001' in ants:
        return True
    else:
        return False

def fixbeam_ST001(H5name):
    
   H5 = h5parm.h5parm(H5name, readonly=False)
   
   ants = H5.getSolset('sol000').getAnt().keys()
   antsrs = fnmatch.filter(ants,'RS*')
   ST001 = False
   
   if 'ST001' in ants:
     ST001 = True
     amps    = H5.getSolset('sol000').getSoltab('amplitude000').getValues()
     ampvals = H5.getSolset('sol000').getSoltab('amplitude000').getValues()[0]
     phasevals = H5.getSolset('sol000').getSoltab('phase000').getValues()[0]
     
     idx = np.where(amps[1]['ant'] == 'ST001')[0][0]
     idxrs = np.where(amps[1]['ant'] == antsrs[0])[0][0]
     #idx106 = np.where(amps[1]['ant'] == 'RS106HBA')[0][0]
     #idx305 = np.where(amps[1]['ant'] == 'RS305HBA')[0][0]
     #idx508 = np.where(amps[1]['ant'] == 'RS508HBA')[0][0]
     #idx406 = np.where(amps[1]['ant'] == 'RS406HBA')[0][0]
     #idxnonecheck = np.where(amps[1]['ant'] == 'blaHBA')[0][0]

     #ampvals[:,:, idx, 0,:] = 1.0  # set amplitude to 1.
     #phasevals[:,:, idx, 0,:] = 0.0 # set phase to 0.

     ampvals[:,:, idx, 0,:] = ampvals[:,:, idxrs, 0,:]
     phasevals[:,:, idx, 0,:] = 0.0 
     
     H5.getSolset('sol000').getSoltab('amplitude000').setValues(ampvals)
     H5.getSolset('sol000').getSoltab('phase000').setValues(phasevals)
   
   H5.close()
    
   return ST001

def circular(ms, linear=False):
    """
    convert to circular correlations
    """
    taql = 'taql'
    scriptn = 'python lin2circ.py'
    if linear:
      cmdlin2circ = scriptn + ' -i ' + ms + ' --column=DATA --lincol=CORRECTED_DATA --back' 
    else:
      cmdlin2circ = scriptn + ' -i ' + ms + ' --column=DATA --outcol=CORRECTED_DATA' 
    print(cmdlin2circ)
    os.system(cmdlin2circ)
    os.system(taql + " 'update " + ms + " set DATA=CORRECTED_DATA'")
    return

def beamcor(ms, usedppp=True):
    """
    correct a ms for the beam in the phase center (array_factor only)
    """
    losoto = 'losoto'
    taql = 'taql'
    H5name = create_beamcortemplate(ms)
    parset = create_losoto_beamcorparset(ms)

    losotolofarbeam(H5name, 'phase000', ms, useElementResponse=False, useArrayFactor=True, useChanFreq=True)
    losotolofarbeam(H5name, 'amplitude000', ms, useElementResponse=False, useArrayFactor=True, useChanFreq=True)   

    phasedup = fixbeam_ST001(H5name)

    if usedppp and not phasedup :
        cmddppp = 'DPPP numthreads='+str(multiprocessing.cpu_count())+ ' msin=' + ms + ' msin.datacolumn=DATA msout=. '
        cmddppp += 'msin.weightcolumn=WEIGHT_SPECTRUM '
        cmddppp += 'msout.datacolumn=CORRECTED_DATA steps=[beam] msout.storagemanager=dysco '
        cmddppp += 'beam.type=applybeam beam.updateweights=True ' 
        # weights True ???
        cmddppp += 'beam.direction=[] ' # correction for the current phase center
        #cmddppp += 'beam.beammode= ' default is full, will undo element as well(!)
        print('DPPP applybeam:', cmddppp)
        os.system(cmddppp)
        os.system(taql + " 'update " + ms + " set DATA=CORRECTED_DATA'")
    else:
        #print('Phase up dataset, cannot use DPPP beam, do manual correction')
        cmdlosoto = losoto + ' ' + H5name + ' ' + parset
        print(cmdlosoto)
        os.system(cmdlosoto)
    
        cmd = 'DPPP numthreads='+str(multiprocessing.cpu_count())+ ' msin=' + ms + ' msin.datacolumn=DATA msout=. '
        cmd += 'msin.weightcolumn=WEIGHT_SPECTRUM '
        cmd += 'msout.datacolumn=CORRECTED_DATA steps=[ac1,ac2] msout.storagemanager=dysco '
        cmd += 'ac1.parmdb='+H5name + ' ac2.parmdb='+H5name + ' '
        cmd += 'ac1.type=applycal ac2.type=applycal '
        cmd += 'ac1.correction=phase000 ac2.correction=amplitude000 ac2.updateweights=True ' 
        print('DPPP applycal:', cmd)
        os.system(cmd)
        os.system(taql + " 'update " + ms + " set DATA=CORRECTED_DATA'")
 
        # Add beam correction keyword here.
        # This code only applies the array factor and assumes the element beam was corrected already.
        # Valid values are Element, ArrayFactor or Full.
        try:
           t = pt.table(ms, readonly=False)
           t.putcolkeywords('DATA', {'LOFAR_APPLIED_BEAM_MODE': 'Full'})
           t2 = pt.table(ms + '::FIELD')
           phasedir = t2.getcol('PHASE_DIR').squeeze()
           t2.close()
           beamdir = t.getcolkeyword('DATA', 'LOFAR_APPLIED_BEAM_DIR')
           # Right ascension in radians is set in m0
           # Declination in radians is set in m1
           beamdir['m0']['value'] = phasedir[0]
           beamdir['m1']['value'] = phasedir[1]
           t.putcolkeywords('DATA', {'LOFAR_APPLIED_BEAM_DIR': beamdir})
           t.close()
        except:
           print('Warning could not update LOFAR BEAM keywords in ms, it seems this data was preprocessed with a very old DPPP version')  
    
    return

def beamcormodel(ms):
    """
    create MODEL_DATA_BEAMCOR where we store beam corrupted model data
    """   
    H5name = ms + '_templatejones.h5'   
    
    cmd = 'DPPP numthreads='+str(multiprocessing.cpu_count())+' msin=' + ms + ' msin.datacolumn=MODEL_DATA msout=. '
    cmd += 'msout.datacolumn=MODEL_DATA_BEAMCOR steps=[ac1,ac2] msout.storagemanager=dysco '
    cmd += 'ac1.parmdb='+H5name + ' ac2.parmdb='+H5name + ' '
    cmd += 'ac1.type=applycal ac2.type=applycal '
    cmd += 'ac1.correction=phase000 ac2.correction=amplitude000 ac2.updateweights=False '
    cmd += 'ac1.invert=False ac2.invert=False ' # Here we corrupt with the beam !
    print('DPPP applycal:', cmd)
    os.system(cmd)
   
    return

def findrms(mIn,maskSup=1e-7):
    """
    find the rms of an array, from Cycil Tasse/kMS
    """
    m=mIn[np.abs(mIn)>maskSup]
    rmsold=np.std(m)
    diff=1e-1
    cut=3.
    bins=np.arange(np.min(m),np.max(m),(np.max(m)-np.min(m))/30.)
    med=np.median(m)
    for i in range(10):
        ind=np.where(np.abs(m-med)<rmsold*cut)[0]
        rms=np.std(m[ind])
        if np.abs((rms-rmsold)/rmsold)<diff: break
        rmsold=rms
    return rms


def findamplitudenoise(parmdb):
      """
      find the 'amplitude noise' in a parmdb, return non-clipped rms value
      """
      H5 = h5parm.h5parm(parmdb, readonly=True) 
      amps =H5.getSolset('sol000').getSoltab('amplitude000').getValues()[0]
      weights = H5.getSolset('sol000').getSoltab('amplitude000').getValues(weight=True)[0]
      H5.close()

      idx = np.where(weights != 0.0)
      
      amps = amps[idx]
      amps = amps[np.isfinite(amps)]
      amps = np.log10(np.ndarray.flatten(amps))
      
      
      noise = findrms(amps)
      
      logger.info('Noise and clipped noise' + str(parmdb) + ' ' + str(np.std(amps)) + ' ' + str(noise))

      # do not return clipped noise, we are intersted in finding data with high outliers
      return np.std(amps)


def getimsize(boxfile, cellsize=1.5):
   """
   find imsize need to image a DS9 boxfile region
   """
   r = pyregion.open(boxfile)
   
   xs = np.ceil((r[0].coord_list[2])*1.2*3600./cellsize)
   ys = np.ceil((r[0].coord_list[3])*1.2*3600./cellsize)

   imsize = np.ceil(xs) # // Round up decimals to an integer
   if(imsize % 2 == 1): 
       imsize = imsize + 1
   
   #if np.int(imsize) < 512:
   #    imsize = 512
   return np.int(imsize)


def smoothsols(parmdb, ms, longbaseline, includesphase=True):
    
    losoto = 'losoto'    
    
    cmdlosoto = losoto + ' ' + parmdb + ' '
    noise = findamplitudenoise(parmdb)
    smooth = False
    if noise >= 0.1:
      cmdlosoto += create_losoto_mediumsmoothparset(ms, '9', longbaseline, includesphase=includesphase)
      smooth = True      
    if noise < 0.1 and noise >= 0.08:
      cmdlosoto += create_losoto_mediumsmoothparset(ms, '7', longbaseline, includesphase=includesphase)
      smooth = True    
    if noise < 0.08 and noise >= 0.07:
      cmdlosoto += create_losoto_mediumsmoothparset(ms, '5', longbaseline, includesphase=includesphase)
      smooth = True
    if noise < 0.07 and noise >= 0.04:
      cmdlosoto += create_losoto_mediumsmoothparset(ms, '3', longbaseline, includesphase=includesphase)
      smooth = True
    print(cmdlosoto)
    if smooth:
       os.system(cmdlosoto)
    return



def change_refant(parmdb, soltab):
    '''
    Changes the reference antenna, if needed, for phase
    '''
    H5     = h5parm.h5parm(parmdb, readonly=False) 
    phases = H5.getSolset('sol000').getSoltab(soltab).getValues()[0]
    weights= H5.getSolset('sol000').getSoltab(soltab).getValues(weight=True)[0]
    axesnames = H5.getSolset('sol000').getSoltab(soltab).getAxesNames() 
    print('axesname', axesnames)
    #print 'SHAPE', np.shape(weights)#, np.size(weights[:,:,0,:,:])

    
    antennas = list(H5.getSolset('sol000').getSoltab(soltab).getValues()[1]['ant'])
    #print antennas
    
    if 'pol' in axesnames:
      idx0    = np.where((weights[:,:,0,:,:] == 0.0))[0]
      idxnan  = np.where((~np.isfinite(phases[:,:,0,:,:])))[0]
    
      refant = ' '
      tmpvar = np.float(np.size(weights[:,:,0,:,:]))
    else:
      idx0    = np.where((weights[:,0,:,:] == 0.0))[0]
      idxnan  = np.where((~np.isfinite(phases[:,0,:,:])))[0]
    
      refant = ' '
      tmpvar = np.float(np.size(weights[:,0,:,:]))   
    
    if ((np.float(len(idx0))/tmpvar) > 0.5) or ((np.float(len(idxnan))/tmpvar) > 0.5):
      logger.info('Trying to changing reference anntena')
    

      for antennaid,antenna in enumerate(antennas[1::]):
            print(antenna)
            if 'pol' in axesnames:
              idx0    = np.where((weights[:,:,antennaid+1,:,:] == 0.0))[0]
              idxnan  = np.where((~np.isfinite(phases[:,:,antennaid+1,:,:])))[0]
              tmpvar = np.float(np.size(weights[:,:,antennaid+1,:,:]))
            else:
              idx0    = np.where((weights[:,antennaid+1,:,:] == 0.0))[0]
              idxnan  = np.where((~np.isfinite(phases[:,antennaid+1,:,:])))[0]
              tmpvar = np.float(np.size(weights[:,antennaid+1,:,:]))
            
            print(idx0, idxnan, ((np.float(len(idx0))/tmpvar)))
            if  ((np.float(len(idx0))/tmpvar) < 0.5) and ((np.float(len(idxnan))/tmpvar) < 0.5):
              logger.info('Found new reference anntena,' + str(antenna))
              refant = antenna
              break
    
    
    if refant != ' ':
        for antennaid,antenna in enumerate(antennas):
            if 'pol' in axesnames:
              phases[:,:,antennaid,:,:] = phases[:,:,antennaid,:,:] - phases[:,:,antennas.index(refant),:,:]
            else:
              #phases[:,antennaid,:,:] = phases[:,antennaid,:,:] - phases[:,antennas.index(refant),:,:]
              phases[:,:,antennaid,:] = phases[:,:,antennaid,:] - phases[:,:,antennas.index(refant),:]   
        H5.getSolset('sol000').getSoltab(soltab).setValues(phases)     

    H5.close()
    return


def calculate_solintnchan(compactflux):
    
    if compactflux >= 3.5:
        nchan = 5.
        solint_phase = 1.
        
    if compactflux <= 3.5:
        nchan = 5.
        solint_phase = 1.
  
    if compactflux <= 1.0:
        nchan= 10.
        solint_phase = 2
 
    if compactflux <= 0.75:
        nchan= 15.
        solint_phase = 3.

 
    #solint_ap = 100. / np.sqrt(compactflux)
    solint_ap = 120. /(compactflux**(1./3.)) # do third power-scaling
    #print solint_ap
    if solint_ap < 60.:
        solint_ap = 60.  # shortest solint_ap allowed
    if solint_ap > 180.:
        solint_ap = 180.  # longest solint_ap allowed
 
    if compactflux <= 0.4:
        nchan= 15.
        solint_ap = 180.
 
    return np.int(nchan), np.int(solint_phase), np.int(solint_ap)




def determine_compactsource_flux(fitsimage):
    
    hdul = fits.open(fitsimage)
    bmaj = hdul[0].header['BMAJ']
    bmin = hdul[0].header['BMIN']
    avgbeam = 3600.*0.5*(bmaj + bmin)
    pixsize = 3600.*(hdul[0].header['CDELT2'])
    rmsbox1 = np.int(7.*avgbeam/pixsize)
    rmsbox2 = np.int((rmsbox1/10.) + 1.)
    
    img = bdsf.process_image(fitsimage,mean_map='zero', rms_map=True, rms_box = (rmsbox1,rmsbox2))
    total_flux_gaus = np.copy(img.total_flux_gaus)
    hdul.close()
    # trying to reset.....
    del img
    
    return total_flux_gaus


def getdeclinationms(ms):
    '''
    return approximate declination of pointing center of the ms
    input: a ms
    output: declination in degrees
    '''
    t = pt.table(ms +'/FIELD', readonly=True)
    direction = np.squeeze ( t.getcol('PHASE_DIR') )
    t.close()
    return 360.*direction[1]/(2.*np.pi)

#print getdeclinationms('1E216.dysco.sub.shift.avg.weights.set0.ms')
#sys.exit()

def declination_sensivity_factor(declination):
    '''
    compute sensitivy factor lofar data, reduced by delclination, eq. from G. Heald.
    input declination is units of degrees
    '''
    factor = 1./(np.cos(2.*np.pi*(declination - 52.9)/360.)**2)

    return factor

#print declination_sensivity_factor(-3.7)
#sys.exit()

def flaglowamps(parmdb, lowampval=0.1, flagging=True, setweightsphases=True):
    '''
    flag bad amplitudes in H5 parmdb, those with values < lowampval
    '''
    H5 = h5parm.h5parm(parmdb, readonly=False) 
    amps =H5.getSolset('sol000').getSoltab('amplitude000').getValues()[0]
    idx = np.where(amps < lowampval)
    weights = H5.getSolset('sol000').getSoltab('amplitude000').getValues(weight=True)[0]
    
    if flagging: # no flagging
      weights[idx] = 0.0
      print('Settting some weights to zero in flaglowamps')
    amps[idx] = 1.0
    H5.getSolset('sol000').getSoltab('amplitude000').setValues(weights,weight=True)
    H5.getSolset('sol000').getSoltab('amplitude000').setValues(amps)

    #also put phases weights and phases to zero
    if setweightsphases:
        phases =H5.getSolset('sol000').getSoltab('phase000').getValues()[0]
        weights_p = H5.getSolset('sol000').getSoltab('phase000').getValues(weight=True)[0]
        if flagging: # no flagging
            weights_p[idx] = 0.0
            phases[idx] = 0.0
            H5.getSolset('sol000').getSoltab('phase000').setValues(weights_p,weight=True)
            H5.getSolset('sol000').getSoltab('phase000').setValues(phases)
    
    H5.close()
    return



def flagbadamps(parmdb, setweightsphases=True):
    '''
    flag bad amplitudes in H5 parmdb, those with amplitude==1.0
    '''
    H5 = h5parm.h5parm(parmdb, readonly=False) 
    amps =H5.getSolset('sol000').getSoltab('amplitude000').getValues()[0]
    idx = np.where(amps == 1.0)
    weights = H5.getSolset('sol000').getSoltab('amplitude000').getValues(weight=True)[0]
    
    weights[idx] = 0.0
    H5.getSolset('sol000').getSoltab('amplitude000').setValues(weights,weight=True)

    #also put phases weights and phases to zero
    if setweightsphases:
        phases =H5.getSolset('sol000').getSoltab('phase000').getValues()[0]
        weights_p = H5.getSolset('sol000').getSoltab('phase000').getValues(weight=True)[0]
        weights_p[idx] = 0.0
        phases[idx] = 0.0

        H5.getSolset('sol000').getSoltab('phase000').setValues(weights_p,weight=True)
        H5.getSolset('sol000').getSoltab('phase000').setValues(phases)
    
    H5.close()
    return


def medianamp(parmdb):
    H5 = h5parm.h5parm(parmdb, readonly=True) 
    amps =H5.getSolset('sol000').getSoltab('amplitude000').getValues()[0]
    weights = H5.getSolset('sol000').getSoltab('amplitude000').getValues(weight=True)[0]
    idx = np.where(weights != 0.0)
    medamps = 10**(np.nanmedian(np.log10(amps[idx])))
    H5.close()
    print('Median amplitude of ', parmdb, ':', medamps)
    return medamps

def normamplitudes(parmdb):
    '''
    normalize amplitude solutions to one
    '''
    
    if len(parmdb) == 1:
      H5 = h5parm.h5parm(parmdb[0], readonly=False) 
      amps =H5.getSolset('sol000').getSoltab('amplitude000').getValues()[0]
      weights = H5.getSolset('sol000').getSoltab('amplitude000').getValues(weight=True)[0]
      idx = np.where(weights != 0.0)
    
      amps = np.log10(amps)
      logger.info('Mean amplitudes before normalization: ' + str(10**(np.nanmean(amps[idx]))))
      amps = amps - (np.nanmean(amps[idx]))
      logger.info('Mean amplitudes after normalization: ' + str(10**(np.nanmean(amps[idx]))))
      amps = 10**(amps)

      H5.getSolset('sol000').getSoltab('amplitude000').setValues(amps) 
      H5.close()

    else:
      #amps = []  
      for i, parmdbi in enumerate(parmdb):
          H5 = h5parm.h5parm(parmdbi, readonly=True) 
          ampsi = np.copy(H5.getSolset('sol000').getSoltab('amplitude000').getValues()[0])
          weights = H5.getSolset('sol000').getSoltab('amplitude000').getValues(weight=True)[0]
          idx = np.where(weights != 0.0)
          logger.info(parmdbi + '  Normfactor: '+ str(10**(np.nanmean(np.log10(ampsi[idx])))))
          if i == 0:
            amps = np.ndarray.flatten(ampsi[idx])
          else:
            amps = np.concatenate((amps, np.ndarray.flatten(ampsi[idx])),axis=0)

          #print np.shape(amps), parmdbi
          H5.close()
      normmin = (np.nanmean(np.log10(amps))) 
      logger.info('Global normfactor: ' + str(10**normmin))
      # now write the new H5 files
      for parmdbi in parmdb:  
         H5   = h5parm.h5parm(parmdbi, readonly=False) 
         ampsi = np.copy(H5.getSolset('sol000').getSoltab('amplitude000').getValues()[0])
         ampsi = (np.log10(ampsi)) - normmin
         ampsi = 10**ampsi
         H5.getSolset('sol000').getSoltab('amplitude000').setValues(ampsi) 
         H5.close()
    return






def flaghighgamps(parmdb, highampval=10.,flagging=True, setweightsphases=True):
    '''
    flag bad amplitudes in H5 parmdb, those with values > highampval
    '''
    H5 = h5parm.h5parm(parmdb, readonly=False) 
    amps =H5.getSolset('sol000').getSoltab('amplitude000').getValues()[0]
    idx = np.where(amps > highampval)
    weights = H5.getSolset('sol000').getSoltab('amplitude000').getValues(weight=True)[0]
    
    
    if flagging: 
      weights[idx] = 0.0
      print('Settting some weights to zero in flaghighgamps')
    amps[idx] = 1.0
    H5.getSolset('sol000').getSoltab('amplitude000').setValues(weights,weight=True)
    H5.getSolset('sol000').getSoltab('amplitude000').setValues(amps)

    #also put phases weights and phases to zero
    if setweightsphases:
        phases =H5.getSolset('sol000').getSoltab('phase000').getValues()[0]
        weights_p = H5.getSolset('sol000').getSoltab('phase000').getValues(weight=True)[0]
        if flagging: 
            weights_p[idx] = 0.0
            phases[idx] = 0.0
            #print(idx)
            H5.getSolset('sol000').getSoltab('phase000').setValues(weights_p,weight=True)
            H5.getSolset('sol000').getSoltab('phase000').setValues(phases)
    
    #H5.getSolset('sol000').getSoltab('phase000').flush()
    #H5.getSolset('sol000').getSoltab('amplitude000').flush()
    H5.close()
    return


def removenegativefrommodel(imagenames):
    '''
    replace negative pixel values in WSCLEAN model images with zeros
    '''
    perseus = False
    A1795   = False
    A1795imlist = sorted(glob.glob('/net/nieuwerijn/data2/rtimmerman/A1795_HBA/A1795/selfcal/selfcal_pix0.15_wide-????-model.fits'))
    
    for image_id, image in enumerate(imagenames):
        print('remove negatives from model: ', image)
        hdul = fits.open(image)
        data = hdul[0].data
        
        data[np.where(data < 0.0)] = 0.0
        hdul[0].data = data
        hdul.writeto(image, overwrite=True)
        hdul.close()
    
        if perseus:
          os.system('python /net/rijn/data2/rvweeren/LoTSS_ClusterCAL/editmodel.py {} /net/ouderijn/data2/rvweeren/PerseusHBA/inner_ring_j2000.reg /net/ouderijn/data2/rvweeren/PerseusHBA/outer_ring_j2000.reg'.format(image))
         #os.system('python /net/rijn/data2/rvweeren/LoTSS_ClusterCAL/editmodel.py ' + image)
        if A1795: 
          cmdA1795 = 'python /net/rijn/data2/rvweeren/LoTSS_ClusterCAL/insert_highres.py '
          cmdA1795 +=  image + ' '
          cmdA1795 +=  A1795imlist[image_id] + ' '
          cmdA1795 += '/net/rijn/data2/rvweeren/LoTSS_ClusterCAL/A1795core.reg '
          print(cmdA1795)
          os.system(cmdA1795)
    
    return

def makeimage(mslist, imageout, pixsize, imsize, channelsout, niter, robust, \
              uvtaper=False, multiscale=False, predict=True, onlypredict=False, fitsmask=None, \
              idg=False, deepmultiscale=False, uvminim=80, fitspectralpol=True, \
              fitspectralpolorder=3, imager='WSCLEAN', restoringbeam=15, automask=2.5, \
              removenegativecc=True, usewgridder=False, paralleldeconvolution=0, \
              deconvolutionchannels=0, parallelgridding=1, multiscalescalebias=0.8):
    fitspectrallogpol = False # for testing Perseus
    msliststring = ' '.join(map(str, mslist))
    
    #  --- predict only when starting from external model images ---
    if onlypredict:
      if predict:
        cmd = 'wsclean -channels-out ' + str(channelsout) + ' -padding 1.8 -predict ' 
        if idg:
          cmd += '-use-idg -grid-with-beam -use-differential-lofar-beam -idg-mode cpu '
          cmd += '-beam-aterm-update 800 '
          cmd += '-pol iquv '
        else:
          if usewgridder:
            cmd +='-use-wgridder '  
          if parallelgridding > 1:
            cmd += '-parallel-gridding ' + str(parallelgridding) + ' ' 
        cmd += '-name ' + imageout + ' ' + msliststring
        print('PREDICT STEP: ', cmd)
        os.system(cmd)    
      return    
    #  --- end predict only ---
    
    os.system('rm -f ' + imageout + '-*.fits')
    imcol = 'CORRECTED_DATA'
    t = pt.table(mslist[0],readonly=True) # just test for first ms in mslist
    colnames =t.colnames()
    if 'CORRECTED_DATA' not in colnames: # for first imaging run
      imcol = 'DATA'
    t.close()
    baselineav = str (1.5e3*60000.*2.*np.pi *np.float(1.5)/(24.*60.*60*np.float(imsize)) )
   
    if imager == 'WSCLEAN':
      cmd = 'wsclean '
      cmd += '-no-update-model-required -minuv-l ' + str(uvminim) + ' '
      cmd += '-size ' + str(np.int(imsize)) + ' ' + str(np.int(imsize)) + ' -reorder '
      cmd += '-weight briggs ' + str(robust) + ' -weighting-rank-filter 3 -clean-border 1 -parallel-reordering 4 '
      cmd += '-mgain 0.8 -fit-beam -data-column ' + imcol +' -join-channels -channels-out '
      cmd += str(channelsout) + ' -padding 1.4 '
      if paralleldeconvolution > 0:
        cmd += '-parallel-deconvolution ' +  str(paralleldeconvolution) + ' '
      if parallelgridding > 1:
        cmd += '-parallel-gridding ' + str(parallelgridding) + ' '  
      if deconvolutionchannels > 0:
        cmd += '-deconvolution-channels ' +  str(deconvolutionchannels) + ' '
      if automask > 0.5:
        cmd += '-auto-mask '+ str(automask)  + ' -auto-threshold 0.5 ' # to avoid automask 0
      
      if multiscale:
         #cmd += '-multiscale '+' -multiscale-scales 0,4,8,16,32,64 -multiscale-scale-bias 0.6 '
         #cmd += '-multiscale '+' -multiscale-scales 0,6,12,16,24,32,42,64,72,128,180,256,380,512,650 '
         cmd += '-multiscale '
         cmd += '-multiscale-scale-bias ' + str(multiscalescalebias) + ' '
         cmd += '-multiscale-max-scales ' + str(np.int(np.rint(np.log2(np.float(imsize)) -3))) + ' '
      if fitsmask != None:
        if os.path.isfile(fitsmask): 
          cmd += '-fits-mask '+ fitsmask + ' '
        else:
          print('fitsmask: ', fitsmask, 'does not exist')
          sys.exit(1)
      if uvtaper:
         cmd += '-taper-gaussian 15arcsec '

      if idg:
        cmd += '-use-idg -grid-with-beam -use-differential-lofar-beam -idg-mode cpu '
        cmd += '-beam-aterm-update 800 '
        cmd += '-pol iquv -link-polarizations i '
      else:
        if fitspectralpol:
           if fitspectrallogpol: 
             cmd += '-fit-spectral-log-pol ' + str(fitspectralpolorder) + ' '   
           else:
             cmd += '-fit-spectral-pol ' + str(fitspectralpolorder) + ' '        
        cmd += '-pol i '
        cmdbtmp = '-baseline-averaging ' + baselineav + ' '
        cmd += '-baseline-averaging ' + baselineav + ' '
        if usewgridder:
          cmd +='-use-wgridder '  
          #cmd +='-wgridder-accuracy 1e-4 '
    
      cmd += '-name ' + imageout + ' -scale ' + str(pixsize) + 'arcsec ' 
      print('WSCLEAN: ', cmd + '-niter ' + str(niter) + ' ' + msliststring)
      logger.info(cmd + '-niter ' + str(niter) + ' ' + msliststring)
      os.system(cmd + '-niter ' + str(niter) + ' ' + msliststring)        
        

      if deepmultiscale:
        
        # predict first to fill MODEL_DATA so we can continue with clean
        cmdp = 'wsclean -size ' 
        cmdp += str(np.int(imsize)) + ' ' + str(np.int(imsize)) + ' -channels-out ' + str(channelsout) + ' -padding 1.8 -predict ' 
        if idg:
          cmdp += '-use-idg -grid-with-beam -use-differential-lofar-beam -idg-mode cpu '
          cmdp += '-beam-aterm-update 800 '
          cmdp += '-pol iquv '
        else:
          if usewgridder:    
            cmd +='-use-wgridder '  
            #cmd +='-wgridder-accuracy 1e-4 '
          
        cmdp += '-name ' + imageout + ' -scale ' + str(pixsize) + 'arcsec ' + msliststring
        print('PREDICT STEP for continue: ', cmdp)
        os.system(cmdp)
       
        # NOW continue cleaning
        if not multiscale: # if multiscale is true then this is already set above
          #cmd += '-multiscale '+' -multiscale-scales 0,4,8,16,32,64 '
          cmd += '-multiscale '
          cmd += '-multiscale-scale-bias ' + str(multiscalescalebias) + ' '
          cmd += '-multiscale-max-scales ' + str(np.int(np.rint(np.log2(np.float(imsize)) -3))) + ' '
        cmd += '-niter ' + str(np.int(niter/5)) + ' -continue ' + msliststring

        # Remove baselinedepedent averaging because of -continue from MODEL_DATA  
        if not idg:
          cmd = cmd.replace(cmdbtmp,'')
        print('WSCLEAN continue: ', cmd)
        os.system(cmd)

      # REMOVE nagetive model components, these are artifacts (only for Stokes I)
      if removenegativecc:
        if idg:
            removenegativefrommodel(sorted(glob.glob(imageout +'-????-I-model*.fits')))  # only Stokes I
        else:    
            removenegativefrommodel(sorted(glob.glob(imageout + '-????-model.fits')))


      if predict:
        cmd = 'wsclean -size ' 
        cmd += str(np.int(imsize)) + ' ' + str(np.int(imsize)) + ' -channels-out ' + str(channelsout) + ' -padding 1.8 -predict ' 
        if idg:
          cmd += '-use-idg -grid-with-beam -use-differential-lofar-beam -idg-mode cpu '
          cmd += '-beam-aterm-update 800 '
          cmd += '-pol iquv '
        else:
          if usewgridder:
            cmd +='-use-wgridder '  
            #cmd +='-wgridder-accuracy 1e-4 '    
          if parallelgridding > 1:
            cmd += '-parallel-gridding ' + str(parallelgridding) + ' ' 

      
        cmd += '-name ' + imageout + ' -scale ' + str(pixsize) + 'arcsec ' + msliststring
        print('PREDICT STEP: ', cmd)
        os.system(cmd)
        
        
    if imager == 'DDFACET':
        makemslist(mslist)
        #restoringbeam = '15'
        cmd = 'DDF.py --Data-MS=mslist.txt --Deconv-PeakFactor=0.001 --Data-ColName=' + imcol + ' ' + \
              '--Parallel-NCPU=32 --Output-Mode=Clean --Deconv-CycleFactor=0 ' + \
              '--Deconv-MaxMinorIter=' + str(niter) + ' --Deconv-MaxMajorIter=5 ' + \
              '--Deconv-Mode=SSD --Weight-Robust=' + str(robust) + ' --Image-NPix=' + str(np.int(imsize)) + ' ' + \
              '--CF-wmax=50000 --CF-Nw=100 --Beam-Model=None --Beam-LOFARBeamMode=A --Beam-NBand=1 ' + \
              '--Output-Also=onNeds --Image-Cell=' + str(pixsize) + ' --Facets-NFacets=1 --Freq-NDegridBand=1 ' + \
              '--Deconv-RMSFactor=3.0 --Deconv-FluxThreshold=0.0 --Data-Sort=1 --Cache-Dir=. --Freq-NBand=2 ' + \
              '--GAClean-MinSizeInit=10 --Facets-DiamMax=1.5 --Facets-DiamMin=0.1 ' + \
              '--Cache-Dirty=auto --Weight-ColName=WEIGHT_SPECTRUM --Output-Name=' + imageout + ' ' + \
              '--Comp-BDAMode=1 --DDESolutions-DDModeGrid=AP --DDESolutions-DDModeDeGrid=AP --Cache-Reset=1 ' + \
              '--RIME-ForwardMode=BDA-degrid --Predict-ColName=MODEL_DATA --Selection-UVRange=[0.1,2000.] ' + \
              '--Output-RestoringBeam=' + str(restoringbeam) + ' --Mask-SigTh=5.0 '
        if fitsmask != None:
           cmd += '--Mask-External=' + fitsmask + ' --Mask-Auto=0 ' 
        else:
           cmd += '--Mask-Auto=1 '
        
        print(cmd)
        os.system(cmd)


def calibrateandapplycal(mslist, selfcalcycle, args, solint_list, nchan_list, soltype_list, soltypecycles_list, \
              smoothnessconstraint_list, smoothnessreffrequency_list, antennaconstraint_list, uvmin=0, normamps=False, skymodel=None, predictskywithbeam=False, restoreflags=False, \
              flagging=False, longbaseline=False, BLsmooth=False, flagslowphases=True, flagslowamprms=7.0, flagslowphaserms=7.0, skymodelsource=None, skymodelpointsource=None, wscleanskymodel=None):

   soltypecycles_list_array = np.array(soltypecycles_list) # needed to slice (slicing does not work in nested l
   incol = [] # len(mslist)
   pertubation = [] # len(mslist)
   for ms in mslist:
     incol.append('DATA') # start here, will be updated at applycal step for next solve if needed
     pertubation.append(False) 
   
   parmdbmergelist =  [[] for x in range(len(mslist))]   #  [[],[],[],[]] nested list length mslist used for Jurjen's h5_merge
   # LOOP OVER THE ENTIRE SOLTYPE LIST (so includes pertubations via a pre-applycal)
   for soltypenumber, soltype in enumerate(soltype_list):
     # SOLVE LOOP OVER MS
     parmdbmslist = []
     for msnumber, ms in enumerate(mslist):
       # check we are above far enough in the selfcal to solve for the extra pertubation
       if selfcalcycle >= soltypecycles_list[soltypenumber][msnumber]: 
         print('selfcalcycle, soltypenumber',selfcalcycle, soltypenumber)
         if (soltypenumber < len(soltype_list)-1):
             
           print(selfcalcycle,soltypecycles_list[soltypenumber+1][msnumber])
           print('Array soltypecycles_list ahead',soltypecycles_list_array[soltypenumber+1:len(soltypecycles_list_array[:,0]),msnumber])
           #if (selfcalcycle >= soltypecycles_list[soltypenumber+1][msnumber]): # this looks one soltpype ahead...hmmm, not good 
           if selfcalcycle >= np.min(soltypecycles_list_array[soltypenumber+1:len(soltypecycles_list_array[:,0]),msnumber]): # this looks all soltype ahead   
             pertubation[msnumber] = True
           else:
             pertubation[msnumber] = False   
         else:
           pertubation[msnumber] = False     
             
         if skymodel != None and selfcalcycle == 0:  
           parmdb = soltype + str(soltypenumber) + '_skyselfcalcyle' + str(selfcalcycle).zfill(3) + '_' + ms + '.h5'
         else:
           parmdb = soltype + str(soltypenumber) + '_selfcalcyle' + str(selfcalcycle).zfill(3) + '_' + ms + '.h5'
          
         runDPPPbase(ms, solint_list[soltypenumber][msnumber], nchan_list[soltypenumber][msnumber], parmdb, soltype, \
                     longbaseline=longbaseline, uvmin=uvmin, \
                     SMconstraint=smoothnessconstraint_list[soltypenumber][msnumber], \
                     SMconstraintreffreq=smoothnessreffrequency_list[soltypenumber][msnumber],\
                     antennaconstraint=antennaconstraint_list[soltypenumber][msnumber], \
                     restoreflags=restoreflags, maxiter=100, flagging=flagging, skymodel=skymodel, \
                     flagslowphases=flagslowphases, flagslowamprms=flagslowamprms, \
                     flagslowphaserms=flagslowphaserms, incol=incol[msnumber], \
                     predictskywithbeam=predictskywithbeam, BLsmooth=BLsmooth, skymodelsource=skymodelsource, \
                     skymodelpointsource=skymodelpointsource, wscleanskymodel=wscleanskymodel)
         parmdbmslist.append(parmdb)
         parmdbmergelist[msnumber].append(parmdb) # for h5_merge
       
     # NORMALIZE amplitudes
     if normamps and (soltype in ['complexgain','scalarcomplexgain','rotation+diagonal',\
                                  'amplitudeonly','scalaramplitude']) and len(parmdbmslist) > 0:
       print('Doing global gain normalization')
       normamplitudes(parmdbmslist) # list of h5 for different ms, all same soltype

     # APPLYCAL or PRE-APPLYCAL
     count = 0
     for msnumber, ms in enumerate(mslist):
       if selfcalcycle >= soltypecycles_list[soltypenumber][msnumber]: #
         print(pertubation[msnumber], parmdbmslist[count], msnumber, count)
         if pertubation[msnumber]: # so another solve follows after this
           if soltypenumber == 0:  
             applycal(ms, parmdbmslist[count], msincol='DATA',msoutcol='CORRECTED_PREAPPLY' + str(soltypenumber))
           else:
#             applycal(ms, parmdbmslist[count], msincol='CORRECTED_PREAPPLY' + str(soltypenumber-1),\
#                      msoutcol='CORRECTED_PREAPPLY' + str(soltypenumber))   
             applycal(ms, parmdbmslist[count], msincol=incol[msnumber], msoutcol='CORRECTED_PREAPPLY' + str(soltypenumber)) # msincol gets incol from previous solve 
           incol[msnumber] = 'CORRECTED_PREAPPLY' + str(soltypenumber) # SET NEW incol for next solve
         else: # so this is the last solve, no other pertubation
           if soltypenumber == 0:  
             applycal(ms, parmdbmslist[count], msincol='DATA',msoutcol='CORRECTED_DATA')
           else:
             #applycal(ms, parmdbmslist[count], msincol='CORRECTED_PREAPPLY' + str(soltypenumber-1),msoutcol='CORRECTED_DATA')
             applycal(ms, parmdbmslist[count], msincol=incol[msnumber], msoutcol='CORRECTED_DATA') # msincol gets incol from previous solve
         count = count + 1 # extra counter because parmdbmslist can have less length than mslist as soltypecycles_list goes per ms
   

   # merge all solutions
   print(parmdbmergelist)
   #try:
   if True:
     import h5_merger
     for msnumber, ms in enumerate(mslist):
       if skymodel != None and selfcalcycle == 0: 
         parmdbmergename = 'merged_skyselfcalcyle' + str(selfcalcycle).zfill(3) + '_' + ms + '.h5'
       else:
         parmdbmergename = 'merged_selfcalcyle' + str(selfcalcycle).zfill(3) + '_' + ms + '.h5'    
       if os.path.isfile(parmdbmergename):
         os.system('rm -f ' + parmdbmergename)
       
       # add extra from preapplyH5_list
       if args['preapplyH5_list'][0] != None:
         preapplyh5parm = time_match_mstoH5(args['preapplyH5_list'], ms)  
         # replace the source direction coordinates so that the merge goes correctly
         #copy_over_sourcedirection_h5(parmdbmergelist[msnumber][0], preapplyh5parm)
         parmdbmergelist[msnumber].append(preapplyh5parm)
       
       print(parmdbmergename,parmdbmergelist[msnumber],ms)
       h5_merger.merge_h5(h5_out=parmdbmergename,h5_tables=parmdbmergelist[msnumber],ms_files=ms,\
                          convert_tec=True, make_new_direction=False)
       
       #testing only
       #applycal(ms, parmdbmergename, msincol='DATA',msoutcol='CORRECTED_DATA')
   #except:
   #  pass 
 
   return 


def predictsky(ms, skymodel, modeldata='MODEL_DATA', predictskywithbeam=False, sources=None):
   
   if skymodel.split('.')[-1] != 'sourcedb':
      #make sourcedb
      sourcedb = skymodel + 'sourcedb'
      if os.path.isfile(sourcedb):
         os.system('rm -rf ' + sourcedb)
      cmdmsdb = "makesourcedb in=" + skymodel + " "
      cmdmsdb += "out=" + sourcedb + " outtype='blob' format='<' append=False"
      print(cmdmsdb)
      os.system(cmdmsdb)
   else:
      sourcedb = skymodel    
   
   
   cmd = 'DPPP numthreads='+str(multiprocessing.cpu_count())+ ' msin=' + ms + ' msout=. ' 
   cmd += 'p.sourcedb=' + sourcedb + ' steps=[p] p.type=predict msout.datacolumn=' + modeldata + ' '
   if sources != None:
      cmd += 'p.sources=[' + str(sources) + '] '    
   if predictskywithbeam:
      cmd += 'p.usebeammodel=True p.usechannelfreq=True p.beammode=array_factor ' 
   print(cmd)
   os.system(cmd)
   return    

def runDPPPbase(ms, solint, nchan, parmdb, soltype, longbaseline=False, uvmin=0, \
                SMconstraint=0.0, SMconstraintreffreq=0.0, antennaconstraint=None, restoreflags=False, \
                maxiter=100, flagging=False, skymodel=None, flagslowphases=True, \
                flagslowamprms=7.0, flagslowphaserms=7.0, incol='DATA', \
                predictskywithbeam=False, BLsmooth=False, skymodelsource=None, \
                skymodelpointsource=None, wscleanskymodel=None):
    
    soltypein = soltype # save the input soltype is as soltype could be modified (for example by scalarphasediff)
    
    modeldata = 'MODEL_DATA' # the default, update if needed for scalarphasediff and phmin solves
    if BLsmooth:
      os.system('python BLsmooth.py -n 8 -i '+ incol + ' -o SMOOTHED_DATA ' + ms)        
      incol = 'SMOOTHED_DATA'    

    if soltype == 'scalarphasediff':
      create_phasediff_column(ms, incol=incol)
      soltype = 'phaseonly' # do this type of solve, maybe scalarphase is fine? 'scalarphase' #
      incol='DATA_CIRCULAR_PHASEDIFF'
      skymodel = None # solve out of MODEL_DATA complex(1,0)
      create_MODEL_DATA_PDIFF(ms)
      modeldata = 'MODEL_DATA_PDIFF'

    if skymodel !=None and soltypein != 'scalarphasediff':
        predictsky(ms, skymodel, modeldata='MODEL_DATA', predictskywithbeam=predictskywithbeam, sources=skymodelsource)

    if wscleanskymodel !=None and soltypein != 'scalarphasediff':
        makeimage([ms], wscleanskymodel, 1., 1., len(glob.glob(wscleanskymodel + '-????-model.fits')), 0, 0.0, \
               onlypredict=True, idg=False, usewgridder=True)



    if skymodelpointsource !=None and soltypein != 'scalarphasediff':
        # create MODEL_DATA (no dysco!)
        os.system('DPPP msin=' + ms + ' msout=. msout.datacolumn=MODEL_DATA steps=[]')
        # do the predict with taql
        os.system("taql" + " 'update " + ms + " set MODEL_DATA[,0]=(" + str(skymodelpointsource)+ "+0i)'")
        os.system("taql" + " 'update " + ms + " set MODEL_DATA[,3]=(" + str(skymodelpointsource)+ "+0i)'")
        os.system("taql" + " 'update " + ms + " set MODEL_DATA[,1]=(0+0i)'")
        os.system("taql" + " 'update " + ms + " set MODEL_DATA[,2]=(0+0i)'")
        

    if soltype in ['phaseonly_phmin', 'rotation_phmin', 'tec_phmin', 'tecandphase_phmin','scalarphase_phmin']:
      create_phase_column(ms, incol=incol, outcol='DATA_PHASEONLY')
      create_phase_column(ms, incol='MODEL_DATA', outcol='MODEL_DATA_PHASEONLY')
      soltype = soltype.split('_phmin')[0]
      incol = 'DATA_PHASEONLY'
      modeldata = 'MODEL_DATA_PHASEONLY'

    if soltype in ['phaseonly_slope', 'scalarphase_slope']:
      create_phase_slope(ms, incol=incol, outcol='DATA_PHASE_SLOPE', ampnorm=False)
      create_phase_slope(ms, incol='MODEL_DATA', outcol='MODEL_DATA_PHASE_SLOPE', ampnorm=False)
      soltype = soltype.split('_slope')[0]
      incol = 'DATA_PHASE_SLOPE'
      modeldata = 'MODEL_DATA_PHASE_SLOPE'      


    if soltype in ['phaseonly','complexgain','fulljones','rotation+diagonal','amplitudeonly']: # for 1D plotting
      onepol = False
    if soltype in ['scalarphase','tecandphase','tec','scalaramplitude','scalarcomplexgain','rotation']:
      onepol = True

    if restoreflags:
      cmdtaql = "'update " + ms + " set FLAG=FLAG_BACKUP'"
      print("Restore flagging column: " + "taql " + cmdtaql)
      os.system("taql " + cmdtaql)  


    t    = pt.table(ms + '/SPECTRAL_WINDOW',ack=False)
    freq = np.median(t.getcol('CHAN_FREQ')[0])
    t.close()
    
    t          = pt.table(ms + '/ANTENNA',ack=False)
    antennasms = t.getcol('NAME')
    t.close()
    if freq > 100e6:
      HBAorLBA = 'HBA'
    else:
      HBAorLBA = 'LBA'   
    print('This is', HBAorLBA, 'data')
    print('This ms contains', antennasms)

    # determine if phases needs to be included, important if slowgains do not contain phase solutions    
    includesphase = True
    if soltype == 'scalaramplitude' or soltype == 'amplitudeonly':
      includesphase = False

    # figure out which weight_spectrum column to use
    t = pt.table(ms)
    if 'WEIGHT_SPECTRUM_SOLVE' in t.colnames():
       weight_spectrum =  'WEIGHT_SPECTRUM_SOLVE'
    else:
       weight_spectrum =  'WEIGHT_SPECTRUM'
    t.close()   
    
    # check for previous old parmdb and remove them   
    if os.path.isfile(parmdb):
      print('H5 file exists  ', parmdb)
      os.system('rm -f ' + parmdb)
     
    cmd = 'DPPP numthreads='+str(multiprocessing.cpu_count())+ ' msin=' + ms + ' msin.datacolumn=' + incol + ' '
    cmd += 'msout=. ddecal.mode=' + soltype + ' '
    cmd += 'msin.weightcolumn='+weight_spectrum + ' '
    cmd += 'steps=[ddecal] ' + 'msout.storagemanager=dysco ddecal.type=ddecal '
    cmd += 'ddecal.maxiter='+str(np.int(maxiter)) + ' ddecal.propagatesolutions=True '
    cmd += 'ddecal.usemodelcolumn=True '
    cmd += 'msin.modelcolumn=' + modeldata + ' '  
    cmd += 'ddecal.solint=' + str(solint) + ' '
    cmd += 'ddecal.nchan=' + str(nchan) + ' '
    cmd += 'ddecal.h5parm=' + parmdb + ' '
    
    if uvmin != 0:
        cmd += 'ddecal.uvlambdamin=' + str(uvmin) + ' '      

    if antennaconstraint != None:
        cmd += 'ddecal.antennaconstraint=' + antennaconstraintstr(antennaconstraint, antennasms, HBAorLBA) + ' '

    if SMconstraint > 0.0 and nchan != 0:
        cmd += 'ddecal.smoothnessconstraint=' + str(SMconstraint*1e6) + ' ' 
        cmd += 'ddecal.smoothnessreffrequency=' + str(SMconstraintreffreq*1e6) + ' ' 
        
    if soltype in ['phaseonly','scalarphase','tecandphase','tec','rotation']:
       cmd += 'ddecal.tolerance=1.e-4 '
       if soltype in ['tecandphase','tec']:             
          cmd += 'ddecal.approximatetec=True '
          cmd += 'ddecal.stepsize=0.2 '
          cmd += 'ddecal.maxapproxiter=45 '
          cmd += 'ddecal.approxtolerance=6e-3 '
    if soltype in ['complexgain','scalarcomplexgain','scalaramplitude','amplitudeonly','rotation+diagonal','fulljones']:   
       cmd += 'ddecal.tolerance=1.e-4 ' # for now the same as phase soltypes
 
    print('DPPP solve:', cmd)
    os.system(cmd)
    if np.int(maxiter) == 1: # this is a template solve only
      print('Template solve, not going to make plots or do solution flagging')
      return    

    if incol == 'DATA_CIRCULAR_PHASEDIFF':
      print('Manually updating H5 to get the phase difference correct')
      makephasediffh5(parmdb)
    if incol == 'DATA_PHASE_SLOPE':
      print('Manually updating H5 to get the cumulative phase')
      makephaseCDFh5(parmdb)

    if number_freqchan_h5(parmdb) > 1:
      onechannel = False
    else:
      onechannel = True  

    # Roation checking
    if soltype in ['rotation','rotation+diagonal']:
      removenans(parmdb, 'rotation000')

    # Check for bad values  
    if soltype in ['scalarcomplexgain','complexgain','amplitudeonly','scalaramplitude','fulljones','rotation+diagonal']:
      flagbadamps(parmdb, setweightsphases=includesphase)
      if soltype != 'fulljones':
        medamp = medianamp(parmdb) # fu
      else:
        print('Not implemented, medamp of fulljones')
        sys.exit() 

      flaglowamps(parmdb, lowampval=medamp*0.1, flagging=flagging, setweightsphases=includesphase)
      flaghighgamps(parmdb, highampval=medamp*10., flagging=flagging, setweightsphases=includesphase)
      if soltype != 'amplitudeonly' and soltype != 'scalaramplitude':
         try:
           change_refant(parmdb,'phase000')
         except:
           pass
         removenans(parmdb, 'phase000')
      removenans(parmdb, 'amplitude000')

    # makes plots and do LOSOTO flagging
    outplotname = parmdb.split('_' + ms + '.h5')[0]

    if soltype in ['rotation','rotation+diagonal']:

      if check_phaseup(parmdb):
        losotoparset_rotation = create_losoto_rotationparset(ms, onechannel=onechannel, outplotname=outplotname, refant='ST001') # phase matrix plot
      else:
        losotoparset_rotation = create_losoto_rotationparset(ms, onechannel=onechannel, outplotname=outplotname) # phase matrix plot    

      cmdlosoto = 'losoto ' + parmdb + ' ' + losotoparset_rotation
      os.system(cmdlosoto)
        
    
    if soltype in ['phaseonly','scalarphase']:
      if check_phaseup(parmdb):
        losotoparset_phase = create_losoto_fastphaseparset(ms, onechannel=onechannel, onepol=onepol, outplotname=outplotname, refant='ST001') # phase matrix plot
      else:
        losotoparset_phase = create_losoto_fastphaseparset(ms, onechannel=onechannel, onepol=onepol, outplotname=outplotname) # phase matrix plot    

      cmdlosoto = 'losoto ' + parmdb + ' ' + losotoparset_phase
      os.system(cmdlosoto)


    if soltype in ['tecandphase', 'tec']:
       tecandphaseplotter(parmdb, ms, outplotname=outplotname) # use own plotter because losoto cannot add tec and phase
        
    if soltype in ['tec']:
       if check_phaseup(parmdb): 
         losotoparset_tec = create_losoto_tecparset(ms, outplotname=outplotname, refant='ST001')
       else:
         losotoparset_tec = create_losoto_tecparset(ms, outplotname=outplotname)
       cmdlosoto = 'losoto ' + parmdb + ' ' + losotoparset_tec
       os.system(cmdlosoto)    


    if soltype in ['scalarcomplexgain','complexgain','amplitudeonly','scalaramplitude', \
                   'fulljones','rotation+diagonal'] and (ntimesH5(parmdb) > 1): # plotting/flagging fails if only 1 timeslot
       print('Do flagging?:', flagging)
       if flagging and not onechannel:
          if soltype == 'fulljones':
            print('Fulljones and flagging not implemtened')
            sys.exit()
          else:    
            losotoparset = create_losoto_flag_apgridparset(ms, flagging=True, maxrms=flagslowamprms, \
                                                           maxrmsphase=flagslowphaserms, \
                                                           includesphase=includesphase, onechannel=onechannel, \
                                                           medamp=medamp, flagphases=flagslowphases, onepol=onepol,\
                                                           outplotname=outplotname)
       else:
          losotoparset = create_losoto_flag_apgridparset(ms, flagging=False, includesphase=includesphase, \
                         onechannel=onechannel, medamp=medamp, onepol=onepol, outplotname=outplotname)  

       # MAKE losoto command    
       if flagging:
         os.system('cp -f ' + parmdb + ' ' + parmdb + '.backup')
       cmdlosoto = 'losoto ' + parmdb + ' ' + losotoparset
       os.system(cmdlosoto)
    return 




# to remove H5/h5 and other files out of a wildcard selection if needed
def removenonms(mslist):
  newmslist = []  
  for ms in mslist:      
   if ms.lower().endswith(('.h5', '.png', '.parset', '.fits', '.backup', '.obj', '.log', '.p', '.reg', '.gz', '.tar', '.tmp', '.ddfcache')) or \
      ms.lower().startswith(('plotlosoto','solintimage')):
     print('WARNING, removing ', ms, 'not a ms-type? Removed it!') 
   else:
     newmslist.append(ms)  
  return newmslist     

# check is there are enough timesteps in the ms
# for example this will remove an observations of length 600s
# in that case a lot of assumptions break, for example amplitude flagging in losoto
def select_valid_ms(mslist):
  newmslist = []  
  
  for ms in mslist:
    if not os.path.isdir(ms):
      print(ms, ' does not exist')
      sys.exit(1)
  
  for ms in mslist:  
    t = pt.table(ms, ack=False)
    times = np.unique(t.getcol('TIME'))
    
    if len(times) > 30:
      newmslist.append(ms)
    else:
      print('---------------------------------------------------------------------------')  
      print('WARNING, removing ', ms, 'not enough timesteps in ms/too short observation')  
      print('---------------------------------------------------------------------------')   
    t.close()
    
    
  return newmslist    

def arg_as_list(s):                                                            
    v = ast.literal_eval(s)                                                    
    if type(v) is not list:                                                    
        raise argparse.ArgumentTypeError("Argument \"%s\" is not a list" % (s))
    return v  

def makemaskthresholdlist(maskthresholdlist, stop):
   maskthresholdselfcalcycle = []
   for mm in range(stop):
      try: 
        maskthresholdselfcalcycle.append(maskthresholdlist[mm])
      except:
        maskthresholdselfcalcycle.append(maskthresholdlist[-1]) # add last value
   return maskthresholdselfcalcycle


###############################
############## MAIN ###########
###############################

def main():
   
   #flagms_startend('P217+57_object.dysco.sub.shift.avg.weights.ms.archive0','tecandphase0_selfcalcyle1_P217+57_object.dysco.sub.shift.avg.weights.ms.archive0.h5',1)
   #sys.exit()
   
   parser = argparse.ArgumentParser(description='Self-Calibrate a facet from a LOFAR observation')
   parser.add_argument('-b','--boxfile', help='boxfile', type=str)
   parser.add_argument('--imsize', help='image size, required if boxfile is not used', type=int)
   parser.add_argument('--pixelscale','--pixelsize', help='pixels size in arcsec, default=3.0/1.5 (LBA/HBA)', type=float)
   parser.add_argument('-i','--imagename', help='imagename, default=image', default='image', type=str)
   parser.add_argument('--fitsmask', help='fitsmask for deconvolution (needs to match image size), if not provided use automasking', type=str)
   parser.add_argument('-n', '--niter', help='niter, default=15000', default=15000, type=int)
   parser.add_argument('--robust', help='Briggs robust paramter, default=-0.5', default=-0.5, type=float)
   parser.add_argument('--channelsout', help='channelsout, default=6', default=6, type=int)
   parser.add_argument('--multiscale', help='use multiscale deconvolution, not recommended/unstable', action='store_true')
   parser.add_argument('--multiscale-start', help='start multiscale deconvolution at this selfcal cycle (default=1)', default=1, type=int)
   parser.add_argument('--multiscalescalebias', help='multiscalescale bias scale paramter (see WSClean documentation), default=0.8', default=0.8, type=float)
   parser.add_argument('--deepmultiscale', help='do extra multiscale deconvolution on the residual', action='store_true')
   parser.add_argument('--uvminim', help='inner uv-cut for imaging in lambda, default=80', default=80., type=float)
   parser.add_argument('--usewgridder', help='use wgridder in WSClean, mainly useful for very large images (True/False, default=True)', type=ast.literal_eval, default=True)

   parser.add_argument('--phaseupstations', help='phase up to a superstation (core or superterp, default None)', default=None, type=str)
   parser.add_argument('--paralleldeconvolution', help='parallel-deconvolution size for wsclean, default=0 (means no parallel deconvolution, suggested value in about 2000, only use for very large images)', default=0, type=int)
   parser.add_argument('--parallelgridding', help='parallel-gridding for wsclean, default=1 (means no parallel gridding)', default=1, type=int)
   parser.add_argument('--deconvolutionchannels', help='deconvolution-channels value for wsclean, default=0 (means deconvolution-channels equals channels-out)', default=0, type=int)


   parser.add_argument('--idg', help='use the Image Domain gridder', action='store_true')
   parser.add_argument('--maskthreshold', help='Maskthresholds used from image1 onwards made by MakeMask.py, default= default=[5.0,4.5,4.0,4.0,3.5,3.5,3.5,...]', default=[5.0,4.5,4.0,4.0,3.5], type=arg_as_list)
   #parser.add_argument('--maskthreshold-selfcalcycle', type=arg_as_list, default=[5.0],help='Maskthresholds from image1 onwards')
   parser.add_argument('--imager', help='Imager to use WSClean or DDFACET, default WSCLEAN', default='WSCLEAN', type=str)
   parser.add_argument('--fitspectralpol', help='use fit-spectral-pol in WSClean (True/False, default=True)', type=ast.literal_eval, default=True)
   parser.add_argument('--fitspectralpolorder', help='fit-spectral-pol order for WSClean, default=3', default=3, type=int)
   parser.add_argument('--removenegativefrommodel', help='remove negative clean components in model predict (True/False, default=True)', type=ast.literal_eval, default=True)
   parser.add_argument('--autofrequencyaverage', help='By default frequency averaging is tried if it does not result in bandwidth  smearing (True/False, default=False)', type=ast.literal_eval, default=True)
   parser.add_argument('--autofrequencyaverage-calspeedup', help='Extra avearging during some selfcalcycles to speed up calibration (True/False, default=False)', type=ast.literal_eval, default=False)
   
   parser.add_argument('--avgfreqstep', help='Extra DPPP frequnecy averaging to speed up a solve, this is done before any other correction, could be useful for long baseline infield calibrators', type=int, default=None)
   parser.add_argument('--avgtimestep', help='Extra DPPP time averaging to speed up a solve, this is done before any other correction, could be useful for long baseline infield calibrators', type=int, default=None)
   parser.add_argument('--msinnchan', help='Before averarging, only take this number input channels', type=int, default=None)

   # calibration options
   parser.add_argument('-u', '--uvmin', help='inner uv-cut for calibration in lambda, default=80/350 (LBA/HBA)', type=float)
   parser.add_argument("--update-uvmin", help='Update uvmin automatically for the Dutch array', action='store_true')
   parser.add_argument("--soltype-list", type=arg_as_list, default=['tecandphase','tecandphase','scalarcomplexgain'],help="List of complexgain,scalarcomplexgain,scalaramplitude,amplitudeonly,phaseonly,fulljones,rotation,rotation+diagonal,tec,tecandphase,scalarphase,scalarphasediff,phaseonly_phmin,rotation_phmin,tec_phmin,tecandphase_phmin,scalarphase_phmin,scalarphase_slope,phaseonly_slope")
   parser.add_argument("--solint-list", type=arg_as_list, default=[1,1,120],help="List of values")
   parser.add_argument("--nchan-list", type=arg_as_list, default=[1,1,10],help="List of values")
   parser.add_argument("--smoothnessconstraint-list", type=arg_as_list, default=[0.,0.,5.],help="List of values")
   parser.add_argument("--smoothnessreffrequency-list", type=arg_as_list, default=[0.,0.,0.],help="An optional reference frequency (in MHz) for the smoothness constraint. When unequal to 0, the size of the smoothing kernel will vary over frequency by a factor of frequency/smoothnessreffrequency, i.e., the kernel will be smaller for lower frequencies, default is [0.0]")
   parser.add_argument("--antennaconstraint-list", type=arg_as_list, default=[None,None,None],help="List of values")
   parser.add_argument("--soltypecycles-list", type=arg_as_list, default=[0,999,3],help="List of values, first entry is required to be 0")
   parser.add_argument("--BLsmooth", help='Employ BLsmooth for low S/N data', action='store_true')
   parser.add_argument('--usemodeldataforsolints', help='Determine solints from MODEL_DATA', action='store_true')
   parser.add_argument('--tecfactorsolint', help='Experts only', type=float, default=1.0)
   parser.add_argument("--preapplyH5-list", type=arg_as_list, default=[None],help="List of H5 files, one per ms")


   # general options
   parser.add_argument('--skymodel', help='skymodel for first selfcalcycle', type=str)
   parser.add_argument('--skymodelsource', help='source name (string) in skymodel, default=None (means the skymodel only contains one source/patch', type=str)
   parser.add_argument('--skymodelpointsource', help='If set, start from a point source in the phase center with the flux density given by this parameter, default=None (None means do not use this option)', type=float, default=None)
      # general options
   parser.add_argument('--wscleanskymodel', help='WSclean basename for model images (for a WSClean predict)', type=str, default=None)
   parser.add_argument('--predictskywithbeam', help='predict the skymodel with the beam array factor', action='store_true')
   parser.add_argument('--startfromtgss', help='Start from TGSS skymodel for positions (boxfile required)', action='store_true')
   parser.add_argument('--tgssfitsimage', help='Start TGSS fits image for model (if not provided use SkyView', type=str)
   parser.add_argument('--no-beamcor', help='Do not correct the visilbities for the array factor', action='store_true')
   parser.add_argument('--use-dpppbeamcor', help='Use DP3 for beam correction, requires recent DP3 version and no phased-up stations', action='store_true')
   parser.add_argument('--docircular', help='Convert linear to circular correlations', action='store_true')
   parser.add_argument('--dolinear', help='Convert circular to linear correlations', action='store_true')
   parser.add_argument('--forwidefield', help='Keep solutions such that they can be used for widefield imaging/screens', action='store_true')
   parser.add_argument('--doflagging', help='Flag on complexgain solutions (True/False, default=True)', type=ast.literal_eval, default=True)
   parser.add_argument('--restoreflags', help='Restore flagging column after each selfcal cycle, only relevant if --doflagging=True', action='store_true')
   parser.add_argument('--flagslowamprms', help='RMS outlier value to flag on slow amplitudes (default=7.0)', default=7.0, type=float)
   parser.add_argument('--flagslowphaserms', help='RMS outlier value to flag on slow phases (default=7.0)', default=7.0, type=float)
   parser.add_argument('--doflagslowphases', help='If solution flagging is done also flag outliers phases in the slow phase solutions (True/False, default=True)', type=ast.literal_eval, default=True)
   parser.add_argument('--useaoflagger', help='Run AOflagger on input data', action='store_true')
   parser.add_argument('--useaoflaggerbeforeavg', help='run before (True) or after averaging (False), default=True', type=ast.literal_eval, default=True)
   parser.add_argument('--normamps', help='Normalize global amplitudes to 1.0 (True/False, default=True, turned off if fulljones is used)', type=ast.literal_eval, default=True)
   parser.add_argument('--normampsskymodel', help='Normalize global amplitudes to 1.0 when solving against an external skymodel (True/False, default=False, turned off if fulljones is used)', type=ast.literal_eval, default=False)
   parser.add_argument('--resetweights', help='If you want to ignore weight_spectrum_solve', action='store_true')
   parser.add_argument('--start', help='Start selfcal cycle at this iteration, default=0', default=0, type=int)
   parser.add_argument('--stop', help='Stop selfcal cycle at this iteration, default=10', default=10, type=int)
   parser.add_argument('--stopafterskysolve', help='Stop calibration after solving against external skymodel', action='store_true')
   parser.add_argument('--noarchive', help='Do not archive the data', action='store_true')
   parser.add_argument('--helperscriptspath', help='location were additional helper scripts are located', default='/net/rijn/data2/rvweeren/LoTSS_ClusterCAL/', type=str)
   parser.add_argument('--auto', help='Trigger fully automated processing (still under construction, HBA-dutch only)', action='store_true')

   parser.add_argument('ms', nargs='+', help='msfile(s)')  

   args = vars(parser.parse_args())
   options = parser.parse_args() # start of replacing args dictionary with objects options
   #print (options.preapplyH5_list)

   version = '2.9.3'
   print_title(version)

   os.system('cp ' + args['helperscriptspath'] + '/lib_multiproc.py .')
   os.system('cp ' + args['helperscriptspath'] + '/h5_merger.py .')
   os.system('cp ' + args['helperscriptspath'] + '/plot_tecandphase.py .')
   os.system('cp ' + args['helperscriptspath'] + '/lin2circ.py .')
   os.system('cp ' + args['helperscriptspath'] + '/BLsmooth.py .')

   inputchecker(args)
   maskthreshold_selfcalcycle = makemaskthresholdlist(args['maskthreshold'], args['stop'])
   

   for h5parm_id, h5parm in enumerate(args['preapplyH5_list']):
     if h5parm != None:
       os.system('cp ' + h5parm +  ' .') # make them local because we are going to update the source direction for merging    
       args['preapplyH5_list'][h5parm_id] = h5parm.split('/')[-1] # update input list to local location

   mslist = sorted(args['ms'])
   for ms_id, ms in enumerate(mslist):
      mslist[ms_id] = ms.replace('/', '') # remove possible / at end of ms

   # remove non-ms that ended up in mslist
   mslist = removenonms(mslist)

   # remove ms which are too short (to catch Elais-N1 case of 600s of data)
   mslist = sorted(select_valid_ms(mslist))

   # extra flagging if requested
   if args['start'] == 0 and args['useaoflagger'] and args['useaoflaggerbeforeavg']:  
     runaoflagger(mslist) 

   # reset weights if requested
   if args['resetweights']:
     for ms in mslist:
       cmd = "'update " + ms + " set WEIGHT_SPECTRUM=WEIGHT_SPECTRUM_SOLVE'"
       os.system("taql " + cmd)
   
   longbaseline =  checklongbaseline(mslist[0])

   # Determine HBA or LBA
   t    = pt.table(mslist[0] + '/SPECTRAL_WINDOW',ack=False)
   freq = np.median(t.getcol('CHAN_FREQ')[0])
   t.close()

   if freq < 100e6:
     LBA = True
   else: 
     LBA = False      



   # set some default values if not provided
   if args['uvmin'] == None:
     if LBA:
         if freq >= 40e6:
           args['uvmin'] = 80.
         if freq < 40e6:
           args['uvmin'] = 60.  
     else:
         args['uvmin'] = 350.

   if args['pixelscale'] == None:  
     if LBA:
       if longbaseline:
         args['pixelscale'] = 0.08 
       else:
         args['pixelscale'] = np.rint(3.0*54e6/freq)  
     else:
       if longbaseline:
         args['pixelscale'] = 0.03
       else:
         args['pixelscale'] = 1.5

   if args['boxfile'] != None:
     args['imsize']   = getimsize(args['boxfile'], args['pixelscale'])
     
   # check if we could average more
   avgfreqstep = []  # vector of len(mslist) with average values, 0 means no averaging
   for ms in mslist:
      if not longbaseline and args['avgfreqstep'] == None and args['autofrequencyaverage'] and not LBA: # autoaverage
        avgfreqstep.append(findfreqavg(ms,np.float(args['imsize'])))
      else:
        if args['avgfreqstep'] != None:
           avgfreqstep.append(args['avgfreqstep']) # take over handpicked average value
        else:
           avgfreqstep.append(0) # put to zero, zero means no average

   if args['auto']:
     args['update-uvmin'] = True
     args['usemodeldataforsolints'] = True
     args['forwidefield'] = True
     args['autofrequencyaverage'] = True
     if LBA:
       args['BLsmooth'] = True    

   # reset tecandphase -> tec for LBA 
   if LBA and args['usemodeldataforsolints']:
     args['soltype_list'][1] = 'tec'  
     #args['soltype_list'][0] = 'tec'  
     if freq < 30e6:
       args['soltype_list'] = ['tecandphase','tec']    # no scalarcomplexgain in the list

   if args['forwidefield']:
      args['doflagging'] = False
   #print( args['doflagging'])

   # average if requested
   mslist = average(mslist, freqstep=avgfreqstep, timestep=args['avgtimestep'], \
                    start=args['start'], msinnchan=args['msinnchan'])


   # extra flagging if requested
   if args['start'] == 0 and args['useaoflagger'] and not args['useaoflaggerbeforeavg']:
     runaoflagger(mslist) 


   t    = pt.table(mslist[0] + '/SPECTRAL_WINDOW',ack=False)
   bwsmear = bandwidthsmearing(np.median(t.getcol('CHAN_WIDTH')), np.min(t.getcol('CHAN_FREQ')[0]), np.float(args['imsize']))
   t.close() # close pt.table(ms + '/SPECTRAL_WINDOW',ack=False) here



   # backup flagging column for option --restoreflags if needed
   if args['restoreflags']:
     for ms in mslist:
       create_backup_flag_col(ms)
    


   automask = 2.5
   if args['maskthreshold'][-1] < automask:
     automask = args['maskthreshold'][-1] # in case we use a low value for maskthreshold, like Herc A    



   args['imagename']  = args['imagename'] + '_'
   if args['fitsmask'] != None:
     fitsmask = args['fitsmask']
   else:
     fitsmask = None



   if args['boxfile'] != None:
     outtarname = (args['boxfile'].split('/')[-1]).split('.reg')[0] + '.tar.gz'
   else:
     outtarname = 'calibrateddata' + '.tar.gz' 

   # LOG INPUT SETTINGS
   logbasicinfo(args, fitsmask, mslist, version, sys.argv)



   if args['startfromtgss'] and args['start'] == 0:
     if args['boxfile'] != None and args['skymodel'] == None:
       args['skymodel'] = makeBBSmodelforTGSS(args['boxfile'],fitsimage = args['tgssfitsimage'], \
                                              pixelscale=args['pixelscale'], imsize=args['imsize'])
     else:
       print('You need to provide a boxfile to use --startfromtgss')
       print('And you cannot provide a skymodel file manually')
       sys.exit(1)



   if args['start'] == 0:
     os.system('rm -f nchan.p solint.p smoothnessconstraint.p smoothnessreffrequency.p antennaconstraint.p soltypecycles.p') 



   nchan_list,solint_list,smoothnessconstraint_list, smoothnessreffrequency_list,  antennaconstraint_list, soltypecycles_list = \
                                              setinitial_solint(mslist, args['soltype_list'],longbaseline, LBA, \
                                              args['nchan_list'], args['solint_list'], \
                                              args['smoothnessconstraint_list'], \
                                              args['smoothnessreffrequency_list'], args['antennaconstraint_list'],\
                                              args['soltypecycles_list'])


   # Get restoring beam for DDFACET in case it is needed
   restoringbeam = calculate_restoringbeam(mslist, LBA)



   # ----- START SELFCAL LOOP -----
   for i in range(args['start'],args['stop']):

     # AUTOMATICALLY PICKUP PREVIOUS MASK (in case of a restart)
     if (i > 0) and (args['fitsmask'] == None):
       if args['idg']:  
         if os.path.isfile(args['imagename'] + str(i-1).zfill(3) + '-MFS-I-image.fits.mask.fits'):
             fitsmask = args['imagename'] + str(i-1).zfill(3) + '-MFS-I-image.fits.mask.fits'
       else:
         if args['imager'] == 'WSCLEAN':
           if os.path.isfile(args['imagename'] + str(i-1).zfill(3) + '-MFS-image.fits.mask.fits'):
               fitsmask = args['imagename'] + str(i-1).zfill(3) + '-MFS-image.fits.mask.fits'
         if args['imager'] == 'DDFACET':
           if os.path.isfile(args['imagename'] + str(i-1).zfill(3) + '.app.restored.fits'):
               fitsmask = args['imagename'] + str(i-1).zfill(3) + '.app.restored.fits.mask.fits'

       
     # BEAM CORRECTION
     if not args['no_beamcor'] and i == 0:
         for ms in mslist:
           beamcor(ms, usedppp=args['use_dpppbeamcor'])

     # CONVERT TO CIRCULAR/LINEAR CORRELATIONS      
     if (args['docircular'] or args['dolinear']) and i == 0:
         for ms in mslist:
           circular(ms, linear=args['dolinear'])

     # PRE-APPLY SOLUTIONS (from a nearby direction for example)
     if (args['preapplyH5_list'][0]) != None and i == 0:
         preapply(args['preapplyH5_list'], mslist)

     # TMP AVERAGE TO SPEED UP CALIBRATION
     if args['autofrequencyaverage_calspeedup'] and i == 0:
         avgfreqstep = []
         mslist_backup = mslist[:] # make a backup list, note copy by slicing otherwise list refers to original
         for ms in mslist:
            avgfreqstep.append(findfreqavg(ms,np.float(args['imsize']),bwsmearlimit=3.0))
         mslist = average(mslist, freqstep=avgfreqstep, timestep=4)
     if args['autofrequencyaverage_calspeedup'] and i == args['stop'] - 3:
         mslist = mslist_backup[:]  # reset back, note copy by slicing otherwise list refers to original 
         preapply(create_mergeparmdbname(mslist, i-1), mslist)

     # PHASE-UP if requested
     if args['phaseupstations'] != None and i== 0:
         mslist = phaseup(mslist,datacolumn='DATA',superstation=args['phaseupstations'])


  
     # CALIBRATE AGAINST SKYMODEL
     if (args['skymodel'] != None or args['skymodelpointsource'] != None or args['wscleanskymodel'] != None) and (i ==0):
        calibrateandapplycal(mslist, i, args, solint_list, nchan_list, args['soltype_list'], \
                             soltypecycles_list, smoothnessconstraint_list, smoothnessreffrequency_list, \
                             antennaconstraint_list, uvmin=args['uvmin'], normamps=args['normampsskymodel'], \
                             skymodel=args['skymodel'], \
                             predictskywithbeam=args['predictskywithbeam'], \
                             restoreflags=args['restoreflags'], flagging=args['doflagging'], \
                             longbaseline=longbaseline, BLsmooth=args['BLsmooth'], \
                             flagslowphases=args['doflagslowphases'], \
                             flagslowamprms=args['flagslowamprms'], flagslowphaserms=args['flagslowphaserms'],\
                             skymodelsource=args['skymodelsource'], skymodelpointsource=args['skymodelpointsource'],\
                             wscleanskymodel=args['wscleanskymodel']) 


  
     # TRIGGER MULTISCALE
     if args['multiscale'] and i >= args['multiscale_start']:
       multiscale = True
     else:
       multiscale = False  

     # IMAGE WITH WSCLEAN OR DDF.py
     makeimage(mslist, args['imagename'] + str(i).zfill(3), args['pixelscale'], args['imsize'], args['channelsout'], \
               args['niter'], args['robust'], \
               uvtaper=False, multiscale=multiscale, idg=args['idg'], fitsmask=fitsmask, \
               deepmultiscale=args['deepmultiscale'], uvminim=args['uvminim'], fitspectralpol=args['fitspectralpol'], \
               imager=args['imager'], restoringbeam=restoringbeam, automask=automask, \
               removenegativecc=args['removenegativefrommodel'], fitspectralpolorder=args['fitspectralpolorder'], \
               usewgridder=args['usewgridder'], paralleldeconvolution=args['paralleldeconvolution'],\
               deconvolutionchannels=args['deconvolutionchannels'], parallelgridding=args['parallelgridding'],\
               multiscalescalebias=args['multiscalescalebias'])
  
     # MAKE FIGURE WITH APLPY
     if args['imager'] == 'WSCLEAN':
       if args['idg']:
         plotimage(args['imagename'] + str(i).zfill(3) +'-MFS-I-image.fits',args['imagename'] + str(i).zfill(3) + '.png' , \
                   mask=fitsmask, rmsnoiseimage=args['imagename'] + str(0).zfill(3) +'-MFS-I-image.fits')
       else:
         plotimage(args['imagename'] + str(i).zfill(3) +'-MFS-image.fits',args['imagename'] + str(i) + '.png' , \
                   mask=fitsmask, rmsnoiseimage=args['imagename'] + str(0).zfill(3) +'-MFS-image.fits')
     if args['imager'] == 'DDFACET':
       plotimage(args['imagename'] + str(i).zfill(3) +'.app.restored.fits',args['imagename'] + str(i) + '.png' , \
                   mask=fitsmask, rmsnoiseimage=args['imagename'] + str(0).zfill(3) +'.app.restored.fits')

     if args['stopafterskysolve']:
       print('Stopping as requested via --stopafterskysolve')
       sys.exit(0)
     
     # REDETERMINE SOLINTS IF REQUESTED
     if (i >= 0) and (args['usemodeldataforsolints']):
       print('Recomputing solints .... ')
       nchan_list,solint_list,smoothnessconstraint_list,smoothnessreffrequency_list, \
       antennaconstraint_list,soltypecycles_list  = \
                              auto_determinesolints(mslist, args['soltype_list'], \
                              longbaseline, LBA, \
                              innchan_list=nchan_list, insolint_list=solint_list, \
                              insmoothnessconstraint_list=smoothnessconstraint_list, \
                              insmoothnessreffrequency_list=smoothnessreffrequency_list,\
                              inantennaconstraint_list=antennaconstraint_list, \
                              insoltypecycles_list=soltypecycles_list, redo=True, \
                              tecfactorsolint=args['tecfactorsolint'])  

     # CALIBRATE AND APPLYCAL
     calibrateandapplycal(mslist, i, args, solint_list, nchan_list, args['soltype_list'], soltypecycles_list,\
                           smoothnessconstraint_list, smoothnessreffrequency_list, antennaconstraint_list, uvmin=args['uvmin'], \
                           normamps=args['normamps'], restoreflags=args['restoreflags'], \
                           flagging=args['doflagging'], longbaseline=longbaseline, \
                           BLsmooth=args['BLsmooth'], flagslowphases=args['doflagslowphases'], \
                           flagslowamprms=args['flagslowamprms'], flagslowphaserms=args['flagslowphaserms'])


 
     # MAKE MASK AND UPDATE UVMIN IF REQUESTED
     if args['fitsmask'] == None:
       if args['imager'] == 'WSCLEAN':   
         if args['idg']:  
           imagename  = args['imagename'] + str(i).zfill(3) + '-MFS-I-image.fits'
         else:
           imagename  = args['imagename'] + str(i).zfill(3) + '-MFS-image.fits'
       if args['imager'] == 'DDFACET':
         imagename  = args['imagename'] + str(i).zfill(3) +'.app.restored.fits'

       if maskthreshold_selfcalcycle[i] > 0.0:    
         cmdm  = 'MakeMask.py --Th='+ str(maskthreshold_selfcalcycle[i]) + ' --RestoredIm=' + imagename
         if fitsmask != None:
           if os.path.isfile(imagename + '.mask.fits'):
             os.system('rm -f ' + imagename + '.mask.fits')
         os.system(cmdm)
         fitsmask = imagename + '.mask.fits'
      
         # update uvmin if allowed/requested
         if not longbaseline and args['update_uvmin']:
           if getlargestislandsize(fitsmask) > 1000:
             if not LBA:
               print('Extended emission found, setting uvmin to 750 klambda')
               args['uvmin'] = 750
             else:
               print('Extended emission found, setting uvmin to 250 klambda')
               args['uvmin'] = 250   
       else:
         fitsmask = None # no masking requested as args['maskthreshold'] less/equal 0
          
     # CUT FLAGGED DATA FROM MS AT START&END to win some compute time if possible
     #if TEC and not args['forwidefield']: # does not work for phaseonly sols
     #  if (i == 0) or (i == args['phasecycles']) or (i == args['phasecycles'] + 1) or (i == args['phasecycles'] + 2) \
     #    or (i == args['phasecycles'] + 3) or (i == args['phasecycles'] + 4):
     #     for msnumber, ms in enumerate(mslist): 
     #         flagms_startend(ms, 'phaseonly' + ms + parmdb + str(i) + '.h5', np.int(solint_phase[msnumber]))
  
  
   if not longbaseline and not args['noarchive'] :
     if not LBA:   
      archive(mslist, outtarname, args['boxfile'], fitsmask, imagename)    
      cleanup(mslist)

if __name__ == "__main__":
   main()
