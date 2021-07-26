#!/usr/bin/env python

# to do, more logging fixes
# auto solints tune, HBA slow
# auto solints tune, LBA slow
# auto solints tune, HBA-international slow
# second slow gain pertubation
# solnorm fulljones fix?
# pre-apply close direction?
# IDG and beamcors

# put this first to avoid issues with the logging
import logging
logging.basicConfig(filename='selfcal.log', format='%(levelname)s:%(asctime)s ---- %(message)s', datefmt='%m/%d/%Y %H:%M:%S', level=logging.DEBUG)

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



from lofar.stationresponse import stationresponse

def reset_gains_noncore(h5parm):
   hasphase = True
   hasamps  = True
   H=tables.open_file(H5file, mode='r')
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
   
   antennaxis = axisn.index('ant')
   #for antennaid,antenna in enumerate(antennas):
   
   H.close()    
   return

def phaseup(msinlist,datacolumn='DATA',superstation='core'):
  msoutlist = []
  for ms in msinlist:
    msout=ms + '.phaseup'
    msoutlist.append(msout)
    cmd = "DPPP msin=" + ms + " msout.storagemanager=dysco steps=[add,filter] "
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

def findfreqavg(ms, imsize):
    
  t = pt.table(ms + '/SPECTRAL_WINDOW',ack=False)
  bwsmear = bandwidthsmearing(np.median(t.getcol('CHAN_WIDTH')), \
            np.min(t.getcol('CHAN_FREQ')[0]), np.float(imsize), verbose=False)
  nfreq = len(t.getcol('CHAN_FREQ')[0])
  avgfactor = 0
  t.close()
  #bwsmear = 0.4
  #print(bwsmear, nfreq)  
    
  if bwsmear  < 0.5: # factor 2 avg
    if nfreq % 2 == 0:
      avgfactor = 2

  if bwsmear  < 0.333: # factor 3 avg
    if nfreq % 3 == 0:
      avgfactor = 3

  if bwsmear  < 0.25: # factor 4 avg
    if nfreq % 4 == 0:
      avgfactor = 4

  if bwsmear  < 0.2: # factor 5 avg
    if nfreq % 5 == 0:
      avgfactor = 5

  if bwsmear  < 0.0166: # factor 6 avg
    if nfreq % 6 == 0:
      avgfactor = 6            

  if bwsmear  < 0.125: # factor 8 avg
    if nfreq % 8 == 0:
      avgfactor = 8

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

def average(mslist, freqstep=None, timestep=None, start=0):
    if freqstep == None and timestep == None:
      return mslist
    
    outmslist = []
    for ms in mslist:
      msout = ms + '.avg'  
      cmd = 'DPPP msin=' + ms + ' msout.storagemanager=dysco steps=[av] av.type=averager '
      cmd+= 'msout='+ msout + ' msin.weightcolumn=WEIGHT_SPECTRUM  '
      if freqstep != None:
        cmd+='av.freqstep=' + str(freqstep) + ' '
      if timestep != None:  
        cmd+='av.timestep=' + str(timestep) + ' '
      if start == 0:
        print('Average with default WEIGHT_SPECTRUM:', cmd)  
        os.system(cmd)

      msouttmp = ms + '.avgtmp'  
      cmd = 'DPPP msin=' + ms + ' msout.storagemanager=dysco steps=[av] av.type=averager '
      cmd+= 'msout='+ msouttmp + ' msin.weightcolumn=WEIGHT_SPECTRUM_SOLVE  '
      if freqstep != None:
        cmd+='av.freqstep=' + str(freqstep) + ' '
      if timestep != None:  
        cmd+='av.timestep=' + str(timestep) + ' '
      if start == 0:
        t = pt.table(ms)
        if 'WEIGHT_SPECTRUM_SOLVE' in t.colnames(): # check if present otherwise this is not needed
          t.close()   
          print('Average with default WEIGHT_SPECTRUM_SOLVE:', cmd)  
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
    return outmslist

def makeh5templates(mslist, parmdb, TEC, puretec, solint_phase, solint_ap, nchan_phase, nchan_ap):
  for msnumber, ms in enumerate(mslist):
    if args['phase_soltype'] == 'scalarphase' and TEC == False:
       runDPPP(ms, np.int(solint_ap[msnumber]), np.int(solint_phase[msnumber]), \
               np.int(nchan_phase[msnumber]), np.int(nchan_ap[msnumber]), \
               ms + parmdb + str(0) + '_polversion.h5' ,args['phase_soltype'], \
               preapplyphase=False, TEC=TEC, puretec=args['pure_tec'], maxiter=1)
    if args['slow_soltype'] != 'fulljones':
       runDPPP(ms, np.int(solint_ap[msnumber]), np.int(solint_phase[msnumber]), \
               np.int(nchan_phase[msnumber]), np.int(nchan_ap[msnumber]), \
               ms + parmdb + str(0) + '_slowgainversion.h5' ,'complexgain', \
               preapplyphase=False, TEC=False, puretec=False, maxiter=1)
    else:  
       runDPPP(ms, np.int(solint_ap[msnumber]), np.int(solint_phase[msnumber]), \
               np.int(nchan_phase[msnumber]), np.int(nchan_ap[msnumber]), \
               ms + parmdb + str(0) + '_slowgainversion.h5' ,'fulljones', \
               preapplyphase=False, TEC=False, puretec=False, maxiter=1)
    resetgains(ms + parmdb + str(0) + '_slowgainversion.h5')
  return


def tecandphaseplotter(h5, ms):
    if not os.path.isdir('plotlosoto%s'  % ms): # needed because if this is the first plot this directory does not yet exist
      os.system('mkdir plotlosoto%s'  % ms)
    cmd = 'python /net/rijn/data2/rvweeren/LoTSS_ClusterCAL/plot_tecandphase.py  '
    cmd += '--H5file=' + h5 + ' --outfile=plotlosoto%s/fasttecandphase_nolosoto.png' % ms
    print(cmd)
    os.system(cmd)
    return

def runaoflagger(mslist):
    for ms in mslist:
       cmd = 'aoflagger ' + ms
       os.system(cmd)
    return

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

def applycal(ms, parmdblist, msoutcol='CORRECTED_DATA'):
    cmd = 'DPPP numthreads='+ str(multiprocessing.cpu_count()) + ' msin=' + ms + ' msin.datacolumn=DATA msout=. '
    cmd += 'msout.datacolumn=' + msoutcol + ' msout.storagemanager=dysco '
    count = 0
    for parmdb in parmdblist:
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
  if args['antennaconstraint_phase'] != 'superterp' and  args['antennaconstraint_phase'] != 'coreandfirstremotes' and \
     args['antennaconstraint_phase'] != 'core' and args['antennaconstraint_phase'] != 'remote' and \
     args['antennaconstraint_phase'] != 'all' and args['antennaconstraint_phase'] != 'international' and \
     args['antennaconstraint_phase'] != 'alldutch' and args['antennaconstraint_phase'] != 'core-remote' and \
     args['antennaconstraint_phase'] != None:
       print('Invalid input, --antennaconstraint_phase can only be core, superterp, coreandfirstremotes, remote, alldutch, international, or all')
       sys.exit(1)

  if args['antennaconstraint_slow'] != 'superterp' and  args['antennaconstraint_slow'] != 'coreandfirstremotes' and \
     args['antennaconstraint_slow'] != 'core' and args['antennaconstraint_slow'] != 'remote' and \
     args['antennaconstraint_slow'] != 'all' and args['antennaconstraint_slow'] != 'international' and \
     args['antennaconstraint_slow'] != 'alldutch' and args['antennaconstraint_slow'] != 'core-remote' and \
     args['antennaconstraint_slow'] != None:
       print('Invalid input, --antennaconstraint_slow can only be core, superterp, coreandfirstremotes, remote, alldutch, international, or all')
       sys.exit(1)

  if args['slow_soltype'] not in ['complexgain','scalarcomplexgain','scalaramplitude','amplitudeonly','fulljones']:
       print('Invalid input, --slow-soltype can only be complexgain, scalarcomplexgain, scalaramplitude, or amplitudeonly')
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
       #print('Your solutions contain do NOT contain TEC values')
       freq = H.root.sol000.amplitude000.freq[:] # apparently we only have amplitudes
    
    #if soltype != 'phase':
    #    freq = H.root.sol000.amplitude000.freq[:]
    #else:
    #    freq = H.root.sol000.phase000.freq[:]
    H.close()
    print('Number of frequency channels in this solutions file is:', len(freq))
    return len(freq)


def calculate_restoringbeam(mslist, LBA):
    
    if LBA: # so we have LBA
      restoringbeam = 15.
    else : # so we have HBA
      restoringbeam = 6.  
    
    return restoringbeam



def print_title():
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
                                                                              
    
                      Reinout van Weeren (2020, A&A, submitted)

                              Starting.........
          """)


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
    if ctype != 'superterp' and ctype != 'core' and ctype != 'coreandfirstremotes' and \
       ctype != 'remote' and ctype != 'alldutch' and ctype != 'all' and \
       ctype != 'international' and ctype != 'core-remote' :
        print('Invalid input, ctype can only be "superterp" or "core"')
        sys.exit(1)
    if HBAorLBA == 'LBA':  
      if ctype == 'superterp':  
        antstr=['CS002LBA','CS003LBA','CS004LBA','CS005LBA','CS006LBA','CS007LBA']
      if ctype == 'core':
        antstr=['CS001LBA','CS002LBA','CS003LBA','CS004LBA','CS005LBA','CS006LBA','CS007LBA', \
                'CS011LBA','CS013LBA','CS017LBA','CS021LBA','CS024LBA','CS026LBA','CS028LBA', \
                'CS030LBA','CS031LBA','CS032LBA','CS101LBA','CS103LBA','CS201LBA','CS301LBA', \
                'CS302LBA','CS401LBA','CS501LBA']
      if ctype == 'coreandfirstremotes':
        antstr=['CS001LBA','CS002LBA','CS003LBA','CS004LBA','CS005LBA','CS006LBA','CS007LBA', \
                'CS011LBA','CS013LBA','CS017LBA','CS021LBA','CS024LBA','CS026LBA','CS028LBA', \
                'CS030LBA','CS031LBA','CS032LBA','CS101LBA','CS103LBA','CS201LBA','CS301LBA', \
                'CS302LBA','CS401LBA','CS501LBA','RS503LBA','RS305LBA','RS205LBA','RS306LBA', \
                'RS106LBA']
      if ctype == 'remote':
        antstr=['RS503LBA','RS305LBA','RS205LBA','RS306LBA', 'RS310LBA','RS406LBA','RS407LBA',\
                'RS106LBA','RS307LBA','RS208LBA','RS210LBA', 'RS409LBA','RS508LBA','RS509LBA']
      if ctype == 'alldutch':
        antstr=['CS001LBA','CS002LBA','CS003LBA','CS004LBA','CS005LBA','CS006LBA','CS007LBA', \
                'CS011LBA','CS013LBA','CS017LBA','CS021LBA','CS024LBA','CS026LBA','CS028LBA', \
                'CS030LBA','CS031LBA','CS032LBA','CS101LBA','CS103LBA','CS201LBA','CS301LBA', \
                'CS302LBA','CS401LBA','CS501LBA','RS503LBA','RS305LBA','RS205LBA','RS306LBA', \
                'RS310LBA','RS406LBA','RS407LBA','RS106LBA','RS307LBA','RS208LBA','RS210LBA', \
                'RS409LBA','RS508LBA','RS509LBA']
      if ctype == 'all':
        antstr=['CS001LBA','CS002LBA','CS003LBA','CS004LBA','CS005LBA','CS006LBA','CS007LBA', \
                'CS011LBA','CS013LBA','CS017LBA','CS021LBA','CS024LBA','CS026LBA','CS028LBA', \
                'CS030LBA','CS031LBA','CS032LBA','CS101LBA','CS103LBA','CS201LBA','CS301LBA', \
                'CS302LBA','CS401LBA','CS501LBA','RS503LBA','RS305LBA','RS205LBA','RS306LBA', \
                'RS310LBA','RS406LBA','RS407LBA','RS106LBA','RS307LBA','RS208LBA','RS210LBA', \
                'RS409LBA','RS508LBA','RS509LBA', \
                'DE601LBA','DE602LBA','DE603LBA','DE604LBA', 'DE605LBA','DE609LBA','FR606LBA', \
                'SE607LBA','UK608LBA','PL610LBA','PL611LBA', 'PL612LBA','IE613LBA','LV614LBA']          
      if ctype == 'international':
        antstr=['DE601LBA','DE602LBA','DE603LBA','DE604LBA', 'DE605LBA','DE609LBA','FR606LBA', \
                'SE607LBA','UK608LBA','PL610LBA','PL611LBA', 'PL612LBA','IE613LBA','LV614LBA']    
      if ctype == 'core-remote':
        antstr1=['CS001LBA','CS002LBA','CS003LBA','CS004LBA','CS005LBA','CS006LBA','CS007LBA', \
                'CS011LBA','CS013LBA','CS017LBA','CS021LBA','CS024LBA','CS026LBA','CS028LBA', \
                'CS030LBA','CS031LBA','CS032LBA','CS101LBA','CS103LBA','CS201LBA','CS301LBA', \
                'CS302LBA','CS401LBA','CS501LBA']
        antstr2=['RS503LBA','RS305LBA','RS205LBA','RS306LBA', 'RS310LBA','RS406LBA','RS407LBA',\
                'RS106LBA','RS307LBA','RS208LBA','RS210LBA', 'RS409LBA','RS508LBA','RS509LBA']
          

    if HBAorLBA == 'HBA':    
      if ctype == 'superterp': 
         antstr=['CS002HBA0','CS003HBA0','CS004HBA0','CS005HBA0','CS006HBA0','CS007HBA0', \
                 'CS002HBA1','CS003HBA1','CS004HBA1','CS005HBA1','CS006HBA1','CS007HBA1']
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
                'CS302HBA1','CS401HBA1','CS501HBA1']
      if ctype == 'coreandfirstremotes':
        antstr=['CS001HBA0','CS002HBA0','CS003HBA0','CS004HBA0','CS005HBA0','CS006HBA0','CS007HBA0', \
                'CS011HBA0','CS013HBA0','CS017HBA0','CS021HBA0','CS024HBA0','CS026HBA0','CS028HBA0', \
                'CS030HBA0','CS031HBA0','CS032HBA0','CS101HBA0','CS103HBA0','CS201HBA0','CS301HBA0', \
                'CS302HBA0','CS401HBA0','CS501HBA0',\
                'CS001HBA1','CS002HBA1','CS003HBA1','CS004HBA1','CS005HBA1','CS006HBA1','CS007HBA1', \
                'CS011HBA1','CS013HBA1','CS017HBA1','CS021HBA1','CS024HBA1','CS026HBA1','CS028HBA1', \
                'CS030HBA1','CS031HBA1','CS032HBA1','CS101HBA1','CS103HBA1','CS201HBA1','CS301HBA1', \
                'CS302HBA1','CS401HBA1','CS501HBA1','RS503HBA' ,'RS305HBA' ,'RS205HBA' ,'RS306HBA',  \
                'RS106HBA']
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
                'RS106HBA','RS307HBA','RS208HBA','RS210HBA', 'RS409HBA','RS508HBA','RS509HBA']
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
                'SE607HBA','UK608HBA','PL610HBA','PL611HBA', 'PL612HBA','IE613HBA','LV614HBA']
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
                'CS302HBA1','CS401HBA1','CS501HBA1']
        antstr2=['RS503HBA','RS305HBA','RS205HBA','RS306HBA', 'RS310HBA','RS406HBA','RS407HBA', \
                'RS106HBA','RS307HBA','RS208HBA','RS210HBA', 'RS409HBA','RS508HBA','RS509HBA']

    if ctype != 'core-remote':
        antstrtmp = antstr
        for ant in antstr:
            if ant not in antennasms:
                antstrtmp.remove(ant)    
            
        antstr = ','.join(map(str, antstrtmp))
        antstr = '[[' + antstr + ']]'
    else:
        antstrtmp1 = antstr1
        for ant in antstr1:
            if ant not in antennasms:
                antstrtmp1.remove(ant)    
        antstr1 = ','.join(map(str, antstrtmp1))

        antstrtmp2 = antstr2
        for ant in antstr2:
            if ant not in antennasms:
                antstrtmp2.remove(ant)    
        antstr2 = ','.join(map(str, antstrtmp2))        
        
        antstr =  '[[' + antstr1 + '],[' + antstr2 + ']]'

    return antstr    

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
            logging.debug('Working on station number %i' % stationnum)
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
                    logging.error('Beam prediction works only for amplitude/phase solution tables.')
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
    cmd+=  'msin.weightcolumn=WEIGHT_SPECTRUM steps=[] ' 
    if starttime is not None:
      cmd+= 'msin.starttime=' + starttime + ' '
    if endtime is not None:  
      cmd+= 'msin.endtime=' + endtime   + ' '   
    print(cmd)  
    os.system(cmd)
    
    cmd = 'DPPP msin=' + ms + ' ' + 'msout.storagemanager=dysco msout=' + ms + '.cuttmp '
    cmd+= 'msin.weightcolumn=WEIGHT_SPECTRUM_SOLVE steps=[] '  
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
  
  #logging.basicConfig(level=logging.ERROR)   # to block astropy/aplpy warnings about fits headers
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
  try:
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
  #logging.basicConfig(level=logging.DEBUG)
  logging.info(fitsimagename + ' RMS noise: ' + str(imagenoiseinfo))
  return



def archive(mslist, outtarname, regionfile, fitsmask, imagename):
  path = '/disks/ftphome/pub/vanweeren'
  for ms in mslist:
    msout = ms + '.calibrated'
    if os.path.isdir(msout):
      os.system('rm -rf ' + msout)
    cmd  ='DPPP numthreads='+ str(multiprocessing.cpu_count()) +' msin=' + ms + ' msout=' + msout + ' '
    cmd +='msin.datacolumn=CORRECTED_DATA msout.storagemanager=dysco steps=[]'
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
  logging.info('Creating archived calibrated tarball: ' + outtarname)    
  os.system(cmd)
  
  for ms in mslist:
    msout = ms + '.calibrated'   
    os.system('rm -rf ' + msout)
  return


def reweight(mslist, pixsize, imsize, channelsout, niter, robust, multiscale=False, fitsmask=None):
   """
   determine the solution time and frequency intervals based on the amount of compact source flux
   """
   
   rmslist = []

   logging.info('Adjusting weights')

   for ms in mslist:
          imageout =  'rmsimage' + ms.split('.ms')[0] 
          makeimage([ms], imageout, pixsize, imsize, channelsout, np.int(niter/(len(mslist)**(1./3.))), robust, multiscale=multiscale, predict=False,fitsmask=fitsmask)
          
          hdulist = fits.open(imageout + '-MFS-image.fits')
          imagenoise = findrms(np.ndarray.flatten(hdulist[0].data))
          hdulist.close() 
          rmslist.append(imagenoise)
          
   weightslist = []       
   return 


def setinitial_solint(mslist, TEC, longbaseline, LBA,\
                      innchan_phase=None,innchan_ap=None,insolint_phase=None,insolint_ap=None):
   """
   determine the solution time and frequency intervals based on the amount of compact source flux
   Use some sensible defaults if values are not set
   """
   if os.path.isfile('nchan_phase.p') and os.path.isfile('solint_phase.p') and \
      os.path.isfile('solint_ap.p') and os.path.isfile('nchan_ap.p'):
    
      f = open('nchan_phase.p', 'r') 
      nchan_phase_F = pickle.load(f)        
      f.close()   
  
      f = open('solint_phase.p', 'r') 
      solint_phase_F = pickle.load(f)        
      f.close()   

      f = open('solint_ap.p', 'r') 
      solint_ap_F = pickle.load(f)        
      f.close()         
  
      f = open('nchan_ap.p', 'r') 
      nchan_ap_F = pickle.load(f)        
      f.close()   
  
   else:
      
      nchan_phase_F  = []
      solint_phase_F = []
      solint_ap_F    = []
      nchan_ap_F     = [] 

      for ms in mslist:
          
          #solint phase
          if insolint_phase == None:
            solint_phase = 1
          else:
            solint_phase = insolint_phase

          # solint ap
          if insolint_ap == None:
            solint_ap = 120
          else:
            solint_ap = insolint_ap

          # nchan phase
          if innchan_phase == None:
            nchan_phase = 5
          else:
            nchan_phase = innchan_phase

          # nchan ap
          if innchan_ap == None:
            nchan_ap = 20
          else:  
            nchan_ap = innchan_ap
          
          
          if TEC:
            nchan_phase  = 1

          logging.info('MS, NCHAN_PHASE - SOLINT_PHASE ||| NCHAN_SLOW - SOLINT_SLOW: ' + str(ms) + ':: ' + str(nchan_phase) + \
                       ' - ' + str (solint_phase) + ' ||| ' + str(nchan_ap) + ' - ' + str(solint_ap))
        
          nchan_phase_F.append(nchan_phase)
          solint_phase_F.append(solint_phase)
          solint_ap_F.append(solint_ap)
          nchan_ap_F.append(nchan_ap)

      f = open('nchan_phase.p', 'wb') 
      pickle.dump(nchan_phase_F,f)        
      f.close()   
  
      f = open('solint_phase.p', 'wb') 
      pickle.dump(solint_phase_F,f)        
      f.close()   

      f = open('solint_ap.p', 'wb') 
      pickle.dump(solint_ap_F,f)        
      f.close()         
  
      f = open('nchan_ap.p', 'wb') 
      pickle.dump(nchan_ap_F,f)        
      f.close()     

   return nchan_phase_F, solint_phase_F, solint_ap_F, nchan_ap_F



def auto_determinesolints(mslist, TEC, puretec, longbaseline, LBA, \
                          smoothnessconstraint_phase=0.0, smoothnessconstraint_slow=0.0,\
                          innchan_phase=5,innchan_ap=10,insolint_phase=1,\
                          insolint_ap=None, uvdismod=None, modelcolumn='MODEL_DATA', redo=False,
                          antennaconstraint=None,antennaconstraintslow=None):
   """
   determine the solution time and frequency intervals based on the amount of compact source flux and noise
   """
   

   trigger_antennaconstraint = [] # for HBA LoTSS type data
   trigger_phasecycles       = []
   
   if os.path.isfile('nchan_phase.p') and os.path.isfile('solint_phase.p') and \
      os.path.isfile('solint_ap.p') and os.path.isfile('nchan_ap.p') and not redo:
    
      f = open('nchan_phase.p', 'r') 
      nchan_phase_F = pickle.load(f)        
      f.close()   
  
      f = open('solint_phase.p', 'r') 
      solint_phase_F = pickle.load(f)        
      f.close()   

      f = open('solint_ap.p', 'r') 
      solint_ap_F = pickle.load(f)        
      f.close()         
  
      f = open('nchan_ap.p', 'r') 
      nchan_ap_F = pickle.load(f)        
      f.close()   
  
   else:

      nchan_phase_F  = []
      solint_phase_F = []
      solint_ap_F    = []
      nchan_ap_F     = [] 

      for ms in mslist:
          
          uvdismod = get_uvwmax(ms)*0.333 # take range [0.333uvmax - 1.0uvmax]
          
          t = pt.taql('SELECT ' + modelcolumn + ',DATA,UVW,TIME,FLAG FROM ' + ms + ' WHERE SQRT(SUMSQR(UVW[:2])) > '+ str(uvdismod) )
          model = np.abs(t.getcol(modelcolumn))
          flags = t.getcol('FLAG')
          data  = t.getcol('DATA')
          datashape = data.shape
          print('Compute visibility noise of the dataset with robust sigma clipping')
          #noise = astropy.stats.sigma_clipping.sigma_clipped_stats(data[:,:,1:3],mask=flags[:,:,1:3])[2] # use XY and YX
          # take only every fifth element of the array to speed up the computation
          noise = astropy.stats.sigma_clipping.sigma_clipped_stats(data[0:data.shape[0]:5,:,1:3],\
                  mask=flags[0:data.shape[0]:5,:,1:3])[2] # use XY and YX

          
          model = np.ma.masked_array(model, flags)
          flux  = np.ma.mean((model[:,:,0] + model[:,:,3])*0.5) # average XX and YY (ignore XY and YX, they are zero, or nan)
          time  = np.unique(t.getcol('TIME'))
          tint  = np.abs(time[1]-time[0])
          print('Integration time visibilities', tint)
          t.close()
          t = pt.table(ms + '/SPECTRAL_WINDOW')
          chanw = np.median(t.getcol('CHAN_WIDTH'))
          freq = np.median(t.getcol('CHAN_FREQ'))
          t.close()
  
          # -----------------SOLINT-PHASE (tested for solves that include TEC/TECANDPHASE)
          if LBA: 
            if longbaseline:
              solintphase_sf = 0.5e-3 # untested
            else: #for -- LBA dutch --
              solintphase_sf = 0.5e-3 # for tecandphase and coreconstraint
          
          else: # for -- HBA --
            if longbaseline:
              solintphase_sf = 3.0e-3 # for tecandphase, no coreconstraint          
            else: #for -- HBA dutch --
              solintphase_sf = 4.0e-2 # for tecandphase, no coreconstraint          
          
          if puretec:
            solintphase_sf = solintphase_sf/np.sqrt(2.) # tec and coreconstraint
          

          # trigger antennaconstraint_phase core if solint > 6 min
          if not LBA and not longbaseline and (solintphase_sf* ((noise/flux)**2) * (chanw/390.625e3) > 0.1):
            solintphase_sf = solintphase_sf/30. 
            trigger_antennaconstraint.append(True) # list
          else:  
            trigger_antennaconstraint.append(False) # list
            
          # round to nearest integer  
          solint_phase = np.rint(solintphase_sf* ((noise/flux)**2) * (chanw/390.625e3) )
          # frequency scaling is need because if we avearge in freqeuncy the solint should not change for a tec(andphasse) solve
          
          print('Noise visibilities:', noise)
          print('Flux in model', flux, 'Jy')
          print('UV-selection to compute model flux', str(uvdismod/1e3), 'km')
          print(solintphase_sf*((noise/flux)**2)*(chanw/390.625e3), 'want to use solint:', solint_phase)
          if solint_phase < 1:
            solint_phase = 1
          if (np.float(solint_phase)*tint/3600.) > 0.5: # so check if larger than 30 min
            print('Warning, it seems there is not enough flux density on the longer baselines for solving')
            solint_phase = np.rint(0.5*3600./tint) # max is 30 min 
          print('Conclusion: Using phase solint [s]:', np.float(solint_phase)*tint)
          

          
          nchan_ap = innchan_ap
          solint_ap = insolint_ap
          
          
          #  ----------------- SOLINT-AP ------------------- 
          if LBA: 
            if longbaseline:
              solintap_sf = 100.*0.5e-3 # untested
            else: #for -- LBA dutch --
              solintap_sf = 100.*0.5e-3 # for tecandphase and coreconstraint
          
          else: # for -- HBA --
            if longbaseline:
              solintap_sf = 0.4 # for tecandphase, no coreconstraint    
            else: #for -- HBA dutch --
              solintap_sf = 10.0 # for tecandphase, no coreconstraint          

          
          # round to nearest integer  
          solint_ap = np.rint(solintap_sf*(noise/flux)**2)      
          print(solintap_sf*(noise/flux)**2, 'want to use slow solint:', solint_ap)



          # do not allow very short ap solves
          if (np.float(solint_ap)*tint/3600.) < 0.3333: #  check if less than 20 min
            solint_ap = np.rint(0.3333*3600./tint) # min is 20 min  
          
          # do not allow ap solves that are more than 4 hrs
          if (np.float(solint_ap)*tint/3600.) > 4.0: # so check if larger than 30 min
            print('Warning, it seems there is not enough flux density on the longer baselines for slow solving')
            solint_ap = np.rint(4.0*3600./tint) # max is 4 hrs  
          print('Conlusion: Using slow gain solint [s]:', np.float(solint_ap)*tint)  
  
          # disable slow solve if teh solints get too long, target is too faint
          if not LBA and not longbaseline and ((np.float(solint_ap)*tint/3600.) > 16):
            trigger_phasecycles.append(True) # will disable slow solve
          else:  
            trigger_phasecycles.append(False)  
  
  
          # --------------- NCHAN-AP ---------------------
          
          if LBA: 
            if longbaseline:
              print('Not supported')
              sys.exit()
            else: #for -- LBA dutch, untested --
              nchanap_sf = 0.75 # for tecandphase and coreconstraint
          
          else: # for -- HBA --
            if longbaseline:
              nchanap_sf = 0.0075 #   
            else: #for -- HBA dutch --
              nchanap_sf = 0.75 #          
          
          if smoothnessconstraint_slow == 0.0:
             # round to nearest integer  
             nchan_ap = np.rint(nchanap_sf*(noise/flux)**2)
             print(nchanap_sf*(noise/flux)**2, 'so using nchan:', nchan_ap)
             print('Using nchan slow  [MHz]:', np.float(nchan_ap)*chanw/1e6) 

             # do not allow very short ap solves
             if (np.float(nchan_ap)*chanw/1e6) < 2.0: #  check if less than 2 MHz
               nchan_ap = np.rint(2.0*1e6/chanw) # 2 MHz
          
             # do not allow nchan-ap solves that are more than 15 MHz
             if (np.float(nchan_ap)*chanw/1e6) > 4.0: # so check if larger than 30 min
               print('Warning, it seems there is not enough flux density on the longer baselines for solving')
               nchan_ap = np.rint(15*1e6/chanw) # 15 MHz 
          else:
             if innchan_ap == None: # so smoothnessconstraint_slow is used
               nchan_ap     = 2 # just set if not provided, 2 should be ok?
   
  
          #if longbaseline:

              #if (flux/declf) >= 5.0:
                #nchan_ap = np.rint (10.*(195312.5)/chanw) # means 24 solutions across the freq band, assuming 24 blocks
              #if ((flux/declf) >= 1.5) and  ((flux/declf) < 5.0):
                #nchan_ap = np.rint (20.*(195312.5)/chanw) # means 12 solutions across the freq band, assuming 24 blocks              
              #if (flux/declf) < 1.5 and (flux/declf) >= 0.01:
                #nchan_ap = np.rint (40.*(195312.5)/chanw) # means 6 solutions across the freq band, assuming 24 blocks         
              #if (flux/declf) < 0.01:
                #nchan_ap = np.rint (80.*(195312.5)/chanw) # means 3 solutions across the freq band, assuming 24 blocks    
          
          if TEC:
            nchan_phase  = 1 
          
          # make integers
          solint_phase = np.int(solint_phase)
          solint_ap    = np.int(solint_ap)
          nchan_phase  = np.int(nchan_phase)
          nchan_ap     = np.int(nchan_ap)

          print('MS, NCHAN_PHASE, SOLINT_PHASE, SOLINT_AP, NCHAN_AP: ' + str(ms) + ' ' + str(nchan_phase) + \
                ' ' + str (solint_phase) + ' ' + str(solint_ap) + ' ' + str(nchan_ap))
          logging.info('MS, NCHAN_PHASE, SOLINT_PHASE, SOLINT_AP, NCHAN_AP: ' + str(ms) + ' ' + str(nchan_phase) + \
                       ' ' + str (solint_phase) + ' ' + str(solint_ap) + ' ' + str(nchan_ap))
        
          nchan_phase_F.append(nchan_phase)
          solint_phase_F.append(solint_phase)
          solint_ap_F.append(solint_ap)
          nchan_ap_F.append(nchan_ap)

      f = open('nchan_phase.p', 'wb') 
      pickle.dump(nchan_phase_F,f)        
      f.close()   
  
      f = open('solint_phase.p', 'wb') 
      pickle.dump(solint_phase_F,f)        
      f.close()   

      f = open('solint_ap.p', 'wb') 
      pickle.dump(solint_ap_F,f)        
      f.close()         
  
      f = open('nchan_ap.p', 'wb') 
      pickle.dump(nchan_ap_F,f)        
      f.close()     

   return nchan_phase_F, solint_phase_F, solint_ap_F, nchan_ap_F, trigger_antennaconstraint



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

def create_losoto_tecandphaseparset(ms, refant='CS003HBA0'):
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

def create_losoto_tecparset(ms, refant='ST001'):
    parset = 'losoto_plotfasttec.parset'
    os.system('rm -f ' + parset)
    f=open(parset, 'w')
  
    f.write('pol = []\n')
    f.write('Ncpu = 0\n\n\n')

    f.write('[plottecandphase]\n')
    f.write('operation = PLOT\n')
    f.write('soltab = [sol000/tec000]\n')
    f.write('axesInPlot = [time]\n')
    f.write('axisInTable = ant\n')
    f.write('minmax = [-0.2,0.2]\n')
    f.write('figSize=[120,20]\n')
    f.write('prefix = plotlosoto%s/fasttec\n' % ms)
    f.write('refAnt = %s\n' % refant)
  
    f.close()
    return parset

def create_losoto_fastphaseparset(ms, refant='CS003HBA0', onechannel=False, twopol=False):
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
      if twopol:
        f.write('axisInCol = pol\n')
      
    else:
      f.write('axesInPlot = [time,freq]\n')
    f.write('axisInTable = ant\n')
    f.write('minmax = [-3.14,3.14]\n')
    f.write('figSize=[120,20]\n')
    f.write('prefix = plotlosoto%s/fastphase\n' % ms)
    f.write('refAnt = %s\n' % refant)

    if twopol:
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
      f.write('prefix = plotlosoto%s/fastphasepoldiff\n' % ms)
      f.write('refAnt = %s\n' % refant)  
      f.write('axisDiff=pol\n')
 
        
    f.close()
    return parset


def create_losoto_flag_apgridparset(ms, flagging=True, maxrms=7.0, maxrmsphase=7.0, includesphase=True, \
                                    refant='CS003HBA0', onechannel=False, medamp=2.5, flagphases=True, onepol=False):

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
    f.write('prefix = plotlosoto%s/slowamp\n\n\n' % ms)

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
        f.write('prefix = plotlosoto%s/slowphase\n' % ms)
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
        f.write('prefix = plotlosoto%s/slowampfl\n\n\n' % ms)

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
            f.write('prefix = plotlosoto%s/slowphasefl\n' % ms)
            f.write('refAnt = %s\n' % refant)
  
  
    f.close()
    return parset

def create_losoto_mediumsmoothparset(ms, boxsize, longbaseline, includesphase=True, refant='CS003HBA0', onechannel=False):
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
    H5 = h5parm.h5parm(H5name, readonly=False)
    ants = H5.getSolset('sol000').getAnt().keys()
    H5.close()
    if 'ST001' in ants:
        return True
    else:
        return False



def calibrateandapplycal_ov(mslist, selfcalcycle, solint_list, nchan_list, soltype_list, soltypecycles_list, \
              smoothnessconstraint_list,antennaconstraint_list, uvmin=0, normamps=False, skymodel=None, \
              addslowsolve_skymodel=False, predictskywithbeam=False, restoreflags=False, \
              flagging=False, longbaseline=False, BLsmooth=False, flagslowphases=True, flagslowamprms=7.0, flagslowphaserms=7.0):

   incol = 'DATA' # start here, will be updated at applycal step for next solve if needed
   pertubation = False   
   # LOOP OVER THE ENTIRE SOLTYPE LIST (so includes pertubations via a pre-applycal)
   for soltypenumber, soltype in enumerate(soltype_list):
     
     # check we are above far enough in the selfcal to solve for the extra pertubation
     if selfcalcycle >= soltypecycles_list[soltypenumber]: #
       print(selfcalcycle, soltypenumber)
       if (soltypenumber < len(soltype_list)-1):
         if (selfcalcycle >= soltypecycles_list[soltypenumber+1]):  
           pertubation = True # so we are doing a pertubation this round
         else:
           pertubation = False  
       else:
         pertubation = False
       

       # SOLVE LOOP OVER MS
       parmdbmslist = []
       for msnumber, ms in enumerate(mslist):
         if skymodel != None and selfcalcycle == 0:  
           parmdb = soltype + str(soltypenumber) + '_skyselfcalcyle' + str(selfcalcycle) + '_' + ms + '.h5'
         else:
           parmdb = soltype + str(soltypenumber) + '_selfcalcyle' + str(selfcalcycle) + '_' + ms + '.h5'
          
         runDPPPbase(ms, solint_list[soltypenumber][msnumber], nchan_list[soltypenumber][msnumber], parmdb, soltype, \
                     longbaseline=longbaseline, uvmin=uvmin, SMconstraint=smoothnessconstraint_list[soltypenumber][msnumber], \
                     antennaconstraint=antennaconstraint_list[soltypenumber][msnumber], \
                     restoreflags=restoreflags, maxiter=100, flagging=flagging, skymodel=skymodel, \
                     flagslowphases=flagslowphases, flagslowamprms=flagslowamprms, \
                     flagslowphaserms=flagslowphaserms, incol=incol, \
                     predictskywithbeam=predictskywithbeam, BLsmooth=BLsmooth)
         parmdbmslist.append(parmdb)
       
       # NORMALIZE amplitudes
       if normamps and (soltype in ['complexgain','scalarcomplexgain','rotation+diagonal',\
                                    'amplitudeonly','scalaramplitude']):
         print('Doing global gain normalization')
         normamplitudes(parmdbmslist) # list of h5 for different ms, all same soltype

       # APPLYCAL or PRE-APPLYCAL
       for msnumber, ms in enumerate(mslist):
         print(pertubation, parmdbmslist, msnumber)
         if pertubation:
           if soltypenumber == 0:  
             applycal(ms, parmdbmslist[msnumber], msincol='DATA',msoutcol='CORRECTED_PREAPPLY' + str(soltypenumber))
           else:
             applycal(ms, parmdbmslist[msnumber], msincol='CORRECTED_PREAPPLY' + str(soltypenumber-1),\
                      msoutcol='CORRECTED_PREAPPLY' + str(soltypenumber))   
           incol = 'CORRECTED_PREAPPLY' + str(soltypenumber) # SET NEW incol for next solve
         else:
           if soltypenumber == 0:  
             applycal(ms, parmdbmslist[msnumber], msincol='DATA',msoutcol='CORRECTED_DATA')
           else:
             applycal(ms, parmdbmslist[msnumber], msincol='CORRECTED_PREAPPLY' + str(soltypenumber-1),\
                      msoutcol='CORRECTED_DATA')
   return



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
    scriptn = 'python /net/rijn/data2/rvweeren/LoTSS_ClusterCAL/lin2circ.py'
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
      
      logging.info('Noise and clipped noise' + str(parmdb) + ' ' + str(np.std(amps)) + ' ' + str(noise))

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
      logging.info('Trying to changing reference anntena')
    

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
              logging.info('Found new reference anntena,' + str(antenna))
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
    import logging
    logging.basicConfig(filename='selfcal.log', format='%(levelname)s:%(asctime)s ---- %(message)s', datefmt='%m/%d/%Y %I:%M:%S', level=logging.DEBUG)

    
    return total_flux_gaus
#print determine_compactsource_flux('imtry1_0-MFS-image.fits')
#sys.exit()


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

def normamps(parmdb):
    '''
    normalize amplitude solutions to one
    '''
    
    if len(parmdb) == 1:
      H5 = h5parm.h5parm(parmdb[0], readonly=False) 
      amps =H5.getSolset('sol000').getSoltab('amplitude000').getValues()[0]
      weights = H5.getSolset('sol000').getSoltab('amplitude000').getValues(weight=True)[0]
      idx = np.where(weights != 0.0)
    
      amps = np.log10(amps)
      logging.info('Mean amplitudes before normalization: ' + str(10**(np.nanmean(amps[idx]))))
      amps = amps - (np.nanmean(amps[idx]))
      logging.info('Mean amplitudes after normalization: ' + str(10**(np.nanmean(amps[idx]))))
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
          logging.info(parmdbi + '  Normfactor: '+ str(10**(np.nanmean(np.log10(ampsi[idx])))))
          if i == 0:
            amps = np.ndarray.flatten(ampsi[idx])
          else:
            amps = np.concatenate((amps, np.ndarray.flatten(ampsi[idx])),axis=0)

          #print np.shape(amps), parmdbi
          H5.close()
      normmin = (np.nanmean(np.log10(amps))) 
      logging.info('Global normfactor: ' + str(10**normmin))
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
              uvtaper=False, multiscale=True, predict=True, uvmin=' ', fitsmask=None, \
              idg=False, deepmultiscale=False, uvminim=80, fitspectralpol=True, \
              fitspectralpolorder=3, imager='WSCLEAN', restoringbeam=15, automask=2.5, \
              removenegativecc=True, usewgridder=False, paralleldeconvolution=0, \
              deconvolutionchannels=0, parallelgridding=1):
    fitspectrallogpol = False # for testing Perseus

    msliststring = ' '.join(map(str, mslist))
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
      cmd += '-size ' + imsize + ' ' + imsize + ' -reorder '
      cmd += '-weight briggs ' + robust + ' -weighting-rank-filter 3 -clean-border 1 -parallel-reordering 4 '
      cmd += '-mgain 0.8 -fit-beam -data-column ' + imcol +' -join-channels -channels-out '
      cmd += channelsout + ' -padding 1.8 '
      if paralleldeconvolution > 0:
        cmd += '-parallel-deconvolution ' +  str(paralleldeconvolution) + ' '
      if parallelgridding > 1:
        cmd += '-parallel-gridding ' + str(parallelgridding) + ' '  
      if deconvolutionchannels > 0:
        cmd += '-deconvolution-channels ' +  str(deconvolutionchannels) + ' '
      if automask > 0.5:
        cmd += '-auto-mask '+ str(automask)  + ' -auto-threshold 0.5 ' # to avoid automask 0
      
      if multiscale:
         #cmd += '-multiscale '+' -multiscale-scales 0,4,8,16,32,64 '
         cmd += '-multiscale '+' -multiscale-scales 0,6,12,16,24,32,42,64,72,128,180,256,380,512,650 '
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
    
      cmd += '-name ' + imageout + ' -scale ' + pixsize + 'arcsec ' 
      print('WSCLEAN: ', cmd + '-niter ' + str(niter) + ' ' + msliststring)
      logging.info(cmd + '-niter ' + str(niter) + ' ' + msliststring)
      os.system(cmd + '-niter ' + str(niter) + ' ' + msliststring)        
        

      if deepmultiscale:
        
        # predict first to fill MODEL_DATA so we can continue with clean
        cmdp = 'wsclean -size ' 
        cmdp += imsize + ' ' + imsize + ' -channels-out ' + channelsout + ' -padding 1.8 -predict ' 
        if idg:
          cmdp += '-use-idg -grid-with-beam -use-differential-lofar-beam -idg-mode cpu '
          cmdp += '-beam-aterm-update 800 '
          cmdp += '-pol iquv '
        else:
          if usewgridder:    
            cmd +='-use-wgridder '  
            #cmd +='-wgridder-accuracy 1e-4 '
          
        cmdp += '-name ' + imageout + ' -scale ' + pixsize + 'arcsec ' + msliststring
        print('PREDICT STEP for continue: ', cmdp)
        os.system(cmdp)
       
        # NOW continue cleaning
        if not multiscale: # if multiscale is true then this is already set above
          cmd += '-multiscale '+' -multiscale-scales 0,4,8,16,32,64 '
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
        cmd += imsize + ' ' + imsize + ' -channels-out ' + channelsout + ' -padding 1.8 -predict ' 
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

      
        cmd += '-name ' + imageout + ' -scale ' + pixsize + 'arcsec ' + msliststring
        print('PREDICT STEP: ', cmd)
        os.system(cmd)
        
        
    if imager == 'DDFACET':
        makemslist(mslist)
        #restoringbeam = '15'
        cmd = 'DDF.py --Data-MS=mslist.txt --Deconv-PeakFactor=0.001 --Data-ColName=' + imcol + ' ' + \
              '--Parallel-NCPU=32 --Output-Mode=Clean --Deconv-CycleFactor=0 ' + \
              '--Deconv-MaxMinorIter=' + str(niter) + ' --Deconv-MaxMajorIter=5 ' + \
              '--Deconv-Mode=SSD --Weight-Robust=' + robust + ' --Image-NPix=' + imsize + ' ' + \
              '--CF-wmax=50000 --CF-Nw=100 --Beam-Model=None --Beam-LOFARBeamMode=A --Beam-NBand=1 ' + \
              '--Output-Also=onNeds --Image-Cell=' + pixsize + ' --Facets-NFacets=1 --Freq-NDegridBand=1 ' + \
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


 
def runDPPP(ms, solint_ap, solint_phaseonly, nchan_phase, nchan_ap, parmdb, soltype, soltypepert=None, 
             solint_ap_pert = None, nchan_ap_pert=None, \
             preapplyphase=True, preapplyslow=False, longbaseline=False, uvmin=0, TEC=False,\
             SMconstraintslow=0.0, SMconstraintslowpert=0.0, SMconstraintphase=0.0, smoothcal=True, puretec=False, restoreflags=False, soltype_fastphase='scalarphase', maxiter=100, flagging=False, \
             antennaconstraint=None, antennaconstraintslow=None, \
             antennaconstraintslowpert=None, usesourcedb=False, skymodel=None, \
             flagslowphases=True, flagslowamprms=7.0, flagslowphaserms=7.0,incol='DATA', \
             predictskywithbeam=False, BLsmooth=False):
    
    if BLsmooth:
      incol = 'SMOOTHED_DATA'    
    
    # create sourcedb if requested
    if usesourcedb:
      if skymodel.split('.')[-1] != 'sourcedb':
        #make sourcedb
        sourcedb = skymodel + 'sourcedb'
        if os.path.isdir(sourcedb):
          os.system('rm -rf ' + sourcedb)
        cmdmsdb = "makesourcedb in=" + skymodel + " "
        cmdmsdb += "out=" + sourcedb + " outtype='blob' format='<'"
        print(cmdmsdb)
        os.system(cmdmsdb)
      else:
        sourcedb = skymodel


    if soltype_fastphase == 'scalarphase': # for 1D plotting
      twopol=False
    if soltype_fastphase == 'phaseonly': # for 1D plotting
      twopol=True

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

    # determine if phases needs to be updated, important if slowgains do not contain phase solutions    
    includesphase = True
    if soltype == 'scalaramplitude' or soltype == 'amplitudeonly':
      includesphase = False



    losoto = 'losoto'
    # figure out which weight_spectrum column to use
    t = pt.table(ms)
    if 'WEIGHT_SPECTRUM_SOLVE' in t.colnames():
       weight_spectrum =  'WEIGHT_SPECTRUM_SOLVE'
    else:
       weight_spectrum =  'WEIGHT_SPECTRUM'
    t.close()   
    
    
    losotoparset_tecandphase = create_losoto_tecandphaseparset(ms)
    losotoparset_tec = create_losoto_tecparset(ms)



    # check for previous old parmdb and remove them   
    if os.path.isfile(parmdb):
      print('H5 file exists  ', parmdb)
      os.system('rm -f ' + parmdb)
    if os.path.isfile('phaseonly' + parmdb): # and preapplyphase == True:
      print('H5 file exists  ', 'phaseonly' + parmdb)
      os.system('rm -f ' + parmdb)
    
    if (soltype == 'scalarphase' or soltype == 'phaseonly') and preapplyphase == True:
        print('Not supported')
        sys.exit(1)
    
 
    cmd = 'DPPP numthreads='+str(multiprocessing.cpu_count())+ ' msin=' + ms + ' msin.datacolumn=' + incol + ' '
    cmd += 'msout=. msin.modelcolumn=MODEL_DATA '
    cmd += 'msin.weightcolumn='+weight_spectrum + ' '
    cmd += 'steps=[ddecal] ' + 'msout.storagemanager=dysco ddecal.type=ddecal '
    cmd += 'ddecal.maxiter='+str(np.int(maxiter)) + ' ddecal.propagatesolutions=True '
    if usesourcedb:
      cmd += 'ddecal.sourcedb=' + sourcedb + ' '
      if predictskywithbeam:
        cmd += 'ddecal.usebeammodel=True ddecal.usechannelfreq=True ddecal.beammode=array_factor '   
    else:
      cmd += 'ddecal.usemodelcolumn=True '
    
    if uvmin != 0:
        cmd += 'ddecal.uvlambdamin=' + str(uvmin) + ' '      
   
    # CASE I   
    if (soltype == 'complexgain' or soltype == 'rotation+diagonal' or soltype == 'scalaramplitude' or \
        soltype == 'scalarcomplexgain' or soltype == 'amplitudeonly' or soltype == 'fulljones') \
        and preapplyphase == True:
        
        if antennaconstraint != None:
           cmd += 'ddecal.antennaconstraint=' + antennaconstraintstr(antennaconstraint, antennasms, HBAorLBA) + ' '
        
        if TEC == False:
           cmd += 'ddecal.mode=' + soltype_fastphase + ' ' # scalarphase, phaseonly, complexgain, tec, tecandphase
           cmd += 'ddecal.tolerance=1.e-4 '
           if SMconstraintphase > 0.0:
             cmd += 'ddecal.smoothnessconstraint=' + str(SMconstraintphase*1e6) + ' '      
        if TEC == True:
           if puretec: 
             cmd += 'ddecal.mode=tec ' 
           else:
             cmd += 'ddecal.mode=tecandphase '                
           cmd += 'ddecal.approximatetec=True '
           cmd += 'ddecal.stepsize=0.2 '
           cmd += 'ddecal.maxapproxiter=45 '
           cmd += 'ddecal.tolerance=1e-4 '
           cmd += 'ddecal.approxtolerance=6e-3 '
        
        cmd += 'ddecal.solint=' + str(solint_phaseonly) + ' '
        cmd += 'ddecal.nchan=' + str(nchan_phase) + ' '
        cmd += 'ddecal.h5parm=phaseonly' + parmdb + ' '

    # CASE II (rare)  
    if (soltype == 'complexgain' or soltype == 'rotation+diagonal' or soltype == 'scalaramplitude' \
        or soltype == 'scalarcomplexgain' or soltype == 'amplitudeonly' or soltype == 'fulljones') \
        and preapplyphase == False:
        
        cmd += 'ddecal.mode='+ soltype + ' ' # scalarphase, phaseonly
        cmd += 'ddecal.solint=' + str(solint_ap) + ' '
        cmd += 'ddecal.nchan=' + str(nchan_ap) + ' '
        cmd += 'ddecal.h5parm=' + parmdb + ' '
        cmd += 'ddecal.tolerance=1e-4 '
        if SMconstraintslow > 0.0:
          cmd += 'ddecal.smoothnessconstraint=' + str(SMconstraintslow*1e6) + ' ' 
        if antennaconstraintslow != None:
           cmd += 'ddecal.antennaconstraint=' + antennaconstraintstr(antennaconstraintslow, antennasms, HBAorLBA) + ' '
 
    # CASE III      
    if (soltype == 'scalarphase' or soltype == 'phaseonly') and preapplyphase == False:   
        if antennaconstraint != None:
           cmd += 'ddecal.antennaconstraint=' + antennaconstraintstr(antennaconstraint, antennasms, HBAorLBA) + ' '
        
        if TEC == False:
           if np.int(maxiter) == 1: # this is a template solve only for pol version
             cmd += 'ddecal.mode=phaseonly ' # to get XX and YY pol
           else:
             cmd += 'ddecal.mode=' + soltype + ' ' # scalarphase, phaseonly
           cmd += 'ddecal.tolerance=1e-4 '
           if SMconstraintphase > 0.0:
            cmd += 'ddecal.smoothnessconstraint=' + str(SMconstraintphase*1e6) + ' '            
        if TEC == True:
           if puretec:
             cmd += 'ddecal.mode=tec '  
           else:    
             cmd += 'ddecal.mode=tecandphase '
           cmd += 'ddecal.approximatetec=True '
           cmd += 'ddecal.stepsize=0.2 '
           cmd += 'ddecal.maxapproxiter=45 '
           cmd += 'ddecal.tolerance=1e-4 '
           cmd += 'ddecal.approxtolerance=6e-3 '
        
        cmd += 'ddecal.solint=' + str(solint_phaseonly) + ' '
        cmd += 'ddecal.nchan=' + str(nchan_phase) + ' '
        cmd += 'ddecal.h5parm=phaseonly' + parmdb + ' ' 
        
 
    print('DPPP solve:', cmd)
    os.system(cmd)
    if np.int(maxiter) == 1: # this is a template solve only
      return    

    if preapplyphase: # APPLY FIRST 
        cmd = 'DPPP numthreads=' + str(multiprocessing.cpu_count()) + ' msin=' + ms + ' msin.datacolumn=DATA msout=. '
        cmd += 'msin.weightcolumn='+weight_spectrum + ' msout.storagemanager=dysco '
        if TEC == False:
          cmd += 'msout.datacolumn=CORRECTED_DATA_PHASE steps=[ac1] '
          cmd += 'ac1.parmdb=phaseonly'+parmdb + ' ac1.type=applycal '
          cmd += 'ac1.correction=phase000 '
          print('DPPP PRE-APPLY PHASE-ONLY:', cmd)
          os.system(cmd)
          if number_freqchan_h5('phaseonly' + parmdb) > 1:
            if check_phaseup('phaseonly' + parmdb):
              losotoparset_phase = create_losoto_fastphaseparset(ms, refant='ST001', twopol=twopol)
            else:
              losotoparset_phase = create_losoto_fastphaseparset(ms, twopol=twopol)
          else: 
            losotoparset_phase = create_losoto_fastphaseparset(ms, onechannel=True, twopol=twopol)  
          cmdlosotophase = losoto + ' '
          cmdlosotophase += 'phaseonly' + parmdb + ' ' + losotoparset_phase
          os.system(cmdlosotophase)
        if TEC == True:
          if puretec:  
            cmd += 'msout.datacolumn=CORRECTED_DATA_PHASE steps=[ac2] '
          else:
            cmd += 'msout.datacolumn=CORRECTED_DATA_PHASE steps=[ac1,ac2] '
          cmd += 'ac1.parmdb=phaseonly'+parmdb + ' ac1.type=applycal '
          cmd += 'ac1.correction=phase000 '
          cmd += 'ac2.parmdb=phaseonly'+parmdb + ' ac2.type=applycal '
          cmd += 'ac2.correction=tec000 '
          print('DPPP PRE-APPLY TECANDPHASE:', cmd)
          os.system(cmd)
          cmdlosotophase = losoto + ' '
          if puretec:            
            cmdlosotophase += 'phaseonly' + parmdb + ' ' + losotoparset_tec
          else:    
            cmdlosotophase += 'phaseonly' + parmdb + ' ' + losotoparset_tecandphase
            tecandphaseplotter('phaseonly' + parmdb, ms) # own plotter because losoto cannot add tec and phase
          os.system(cmdlosotophase)

        # RUN DPPP again
        cmd = 'DPPP numthreads=' + str(multiprocessing.cpu_count()) + ' msin=' + ms + ' msin.datacolumn=CORRECTED_DATA_PHASE msout=. '
        cmd += 'msin.weightcolumn='+weight_spectrum + ' '
        cmd += 'msin.modelcolumn=MODEL_DATA '
        cmd += 'steps=[ddecal] ' + 'msout.storagemanager=dysco ddecal.type=ddecal '
        cmd += 'ddecal.maxiter='+str(np.int(maxiter)) +' ddecal.tolerance=1.e-4 ddecal.propagatesolutions=True '

        if usesourcedb:
          cmd += 'ddecal.sourcedb=' + sourcedb + ' '
          if predictskywithbeam:
            cmd += 'ddecal.usebeammodel=True ddecal.usechannelfreq=True ddecal.beammode=array_factor '  
        else:
          cmd += 'ddecal.usemodelcolumn=True '

        cmd += 'ddecal.nchan=' + str(nchan_ap) + ' '
        cmd += 'ddecal.mode=' + soltype + ' ' # complexgain
        cmd += 'ddecal.h5parm=' + parmdb + ' ' 
        cmd += 'ddecal.solint=' + str(solint_ap) + ' ' 
        if SMconstraintslow > 0.0:
          cmd += 'ddecal.smoothnessconstraint=' + str(SMconstraintslow*1e6) + ' ' 
        if antennaconstraintslow != None:
           cmd += 'ddecal.antennaconstraint=' + antennaconstraintstr(antennaconstraintslow, antennasms, HBAorLBA) + ' '

        if uvmin != 0:
           cmd += 'ddecal.uvlambdamin=' + str(uvmin) + ' '
        print('DPPP SLOW GAIN solve:', cmd)
        os.system(cmd)
        os.system('cp ' + parmdb + ' ' + parmdb + '.backup')

        #  ---- PROCESS SOLUTIONS, FILTERING AND PLOTTING ----    
        if number_freqchan_h5(parmdb) > 1:
          onechannel = False
        else:
          onechannel = True  
        
        medamp = medianamp(parmdb)
        
        if soltype =='scalarcomplexgain' or soltype == 'scalaramplitude':
          onepolslow = True
        else:
          onepolslow = False    
                
        if flagging and (ntimesH5(parmdb) > 1):
          losotoparset = create_losoto_flag_apgridparset(ms, flagging=True, maxrms=flagslowamprms, maxrmsphase=flagslowphaserms, includesphase=includesphase, onechannel=onechannel, \
                                                         medamp=medamp, flagphases=flagslowphases, onepol=onepolslow)
        else:
          losotoparset = create_losoto_flag_apgridparset(ms, flagging=False, includesphase=includesphase, onechannel=onechannel, medamp=medamp, onepol=onepolslow)  

        # MAKE losoto command    
        cmdlosoto = losoto + ' ' + parmdb + ' ' + losotoparset


      
    #  CASE I FLAGGING (rare) 
    if (soltype == 'complexgain' or soltype == 'rotation+diagonal' or soltype == 'scalaramplitude' \
        or soltype == 'scalarcomplexgain' or soltype == 'amplitudeonly' or soltype == 'fulljones') \
        and (preapplyphase == False):
      flagbadamps(parmdb, setweightsphases=includesphase)
      medamp = medianamp(parmdb)
      flaglowamps(parmdb, lowampval=medamp*0.1, flagging=flagging, setweightsphases=includesphase)
      flaghighgamps(parmdb, highampval=medamp*10., flagging=flagging, setweightsphases=includesphase)
      if soltype != 'amplitudeonly' and soltype != 'scalaramplitude':
        try:
          change_refant(parmdb,'phase000')
        except:
          pass    
        removenans(parmdb, 'phase000')
      removenans(parmdb, 'amplitude000')
      if rotation:
        removenans(parmdb, 'rotation000')
      
      if ntimesH5(parmdb) > 1: # plotting does not work if time axis has length 1
      # FLAG/SMOOTH solutions
        os.system(cmdlosoto)
      
      if smoothcal:
        smoothsols(parmdb, ms, longbaseline, includesphase=includesphase)
      

    #  CASE II FLAGGING
    if (soltype == 'scalarphase' or soltype == 'phaseonly') and TEC == False:
      removenans('phaseonly' + parmdb, 'phase000')
      if number_freqchan_h5('phaseonly' + parmdb) > 1:
        if check_phaseup('phaseonly' + parmdb):  
          losotoparset_phase = create_losoto_fastphaseparset(ms, refant='ST001', twopol=twopol)
        else:
          losotoparset_phase = create_losoto_fastphaseparset(ms, twopol=twopol)
      else: 
        losotoparset_phase = create_losoto_fastphaseparset(ms, onechannel=True, twopol=twopol)  
      cmdlosotophase = losoto + ' '
      cmdlosotophase += 'phaseonly' + parmdb + ' ' + losotoparset_phase
      os.system(cmdlosotophase)
    
    # CASE III FLAGGING 
    if soltype == 'scalarphase' and TEC == True:
      cmdlosotophase = losoto + ' '
      if puretec:
        cmdlosotophase += 'phaseonly' + parmdb + ' ' + losotoparset_tec
      else:    
        cmdlosotophase += 'phaseonly' + parmdb + ' ' + losotoparset_tecandphase
        tecandphaseplotter('phaseonly' + parmdb, ms)
      os.system(cmdlosotophase)

    #  CASE IV FLAGGING      
    if (soltype == 'complexgain' or soltype == 'rotation+diagonal' or soltype == 'scalaramplitude' or \
        soltype == 'scalarcomplexgain' or soltype == 'amplitudeonly' or soltype == 'fulljones') \
        and (preapplyphase == True):
      
      flagbadamps(parmdb, setweightsphases=includesphase)
      medamp = medianamp(parmdb)
      if longbaseline:
        if soltype != 'fulljones':  
          flaglowamps(parmdb,lowampval=medamp*0.1,flagging=flagging, setweightsphases=includesphase)
          flaghighgamps(parmdb, highampval=medamp*10.0,flagging=flagging, setweightsphases=includesphase)   
      else:
        if soltype != 'fulljones':
          flaglowamps(parmdb,lowampval=medamp*0.1,flagging=flagging, setweightsphases=includesphase)
          flaghighgamps(parmdb, highampval=medamp*10.0,flagging=flagging, setweightsphases=includesphase)   
      
      # FLAG/SMOOTH solutions
      removenans(parmdb, 'amplitude000')
      if soltype != 'amplitudeonly' and soltype != 'scalaramplitude':
        removenans(parmdb, 'phase000')      
        try:
          change_refant(parmdb,'phase000')
        except:
          pass  
      os.system(cmdlosoto)
      
      if smoothcal:
         smoothsols(parmdb, ms, longbaseline, includesphase=includesphase)

    if preapplyslow:
      cmd = 'DPPP numthreads='+str(multiprocessing.cpu_count()) + ' msin=' + ms + ' msin.datacolumn=CORRECTED_DATA_PHASE msout=. '
      cmd += 'msin.weightcolumn='+weight_spectrum + ' msout.storagemanager=dysco '
      cmd += 'msout.datacolumn=CORRECTED_DATA_SLOW '
      cmd += 'ac1.parmdb='+parmdb + ' ac2.parmdb='+parmdb + ' ac3.parmdb='+parmdb + ' '
      cmd += 'ac1.type=applycal ac2.type=applycal ac3.type=applycal'
      if soltype == 'complexgain' or soltype == 'scalarcomplexgain':    
        cmd += 'steps=[ac1, ac2] '
        cmd += 'ac1.correction=phase000 ac2.correction=amplitude000 '
      if soltype == 'scalaramplitude' or soltype == 'amplitudeonly': 
        cmd += 'steps=[ac2] '
        cmd += 'ac2.correction=amplitude000 ' 
      if soltype == 'fulljones':  
        cmd += 'steps=[ac2] '
        cmd += 'ac2.correction=fulljones '
        cmd += 'ac2.soltab=[phase000,amplitude000] '
      if soltype == 'rotation+diagonal':  
        cmd += 'steps=[ac1, ac2, ac3] '
        cmd += 'ac1.correction=phase000 ac2.correction=amplitude000 ac3.correction=rotation000'

      print('DPPP preapplyslow', cmd)
      os.system(cmd)

      # RUN DPPP again
      cmd = 'DPPP numthreads='+str(multiprocessing.cpu_count()) + ' msin=' + ms + ' msin.datacolumn=CORRECTED_DATA_AMP msout=. '
      cmd += 'msin.weightcolumn='+weight_spectrum + ' '
      cmd += 'msin.modelcolumn=MODEL_DATA '
      cmd += 'steps=[ddecal] ' + 'msout.storagemanager=dysco ddecal.type=ddecal '
      cmd += 'ddecal.maxiter='+str(np.int(maxiter)) +' ddecal.tolerance=1.e-4 ddecal.propagatesolutions=True '

      if usesourcedb:
        cmd += 'ddecal.sourcedb=' + sourcedb + ' '  
        if predictskywithbeam:
          cmd += 'ddecal.usebeammodel=True ddecal.usechannelfreq=True ddecal.beammode=array_factor '      
      else:
        cmd += 'ddecal.usemodelcolumn=True '

      cmd += 'ddecal.nchan=' + str(nchan_ap) + ' '
      cmd += 'ddecal.mode=' + soltype + ' ' # complexgain
      cmd += 'ddecal.h5parm=pertubation' + parmdb + ' ' 
      cmd += 'ddecal.solint=' + str(solint_ap) + ' ' 
      if SMconstraintslow > 0.0:
        cmd += 'ddecal.smoothnessconstraint=' + str(SMconstraintslow*1e6) + ' ' 
      if antennaconstraintslow != None:
         cmd += 'ddecal.antennaconstraint=' + antennaconstraintstr(antennaconstraintslow, antennasms, HBAorLBA) + ' '

      if uvmin != 0:
         cmd += 'ddecal.uvlambdamin=' + str(uvmin) + ' '
      print('DPPP SLOW GAIN SECOND PERTUBATION solve:', cmd)
      os.system(cmd)
      os.system('cp pertubation' + parmdb + ' pertubation' + parmdb + '.backup')





# to remove H5/h5 and other files out of a wildcard selection if needed
def removenonms(mslist):
  newmslist = []  
  for ms in mslist:      
   if ms.lower().endswith(('.h5', '.png', '.parset', '.fits', '.backup', '.obj', '.log', '.p', '.reg', '.gz', '.tar')) or \
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
    
    if len(times) > 99:
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

parser = argparse.ArgumentParser(description='Self-Calibrate a facet from a LOFAR observation')

parser.add_argument('-b','--boxfile', help='boxfile', type=str)
parser.add_argument('--imsize', help='image size, required if boxfile is not used', type=int)
parser.add_argument('--pixelscale','--pixelsize', help='pixels size in arcsec, default=3.0/1.5 (LBA/HBA)', type=float)
parser.add_argument('-i','--imagename', help='imagename, default=image', default='image', type=str)
parser.add_argument('--fitsmask', help='fitsmask for deconvolution, if not provided use automasking', type=str)
parser.add_argument('--H5sols', help='prefix name for H5 solution file, default=solsgrid', default='solsgrid', type=str)
parser.add_argument('-n', '--niter', help='niter, default=15000', default=15000, type=int)
parser.add_argument('--robust', help='Briggs robust paramter, default=-0.5', default=-0.5, type=float)
parser.add_argument('--channelsout', help='channelsout, default=6', default=6, type=int)
parser.add_argument('--multiscale', help='use multiscale deconvolution, not recommended/unstable', action='store_true')
parser.add_argument('--deepmultiscale', help='do extra multiscale deconvolution on the residual', action='store_true')
parser.add_argument('--uvminim', help='inner uv-cut for imaging in lambda, default=80', default=80., type=float)
parser.add_argument('--usewgridder', help='use wgridder in WSClean, mainly useful for very large images', action='store_true')
parser.add_argument('--phaseup', help='phase up to a superstation', action='store_true')
parser.add_argument('--phaseupstations', help='phase up to a superstation (core or superterp, default core)', default='core', type=str)

parser.add_argument('--paralleldeconvolution', help='parallel-deconvolution size for wsclean, default=0 (means no parallel deconvolution)', default=0, type=int)
parser.add_argument('--parallelgridding', help='parallel-gridding for wsclean, default=1 (means no parallel gridding)', default=1, type=int)
parser.add_argument('--deconvolutionchannels', help='deconvolution-channels value for wsclean, default=0 (means deconvolution-channels equals channels-out)', default=0, type=int)


parser.add_argument('--idg', help='use the Image Domain gridder', action='store_true')
parser.add_argument('--maskthreshold', help='threshold for MakeMask.py, default=5', default=5.0, type=float)
parser.add_argument('--imager', help='Imager to use WSClean or DDFACET, default WSCLEAN', default='WSCLEAN', type=str)
parser.add_argument('--fitspectralpol', help='use fit-spectral-pol in WSClean (True/False, default=True)', type=eval, choices=[True, False], default='True')
parser.add_argument('--fitspectralpolorder', help='fit-spectral-pol order for WSClean, default=3', default=3, type=int)
parser.add_argument('--removenegativefrommodel', help='remove negative clean components in model predict (True/False, default=True)', type=eval, choices=[True, False], default='True')
parser.add_argument('--autofrequencyaverage', help='by default frequency averaging is tried if it does not result in bandwidth  smearing (True/False, default=False)', type=eval, choices=[True, False], default='False')
parser.add_argument('--avgfreqstep', help='extra DPPP frequnecy averaging to speed up a solve, this is done before any other correction, could be useful for long baseline infield calibrators', type=int, default=None)
parser.add_argument('--avgtimestep', help='extra DPPP time averaging to speed up a solve, this is done before any other correction, could be useful for long baseline infield calibrators', type=int, default=None)

# calibration options
parser.add_argument('-u', '--uvmin', help='inner uv-cut for calibration in lambda, default=80/350 (LBA/HBA)', type=float)
parser.add_argument('--no-tec', help='do not use TEC fitting', action='store_false')
parser.add_argument('--pure-tec', help='use TEC, mode="tec", instead of mode="tecandphase"', action='store_true')
parser.add_argument('--phase-soltype', help='phase solve type (scalarphase/phaseonly), only relevant if --no-tec is used, default=scalarphase', default='scalarphase', type=str)
parser.add_argument('--slow-soltype', help='slow solve type (complexgain/scalarcomplexgain/scalaramplitude/amplitudeonly/fulljones), default=complexgain', default='scalarcomplexgain', type=str)
#parser.add_argument('--slow-soltype2', help='slow solve type (complexgain/scalarcomplexgain/scalaramplitude/amplitudeonly/fulljones), default=complexgain', default='complexgain', type=str)
parser.add_argument('--solve-rotation', help='slow solve for rotation, mode="rotation+diagonal"', action='store_true')
#parser.add_argument('--solve-rotation2', help='slow solve for rotation, mode="rotation+diagonal"', action='store_true')
parser.add_argument("--soltype-list", type=arg_as_list, default=['tecandphase','scalarcomplexgain'],help="List of values")
parser.add_argument("--solint-list", type=arg_as_list, default=[1,120],help="List of values")
parser.add_argument("--nchan-list", type=arg_as_list, default=[1,10],help="List of values")
parser.add_argument("--smoothnessconstraint-list", type=arg_as_list, default=[0,0],help="List of values")
parser.add_argument("--antennaconstraint-list", type=arg_as_list, default=[None,None],help="List of values")
parser.add_argument("--soltypecycles-list", type=arg_as_list, default=[0,3],help="List of values")
parser.add_argument("--BLsmooth", help='Employ BLsmooth for low S/N data', action='store_true')

parser.add_argument('--phasecycles', help='number of phase-only/tec selfcal cycles before ap solve, default=3', default=3, type=int)
parser.add_argument('--nchan-phase', help='DPPP nchan for fast phase/tec solving', type=int, default=None) #default was 5
parser.add_argument('--nchan-ap', help='DPPP nchan for slow gain solving', type=int, default=None)
#parser.add_argument('--nchan-ap2', help='DPPP nchan for slow gain solving', type=int)
parser.add_argument('--solint-phase', help='DPPP solint for fast phase/tec solving', type=int, default=None) # default was 1
parser.add_argument('--solint-ap', help='DPPP solint, for slow gain solving',type=int,default=None)
#parser.add_argument('--solint-ap2', help='DPPP solint, for slow gain solving',type=int)
parser.add_argument('--usemodeldataforsolints', help='determine solints from MODEL_DATA', action='store_true')
#parser.add_argument('--lb-solintphase-sc', help='multiply solint found for tec/phase from the option --usemodeldataforsolints with this value, default=1.0', default=1.0, type=float)
parser.add_argument('--smoothnessconstraint-slow', help='Kernel size in MHz, default=0. When unequal to 0, will constrain the slow solutions to be smooth over frequency', default=0.0, type=float)
parser.add_argument('--smoothnessconstraint-slow2', help='Kernel size in MHz, default=0. When unequal to 0, will constrain the slow solutions to be smooth over frequency', default=0.0, type=float)
parser.add_argument('--smoothnessconstraint-phase', help='Kernel size in MHz, default=0. When unequal to 0, will constrain the phase solutions (when using --no-tec) to be smooth over frequency', default=0.0, type=float)
parser.add_argument('--antennaconstraint-slow', help='constrain these stations to have the same solutions in the slow solve (options are core, superterp, coreandfirstremotes, remote, alldutch, international, core-remote, and all; default None)', type=str)
parser.add_argument('--antennaconstraint-slow2', help='constrain these stations to have the same solutions in the slow solve (options are core, superterp, coreandfirstremotes, remote, alldutch, international, core-remote, and all; default None)', type=str)
parser.add_argument('--antennaconstraint-phase', help='constrain these stations to have the same solutions (options are core, superterp, coreandfirstremotes, remote, alldutch, international, core-remote, and all; default None)', type=str)

parser.add_argument('--dosecondslow', help='Do a second slow solve after the first one (True/False, default=False, not working at the moment)', type=eval, choices=[True, False], default='False')


# general options
#parser.add_argument('--longbaseline', help='optimisze settings for long baselines', action='store_true')
parser.add_argument('--skymodel', help='skymodel for first phase-only calibration', type=str)
parser.add_argument('--predictskywithbeam', help='predict the skymodel with the beam array factor', action='store_true')
parser.add_argument('--addslowsolve-skymodel', help='add slow solve for the calibration against the skymodel (default=False)', type=eval, choices=[True, False], default='False')
parser.add_argument('--startfromtgss', help='start from TGSS skymodel for positions (boxfile required)', action='store_true')
parser.add_argument('--tgssfitsimage', help='start TGSS fits image for model (if not provided use SkyView', type=str)
parser.add_argument('--no-beamcor', help='do not correct the visilbities for the array factor', action='store_false')
parser.add_argument('--use-dpppbeamcor', help='use DP3 for beam correction, requires recent DP3 version and no phased-up stations', action='store_true')
parser.add_argument('--docircular', help='convert linear to circular correlations', action='store_true')
parser.add_argument('--dolinear', help='convert circular to linear correlations', action='store_true')
parser.add_argument('--forwidefield', help='Keep solutions such that they can be used for widefield imaging/screens', action='store_true')
parser.add_argument('--doflagging', help='flag on complexgain solutions (True/False, default=True)', type=eval, choices=[True, False], default='True')
parser.add_argument('--restoreflags', help='Restore flagging column after each selfcal cycle, only relevant if --doflagging=True', action='store_true')
parser.add_argument('--flagslowamprms', help='RMS outlier value to flag on slow amplitudes (default=7.0)', default=7.0, type=float)
parser.add_argument('--flagslowphaserms', help='RMS outlier value to flag on slow phases (default=7.0)', default=7.0, type=float)
parser.add_argument('--doflagslowphases', help='If solution flagging is done also flag outliers phases in the slow phase solutions (True/False, default=True)', type=eval, choices=[True, False], default='True')
parser.add_argument('--useaoflagger', help='run AOflagger on input data', action='store_true')
parser.add_argument('--useaoflaggerbeforeavg', help='run before (True) or after averaging (False), default=True', type=eval, choices=[True, False], default='True')
parser.add_argument('--normamps', help='Normalize global amplitudes to 1.0 (True/False, default=True, turned off if fulljones is used)', type=eval, choices=[True, False], default='True')
parser.add_argument('--resetweights', help='if you want to ignore weight_spectrum_solve', action='store_true')
parser.add_argument('--start', help='start selfcal cycle at this iteration, default=0', default=0, type=int)
parser.add_argument('--stop', help='stop selfcal cycle at this iteration, default=10', default=10, type=int)
parser.add_argument('--smoothcal', help='median smooth amplitudes (not recommended)', action='store_true')


parser.add_argument('ms', nargs='+', help='msfile(s)')  

#argso1,argso2 = parser.parse_args()

args = vars(parser.parse_args())

#print(args['soltype_list'])
#print(args['solint_list'])

#sys.exit()

print_title()
inputchecker(args)
trigger_antennaconstraint = None # define here, will be overwriten if needed by usemodeldataforsolints HBA LoTSS-type

if not os.path.isfile('lib_multiproc.py'):
    os.system('cp /net/rijn/data2/rvweeren/LoTSS_ClusterCAL/lib_multiproc.py .')


if args['forwidefield']:
   args['doflagging'] = False
   args['smoothcal']  = False



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
      args['uvmin'] = 80.
  else:
      args['uvmin'] = 350.

if args['pixelscale'] == None:  
  if LBA:
    if longbaseline:
      args['pixelscale'] = 0.08 
    else:
      args['pixelscale'] = 3.0  
  else:
    if longbaseline:
      args['pixelscale'] = 0.03
    else:
      args['pixelscale'] = 1.5

if args['boxfile'] != None:
  imsize   = str(getimsize(args['boxfile'], args['pixelscale']))
if args['imsize'] != None:
  imsize = str(args['imsize']) 
# check if we could average more
avgfreqstep = findfreqavg(mslist[0],np.float(imsize))
if (avgfreqstep > 0) and not longbaseline and args['avgfreqstep'] == None \
                     and args['autofrequencyaverage']: 
  args['avgfreqstep'] = avgfreqstep


# average if requested
mslist = average(mslist, freqstep=args['avgfreqstep'], timestep=args['avgtimestep'], start=args['start'])


# extra flagging if requested
if args['start'] == 0 and args['useaoflagger'] and not args['useaoflaggerbeforeavg']:
  runaoflagger(mslist) 


t    = pt.table(mslist[0] + '/SPECTRAL_WINDOW',ack=False)
bwsmear = bandwidthsmearing(np.median(t.getcol('CHAN_WIDTH')), np.min(t.getcol('CHAN_FREQ')[0]), np.float(imsize))
t.close() # close pt.table(ms + '/SPECTRAL_WINDOW',ack=False) here



# backup flagging column for option --restoreflags if needed
if args['restoreflags']:
  for ms in mslist:
    create_backup_flag_col(ms)
    


automask = 2.5
if args['maskthreshold'] < automask:
  automask = args['maskthreshold'] # in case we use a low value for maskthreshold, like Herc A    



TEC = args['no_tec']
if TEC:
    args['phase_soltype'] = 'scalarphase' # because phaseonly is not relevant/not used
idg = args['idg']
multiscale = args['multiscale']
imageout  = args['imagename'] + '_'
dobeamcor = args['no_beamcor']
if args['fitsmask'] != None:
  fitsmask = args['fitsmask']
else:
  fitsmask = None
uvmin = args['uvmin']  
robust = str(args['robust'])
channelsout = str(args['channelsout'])  
niter = args['niter']
parmdb = args['H5sols']  + '_'
pixsize = str(args['pixelscale'])  


if args['boxfile'] != None:
  outtarname = (args['boxfile'].split('/')[-1]).split('.reg')[0] + '.tar.gz'
else:
  outtarname = 'calibrateddata' + '.tar.gz' 

logging.info('Imsize:                    ' + str(imsize))
logging.info('Pixelscale:                ' + str(pixsize))
logging.info('Niter:                     ' + str(niter))
logging.info('Uvmin:                     ' + str(uvmin))
logging.info('Multiscale:                ' + str(multiscale))
logging.info('Beam correction:           ' + str(dobeamcor))
logging.info('IDG:                       ' + str(idg))
logging.info('TEC:                       ' + str(TEC))
logging.info('Widefield:                 ' + str(args['forwidefield']))
logging.info('Flagslowamprms:            ' + str(args['flagslowamprms']))
logging.info('flagslowphaserms:          ' + str(args['flagslowphaserms']))
logging.info('Do linear:                 ' + str(args['dolinear']))
logging.info('Do circular:               ' + str(args['docircular']))
logging.info('smoothnessconstraint_phase:' + str(args['smoothnessconstraint_phase']))
logging.info('smoothnessconstraint_slow: ' + str(args['smoothnessconstraint_slow']))
logging.info('Pure TEC:                  ' + str(args['pure_tec']))
logging.info('phase soltype              ' + str(args['phase_soltype']))
logging.info('slow soltype               ' + str(args['slow_soltype']))
if args['boxfile'] != None:
  logging.info('Bobxfile:                  ' + args['boxfile'])
logging.info('Mslist:                    ' + ' '.join(map(str,mslist)))
logging.info('User specified clean mask: ' + str(fitsmask))
logging.info('Threshold for MakeMask:    ' + str(args['maskthreshold']))
logging.info('Briggs robust:             ' + str(robust))
logging.info('Imagename prefix:          ' + imageout)
logging.info('Solution file prefix:      ' + parmdb)
logging.info('Output file will be:       ' + outtarname)




if dobeamcor and idg:
  print('beamcor=True and IDG=True is not possible')
  sys.exit(1)

if args['startfromtgss'] and args['start'] == 0:
  if args['boxfile'] != None and args['skymodel'] == None:
    args['skymodel'] = makeBBSmodelforTGSS(args['boxfile'],fitsimage = args['tgssfitsimage'], \
                                           pixelscale=args['pixelscale'], imsize=args['imsize'])
  else:
    print('You need to provide a boxfile to use --startfromtgss')
    print('And you cannot provide a skymodel file manually')
    sys.exit(1)


if args['start'] == 0:
  os.system('rm -f nchan_phase.p') 
  os.system('rm -f nchan_ap.p') 
  os.system('rm -f solint_phase.p') 
  os.system('rm -f solint_ap.p') 



nchan_phase,solint_phase,solint_ap,nchan_ap = setinitial_solint(mslist, \
                                              TEC, longbaseline, LBA,\
                                              innchan_phase=args['nchan_phase'], \
                                              innchan_ap=args['nchan_ap'], \
                                              insolint_phase=args['solint_phase'],\
                                              insolint_ap=args['solint_ap'])

# Get restoring beam for DDFACET in case it is needed
restoringbeam = calculate_restoringbeam(mslist, LBA)



# ----- START SELFCAL LOOP -----
for i in range(args['start'],args['stop']):

  # AUTOMATICALLY PICKUP PREVIOUS MASK (in case of a restart)
  if (i > 0) and (args['fitsmask'] == None):
    if idg:  
      if os.path.isfile(imageout + str(i-1) + '-MFS-I-image.fits.mask.fits'):
          fitsmask = imageout + str(i-1) + '-MFS-I-image.fits.mask.fits'
    else:
      if args['imager'] == 'WSCLEAN':
        if os.path.isfile(imageout + str(i-1) + '-MFS-image.fits.mask.fits'):
            fitsmask = imageout + str(i-1) + '-MFS-image.fits.mask.fits'
      if args['imager'] == 'DDFACET':
        if os.path.isfile(imageout + str(i-1) + '.app.restored.fits'):
            fitsmask = imageout + str(i-1) + '.app.restored.fits.mask.fits'

       
  # BEAM CORRECTION
  if dobeamcor and i == 0:
      for ms in mslist:
        beamcor(ms, usedppp=args['use_dpppbeamcor'])
        #sys.exit()
  # CONVERT TO CIRCULAR/LINEAR CORRELATIONS      
  if (args['docircular'] or args['dolinear']) and i == 0:
      for ms in mslist:
        circular(ms, linear=args['dolinear'])

  # PHASE-UP if requested
  if args['phaseup'] and i== 0:
      mslist = phaseup(mslist,datacolumn='DATA',superstation=args['phaseupstations'])

  # RUN BLsmooth if requested
  if args['BLsmooth'] and i == 0:
      for ms in mslist:    
        os.system('python /net/rijn/data2/rvweeren/LoTSS_ClusterCAL/BLsmooth.py -n 8 -i DATA -o SMOOTHED_DATA ' + ms)



  if (args['skymodel'] != None) and (i ==0):
    for msnumber, ms in enumerate(mslist):   
      if args['addslowsolve_skymodel']:
        skysoltype = args['slow_soltype']  
      else:
        skysoltype = args['phase_soltype']
      
      runDPPP(ms, np.int(solint_ap[msnumber]), np.int(solint_phase[msnumber]), \
               np.int(nchan_phase[msnumber]), np.int(nchan_ap[msnumber]), \
               ms + parmdb + 'SKY_' + str(i) + '.h5' ,skysoltype, preapplyphase=args['addslowsolve_skymodel'], longbaseline=longbaseline, uvmin=uvmin, TEC=TEC, \
               SMconstraintslow=args['smoothnessconstraint_slow'], flagging=args['doflagging'], \
               SMconstraintphase=args['smoothnessconstraint_phase'], puretec=args['pure_tec'], \
               antennaconstraint=args['antennaconstraint_phase'],  soltype_fastphase=args['phase_soltype'], \
               antennaconstraintslow=args['antennaconstraint_slow'],\
               usesourcedb=True, skymodel=args['skymodel'], predictskywithbeam=args['predictskywithbeam'], \
               smoothcal=args['smoothcal'], \
               flagslowphases=args['doflagslowphases'], flagslowamprms=args['flagslowamprms'], flagslowphaserms=args['flagslowphaserms'], BLsmooth=args['BLsmooth'])

      if args['addslowsolve_skymodel']:
        applycal(ms, ['phaseonly' + ms + parmdb + 'SKY_' + str(i) + '.h5', ms + parmdb + 'SKY_' + str(i) + '.h5'])
      else:
        applycal(ms, ['phaseonly' + ms + parmdb + 'SKY_' + str(i) + '.h5'])
      

  # IMAGE
  makeimage(mslist, imageout + str(i), pixsize, imsize, channelsout, np.int(niter), robust, \
            uvtaper=False, multiscale=multiscale, idg=idg, fitsmask=fitsmask, \
            deepmultiscale=args['deepmultiscale'], uvminim=args['uvminim'], fitspectralpol=args['fitspectralpol'], \
            imager=args['imager'], restoringbeam=restoringbeam, automask=automask, \
            removenegativecc=args['removenegativefrommodel'], fitspectralpolorder=args['fitspectralpolorder'], \
            usewgridder=args['usewgridder'], paralleldeconvolution=args['paralleldeconvolution'],\
            deconvolutionchannels=args['deconvolutionchannels'], parallelgridding=args['parallelgridding'])
  
  # MAKE FIGURE WITH APLPY
  if args['imager'] == 'WSCLEAN':
    if idg:
      plotimage(imageout + str(i) +'-MFS-I-image.fits',imageout + str(i) + '.png' , \
                mask=fitsmask, rmsnoiseimage=imageout + str(0) +'-MFS-I-image.fits')
    else:
      plotimage(imageout + str(i) +'-MFS-image.fits',imageout + str(i) + '.png' , \
                mask=fitsmask, rmsnoiseimage=imageout + str(0) +'-MFS-image.fits')
  if args['imager'] == 'DDFACET':
    plotimage(imageout + str(i) +'.app.restored.fits',imageout + str(i) + '.png' , \
                mask=fitsmask, rmsnoiseimage=imageout + str(0) +'.app.restored.fits')


    # MAKE template H5 with pol axis for scalarphase solve, so we can use DDFacet
  if i == 0 and not args['usemodeldataforsolints']:
    makeh5templates(mslist, parmdb, TEC, args['pure_tec'], solint_phase, solint_ap, nchan_phase, nchan_ap)
  #  redetermine solints if requested
  if (i >= 0) and (i <= args['phasecycles']) and (args['usemodeldataforsolints']):
    print('Recomputing solints .... ')

    nchan_phase,solint_phase,solint_ap,nchan_ap, trigger_antennaconstraint = auto_determinesolints(mslist, TEC, \
                           args['pure_tec'], longbaseline, LBA,\
                           smoothnessconstraint_phase=args['smoothnessconstraint_phase'],\
                           smoothnessconstraint_slow=args['smoothnessconstraint_slow'], \
                           innchan_phase=args['nchan_phase'], \
                           innchan_ap=args['nchan_ap'], insolint_phase=args['solint_phase'],\
                           insolint_ap=args['solint_ap'], redo=True)  
    # remake templates because solint and nchans changed probably
    print('Remaking solution templates .... ')
    makeh5templates(mslist, parmdb, TEC, args['pure_tec'], solint_phase, solint_ap, nchan_phase, nchan_ap)

  
  # SOLVE
  print('trigger_antennaconstraint:',trigger_antennaconstraint)
  for msnumber, ms in enumerate(mslist):
    
    # auto trigger antennaconstraint_phase in case of usemodeldataforsolints HBA LoTSS
    if trigger_antennaconstraint != None:
      if trigger_antennaconstraint[msnumber]:
        antennaconstraint_phase = 'core'
      else:
        antennaconstraint_phase = None    
    else:
      antennaconstraint_phase = args['antennaconstraint_phase']  
    
    
    if i < args['phasecycles']:      
      runDPPP(ms, np.int(solint_ap[msnumber]), np.int(solint_phase[msnumber]), \
               np.int(nchan_phase[msnumber]), np.int(nchan_ap[msnumber]), \
               ms + parmdb + str(i) + '.h5' ,args['phase_soltype'], preapplyphase=False, \
               longbaseline=longbaseline, uvmin=uvmin, TEC=TEC, \
               SMconstraintslow=args['smoothnessconstraint_slow'], \
               SMconstraintphase=args['smoothnessconstraint_phase'], puretec=args['pure_tec'], \
               antennaconstraint=antennaconstraint_phase, BLsmooth=args['BLsmooth'])
    else:
      if args['solve_rotation'] == True:
        runDPPP(ms, np.int(solint_ap[msnumber]), np.int(solint_phase[msnumber]), \
               np.int(nchan_phase[msnumber]), np.int(nchan_ap[msnumber]), \
               ms + parmdb + str(i) + '.h5'  ,'rotation+diagonal', preapplyphase=True, \
               longbaseline=longbaseline, uvmin=uvmin, TEC=TEC, \
               SMconstraintslow=args['smoothnessconstraint_slow'], \
               SMconstraintphase=args['smoothnessconstraint_phase'], smoothcal=args['smoothcal'], \
               puretec=args['pure_tec'],soltype_fastphase=args['phase_soltype'], flagging=args['doflagging'], restoreflags=args['restoreflags'],\
               antennaconstraint=antennaconstraint_phase, antennaconstraintslow=args['antennaconstraint_slow'],\
               flagslowphases=args['doflagslowphases'], flagslowamprms=args['flagslowamprms'],\
               flagslowphaserms=args['flagslowphaserms'], preapplyslow=args['dosecondslow'], BLsmooth=args['BLsmooth'])          
      else:
        runDPPP(ms, np.int(solint_ap[msnumber]), np.int(solint_phase[msnumber]), \
               np.int(nchan_phase[msnumber]), np.int(nchan_ap[msnumber]), \
               ms + parmdb + str(i) + '.h5'  ,args['slow_soltype'], preapplyphase=True, \
               longbaseline=longbaseline, uvmin=uvmin, TEC=TEC, \
               SMconstraintslow=args['smoothnessconstraint_slow'], \
               SMconstraintphase=args['smoothnessconstraint_phase'], smoothcal=args['smoothcal'], \
               puretec=args['pure_tec'], soltype_fastphase=args['phase_soltype'], flagging=args['doflagging'], restoreflags=args['restoreflags'], \
               antennaconstraint=antennaconstraint_phase, antennaconstraintslow=args['antennaconstraint_slow'],\
               flagslowphases=args['doflagslowphases'], flagslowamprms=args['flagslowamprms'], \
               flagslowphaserms=args['flagslowphaserms'], preapplyslow=args['dosecondslow'], BLsmooth=args['BLsmooth'])

  # NORMALIZE GLOBAL GAIN (done in log-space)
  if i >= args['phasecycles'] and args['normamps'] and args['slow_soltype'] != 'fulljones':
     print('Doing global gain normalization')
     parmdblist = []  
     for msnumber, ms in enumerate(mslist):
       parmdblist.append(ms + parmdb + str(i) + '.h5')
       if args['dosecondslow']:
           parmdblist.append('pertubation' + ms + parmdb + str(i) + '.h5')
     normamps(parmdblist)

  # APPLYCAL
  for msnumber, ms in enumerate(mslist):
    if i < args['phasecycles']:
      applycal(ms, ['phaseonly' + ms + parmdb + str(i) +'.h5'])
    else:
      applycal(ms, ['phaseonly' + ms + parmdb + str(i) +'.h5', ms + parmdb + str(i) +'.h5'])
      # copy gains for screen imaging
      copyovergain(ms + parmdb + str(i) +'.h5', ms + parmdb + str(0) + '_slowgainversion.h5', args['slow_soltype'])
      
    # Copy over to pol-axis template for screen imaging 
    if args['phase_soltype'] == 'scalarphase' and TEC == False:
      print('Copying over scalephase solutions to XX, YY')
      copyoverscalarphase('phaseonly' + ms + parmdb + str(i) + '.h5', \
                          'phaseonly' + ms + parmdb + str(0) + '_polversion.h5')

 
  # MAKE MASK
  if args['fitsmask'] == None:
    if args['imager'] == 'WSCLEAN':   
      if idg:  
        imagename  = imageout + str(i) + '-MFS-I-image.fits'
      else:
        imagename  = imageout + str(i) + '-MFS-image.fits'
    if args['imager'] == 'DDFACET':
      imagename  = imageout + str(i) +'.app.restored.fits'

    if args['maskthreshold'] > 0:    
      cmdm  = 'MakeMask.py --Th='+ str(args['maskthreshold']) + ' --RestoredIm=' + imagename
      if fitsmask != None:
        if os.path.isfile(fitsmask):
          os.system('rm -f ' + fitsmask)
      os.system(cmdm)
      fitsmask = imagename + '.mask.fits'
    else:
      fitsmask = None # no masking requested as args['maskthreshold'] less/equal 0
        
  
  # CUT FLAGGED DATA FROM MS AT START&END to win some compute time if possible
  #if TEC and not args['forwidefield']: # does not work for phaseonly sols
  #  if (i == 0) or (i == args['phasecycles']) or (i == args['phasecycles'] + 1) or (i == args['phasecycles'] + 2) \
  #    or (i == args['phasecycles'] + 3) or (i == args['phasecycles'] + 4):
  #     for msnumber, ms in enumerate(mslist): 
  #         flagms_startend(ms, 'phaseonly' + ms + parmdb + str(i) + '.h5', np.int(solint_phase[msnumber]))
  
  
if not longbaseline:
 if not LBA:   
  archive(mslist, outtarname, args['boxfile'], fitsmask, imagename)    
  cleanup(mslist)
