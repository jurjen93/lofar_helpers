import numpy as np
#import matplotlib.pyplot as plt
#import astroquery
#from astroquery.vizier import Vizier
#from astropy.coordinates import Angle
import astropy
import matplotlib.pyplot as plt
import aplpy
import sys
import math
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
import montage_wrapper as montage
from astropy.visualization import AsymmetricPercentileInterval, ManualInterval
from astropy.visualization import LogStretch
from PIL import Image
Image.MAX_IMAGE_PIXELS = None

def round_to_n(x, n):
    if not x: return 0
    power = -int(math.floor(math.log10(abs(x)))) + (n - 1)
    factor = (10 ** power)
    return round(x * factor) / factor

def truncate(n, decimals=0):
    multiplier = 10 ** decimals
    return np.int(n * multiplier) / multiplier

def round_up(n, decimals=0):
    multiplier = 10 ** decimals
    return np.ceil(n * multiplier) / multiplier

def round_down(n, decimals=0):
    multiplier = 10 ** decimals
    return np.floor(n * multiplier) / multiplier



def poisson_halo_stat():
  
  NLOTSS_PSZ2 = 595
    
  N_PSZ2 = 26 # number of PSZ2 clusters in DR2_sample
  
  N_diffuse = N_PSZ2 - 7
  
  N_PSZ_halo   = 8 # A1758 counts as one here
  N_PSZ_halo_c = 8
  
  N_PSZ_relic   = 4 # cluster counts
  N_PSZ_relic_c = 3 # cluster counts
  
  N_other_halo   = 1 
  N_other_halo_c = 4
  
  N_other_relic   = 1 # cluster counts
  N_other_relic_c = 0 # cluster counts

  
  vals = (np.random.poisson(lam=np.float(N_PSZ_halo), size=100000)).astype(float)
  print ('PSZ halos', 100.*np.float(N_PSZ_halo)/np.float(N_PSZ2), np.std(100.*vals/np.float(N_PSZ2)))
  
  vals = (np.random.poisson(lam=np.float(N_PSZ_halo+N_PSZ_halo_c), size=100000)).astype(float)
  print ('PSZ halos+candidate', 100.*np.float(N_PSZ_halo+N_PSZ_halo_c)/np.float(N_PSZ2), np.std(100.*vals/np.float(N_PSZ2)))

  vals = (np.random.poisson(lam=np.float(N_diffuse), size=100000)).astype(float)
  print ('PSZ diffuse', 100.*np.float(N_diffuse)/np.float(N_PSZ2), np.std(100.*vals/np.float(N_PSZ2)))

  vals = (np.random.poisson(lam=np.float(N_PSZ_halo), size=100000)).astype(float)
  print ('PSZ halos LOTSS', NLOTSS_PSZ2*np.float(N_PSZ_halo)/np.float(N_PSZ2), NLOTSS_PSZ2*np.std(vals/np.float(N_PSZ2)))

  vals = (np.random.poisson(lam=np.float(N_PSZ_halo+N_PSZ_halo_c), size=100000)).astype(float)
  print ('PSZ halos+candidate LOTSS', NLOTSS_PSZ2*np.float(N_PSZ_halo+N_PSZ_halo_c)/np.float(N_PSZ2), NLOTSS_PSZ2*np.std(vals/np.float(N_PSZ2)))

  vals = (np.random.poisson(lam=np.float(N_diffuse), size=100000)).astype(float)
  print ('PSZ diffuse LOTSS', NLOTSS_PSZ2*np.float(N_diffuse)/np.float(N_PSZ2), NLOTSS_PSZ2*np.std(vals/np.float(N_PSZ2)))

  vals = (np.random.poisson(lam=np.float(N_PSZ_relic), size=100000)).astype(float)
  print ('PSZ relic', 100.*np.float(N_PSZ_relic)/np.float(N_PSZ2), np.std(100.*vals/np.float(N_PSZ2)))

  vals = (np.random.poisson(lam=np.float(N_PSZ_relic), size=100000)).astype(float)
  print ('PSZ relic LOTSS', NLOTSS_PSZ2*np.float(N_PSZ_relic)/np.float(N_PSZ2), NLOTSS_PSZ2*np.std(vals/np.float(N_PSZ2)))
  


  return
poisson_halo_stat()


def paddimage(fitsimage, padsize=2048):
    print fitsimage
    padsize = np.int(padsize) # make integer, just in case
    
    hdu = fits.open(fitsimage)
    imdata = hdu[0].data

    (xsize, ysize) = imdata.shape
    #assert(xsize == ysize)
    print 'size is', xsize, ysize

        
    offsetx = (padsize - xsize) / 2
    offsety = (padsize - ysize) / 2
    print 'padding to', padsize
    print 'offsets are', offsetx, offsety

    newdata=np.zeros((1, 1, padsize, padsize))

    newdata[0, 0, offsetx:offsetx+xsize, offsety:offsety+ysize] = imdata
    hdu[0].data = newdata
    hdu[0].header['CRPIX1'] += offsety
    hdu[0].header['CRPIX2'] += offsetx
    hdu.writeto(fitsimage + '.pad', clobber=True)
    return



def header_addbeam(inputimage, refimage):
    hduref = fits.open(refimage)
    bmaj = hduref[0].header['BMAJ']
    bmin = hduref[0].header['BMIN']
    bpa  = hduref[0].header['BPA']
    hduref.close()

    hduout = fits.open(inputimage)
    hduout[0].header['BMAJ'] = bmaj
    hduout[0].header['BMIN'] = bmin
    hduout[0].header['BPA'] = bpa
    
    fits.writeto(inputimage, hduout[0].data, hduout[0].header, overwrite=True)
    hduout.close()
    
    return


def selfcalimagesLB(tfont=22, lfont=24, dpi=36):
    fitsimagename1 = '/net/voorrijn/data2/rvweeren/8C1654+785/image_sc0.18_0-MFS-image.fits'    
    fitsimagename2 = '/net/voorrijn/data2/rvweeren/8C1654+785/image_sc0.18_3-MFS-image.fits'
    fitsimagename3 = '/net/voorrijn/data2/rvweeren/8C1654+785/image_sc0.18_9-MFS-image.fits'   
    
    center = [252.90309,78.48346222]
    width  = [0.022,0.022]

    outname1 = 'LBstep1.pdf'
    outname2 = 'LBstep2.pdf'
    outname3 = 'LBstep3.pdf'

    imagenoise = 39e-6
    rmsfactor  = 100.
    cmap = 'cubehelix_r'

    fig = aplpy.FITSFigure(fitsimagename1,slices=[0,0])
    fig.show_colorscale(cmap=cmap,vmin=0,vmax=rmsfactor*imagenoise,stretch='power',exponent=0.5)
    fig.recenter(center[0],center[1], width=width[0], height=width[1])
    #fig.add_grid()
    #fig.grid.set_color('black')
    #fig.grid.set_alpha(0.2)
    fig.add_beam()
    fig.tick_labels.set_xformat('hh:mm:ss')
    fig.tick_labels.set_yformat('dd:mm:ss')    
    fig.beam.set_frame(True)
    fig.axis_labels.set_font(size=lfont)
    fig.tick_labels.set_font(size=tfont)
    fig.set_title('DI',size=lfont+1)
    fig.savefig(outname1, dpi=dpi)
    fig.close()
    
    fig = aplpy.FITSFigure(fitsimagename2,slices=[0,0])
    fig.show_colorscale(cmap=cmap,vmin=0,vmax=rmsfactor*imagenoise,stretch='power',exponent=0.5)
    fig.recenter(center[0],center[1], width=width[0], height=width[1])
    #fig.add_grid()
    #fig.grid.set_color('black')
    #fig.grid.set_alpha(0.2)
    fig.tick_labels.set_xformat('hh:mm:ss')
    fig.tick_labels.set_yformat('dd:mm:ss')    
    fig.add_beam()
    fig.beam.set_frame(True)
    fig.axis_labels.set_font(size=lfont)
    fig.tick_labels.set_font(size=tfont)
    fig.set_title('tec',size=lfont+1)
    fig.savefig(outname2, dpi=dpi)
    fig.close()
    
    fig = aplpy.FITSFigure(fitsimagename3,slices=[0,0])
    fig.show_colorscale(cmap=cmap,vmin=0,vmax=rmsfactor*imagenoise,stretch='power',exponent=0.5)
    fig.recenter(center[0],center[1], width=width[0], height=width[1])
    #fig.add_grid()
    #fig.grid.set_color('black')
    fig.tick_labels.set_xformat('hh:mm:ss')
    fig.tick_labels.set_yformat('dd:mm:ss')    
    #fig.grid.set_alpha(0.2)
    fig.add_beam()
    fig.beam.set_frame(True)
    fig.axis_labels.set_font(size=lfont)
    fig.tick_labels.set_font(size=tfont)
    fig.set_title('tec + scalarcomplex gain',size=lfont+1)
    fig.savefig(outname3, dpi=dpi)
    fig.close()    

    
    return


def selfcalimagesLBA(tfont=22, lfont=24, dpi=36):
    fitsimagename1 = '/net/tussenrijn/data2/rvweeren/p176+60_selfcal/AP4/imageAP4_0.app.restored.fits'    
    fitsimagename2 = '/net/tussenrijn/data2/rvweeren/p176+60_selfcal/AP4/imageAP4_2.app.restored.fits'
    fitsimagename3 = '/net/tussenrijn/data2/rvweeren/p176+60_selfcal/AP4/imageAP4_4.app.restored.fits'   
    
    center = [180.4325387,58.07346288]
    width  = [0.3,0.3]

    outname1 = 'LBAstep1.pdf'
    outname2 = 'LBAstep2.pdf'
    outname3 = 'LBAstep3.pdf'

    imagenoise = 946e-6
    rmsfactor  = 500.
    cmap = 'cubehelix_r'

    fig = aplpy.FITSFigure(fitsimagename1,slices=[0,0])
    fig.show_colorscale(cmap=cmap,vmin=0,vmax=rmsfactor*imagenoise,stretch='power',exponent=0.5)
    fig.recenter(center[0],center[1], width=width[0], height=width[1])
    #fig.add_grid()
    #fig.grid.set_color('black')
    #fig.grid.set_alpha(0.2)
    fig.add_beam()
    fig.tick_labels.set_xformat('hh:mm:ss')
    fig.tick_labels.set_yformat('dd:mm:ss')    
    fig.beam.set_frame(True)
    fig.axis_labels.set_font(size=lfont)
    fig.tick_labels.set_font(size=tfont)
    fig.set_title('DI',size=lfont+1)
    fig.savefig(outname1, dpi=dpi)
    fig.close()
    
    fig = aplpy.FITSFigure(fitsimagename2,slices=[0,0])
    fig.show_colorscale(cmap=cmap,vmin=0,vmax=rmsfactor*imagenoise,stretch='power',exponent=0.5)
    fig.recenter(center[0],center[1], width=width[0], height=width[1])
    #fig.add_grid()
    #fig.grid.set_color('black')
    #fig.grid.set_alpha(0.2)
    fig.tick_labels.set_xformat('hh:mm:ss')
    fig.tick_labels.set_yformat('dd:mm:ss')    
    fig.add_beam()
    fig.beam.set_frame(True)
    fig.axis_labels.set_font(size=lfont)
    fig.tick_labels.set_font(size=tfont)
    fig.set_title('tec',size=lfont+1)
    fig.savefig(outname2, dpi=dpi)
    fig.close()
    
    fig = aplpy.FITSFigure(fitsimagename3,slices=[0,0])
    fig.show_colorscale(cmap=cmap,vmin=0,vmax=rmsfactor*imagenoise,stretch='power',exponent=0.5)
    fig.recenter(center[0],center[1], width=width[0], height=width[1])
    #fig.add_grid()
    #fig.grid.set_color('black')
    fig.tick_labels.set_xformat('hh:mm:ss')
    fig.tick_labels.set_yformat('dd:mm:ss')    
    #fig.grid.set_alpha(0.2)
    fig.add_beam()
    fig.beam.set_frame(True)
    fig.axis_labels.set_font(size=lfont)
    fig.tick_labels.set_font(size=tfont)
    fig.set_title('tec + scalarcomplex gain',size=lfont+1)
    fig.savefig(outname3, dpi=dpi)
    fig.close()    

    
    return


def selfcalimages(dpi=36):
    fitsimagename1 = '/net/rijn/data2/rvweeren/PSZ2G143.26+65.24_test_largeimpad/imselfcal_0-MFS-image.fits'    
    fitsimagename2 = '/net/rijn/data2/rvweeren/PSZ2G143.26+65.24_test_largeimpad/imselfcal_3-MFS-image.fits'
    fitsimagename3 = '/net/rijn/data2/rvweeren/PSZ2G143.26+65.24_test_largeimpad/imselfcal_9-MFS-image.fits'   
    
    center = [179.934339,49.77317281]
    width  = [0.45,0.3]

    outname1 = 'A1430step1.pdf'
    outname2 = 'A1430step2.pdf'
    outname3 = 'A1430step3.pdf'

    imagenoise = 60e-6
    rmsfactor  = 50.
    cmap = 'cubehelix_r'

    fig = aplpy.FITSFigure(fitsimagename1,slices=[0,0])
    fig.show_colorscale(cmap=cmap,vmin=0,vmax=rmsfactor*imagenoise,stretch='power',exponent=0.5)
    fig.recenter(center[0],center[1], width=width[0], height=width[1])
    #fig.add_grid()
    #fig.grid.set_color('black')
    #fig.grid.set_alpha(0.2)
    fig.add_beam()
    fig.tick_labels.set_xformat('hh:mm:ss')
    fig.tick_labels.set_yformat('dd:mm:ss')
    fig.beam.set_frame(True)
    fig.axis_labels.set_font(size=15)
    fig.tick_labels.set_font(size=14)
    fig.set_title('DR2-DI',size=14)
    fig.savefig(outname1, dpi=dpi)
    fig.close()
    
    fig = aplpy.FITSFigure(fitsimagename2,slices=[0,0])
    fig.show_colorscale(cmap=cmap,vmin=0,vmax=rmsfactor*imagenoise,stretch='power',exponent=0.5)
    fig.recenter(center[0],center[1], width=width[0], height=width[1])
    #fig.add_grid()
    #fig.grid.set_color('black')
    #fig.grid.set_alpha(0.2)
    fig.add_beam()
    fig.tick_labels.set_xformat('hh:mm:ss')
    fig.tick_labels.set_yformat('dd:mm:ss')    
    fig.beam.set_frame(True)
    fig.axis_labels.set_font(size=15)
    fig.tick_labels.set_font(size=14)
    fig.set_title('tecandphase',size=14)
    fig.savefig(outname2, dpi=dpi)
    fig.close()
    
    fig = aplpy.FITSFigure(fitsimagename3,slices=[0,0])
    fig.show_colorscale(cmap=cmap,vmin=0,vmax=rmsfactor*imagenoise,stretch='power',exponent=0.5)
    fig.recenter(center[0],center[1], width=width[0], height=width[1])
    #fig.add_grid()
   # fig.grid.set_color('black')
    fig.tick_labels.set_xformat('hh:mm:ss')
    fig.tick_labels.set_yformat('dd:mm:ss')   
    #fig.grid.set_alpha(0.2)
    fig.add_beam()
    fig.beam.set_frame(True)
    fig.axis_labels.set_font(size=15)
    fig.tick_labels.set_font(size=14)
    fig.set_title('tecandphase + diagonal gain',size=14)
    fig.savefig(outname3, dpi=dpi)
    fig.close()    

 
    return


def makecomparisonimages(dpi=72):

    fitsimagenameDR2 = '/disks/paradata/shimwell/LoTSS-DR2/mosaics/RA13h_field/P18Hetdex03/mosaic.fits'
    fitsimagename    = '/disks/paradata/shimwell/LoTSS-DR2/archive_cluster_images/PSZ2G143.26+65.24/PSZ2G143.26+65.24_maskROBUST-0.5-MFS-image.fits'
    center = [179.934339,49.77317281]
    width  = [0.45,0.3]

    outnameDR2 = 'A1430DR2comp.pdf'
    outname    = 'A1430comp.pdf'
    #hdulist = fits.open(fitsimagename) 
    #imagenoise = findrms(np.ndarray.flatten(hdulist[0].data))
    #hdulist.close()   
    imagenoise = 60e-6
    rmsfactor  = 50.
    cmap = 'cubehelix_r'

    fig = aplpy.FITSFigure(fitsimagename,slices=[0,0])
    fig.show_colorscale(cmap=cmap,vmin=0,vmax=rmsfactor*imagenoise,stretch='power',exponent=0.5)
    fig.recenter(center[0],center[1], width=width[0], height=width[1])
    #fig.add_grid()
    #fig.grid.set_color('black')
    #fig.grid.set_alpha(0.2)
    fig.add_beam()
    fig.beam.set_frame(True)
    fig.axis_labels.set_font(size=15)
    fig.tick_labels.set_font(size=14)
    fig.set_title('DR2 + reprocessing',size=14)
    fig.savefig(outname, dpi=dpi)
    fig.close()
        
    fig = aplpy.FITSFigure(fitsimagenameDR2,slices=[0,0])
    fig.show_colorscale(cmap=cmap,vmin=0,vmax=rmsfactor*imagenoise,stretch='power',exponent=0.5)
    fig.recenter(center[0],center[1], width=width[0], height=width[1])
    #fig.add_grid()
    #fig.grid.set_color('black')
    #fig.grid.set_alpha(0.2)
    fig.add_beam()
    fig.beam.set_frame(True)
    fig.axis_labels.set_font(size=15)
    fig.tick_labels.set_font(size=14)
    fig.set_title('DR2',size=14)
    fig.savefig(outnameDR2, dpi=dpi)
    fig.close()


    #A1550
    fitsimagenameDR2 = '/disks/paradata/shimwell/LoTSS-DR2/mosaics/RA13h_field/P23Hetdex20/mosaic.fits'
    fitsimagename    = '/disks/paradata/shimwell/LoTSS-DR2/archive_cluster_images/PSZ2G133.60+69.04/PSZ2G133.60+69.04_maskROBUST-0.5-MFS-image.fits'  
    center    = [187.2604866,+47.62249881]
    outname   = 'A1550comp.pdf'
    outnameDR2= 'A1550DR2comp.pdf'
    #width     = [0.3,0.3]

    fig = aplpy.FITSFigure(fitsimagename,slices=[0,0])
    fig.show_colorscale(cmap=cmap,vmin=0,vmax=rmsfactor*imagenoise,stretch='power',exponent=0.5)
    fig.recenter(center[0],center[1], width=width[0], height=width[1])
    #fig.add_grid()
    #fig.grid.set_color('black')
    #fig.grid.set_alpha(0.2)
    fig.add_beam()
    fig.beam.set_frame(True)
    fig.axis_labels.set_font(size=15)
    fig.tick_labels.set_font(size=14)
    fig.set_title('DR2 + reprocessing',size=14)
    fig.savefig(outname, dpi=dpi)
    fig.close()
        
    fig = aplpy.FITSFigure(fitsimagenameDR2,slices=[0,0])
    fig.show_colorscale(cmap=cmap,vmin=0,vmax=rmsfactor*imagenoise,stretch='power',exponent=0.5)
    fig.recenter(center[0],center[1], width=width[0], height=width[1])
    #fig.add_grid()
    #fig.grid.set_color('black')
    #fig.grid.set_alpha(0.2)
    fig.add_beam()
    fig.beam.set_frame(True)
    fig.axis_labels.set_font(size=15)
    fig.tick_labels.set_font(size=14)
    fig.set_title('DR2',size=14)
    fig.savefig(outnameDR2, dpi=dpi)
    fig.close()


    #P173+55
    fitsimagenameDR2 = '/disks/paradata/shimwell/LoTSS-DR2/mosaics/RA13h_field/P173+55/mosaic.fits'
    fitsimagename    = '/disks/paradata/shimwell/LoTSS-DR2/archive_cluster_images/PSZ2G145.65+59.30/PSZ2G145.65+59.30_maskROBUST-0.5-MFS-image.fits'  
    outname   = 'A1294comp.pdf'
    outnameDR2= 'A1294DR2comp.pdf'
    center    = [173.22056,54.2209406]
    #width     = [0.2,0.2]

    fig = aplpy.FITSFigure(fitsimagename,slices=[0,0])
    fig.show_colorscale(cmap=cmap,vmin=0,vmax=rmsfactor*imagenoise,stretch='power',exponent=0.5)
    fig.recenter(center[0],center[1], width=0.5*width[0], height=0.5*width[1])
    #fig.add_grid()
    #fig.grid.set_color('black')
    #fig.grid.set_alpha(0.2)
    fig.add_beam()
    fig.beam.set_frame(True)
    fig.axis_labels.set_font(size=15)
    fig.tick_labels.set_font(size=14)
    fig.set_title('DR2 + reprocessing',size=14)
    fig.savefig(outname, dpi=dpi)
    fig.close()
        
    fig = aplpy.FITSFigure(fitsimagenameDR2,slices=[0,0])
    fig.show_colorscale(cmap=cmap,vmin=0,vmax=rmsfactor*imagenoise,stretch='power',exponent=0.5)
    fig.recenter(center[0],center[1], width=0.5*width[0], height=0.5*width[1])
    #fig.add_grid()
    #fig.grid.set_color('black')
    #fig.grid.set_alpha(0.2)
    fig.add_beam()
    fig.beam.set_frame(True)
    fig.axis_labels.set_font(size=15)
    fig.tick_labels.set_font(size=14)
    fig.set_title('DR2',size=14)
    fig.savefig(outnameDR2, dpi=dpi)
    fig.close()

    
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

def makeimage(fitsimagename, center, width, z, name, outname, rmsfactor=50.,plotcode=None,\
              tfont=15, lfont=17, xs=0.25/3, pexponent=0.5, dpi=72, cb=False, rescalemjy=False, \
              contour=False, contourstartlevel=3.0, plotgrid=False):
    
    oneradinmpc = cosmo.angular_diameter_distance(z)/(360./(2.*np.pi))
    scalebarlengthdeg    = 1.0/oneradinmpc.value
    
    if rescalemjy:
      hdulist = fits.open(fitsimagename) 
      imarray = hdulist[0].data
      hdulist[0].data = imarray*1e3 # to mJy
      hdulist.writeto('tmp.fits', overwrite=True)
      fitsimagename = 'tmp.fits'
    
    hdulist = fits.open(fitsimagename) 
    imagenoise = findrms(np.ndarray.flatten(hdulist[0].data))
    hdulist.close()   


    width = [width[0]/oneradinmpc.value,width[1]/oneradinmpc.value]

    fig = aplpy.FITSFigure(fitsimagename,slices=[0,0])
    #fig.show_colorscale(cmap='gist_heat_r',vmin=-2*8.4e-6,vmax=50*8.4e-6,stretch='power',exponent=0.5)
    fig.show_colorscale(cmap='cubehelix_r',vmin=imagenoise*1.0,vmax=rmsfactor*imagenoise,stretch='power',exponent=pexponent)
    fig.recenter(center[0],center[1], width=width[0], height=width[1])
    
    if contour:
      colors = ('black','black','black','black','white','white','white','white','white')  
      levelsr = np.ndarray.tolist (contourstartlevel*imagenoise*np.array([1.,2,4,8,16,32,64,128,256]))
      fig.show_contour(fitsimagename, slices=[0,0],levels=levelsr, colors=colors,overlap=True, linewidths=0.75)
      levelsr = np.ndarray.tolist (-3.*imagenoise*np.array([1.0,]))
      fig.show_contour(fitsimagename, slices=[0,0],levels=levelsr, colors='red',overlap=True, linewidths=0.75, linestyles='solid'), #linestyles='dashed')



    if plotgrid:
       fig.add_grid()
       fig.grid.set_color('black')
       fig.grid.set_alpha(0.3)
    fig.add_beam()
    fig.beam.set_frame(True)
    fig.tick_labels.set_xformat('hh:mm:ss')
    fig.tick_labels.set_yformat('dd:mm:ss')
    fig.axis_labels.set_font(size=lfont)
    fig.tick_labels.set_font(size=tfont)
    fig.set_title(name,fontsize=lfont+1)
    fig.ticks.set_xspacing(xs)
    if cb:
        fig.add_colorbar()
        fig.colorbar.set_axis_label_text('Surface Brightness (mJy beam$^{-1}$)')
        fig.colorbar.set_axis_label_font(size=tfont)
        fig.colorbar.set_font(size=tfont-2)
        tvector = [round_to_n(round_up(imagenoise*1.0,3),3),\
                   round_to_n(round_up(imagenoise*rmsfactor*0.01,3),2),\
                   round_to_n(round_up(imagenoise*rmsfactor*0.05,3),2),\
                   round_to_n(round_up(imagenoise*rmsfactor*0.1,3),2),\
                   round_to_n(round_up(imagenoise*rmsfactor*0.5,3),2),\
                   round_to_n(round_down(0.9*imagenoise*rmsfactor,3),2)]
        print tvector
        fig.colorbar.set_ticks(tvector)

    fig.add_scalebar(scalebarlengthdeg, '1 Mpc',color = 'black', corner='bottom', linewidth=2, fontsize=lfont+1)
    fig.add_label(0.03, 0.95, '$\sigma_{rms}$=' + str(np.int(imagenoise*1e6)) + ' $\mu$Jy beam$^{-1}$',\
                  relative=True, size=tfont,color='black',horizontalalignment='left')


    if plotcode == 'A1430':
      fig.add_label(179.810144,49.81468146,'A',color='black',size=20)
      fig.add_label(179.8979619, 49.81030436,'B',color='black',size=20)
    if plotcode == 'A1550':
      fig.add_label(187.2442737,47.62096424,'A',color='black',size=20)
      fig.add_label(187.2806192, 47.62037709,'H',color='black',size=20)
      
      fig.add_label(187.2913335, 47.67684966,'B',color='black',size=20)
      fig.add_label(187.1781451, 47.59159319,'C',color='black',size=20)
      fig.add_label(187.1864186 ,47.62257925 ,'D',color='black',size=20)
      fig.add_label(187.1485607 , 47.61463872,'E',color='black',size=20)

    if plotcode == 'A1294':
      fig.add_label(173.0990242,54.21856346,'A',color='black',size=20)
      fig.add_label(173.1321849,54.2244528,'B',color='black',size=20)

    if plotcode == 'A1904':
      fig.add_label(215.5745314,48.55523165,'A',color='black',size=20)
      fig.add_label(215.5919001,48.48754991 ,'B',color='black',size=20)
      fig.add_label(215.5070884, 48.6185479,'C',color='black',size=20)
      #fig.add_label(215.5070884, 48.6185479,'C',color='black',size=20)

    if plotcode == 'A1904LF':
      fig.add_label(215.5745314,48.55,'A',color='black',size=20)
      fig.add_label(215.5919001,48.48754991 ,'B',color='black',size=20)
      fig.add_label(215.5070884, 48.62,'C',color='black',size=20)
      fig.add_label(215.8553675, 48.81684833,'D',color='black',size=20)

    if plotcode == 'A1758':
      fig.add_label(203.2833385,50.4027055,'R',color='black',size=20)
      fig.add_label(203.2357502,50.56683517 ,'S1',color='black',size=20)
      fig.add_label(203.3029763, 50.55155483,'S2',color='black',size=20)

    if plotcode == 'MACSJ1115.2+5320':
      fig.add_label(168.7979925,53.32027878,'A',color='black',size=20)
      fig.add_label(168.8355282,53.35308004,'B',color='black',size=20)
      fig.add_label(168.8472406,53.32752072,'C',color='black',size=20)
      fig.add_label(168.8391675,53.30702502,'D',color='black',size=20)

    if plotcode == 'PSZ2G135.17+65.43':      
      fig.add_label(184.7939839,50.92083181,'A',color='black',size=20)
      fig.add_label(184.7745623,50.9048186,'B',color='black',size=20)

      fig.add_label(184.841162,50.89833856,'C',color='black',size=20)
      fig.add_label(184.7292764,50.91313117,'D',color='black',size=20)

      fig.add_label(184.7414669,50.88015664,'E',color='black',size=20)
      fig.add_label(184.8118615,50.94165743,'F',color='black',size=20)

      fig.add_label(184.7899135,50.86525981,'G',color='black',size=20)
      fig.add_label(184.8019774,50.90345936,'H',color='black',size=20)

    if plotcode == 'A1622':
      fig.add_label(192.403572,49.89186307,'A',color='black',size=20)
      fig.add_label(192.486205,49.82468134 ,'B',color='black',size=20)
      
    if plotcode == 'A1703':
      fig.add_label(198.7742329,51.80985756,'A',color='black',size=20)
      fig.add_label(198.7688436,51.84527469,'B',color='black',size=20)
      fig.add_label(198.8207187,51.79401207,'C',color='black',size=20)

    if plotcode == 'A1682':
      fig.add_label(196.7214475,46.55341137,'CH',color='black',size=20)

    if plotcode == 'PSZ2G118.34':
      fig.add_label(195.3697389,48.25020845,'H',color='black',size=20) 
      fig.add_label(195.3246857,48.25016989,'A',color='black',size=20) 
      fig.add_label(195.406665,48.2406436,'B',color='black',size=20)

    if plotcode == 'A1314':
      fig.add_label(173.4418709,49.02699239 ,'IC 708',color='black',size=20)
      fig.add_label(173.7119277,48.93722749 ,'IC 711',color='black',size=20)
      fig.add_label(173.7336122,49.04638591 ,'IC 712',color='black',size=20)
      

    fig.savefig(outname, dpi=dpi)
    fig.close()

def makeimagewrapper(prefixpath, fname, center, width, z, name, rmsfactor=50., \
                     iname=None, plotcode=None, lfont=17, tfont=15, xs=0.25/3, pexponent=0.5, cb=False):
    
    if iname == None:
      iname = fname
    
    makeimage(prefixpath + fname + '/' + iname + '_maskROBUST-1.25-MFS-image.fits', \
            center, width, z, name, fname + 'R-1.25.pdf', rmsfactor=rmsfactor, lfont=lfont, tfont=tfont, xs=xs, pexponent=pexponent, cb=cb)
    makeimage(prefixpath + fname + '/' + iname + '_maskROBUST-0.5-MFS-image.fits', \
            center, width, z, name, fname + '.pdf', rmsfactor=rmsfactor, plotcode=plotcode, lfont=lfont, tfont=tfont, xs=xs, pexponent=pexponent, cb=cb)
    makeimage(prefixpath + fname+ '/' + iname + '_maskROBUST-0.5TAPER10-MFS-image.fits', \
            center, width, z, name, fname + 'T10.pdf', rmsfactor=rmsfactor, lfont=lfont, tfont=tfont, xs=xs, pexponent=pexponent, cb=cb,contour=True)
    makeimage(prefixpath + fname + '/' + iname + '_maskROBUST-0.5TAPER15-MFS-image.fits', \
            center, width, z,  name, fname + 'T15.pdf', rmsfactor=rmsfactor, lfont=lfont, tfont=tfont, xs=xs, pexponent=pexponent, cb=cb,contour=True)
    makeimage(prefixpath + fname + '/' + iname +'_maskROBUST-0.5TAPER30-MFS-image.fits', \
            center, width, z,  name, fname + 'T30.pdf', rmsfactor=rmsfactor, lfont=lfont, tfont=tfont, xs=xs, pexponent=pexponent, cb=cb,contour=True)

    if fname != 'PSZ2G089.52+62.34':
      makeimage(prefixpath + fname + '/' + iname + '_masksubROBUST-0.5TAPER10-MFS-image.fits', \
            center, width, z, name, fname + 'T10SUB.pdf',rmsfactor=rmsfactor, lfont=lfont, tfont=tfont, xs=xs, pexponent=pexponent, cb=cb, contour=True)
      makeimage(prefixpath + fname + '/' + iname + '_masksubROBUST-0.5TAPER15-MFS-image.fits', \
            center, width, z,  name, fname + 'T15UB.pdf',rmsfactor=rmsfactor, lfont=lfont, tfont=tfont, xs=xs, pexponent=pexponent, cb=cb, contour=True)
      makeimage(prefixpath + fname + '/' + iname + '_masksubROBUST-0.5TAPER30-MFS-image.fits', \
            center, width, z,  name, fname + 'T30SUB.pdf',rmsfactor=rmsfactor, lfont=lfont, tfont=tfont, xs=xs, pexponent=pexponent, cb=cb, contour=True)

    if fname == 'PSZ2G089.52+62.34' and plotcode == 'A1904LF':
      makeimage(prefixpath + fname+ '/' + iname + '_maskROBUST-0.5TAPER10-MFS-image.fits', \
                center, width, z, name, fname + 'T10.pdf', rmsfactor=rmsfactor,  plotcode=plotcode, lfont=lfont, tfont=tfont, xs=xs, pexponent=pexponent, cb=cb,contour=True)
      makeimage(prefixpath + fname + '/' + iname + '_maskROBUST-0.5TAPER15-MFS-image.fits', \
                center, width, z,  name, fname + 'T15.pdf', rmsfactor=rmsfactor,  plotcode=plotcode, lfont=lfont, tfont=tfont, xs=xs, pexponent=pexponent, cb=cb,contour=True)
      makeimage(prefixpath + fname + '/' + iname +'_maskROBUST-0.5TAPER30-MFS-image.fits', \
            center, width, z,  name, fname + 'T30.pdf', rmsfactor=rmsfactor, plotcode=plotcode, lfont=lfont, tfont=tfont, xs=xs, pexponent=pexponent, cb=cb,contour=True)
    


def Xray_radio(xrayimage, radioimage, outfigname, z, center, Mpcwidth=[3.0,3.0],\
               titlename='',vmin=1e-9,vmax=2e-7,smooth=5,pad=False, plotgrid=False, \
               contourstartlevel=3.0, tfont=15, lfont=17, xs=0.25/3, dpi=72):

    if pad:
       paddimage(xrayimage) 
       xrayfitsim = xrayimage + '.pad'
    else:
       xrayfitsim =  xrayimage 

    header_addbeam(xrayfitsim, radioimage)

    oneradinmpc = cosmo.angular_diameter_distance(z)/(360./(2.*np.pi))
    scalebarlengthdeg    = 1.0/oneradinmpc.value

    width = [Mpcwidth[0]/oneradinmpc.value,Mpcwidth[1]/oneradinmpc.value]

    fig = aplpy.FITSFigure(xrayfitsim,slices=[0,0])
    #fig.show_colorscale(cmap='gist_heat',vmin=1*110e-6,vmax=100*110e-6,stretch='sqrt')
    #fig.show_colorscale(cmap='magma',vmin=1*9.40042e-07,vmax=2*9.40042e-07,stretch='power',exponent=0.35,smooth=9)
    fig.show_colorscale(cmap='magma',vmin=vmin,vmax=vmax,stretch='log',smooth=smooth)
    fig.recenter(center[0],center[1], width=width[0], height=width[1])
    if plotgrid:
       fig.add_grid()
       fig.grid.set_color('white')
       fig.grid.set_alpha(0.3)

    fig.tick_labels.set_xformat('hh:mm:ss')
    fig.tick_labels.set_yformat('dd:mm:ss')

    fig.axis_labels.set_font(size=lfont)
    fig.tick_labels.set_font(size=tfont)


    
    hdulist = fits.open(radioimage) 
    imagenoise = findrms(np.ndarray.flatten(hdulist[0].data))
    hdulist.close()       

    fig.add_label(0.03, 0.95, '$\sigma_{rms}$=' + str(np.int(imagenoise*1e6)) + ' $\mu$Jy beam$^{-1}$',\
                  relative=True, size=tfont,color='white',horizontalalignment='left')
    
    levelsr = np.ndarray.tolist (contourstartlevel*imagenoise*np.array([1.,2,4,8,16,32,64,128,256]))
    fig.show_contour(radioimage, slices=[0,0],levels=levelsr, colors='white',overlap=True, linewidths=1.5)
    #fig.show_regions('RX42edge.reg')
    levelsr = np.ndarray.tolist (-3.*imagenoise*np.array([1.0,]))
    fig.show_contour(radioimage, slices=[0,0],levels=levelsr, colors='red',overlap=True, linewidths=1.0, linestyles='solid'), #linestyles='dashed')

    fig.ticks.set_xspacing(xs)

    fig.add_scalebar(scalebarlengthdeg, '1 Mpc',color = 'white', corner='bottom', linewidth=2,fontsize=lfont+1)
    fig.add_beam()
    fig.beam.set_frame(True)
    #fig.set_nan_color('red')

    fig.set_title(titlename, fontsize=lfont+1)
    fig.savefig(outfigname, dpi=dpi)
    #fig.savefig('test.pdf', dpi=72)
    fig.close()


def makeopticalfigure(redim,greenim,blueim,radioimage,filters,z,radec,clustername,imagename,\
                      strength_red=1.,strength_green=1.,strength_blue=1.,Mpcwidth=[0.5,0.5], \
                      lfont=17, tfont=15, xs=0.25/3, plotcode=None, dpi=72, plotgrid=False):
  r = fits.open(redim)[0]
  g = fits.open(greenim)[0]
  b = fits.open(blueim)[0]

  oneradinmpc = cosmo.angular_diameter_distance(z)/(360./(2.*np.pi))
  scalebarlengthdeg    = 1.0/oneradinmpc.value
  width = [Mpcwidth[0]/oneradinmpc.value,Mpcwidth[1]/oneradinmpc.value]
  
 
  
  widthc = [3./oneradinmpc.value,3./oneradinmpc.value]

  if Mpcwidth[0] >3.0:

     widthc = [(Mpcwidth[0]+1.)/oneradinmpc.value,(Mpcwidth[1]+1.)/oneradinmpc.value]

  # Image reprojection
  montage.commands.mGetHdr(redim, 'r_header')
  montage.wrappers.reproject(redim, 'r_reprojected.fits', header='r_header')
  montage.wrappers.reproject(greenim, 'g_reprojected.fits', header='r_header') 
  montage.wrappers.reproject(blueim, 'b_reprojected.fits', header='r_header')


  montage.mSubimage('r_reprojected.fits', 'r_reprojected_crop.fits', radec[0], radec[1], widthc[0], widthc[1])
  montage.mSubimage('g_reprojected.fits', 'g_reprojected_crop.fits', radec[0], radec[1], widthc[0], widthc[1])
  montage.mSubimage('b_reprojected.fits', 'b_reprojected_crop.fits', radec[0], radec[1], widthc[0], widthc[1])  

  r_r = fits.open('r_reprojected_crop.fits')[0]
  g_r = fits.open('g_reprojected_crop.fits')[0]
  b_r = fits.open('b_reprojected_crop.fits')[0]
    
  
  r_r.data[np.where(np.isnan(r_r.data))] = 0.0
  g_r.data[np.where(np.isnan(g_r.data))] = 0.0
  b_r.data[np.where(np.isnan(b_r.data))] = 0.0


  oapinterval = AsymmetricPercentileInterval(75., 99.96)
  #redinterval   = ManualInterval(1.e-7,10000)
  #greeninterval = ManualInterval(1.e-7,10000)
  #blueinterval  = ManualInterval(1.e-7,10000)
  

  logstretch = LogStretch()

  r_t = logstretch(oapinterval(r_r.data)) * 255.
  g_t = logstretch(oapinterval(g_r.data)) * 255.
  b_t = logstretch(oapinterval(b_r.data)) * 255.

  red   =  strength_red   * r_t
  green =  strength_green * g_t
  blue  =  strength_blue  * b_t
 
  print red.min(), red.max(), red.mean()
  print green.min(), green.max(), green.mean()
  print blue.min(), blue.max(), blue.mean()

  fits.writeto('red.fits', red, r_r.header, overwrite=True)
  fits.writeto('green.fits', green, r_r.header, overwrite=True)
  fits.writeto('blue.fits', blue, r_r.header, overwrite=True)

  outputfits= 'optical_'+filters
  outputpng = 'optical_'+filters+'.png'
  aplpy.make_rgb_cube(['red.fits', 'green.fits', 'blue.fits'], outputfits+'.fits')
  aplpy.make_rgb_image(outputfits+'.fits',outputpng, stretch_r = 'linear', stretch_g = 'linear', stretch_b = 'linear', 
                     vmin_r = 0., vmin_g = 0., vmin_b = 0., vmax_r = 255., vmax_g = 255., vmax_b = 255., 
                     embed_avm_tags=True) # play with limits

  header_addbeam(outputfits+'_2d.fits', radioimage)


  fig = aplpy.FITSFigure(outputfits+'_2d.fits')
  fig.show_rgb(outputpng)
  
  #fig.axis_labels.set_xtext('Right Ascension (J2000)')
  #fig.axis_labels.set_ytext('Declination (J2000)')
  fig.axis_labels.set_font(size=lfont)
  fig.tick_labels.set_font(size=tfont)
  #fig.axis_labels.hide()
  
  

  fig.ticks.set_color('white')
  fig.tick_labels.set_xformat('hh:mm:ss')
  fig.tick_labels.set_yformat('dd:mm:ss')
  
  if plotgrid:
    fig.add_grid()
    fig.grid.set_color('white')
    fig.grid.set_alpha(0.3)

  
  fig.set_title(clustername, fontsize=lfont+1)

  fig.recenter(radec[0], radec[1], width=width[0], height=width[1])
  
  fig.add_scalebar(0.5*scalebarlengthdeg,"500 kpc",color="white",corner="bottom",linewidth=1.5,fontsize=lfont+1)

  #fig.show_circles(radec[0],radec[1],scalebarlengthdeg,edgecolor="white",linestyle="--",linewidth=0.75)
  fig.show_markers(radec[0],radec[1],marker='x',s=50,facecolor='white')

  hdulist = fits.open(radioimage) 
  imagenoise = findrms(np.ndarray.flatten(hdulist[0].data))
  hdulist.close()     

  if plotcode == 'A1940' or plotcode == 'PSZ2G151':
    imagenoise = imagenoise*2.0


  fig.add_label(0.03, 0.95, '$\sigma_{rms}$=' + str(np.int(imagenoise*1e6)) + ' $\mu$Jy beam$^{-1}$',\
                  relative=True, size=tfont,color='white',horizontalalignment='left')
  
  fig.ticks.set_xspacing(xs)
  lev_factor = 3.
  levs = np.sqrt([1.,4.,16.,64,256,1024])
  levelsr = np.ndarray.tolist(lev_factor*imagenoise*levs)
  fig.show_contour(radioimage,slices=[0,0],levels=levelsr, colors='white', smooth=None, overlap=True, linewidths=0.5)
  
  levelsr = np.ndarray.tolist (-3.*imagenoise*np.array([1.0,]))
  fig.show_contour(radioimage, slices=[0,0],levels=levelsr, colors='red',overlap=True, linewidths=0.5, linestyles='solid'), #linestyles='dashed')

  
  
  
  fig.add_beam()
  fig.beam.set_frame(True)
  
  fig.save(imagename+"_"+filters+"_opt_radiocontours.pdf", dpi=dpi)
  fig.save(imagename+"_"+filters+"_opt_radiocontours.png", dpi=dpi)
  return




prefixpath= '/disks/paradata/shimwell/LoTSS-DR2/archive_cluster_images/'


if True:

  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G151.62+54.78/PSZ2G151.62+54.78_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G151.62+54.78/PSZ2G151.62+54.78_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G151.62+54.78/PSZ2G151.62+54.78_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'PSZ2G151.62+54.78/' + 'PSZ2G151.62+54.78' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri',  0.4864, [163.71686, +55.35387],'PSZ2 G151.62+54.78','figPSZ2G151.62+54.78',Mpcwidth=[1.0,1.0], xs=0.25/6, plotcode='PSZ2G151')

  sys.exit()  
  
  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/Abell1156/Abell1156_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/Abell1156/Abell1156_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/Abell1156/Abell1156_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'Abell1156/' + 'Abell1156' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.209, [166.23149,+47.42078],'Abell 1156','figAbell1156',Mpcwidth=[1.5,1.5], xs=0.25/3,  lfont=24, tfont=22)
     

  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G089.52+62.34/PSZ2G089.52+62.34_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G089.52+62.34/PSZ2G089.52+62.34_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G089.52+62.34/PSZ2G089.52+62.34_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'PSZ2G089.52+62.34/' + 'PSZ2G089.52+62.34' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.0701,  [215.55496, +48.49864],'PSZ2 G089.52+62.34 / Abell 1904','figPSZ2G089.52+62.34',Mpcwidth=[1.,1.], xs=0.25/3)


  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G096.14+56.24/PSZ2G096.14+56.24_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G096.14+56.24/PSZ2G096.14+56.24_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G096.14+56.24/PSZ2G096.14+56.24_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'PSZ2G096.14+56.24/' + 'PSZ2G096.14+56.24' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.139773, [218.86564,+55.12877],'PSZ2 G096.14+56.24 / Abell 1940','figPSZ2G096.14+56.24',Mpcwidth=[0.7,0.7], xs=0.25/6, plotcode='A1940')
       

  
  redim = '/net/lofar7/data1/botteon/DR2_sample/PSZ2G106.61+66.71/PSZ2G106.61+66.71_mosaic_ifilter.fits.smooth'
  greenim = '/net/lofar7/data1/botteon/DR2_sample/PSZ2G106.61+66.71/PSZ2G106.61+66.71_mosaic_rfilter.fits.smooth'
  blueim = '/net/lofar7/data1/botteon/DR2_sample/PSZ2G106.61+66.71/PSZ2G106.61+66.71_mosaic_gfilter.fits.smooth'
  radioimage = '/net/lofar7/data1/botteon/DR2_sample/PSZ2G106.61+66.71/' + 'PSZ2G106.61+66.71' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.331400, [202.62266667, +49.14669444],'PSZ2 G106.61+66.71','figPSZ2G106.61+66.71',Mpcwidth=[2.0,2.0], lfont=24, tfont=22, xs=0.25/3)



  cluster = 'WHLJ122418.6+490549'
  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.1004,[186.07811, +49.09734],'WHL J122418.6+490549','fig' + cluster, Mpcwidth=[0.8,0.8], lfont=24, tfont=22)

 

  


  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G111.75+70.37/PSZ2G111.75+70.37_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G111.75+70.37/PSZ2G111.75+70.37_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G111.75+70.37/PSZ2G111.75+70.37_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'PSZ2G111.75+70.37/' + 'PSZ2G111.75+70.37' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.1830,[198.28288,+46.27711],'PSZ2G111.75+70.37 / Abell 1697','figPSZ2G111.75+70.37',Mpcwidth=[2.0,2.0], xs=0.25/3) # lfont=24, tfont=22, 
 

  cluster = 'Abell1615'
  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster+'/'+cluster+'_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster+'/'+cluster+'_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster+'/'+cluster+'_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri',0.2106, [191.92974,+48.86555],'Abell 1615','fig' + cluster, Mpcwidth=[0.9,0.9], xs=0.25/6)


  redim = '/net/bovenrijn/data1/digennaro/HighRedshiftClusters/OPTICAL/PSZ2G088.98+55.07/PSZ2G088.98+55.07_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HighRedshiftClusters/OPTICAL/PSZ2G088.98+55.07/PSZ2G088.98+55.07_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HighRedshiftClusters/OPTICAL/PSZ2G088.98+55.07/PSZ2G088.98+55.07_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'PSZ2G088.98+55.07/' + 'PSZ2G088.98+55.07' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.702346, [224.7450589,+52.8173732],'PSZ2G088.98+55.07','figPSZ2G088.98+55.07',Mpcwidth=[1.0,1.0], xs=0.25/12)


  redim = '/net/bovenrijn/data1/digennaro/HighRedshiftClusters/OPTICAL/PSZ2G087.39+50.92/PSZ2G087.39+50.92_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HighRedshiftClusters/OPTICAL/PSZ2G087.39+50.92/PSZ2G087.39+50.92_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HighRedshiftClusters/OPTICAL/PSZ2G087.39+50.92/PSZ2G087.39+50.92_mosaic_gfilter.fits.smooth'
  radioimage = '/net/bovenrijn/data1/digennaro/HighRedshiftClusters/LOFAR/' + 'PSZ2G087.39+50.92/' + 'PSZ2G087.39+50.92' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.748, [231.63821,+54.15206],'PSZ2 G087.39+50.92','figPSZ2G087.39+50.92',Mpcwidth=[0.8,0.8], xs=0.25/12)

  
  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G098.44+56.59/PSZ2G098.44+56.59_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G098.44+56.59/PSZ2G098.44+56.59_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G098.44+56.59/PSZ2G098.44+56.59_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'PSZ2G098.44+56.59/' + 'PSZ2G098.44+56.59' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri',  0.1318, [216.852083,+55.750528],'PSZ2 G098.44+56.59 / Abell 1920','figPSZ2G098.44+56.59',Mpcwidth=[2.6,2.6])
   

  
  
  cluster = 'MaxBCGJ173.04772+47.81041'
  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.2261,[173.04772,+47.81041],'MaxBCG J173.04772+47.81041','fig' + cluster, Mpcwidth=[1.2,1.2])
  

  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/Abell1314/Abell1314_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/Abell1314/Abell1314_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/Abell1314/Abell1314_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'Abell1314/' + 'Abell1314' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.0335, [173.70601,+49.07690],'Abell 1314','figAbell1314',Mpcwidth=[1.,1.], lfont=24, tfont=22, xs=0.25)

  cluster = 'WHLJ132226.8+464630'
  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.3718, [200.61109,+46.77494],  'WHL J132226.8+464630','fig' + cluster, Mpcwidth=[1.3,1.3], xs=0.25/6)


  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/Abell1330/Abell1330_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/Abell1330/Abell1330_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/Abell1330/Abell1330_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'Abell1330/' + 'Abell1330' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri',0.2805, [174.56139,+49.54541],'Abell 1330','figAbell1330',Mpcwidth=[1.3,1.3])

  cluster = 'GMBCGJ211.77332+55.09968'
  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri',  0.2506, [211.72963,+55.06747],'GMBCG J211.77332+55.09968','fig' + cluster, Mpcwidth=[2.5,2.5])




  cluster = 'WHLJ124143.1+490510'
  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.3707, [190.45624, +49.07807], 'WHL J124143.1+490510','fig' + cluster, Mpcwidth=[2.5,2.5])


  redim = '/net/bovenrijn/data1/digennaro/HighRedshiftClusters/OPTICAL/PSZ2G086.93+53.18/PSZ2G086.93+53.18_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HighRedshiftClusters/OPTICAL/PSZ2G086.93+53.18/PSZ2G086.93+53.18_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HighRedshiftClusters/OPTICAL/PSZ2G086.93+53.18/PSZ2G086.93+53.18_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'PSZ2G086.93+53.18/' + 'PSZ2G086.93+53.18' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri',  0.6752, [228.5024535, +52.80381756],'PSZ2G086.93+53.18','figPSZ2G086.93+53.18',Mpcwidth=[1.0,1.0],  lfont=24, tfont=22,  xs=0.25/12)

  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G118.34+68.79/PSZ2G118.34+68.79_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G118.34+68.79/PSZ2G118.34+68.79_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G118.34+68.79/PSZ2G118.34+68.79_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'PSZ2G118.34+68.79/' + 'PSZ2G118.34+68.79' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri',  0.254879, [195.36613, +48.25227],'PSZ2 G118.34+68.79 / ZwCl 1259.0+4830','figPSZ2G118.34+68.79',Mpcwidth=[1.5,1.5],lfont=24, tfont=22,  xs=0.25/3)

  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G144.33+62.85/PSZ2G144.33+62.85_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G144.33+62.85/PSZ2G144.33+62.85_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G144.33+62.85/PSZ2G144.33+62.85_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'PSZ2G144.33+62.85/' + 'PSZ2G144.33+62.85' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri',  0.1320, [177.27167,+51.58205],'PSZ2 G144.33+62.85 / Abell 1387','figPSZ2G144.33+62.85',Mpcwidth=[1.5,1.5], xs=0.25/3)
 

  cluster ='WHLJ133936.0+484859'
  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri',0.3265, [204.90010,+48.81644], 'WHL J133936.0+484859','fig' + cluster, Mpcwidth=[1.5,1.5], xs=0.25/6)


  cluster =  'WHLJ125836.8+440111'
  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri',  0.5339, [194.65371, +44.01997],'WHL J125836.8+440111','fig' + cluster, Mpcwidth=[1.5,1.5], xs=0.25/6)
    


  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G136.92+59.46/PSZ2G136.92+59.46_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G136.92+59.46/PSZ2G136.92+59.46_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G136.92+59.46/PSZ2G136.92+59.46_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'PSZ2G136.92+59.46/' + 'PSZ2G136.92+59.46' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri',  0.065, [180.06851,56.26328],'PSZ2 G136.92+59.46 / Abell 1436','figPSZ2G136.92+59.46',Mpcwidth=[0.55,0.55], xs=0.25/6)



  cluster = 'WHLJ132615.8+485229'
  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.2800, [201.56565,+48.87447], 'WHL J132615.8+485229','fig' + cluster, Mpcwidth=[2.,2.])

  cluster = 'WHLJ134746.8+475214'
  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri',0.1695, [206.97006,+47.87677], 'WHL J134746.8+475214','fig' + cluster, Mpcwidth=[1.2,1.2])


  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G123.66+67.25/PSZ2G123.66+67.25_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G123.66+67.25/PSZ2G123.66+67.25_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G123.66+67.25/PSZ2G123.66+67.25_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'PSZ2G123.66+67.25/' + 'PSZ2G123.66+67.25' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri',  0.2838, [192.42231,+49.87176],'PSZ2G123.66+67.25 / Abell 1622','figPSZ2G123.66+67.25',Mpcwidth=[2.5,2.5], lfont=24, tfont=22)




  redim = '/net/bovenrijn/data1/digennaro/HighRedshiftClusters/OPTICAL/PSZ2G084.10+58.72/PSZ2G084.10+58.72_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HighRedshiftClusters/OPTICAL/PSZ2G084.10+58.72/PSZ2G084.10+58.72_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HighRedshiftClusters/OPTICAL/PSZ2G084.10+58.72/PSZ2G084.10+58.72_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'PSZ2G084.10+58.72/' + 'PSZ2G084.10+58.72' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri',  0.731000, [222.25533, +48.55664],'PSZ2 G084.10+58.72','figPSZ2G084.10+58.72',Mpcwidth=[1.0,1.0],  lfont=24, tfont=22, xs=0.25/12)
  
  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G156.26+59.64/PSZ2G156.26+59.64_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G156.26+59.64/PSZ2G156.26+59.64_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G156.26+59.64/PSZ2G156.26+59.64_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'PSZ2G156.26+59.64/' + 'PSZ2G156.26+59.64' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri',  0.5877, [167.12487,50.26757],'PSZ2 G156.26+59.64','figPSZ2G156.26+59.64',Mpcwidth=[1.0,1.0],lfont=24, tfont=22,xs=0.25/12)
    

  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G150.56+58.32/PSZ2G150.56+58.32_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G150.56+58.32/PSZ2G150.56+58.32_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G150.56+58.32/PSZ2G150.56+58.32_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'PSZ2G150.56+58.32/' + 'PSZ2G150.56+58.32' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri',  0.470, [168.8165798,53.33216888],'PSZ2 G150.56+58.32 / MACS J1115.2+5320','figPSZ2G150.56+58.32',Mpcwidth=[1.5,1.5], lfont=24, tfont=22, xs=0.25/6)



  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G114.31+64.89/PSZ2G114.31+64.89_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G114.31+64.89/PSZ2G114.31+64.89_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G114.31+64.89/PSZ2G114.31+64.89_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'PSZ2G114.31+64.89/' + 'PSZ2G114.31+64.89' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.2836,  [198.77156,+51.81702],'PSZ2G114.31+64.89 / Abell 1703','figPSZ2G114.31+64.89',Mpcwidth=[1.2,1.2], lfont=24, tfont=22, xs=0.25/6)

  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G135.17+65.43/PSZ2G135.17+65.43_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G135.17+65.43/PSZ2G135.17+65.43_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G135.17+65.43/PSZ2G135.17+65.43_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'PSZ2G135.17+65.43/' + 'PSZ2G135.17+65.43' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri',  0.5436, [184.7935891,50.90810997],'PSZ2 G135.17+65.43','figPSZ2G135.17+65.43',Mpcwidth=[2.5,2.5], lfont=24, tfont=22)



  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G114.99+70.36/PSZ2G114.99+70.36_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G114.99+70.36/PSZ2G114.99+70.36_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G114.99+70.36/PSZ2G114.99+70.36_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'PSZ2G114.99+70.36/' + 'PSZ2G114.99+70.36' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri',0.2259,   [196.7076514,+46.55760786],'PSZ2G114.99+70.36 / Abell 1682','figPSZ2G114.99+70.36',Mpcwidth=[1.5,1.5],lfont=24, tfont=22)


  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G107.10+65.32/PSZ2G107.10+65.32_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G107.10+65.32/PSZ2G107.10+65.32_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G107.10+65.32/PSZ2G107.10+65.32_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'PSZ2G107.10+65.32/' + 'PSZ2G107.10+65.32' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.2799,  [203.1765617, 50.47508787],'PSZ2 G107.10+65.32 / Abell 1758','figPSZ2G107.10+65.32',Mpcwidth=[3.7,3.7],  lfont=24, tfont=22)




  

  redim = '/net/bovenrijn/data1/digennaro/HighRedshiftClusters/OPTICAL/PSZ2G099.86+58.45/PSZ2G099.86+58.45_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HighRedshiftClusters/OPTICAL/PSZ2G099.86+58.45/PSZ2G099.86+58.45_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HighRedshiftClusters/OPTICAL/PSZ2G099.86+58.45/PSZ2G099.86+58.45_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'PSZ2G099.86+58.45/' + 'PSZ2G099.86+58.45' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.6305,  [213.69655,+54.78442],'PSZ2 G099.86+58.45','figPSZ2G099.86+58.45',Mpcwidth=[0.9,0.9], lfont=24, tfont=22)
  
    
    

  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G145.65+59.30/PSZ2G145.65+59.30_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G145.65+59.30/PSZ2G145.65+59.30_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G145.65+59.30/PSZ2G145.65+59.30_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'PSZ2G145.65+59.30/' + 'PSZ2G145.65+59.30' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.3475,  [173.1767513,54.21999516],'PSZ2 G145.65+59.30 / Abell 1294','figPSZ2G145.65+59.30',Mpcwidth=[2.2,2.2],lfont=24, tfont=22)


  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G095.22+67.41/PSZ2G095.22+67.41_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G095.22+67.41/PSZ2G095.22+67.41_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G095.22+67.41/PSZ2G095.22+67.41_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'PSZ2G095.22+67.41/' + 'PSZ2G095.22+67.41' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri',  0.062500,  [207.92681, +46.36694],'PSZ2 G095.22+67.41','figPSZ2G095.22+67.41',Mpcwidth=[2.,2.])



  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G143.26+65.24/PSZ2G143.26+65.24_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G143.26+65.24/PSZ2G143.26+65.24_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G143.26+65.24/PSZ2G143.26+65.24_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'PSZ2G143.26+65.24/' + 'PSZ2G143.26+65.24' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri',0.3634, [179.84610,+49.79732],'PSZ2 G143.26+65.24 / Abell 1430','figPSZ2G143.26+65.24',Mpcwidth=[2.2,2.2])
 
  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G133.60+69.04/PSZ2G133.60+69.04_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G133.60+69.04/PSZ2G133.60+69.04_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G133.60+69.04/PSZ2G133.60+69.04_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'PSZ2G133.60+69.04/' + 'PSZ2G133.60+69.04' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri',0.2540, [187.2604866,+47.62249881],'PSZ2 G133.60+69.04 / Abell 1550','figPSZ2G133.60+69.04',Mpcwidth=[2.2,2.2])


  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G080.16+57.65/PSZ2G080.16+57.65_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G080.16+57.65/PSZ2G080.16+57.65_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/PSZ2G080.16+57.65/PSZ2G080.16+57.65_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'PSZ2G080.16+57.65/' + 'PSZ2G080.16+57.65' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.0878, [225.29735,+47.27400],'PSZ2 G080.16+57.65 / Abell 2018','figPSZ2G080.16+57.65',Mpcwidth=[0.7,0.7], lfont=24, tfont=22)
    



  

  cluster = 'NSCJ143825+463744'
  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri',  0.03586,[219.69176,+46.66216],'NSC J143825+463744','fig' + cluster, Mpcwidth=[1.0,1.0])
 
  

  cluster =  'WHLJ151310.6+540116'
  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri',  0.4366, [228.29890, +54.02127],'WHL J151310.6+540116','fig' + cluster, Mpcwidth=[1.5,1.5])




  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/Abell1291/Abell1291_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/Abell1291/Abell1291_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/Abell1291/Abell1291_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'Abell1291/' + 'Abell1291' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri',0.0527, [173.09558, +55.96818],'Abell 1291','figAbell1291',Mpcwidth=[0.9,0.9])


  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/RXCJ1332.6+5419/RXCJ1332.6+5419_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/RXCJ1332.6+5419/RXCJ1332.6+5419_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/RXCJ1332.6+5419/RXCJ1332.6+5419_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'RXCJ1332.6+5419/' + 'RXCJ1332.6+5419' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.10655, [203.15183,+54.31651],'RXC J1332.6+5419','figRXCJ1332.6+5419',Mpcwidth=[1.5,1.5])




    
  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/RXCJ1053.7+5452/RXCJ1053.7+5452_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/RXCJ1053.7+5452/RXCJ1053.7+5452_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/RXCJ1053.7+5452/RXCJ1053.7+5452_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + 'RXCJ1053.7+5452/' + 'RXCJ1053.7+5452' + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri',  0.0704, [163.43500, 54.87250],'RXC J1053.7+5452','figRXCJ1053.7+5452',Mpcwidth=[1.65,1.65])

  

    

    





  cluster =  'WHLJ151001.6+535722'
  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.4923, [227.53222,+53.95423],'WHL J151001.6+535722','fig' + cluster, Mpcwidth=[2.,2.])




  sys.exit()




  cluster = 'WHLJ131527.6+484025'
  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.5158,[198.86497,+48.67361], 'WHL J131527.6+484025','fig' + cluster, Mpcwidth=[2.5,2.5])


  cluster =  'WHLJ125746.3+485446'
  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri',0.2636,[194.44273,+48.91299],'WHL J125746.3+485446','fig' + cluster, Mpcwidth=[1.2,1.2])


  cluster =  'WHLJ125334.1+505414'
  redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  makeopticalfigure(redim,greenim,blueim,radioimage,'gri',0.1226, [193.39217,+50.90403],'WHL J125334.1+505414','fig' + cluster, Mpcwidth=[1.2,1.2])



    
  sys.exit()    
  #cluster = 'WHLJ124833.0+554820'
  #redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  #greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  #blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  #radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  #makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.27,  [192.13760,+55.80526],'WHL J124833.0+554820','fig' + cluster, Mpcwidth=[1.75,1.75])



  #cluster = 'WHLJ122816.1+495021'
  #redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  #greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  #blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  #radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  #makeopticalfigure(redim,greenim,blueim,radioimage,'gri',0.2637, [187.06727,+49.83951], 'WHL J122816.1+495021','fig' + cluster, Mpcwidth=[1.5,1.5])


    
  #cluster = 'NSCSJ123842+553825'
  #redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  #greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  #blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  #radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  #makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.27834,[189.71596,+55.64610],'NSCS J123842+553825','fig' + cluster, Mpcwidth=[2.1,2.1])


  #cluster = 'GMBCGJ205.93744+49.75911'
  #redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  #greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  #blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  #radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  #makeopticalfigure(redim,greenim,blueim,radioimage,'gri',  0.271 , [205.91115,+49.77400], 'GMBCG J205.93744+49.75911','fig' + cluster, Mpcwidth=[2.0,2.0])



  #cluster = 'SDSSJ121057.53+553005.8'
  #redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  #greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  #blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  #radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  #makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.3366 ,[182.73972,+55.50167], 'SDSS J121057.53+553005.8','fig' + cluster, Mpcwidth=[1.1,1.1])


  #cluster = 'SDSS-C4-DR33106'
  #redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  #greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  #blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  #radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  #makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.0603 ,[177.54706,+53.72236],'SDSS-C4-DR3 3106','fig' + cluster, Mpcwidth=[1.5,1.5])


  #cluster = 'GMBCGJ181.88181+52.89922'
  #redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  #greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  #blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  #radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  #makeopticalfigure(redim,greenim,blueim,radioimage,'gri',  0.2759,[181.90365,+52.91709],'GMBCG J181.88181+52.89922','fig' + cluster, Mpcwidth=[1.5,1.5])



  #cluster = 'ZwCl1127.5+4804'
  #redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  #greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  #blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  #radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  #makeopticalfigure(redim,greenim,blueim,radioimage,'gri',  0.1267,[172.56418,+47.79317],'ZwCl 1127.5+4804','fig' + cluster, Mpcwidth=[1.5,1.5])

    






  #cluster = 'MCXCJ1221.4+4918'
  #redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  #greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  #blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  #radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  #makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.700,  [185.36489,+49.30730],'MCXC J1221.4+4918','fig' + cluster, Mpcwidth=[1.5,1.5])


  #cluster = 'MCXCJ1217.7+4729'
  #redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  #greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  #blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  #radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  #makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.2700,   [184.42961,+47.48647],'MCXC J1217.7+4729','fig' + cluster, Mpcwidth=[1.5,1.5])


  #cluster = 'MCXCJ1053.7+4929'
  #redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  #greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  #blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  #radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  #makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.1400,  [163.43380,+49.49889],'MCXC J1053.7+4929','fig' + cluster, Mpcwidth=[1.5,1.5])


  #cluster = 'ClGJ120958.9+495352'
  #redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  #greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  #blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  #radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  #makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.9020, [182.49625, +49.89769],'ClG J120958.9+495352','fig' + cluster, Mpcwidth=[1.5,1.5])




  #cluster = 'Abell2011'
  #redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  #greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  #blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  #radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  #makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.1697, [224.94504, +49.76864],'Abell 2011','fig' + cluster, Mpcwidth=[1.5,1.5])


  #cluster = 'Abell2000'
  #redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  #greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  #blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  #radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  #makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.1012, [223.59811,+54.32444],'Abell 2000 / Abell 1999','fig' + cluster, Mpcwidth=[2.,2.])


  #cluster = 'Abell1855'
  #redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  #greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  #blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  #radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  #makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.2379, [211.21535,+47.08447],'Abell 1855','fig' + cluster, Mpcwidth=[2.5,2.5])

  #cluster = 'Abell1804'
  #redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  #greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  #blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  #radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  #makeopticalfigure(redim,greenim,blueim,radioimage,'gri',  0.1665, [207.25215, +49.31121],'Abell 1804','fig' + cluster, Mpcwidth=[1.6,1.6])
  
  #cluster = 'Abell1788'
  #redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  #greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  #blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  #radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  #makeopticalfigure(redim,greenim,blueim,radioimage,'gri',  0.1637, [206.23556,+53.75058],'Abell 1788','fig' + cluster, Mpcwidth=[2.0,2.0])
  


  #cluster = 'Abell1745'
  #redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'  + cluster +'/' + cluster +'_mosaic_ifilter.fits.smooth'
  #greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster  +'/' + cluster +'_mosaic_rfilter.fits.smooth'
  #blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/' +cluster  +'/' + cluster +'_mosaic_gfilter.fits.smooth'
  #radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  #makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.3665,  [201.69882,+53.82430],'Abell 1745','fig' + cluster, Mpcwidth=[2.0,2.0])
  

  #cluster = 'Abell1643'
  #redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster+'/'+cluster+'_mosaic_ifilter.fits.smooth'
  #greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster+'/'+cluster+'_mosaic_rfilter.fits.smooth'
  #blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/'+cluster+'/'+cluster+'_mosaic_gfilter.fits.smooth'
  #radioimage = prefixpath + cluster + '/' + cluster + '_maskROBUST-0.5-MFS-image.fits'
  #makeopticalfigure(redim,greenim,blueim,radioimage,'gri',0.2328, [193.95482,+44.08842],'Abell 1643','fig' + cluster, Mpcwidth=[1.1,1.1])





  #redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/Abell1368/Abell1368_mosaic_ifilter.fits.smooth'
  #greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/Abell1368/Abell1368_mosaic_rfilter.fits.smooth'
  #blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/Abell1368/Abell1368_mosaic_gfilter.fits.smooth'
  #radioimage = prefixpath + 'Abell1368/' + 'Abell1368' + '_maskROBUST-0.5-MFS-image.fits'
  #makeopticalfigure(redim,greenim,blueim,radioimage,'gri',0.1291, [176.22393,+51.26713],'Abell 1368','figAbell1368',Mpcwidth=[1.1,1.1])


  #redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/Abell1542/Abell1542_mosaic_ifilter.fits.smooth'
  #greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/Abell1542/Abell1542_mosaic_rfilter.fits.smooth'
  #blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/Abell1542/Abell1542_mosaic_gfilter.fits.smooth'
  #radioimage = prefixpath + 'Abell1542/' + 'Abell1542' + '_maskROBUST-0.5-MFS-image.fits'
  #makeopticalfigure(redim,greenim,blueim,radioimage,'gri',0.1218, [186.88430,+49.47882],'Abell 1542','figAbell1542',Mpcwidth=[1.5,1.5])


  #redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/Abell1377/Abell1377_mosaic_ifilter.fits.smooth'
  #greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/Abell1377/Abell1377_mosaic_rfilter.fits.smooth'
  #blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/Abell1377/Abell1377_mosaic_gfilter.fits.smooth'
  #radioimage = prefixpath + 'Abell1377/' + 'Abell1377' + '_maskROBUST-0.5-MFS-image.fits'
  #makeopticalfigure(redim,greenim,blueim,radioimage,'gri',0.0514, [176.83938,+55.72973],'Abell 1377','figAbell1377',Mpcwidth=[1.5,1.5])





  #redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/Abell1261/Abell1261_mosaic_ifilter.fits.smooth'
  #greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/Abell1261/Abell1261_mosaic_rfilter.fits.smooth'
  #blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/Abell1261/Abell1261_mosaic_gfilter.fits.smooth'
  #radioimage = prefixpath + 'Abell1261/' + 'Abell1261' + '_maskROBUST-0.5-MFS-image.fits'
  #makeopticalfigure(redim,greenim,blueim,radioimage,'gri',0.1562, [171.83958,+48.29855],'Abell 1261','figAbell1261',Mpcwidth=[1.1,1.1])

  #redim = '/net/bovenrijn/data1/digennaro/HetdexPaper/Abell1202/Abell1202_mosaic_ifilter.fits.smooth'
  #greenim = '/net/bovenrijn/data1/digennaro/HetdexPaper/Abell1202/Abell1202_mosaic_rfilter.fits.smooth'
  #blueim = '/net/bovenrijn/data1/digennaro/HetdexPaper/Abell1202/Abell1202_mosaic_gfilter.fits.smooth'
  #radioimage = prefixpath + 'Abell1202/' + 'Abell1202' + '_maskROBUST-0.5-MFS-image.fits'
  #makeopticalfigure(redim,greenim,blueim,radioimage,'gri', 0.1121, [168.40705,+47.48700],'Abell 1202','figAbell1202',Mpcwidth=[2.0,2.0])






    
   








    





if False:
    ff = 'pdf'

    radioimage = prefixpath + 'PSZ2G089.52+62.34/' + 'PSZ2G089.52+62.34' + '_maskROBUST-0.5TAPER10-MFS-image.fits'
    Xray_radio('Xrays/PSZ2G089.52+62.34_XMM.fits.gz', radioimage, 'PSZ2G089.52+62.34_xray.' + ff,  0.0701, [215.55496, +48.49864], Mpcwidth=[1.5,1.5], \
            titlename = 'PSZ2 G089.52+62.34 / Abell 1904',vmin=6e-6,vmax=3e-5,smooth=3, pad=True)   


    radioimage = prefixpath + 'PSZ2G133.60+69.04/' + 'PSZ2G133.60+69.04' + '_maskROBUST-0.5TAPER10-MFS-image.fits'
    Xray_radio('Xrays/PSZ2G133.60+69.04_CHANDRA.img.gz', radioimage, 'PSZ2G133.60+69.04_xray.' + ff, 0.2540, [187.2604866,+47.62249881], Mpcwidth=[2.5,2.5], \
            titlename = 'PSZ2 G133.60+69.04 / Abell 1550',vmin=0.01,vmax=0.1,smooth=5,pad=True)


    radioimage = prefixpath + 'PSZ2G145.65+59.30/' + 'PSZ2G145.65+59.30' + '_maskROBUST-0.5TAPER15-MFS-image.fits'
    Xray_radio('Xrays/PSZ2G145.65+59.30_XMM.fits.gz', radioimage, 'PSZ2G145.65+59.30_xray.' + ff,  0.3475, [173.1767513,54.21999516], Mpcwidth=[2.0,2.0], \
            titlename = 'PSZ2 G145.65+59.30 / Abell 1294',vmin=2e-6,vmax=30e-5,smooth=1, lfont=24, tfont=22)    
    
    radioimage = prefixpath + 'PSZ2G111.75+70.37/' + 'PSZ2G111.75+70.37' + '_masksubROBUST-0.5TAPER30-MFS-image.fits'
    Xray_radio('Xrays/mosaic_a1697_flux.fits.gz', radioimage, 'PSZ2G111.75+70.37_xray.' + ff, 0.1830,[198.28288,+46.27711], Mpcwidth=[2.9,2.9], \
            titlename = 'PSZ2G111.75+70.37 / Abell 1697',vmin=2e-6,vmax=5e-5,smooth=3)

    radioimage = prefixpath + 'PSZ2G111.75+70.37/' + 'PSZ2G111.75+70.37' + '_masksubROBUST-0.5TAPER30-MFS-image.fits'
    Xray_radio('Xrays/a1697_0.5-2.0_flux.img.gz', radioimage, 'PSZ2G111.75+70.37_xray.' + ff, 0.1830,[198.28288,+46.27711], Mpcwidth=[2.3,2.3], \
            titlename = 'PSZ2G111.75+70.37 / Abell 1697',vmin=5e-9,vmax=1.5e-7,smooth=5)

    radioimage = prefixpath + 'Abell1156/' + 'Abell1156' + '_masksubROBUST-0.5TAPER30-MFS-image.fits'
    Xray_radio('Xrays/a1156_0.5-2.0_flux.img.gz', radioimage, 'Abell1156_xray.' + ff, 0.209, [166.23149,+47.42078], Mpcwidth=[2.,2.], \
            titlename = 'Abell 1156',vmin=2e-9,vmax=1.5e-7,smooth=5, lfont=24, tfont=22)

    radioimage = '/net/lofar7/data1/botteon/DR2_sample/PSZ2G106.61+66.71/' + 'PSZ2G106.61+66.71' + '_masksubROBUST-0.5TAPER10-MFS-image.fits'
    Xray_radio('Xrays/whlj133029_0.5-2.0_flux.img.gz', radioimage, 'PSZ2G106.61+66.71_xray.' + ff, 0.331400,[202.62266667, +49.14669444], Mpcwidth=[3.1,3.1], \
            titlename = 'PSZ2 G106.61+66.71',vmin=2e-9,vmax=1.e-7,smooth=5, lfont=24, tfont=22)

    radioimage = prefixpath + 'WHLJ122418.6+490549/' + 'WHLJ122418.6+490549' + '_maskROBUST-0.5TAPER10-MFS-image.fits'
    Xray_radio('Xrays/whlj122418_0.5-2.0_flux.img.gz', radioimage, 'WHLJ122418.6+490549_xray.' + ff, 0.1004,[186.07811, +49.09734], Mpcwidth=[1.1,1.1], \
            titlename = 'WHL J122418.6+490549',vmin=2e-9,vmax=5.e-8,smooth=5, lfont=24, tfont=22)

    radioimage = prefixpath + 'PSZ2G143.26+65.24/' + 'PSZ2G143.26+65.24' + '_masksubROBUST-0.5TAPER15-MFS-image.fits'
    Xray_radio('Xrays/PSZ2G143.26+65.24_CHANDRA.img.gz', radioimage, 'PSZ2G143.26+65.24_xray.' + ff, 0.3634, [179.84610,+49.79732], Mpcwidth=[3.0,3.0], \
            titlename = 'PSZ2 G143.26+65.24 / Abell 1430',vmin=1e-9,vmax=2e-7,smooth=5)    

    radioimage = prefixpath + 'PSZ2G135.17+65.43/' + 'PSZ2G135.17+65.43' + '_masksubROBUST-0.5TAPER10-MFS-image.fits'
    Xray_radio('Xrays/PSZ2G135.17+65.43_0.5-2.0_flux.img.gz', radioimage, 'PSZ2G135.17+65.43_xray.' + ff, 0.5436, [184.7935891,50.90810997], Mpcwidth=[3.,3.], \
            titlename = 'PSZ2G 135.17+65.43',vmin=5e-9,vmax=1.5e-7,smooth=3, lfont=24, tfont=22)

    radioimage = prefixpath + 'PSZ2G114.31+64.89/' + 'PSZ2G114.31+64.89' + '_masksubROBUST-0.5TAPER15-MFS-image.fits'
    Xray_radio('Xrays/PSZ2G114.31+64.89_CHANDRA.img.gz', radioimage, 'PSZ2G114.31+64.89_xray.' + ff, 0.2836,  [198.77156,+51.81702], Mpcwidth=[2.3,2.3], \
            titlename = 'PSZ2 G114.31+64.89 / Abell 1703',vmin=1e-9,vmax=2e-7,smooth=3, lfont=24, tfont=22)


    radioimage = prefixpath + 'PSZ2G123.66+67.25/' + 'PSZ2G123.66+67.25' + '_maskROBUST-0.5TAPER10-MFS-image.fits'
    Xray_radio('Xrays/PSZ2G123.66+67.25_CHANDRA.img.gz', radioimage, 'PSZ2G123.66+67.25_xray.' + ff,  0.2838, [192.42231,+49.87176], Mpcwidth=[3.0,3.0], \
            titlename = 'PSZ2 G123.66+67.25 / Abell 1622',vmin=1e-9,vmax=2e-7,smooth=5, lfont=24, tfont=22)


    radioimage = prefixpath + 'Abell1314/' + 'Abell1314' + '_maskROBUST-0.5TAPER10-MFS-image.fits'
    Xray_radio('Xrays/Abell1314_XMM.fits.gz', radioimage, 'Abell1314_xray.' + ff, 0.0335, [173.70601,+49.07690], Mpcwidth=[1.3,1.3], \
            titlename = 'Abell 1314',vmin=2.*2e-6,vmax=5e-5,smooth=3,pad=True, xs=0.25,  lfont=24, tfont=22)

    radioimage = prefixpath + 'PSZ2G086.93+53.18/' + 'PSZ2G086.93+53.18' + '_masksubROBUST-0.5TAPER15-MFS-image.fits'
    Xray_radio('Xrays/PSZ2G086.93+53.18_XMM.fits.gz', radioimage, 'PSZ2G086.93+53.18_xray.' + ff, 0.6752, [228.5024535, +52.80381756], Mpcwidth=[2.0,2.0], \
            titlename = 'PSZ2 G086.93+53.18',vmin=2e-6,vmax=5e-5,smooth=1,  lfont=24, tfont=22, xs=0.25/6)

    radioimage = prefixpath + 'PSZ2G084.10+58.72/' + 'PSZ2G084.10+58.72' + '_masksubROBUST-0.5TAPER15-MFS-image.fits'
    Xray_radio('Xrays/PSZ2G084.10+58.72_XMM.fits.gz', radioimage, 'PSZ2G084.10+58.72_xray.' + ff, 0.731000, [222.25533, +48.55664], Mpcwidth=[2.0,2.0], \
            titlename = 'PSZ2 G084.10+58.72',vmin=2e-6,vmax=5e-5,smooth=1,  lfont=24, tfont=22, xs=0.25/6)


    radioimage = prefixpath + 'RXCJ1053.7+5452/' + 'RXCJ1053.7+5452' + '_maskROBUST-0.5TAPER10-MFS-image.fits'
    Xray_radio('Xrays/RXCJ1053.7+5452_CHANDRA.img.gz', radioimage, 'RXCJ1053.7+5452_xray.' + ff,  0.0704, [163.43500, 54.87250], Mpcwidth=[1.65,1.65], \
            titlename = 'RXC J1053.7+5452',vmin=1e-9,vmax=2e-7,smooth=5,pad=True, xs=0.25)    

    
    radioimage = prefixpath + 'PSZ2G099.86+58.45/' + 'PSZ2G099.86+58.45' + '_masksubROBUST-0.5TAPER30-MFS-image.fits'
    Xray_radio('Xrays/PSZ2G099.86+58.45_XMM.fits.gz', radioimage, 'PSZ2G099.86+58.45_xray.' + ff, 0.6305, [213.69655,+54.78442], Mpcwidth=[3.0,3.0], \
            titlename = 'PSZ2 G099.86+58.45',vmin=2e-6,vmax=15e-5,smooth=1,  lfont=24, tfont=22)

    radioimage = prefixpath + 'Abell1291/' + 'Abell1291' + '_maskROBUST-0.5TAPER15-MFS-image.fits'
    Xray_radio('Xrays/Abell1291_XMM.fits.gz', radioimage, 'Abell1291_xray.' + ff, 0.0527, [173.09558, +55.96818], Mpcwidth=[1.1,1.1], \
            titlename = 'Abell 1291',vmin=10e-5,vmax=30e-5,smooth=1)      
    
    radioimage = prefixpath + 'PSZ2G096.14+56.24/' + 'PSZ2G096.14+56.24' + '_maskROBUST-0.5TAPER15-MFS-image.fits'
    Xray_radio('Xrays/PSZ2G096.14+56.24_XMM.fits.gz', radioimage, 'PSZ2G096.14+56.24_xray.' + ff,  0.139773, [218.86564,+55.12877], Mpcwidth=[2.,2.], \
            titlename = 'PSZ2 G096.14+56.24 / Abell 1940',vmin=10e-5,vmax=30e-5,smooth=1)  

    radioimage = prefixpath + 'PSZ2G136.92+59.46/' + 'PSZ2G136.92+59.46' + '_maskROBUST-0.5TAPER15-MFS-image.fits'
    Xray_radio('Xrays/PSZ2G136.92+59.46_XMM.fits.gz', radioimage, 'PSZ2G136.92+59.46_xray.' + ff,  0.0665, [180.06851,56.26328], Mpcwidth=[1.5,1.5], \
            titlename = 'PSZ2 G136.92+59.46 / Abell 1436',vmin=2e-6,vmax=30e-5,smooth=1)   

    radioimage = '/net/bovenrijn/data1/digennaro/HighRedshiftClusters/LOFAR/' + 'PSZ2G087.39+50.92/' + 'PSZ2G087.39+50.92_maskROBUST-0.5TAPER10-MFS-image.fits'
    Xray_radio('Xrays/PSZ2G087.39+50.92_XMM.fits.gz', radioimage, 'PSZ2G087.39+50.92_xray.' + ff,0.748, [231.63821,+54.15206], Mpcwidth=[3.0,3.0], \
            titlename = 'PSZ2 G087.39+50.92',vmin=2e-6,vmax=15e-5,smooth=1)

    #radioimage = prefixpath + 'PSZ2G080.16+57.65/' + 'PSZ2G080.16+57.65' + '_maskROBUST-0.5TAPER30-MFS-image.fits'
    radioimage = prefixpath + 'PSZ2G080.16+57.65/' + 'PSZ2G080.16+57.65_subROBUST-0.5TAPER90multiscale-MFS-image.fits'
    Xray_radio('Xrays/PSZ2G080.16+57.65_XMM.fits.gz', radioimage, 'PSZ2G080.16+57.65_xray.' + ff,0.0878, [225.29735,+47.27400], Mpcwidth=[3.0,3.0], \
            titlename = 'PSZ2 G080.16+57.65 / Abell 2018',vmin=2e-6,vmax=15e-5,smooth=1,contourstartlevel=1.5, lfont=24, tfont=22, xs=0.25)

    radioimage = prefixpath + 'PSZ2G150.56+58.32/' + 'PSZ2G150.56+58.32' + '_masksubROBUST-0.5TAPER10-MFS-image.fits'
    Xray_radio('Xrays/PSZ2G150.56+58.32_CHANDRA.img.gz', radioimage, 'PSZ2G150.56+58.32_xray.' + ff, 0.470,   [168.8165798,53.33216888], Mpcwidth=[2.5,2.5], \
            titlename = 'PSZ2 G150.56+58.32 / MACS J1115.2+5320',vmin=1e-9,vmax=2e-7,smooth=3, lfont=24, tfont=22)
  
    radioimage = prefixpath + 'PSZ2G114.99+70.36/' + 'PSZ2G114.99+70.36' + '_maskROBUST-0.5TAPER10-MFS-image.fits'
    Xray_radio('Xrays/PSZ2G114.99+70.36_CHANDRA.img.gz', radioimage, 'PSZ2G114.99+70.36_xray.' + ff, 0.2259,   [196.7076514,+46.55760786], Mpcwidth=[3.0,3.0], \
            titlename = 'PSZ2 G114.99+70.36 / Abell 1682',vmin=1e-9,vmax=2e-7,smooth=5, lfont=24, tfont=22)

    radioimage = prefixpath + 'PSZ2G107.10+65.32/' + 'PSZ2G107.10+65.32' + '_maskROBUST-0.5TAPER10-MFS-image.fits'
    Xray_radio('Xrays/PSZ2G107.10+65.32_XMM.fits.gz', radioimage, 'PSZ2G107.10+65.32_xray.' + ff,  0.2799,  [203.1765617, 50.48], Mpcwidth=[4.3,4.3], \
            titlename = 'PSZ2 G107.10+65.32 / Abell 1758',vmin=3e-6,vmax=15e-5,smooth=1, pad=True,  lfont=24, tfont=22)   
   

    radioimage = prefixpath + 'PSZ2G095.22+67.41/' + 'PSZ2G095.22+67.41' + '_masksubROBUST-0.5TAPER30-MFS-image.fits'
    Xray_radio('Xrays/PSZ2G095.22+67.41_XMM.fits.gz', radioimage, 'PSZ2G095.22+67.41_xray.' + ff,  0.062500, [207.92681, +46.36694], Mpcwidth=[2.0,2.0], \
            titlename = 'PSZ2 G095.22+67.41',vmin=3e-6,vmax=5e-5,smooth=3, pad=True)   


    radioimage = prefixpath + 'PSZ2G088.98+55.07/' + 'PSZ2G088.98+55.07' + '_maskROBUST-0.5TAPER10-MFS-image.fits'
    Xray_radio('Xrays/PSZ2G088.98+55.07_XMM.fits.gz', radioimage, 'PSZ2G088.98+55.07_xray.' + ff, 0.702346, [224.7450589,+52.8173732], Mpcwidth=[1.5,1.5], \
            titlename = 'PSZ2G 088.98+55.07',vmin=2e-6,vmax=5e-5,smooth=1)

    radioimage = prefixpath + 'MCXCJ1221.4+4918/' + 'MCXCJ1221.4+4918' + '_maskROBUST-0.5TAPER10-MFS-image.fits'
    Xray_radio('Xrays/MCXCJ1221.4+4918_CHANDRA.img.gz', radioimage, 'MCXCJ1221.4+4918_xray.' + ff,  0.700, [185.36489,+49.30730], Mpcwidth=[1.5,1.5], \
            titlename = 'MCXC J1221.4+4918',vmin=1e-9,vmax=2e-7,smooth=3)

    radioimage = prefixpath + 'Abell1377/' + 'Abell1377' + '_maskROBUST-0.5TAPER10-MFS-image.fits'
    Xray_radio('Xrays/Abell1377_CHANDRA.img.gz', radioimage, 'Abell1377_xray.' + ff,  0.0514, [176.83938,+55.72973], Mpcwidth=[1.5,1.5], \
            titlename = 'Abell 1377',vmin=1e-9,vmax=2e-7,smooth=5,pad=True)

    radioimage = prefixpath + 'Abell1804/' + 'Abell1804' + '_maskROBUST-0.5TAPER10-MFS-image.fits'
    Xray_radio('Xrays/Abell1804_CHANDRA.img.gz', radioimage, 'Abell1804_xray.' + ff,  0.1665, [207.25215, +49.31121], Mpcwidth=[1.5,1.5], \
            titlename = 'Abell 1804',vmin=1e-9,vmax=2e-7,smooth=5)

    radioimage = prefixpath + 'GMBCGJ181.88181+52.89922/' + 'GMBCGJ181.88181+52.89922' + '_maskROBUST-0.5TAPER10-MFS-image.fits'
    Xray_radio('Xrays/GMBCGJ181.88181+52.89922_CHANDRA.img.gz', radioimage, 'GMBCGJ181.88181+52.89922_xray.' + ff,0.2759, [181.90365,+52.91709], Mpcwidth=[1.5,1.5], \
            titlename = 'GMBCG J181.88181+52.89922',vmin=2e-9,vmax=2e-7,smooth=3,pad=True)

    radioimage = prefixpath + 'Abell2011/' + 'Abell2011' + '_maskROBUST-0.5TAPER10-MFS-image.fits'
    Xray_radio('Xrays/Abell2011_XMM.fits.gz', radioimage, 'Abell2011_xray.' + ff, 0.1697, [224.94504, +49.76864], Mpcwidth=[3.0,3.0], \
            titlename = 'Abell 2011',vmin=2e-6,vmax=5e-5,smooth=1,pad=True)

    radioimage = prefixpath + 'MCXCJ1217.7+4729/' + 'MCXCJ1217.7+4729' + '_maskROBUST-0.5TAPER10-MFS-image.fits'
    Xray_radio('Xrays/MCXCJ1217.7+4729_XMM.fits.gz', radioimage, 'MCXCJ1217.7+4729_xray.' + ff,  0.2700, [184.42961,+47.48647], Mpcwidth=[1.5,1.5], \
            titlename = 'MCXC J1217.7+4729',vmin=15*2e-6,vmax=5*5e-5,smooth=1)

    radioimage = prefixpath + 'ClGJ120958.9+495352/' + 'ClGJ120958.9+495352' + '_maskROBUST-0.5TAPER10-MFS-image.fits'
    Xray_radio('Xrays/ClGJ120958.9+495352_XMM.fits.gz', radioimage, 'ClGJ120958.9+495352_xray.' + ff, 0.9020, [182.49625, +49.89769], Mpcwidth=[2.0,2.0], \
            titlename = 'ClG J120958.9+495352',vmin=2e-6,vmax=30e-5,smooth=1)

    radioimage = prefixpath + 'WHLJ134746.8+475214/' + 'WHLJ134746.8+475214' + '_maskROBUST-0.5TAPER10-MFS-image.fits'
    Xray_radio('Xrays/WHLJ134746.8+475214_XMM.fits.gz', radioimage, 'WHLJ134746.8+475214_xray.' + ff, 0.1695, [206.97006,+47.87677], Mpcwidth=[2.0,2.0], \
            titlename = 'WHL J134746.8+475214',vmin=2e-6,vmax=0.5*5e-5,smooth=3)


    sys.exit()





if False:
  makecomparisonimages()

  selfcalimages()

  selfcalimagesLBA()
  selfcalimagesLB()
  sys.exit()


center    = [215.55496, +48.49864] 
name      = 'PSZ2G089.52+62.34 / Abell 1904'
z         = 0.0701
width     = [1.5,1.5]
makeimagewrapper(prefixpath, 'PSZ2G089.52+62.34', center, width, z, name, rmsfactor=30., plotcode='A1904', xs=0.25)

#sys.exit()

center    = [215.55496, +48.65]  # [215.55496, +48.49864] 
name      = 'PSZ2G089.52+62.34 / Abell 1904'
z         = 0.0701
width     = [2.5,2.5]
makeimagewrapper(prefixpath, 'PSZ2G089.52+62.34', center, width, z, name, rmsfactor=30., plotcode='A1904LF', xs=0.25)





center    = [163.43500, 54.87250]
name      = 'RXC J1053.7+5452'
z         = 0.0704
width     = [1.65,1.65]
makeimagewrapper(prefixpath, 'RXCJ1053.7+5452', center, width, z, name, rmsfactor=30., xs=0.25)



center    = [194.65371, +44.01997]
name      = 'WHL J125836.8+440111'
z         = 0.5339
width     = [2.5,2.5]
makeimagewrapper(prefixpath, 'WHLJ125836.8+440111', center, width, z, name, rmsfactor=30.)


center    = [200.61109,+46.77494]
name      = 'WHL J132226.8+464630'
z         = 0.3718
width     = [2.0,2.0]
makeimagewrapper(prefixpath, 'WHLJ132226.8+464630', center, width, z, name, rmsfactor=30.)

center    = [204.90010,+48.81644]
name      = 'WHL J133936.0+484859'
z         = 0.3265
width     = [1.5,1.5]
makeimagewrapper(prefixpath, 'WHLJ133936.0+484859', center, width, z, name, rmsfactor=30.)


center    = [190.45624, +49.07807]
name      = 'WHL J124143.1+490510'
z         = 0.3707
width     = [2.5,2.5]
makeimagewrapper(prefixpath, 'WHLJ124143.1+490510', center, width, z, name, rmsfactor=30.)



center    = [219.69176,+46.66216]
name      = 'NSC J143825+463744'
z         = 0.03586
width     = [1.1,1.1]
makeimagewrapper(prefixpath, 'NSCJ143825+463744', center, width, z, name, rmsfactor=30., xs=0.25)


center    = [166.23149,+47.42078]
name      = 'Abell 1156'
z         = 0.209
width     = [2.,2.]
makeimagewrapper(prefixpath, 'Abell1156', center, width, z, name, rmsfactor=150.,  lfont=24, tfont=22, pexponent=0.35)

center    = [167.12487,50.26757]
name      = 'PSZ2G156.26+59.64'
z         = 0.5877
width     = [2.,2.]
makeimagewrapper(prefixpath, 'PSZ2G156.26+59.64', center, width, z, name, rmsfactor=30.,lfont=24, tfont=22)



center    = [187.2604866,+47.62249881]
name      = 'PSZ2 G133.60+69.04 / Abell 1550'
z         = 0.2540
width     = [3.0,3.0]
makeimagewrapper(prefixpath, 'PSZ2G133.60+69.04', center, width, z, name, rmsfactor=150., pexponent=0.3, plotcode='A1550')


center    = [179.84610,+49.79732]
name      = 'PSZ2 G143.26+65.24 / Abell 1430'
z         = 0.3634
width     = [3.0,3.0]
makeimagewrapper(prefixpath, 'PSZ2G143.26+65.24', center, width, z, name, rmsfactor=150., pexponent=0.3, plotcode='A1430')



center    = [195.36613, +48.25227]
name      = 'PSZ2 G118.34+68.79 / ZwCl 1259.0+4830'
z         = 0.254879
width     = [2.0,2.0]
makeimagewrapper('/net/lofar7/data1/botteon/DR2_sample/', 'PSZ2G118.34+68.79', center, width, z, name, rmsfactor=30., plotcode='PSZ2G118.34',lfont=24, tfont=22)


center    = [198.28288,+46.27711]
name      = 'PSZ2 G111.75+70.37 / Abell 1697'
z         = 0.1830
width     = [2.0,2.0]
makeimagewrapper(prefixpath, 'PSZ2G111.75+70.37', center, width, z, name, rmsfactor=150., pexponent=0.35) #,lfont=24, tfont=22)






center    = [225.29735,+47.27400]
name      = 'PSZ2 G080.16+57.65 / Abell 2018'
z         = 0.0878
width     = [3.0,3.0]
makeimagewrapper(prefixpath, 'PSZ2G080.16+57.65', center, width, z, name, rmsfactor=30., lfont=24, tfont=22, xs=0.25)



center    = [192.42231,+49.87176]
name      = 'PSZ2 G123.66+67.25 / Abell 1622'
z         = 0.2838
width     = [3.0,3.0]
makeimagewrapper(prefixpath, 'PSZ2G123.66+67.25', center, width, z, name, rmsfactor=150.,  pexponent=0.3,plotcode='A1622', lfont=24, tfont=22)




center    = [186.07811, +49.09734]
name      = 'WHL J122418.6+490549'
z         = 0.1004
width     = [1.1,1.1]
makeimagewrapper(prefixpath, 'WHLJ122418.6+490549', center, width, z, name, rmsfactor=150.,  pexponent=0.3, lfont=24, tfont=22)



center    = [184.7935891,50.90810997]
name      = 'PSZ2 G135.17+65.43'
z         = 0.5436
width     = [3.0,3.0]
makeimagewrapper(prefixpath, 'PSZ2G135.17+65.43', center, width, z, name, rmsfactor=150., pexponent=0.3, plotcode='PSZ2G135.17+65.43', lfont=24, tfont=22)






center    =  [196.7076514,+46.55760786]
name      = 'PSZ2 G114.99+70.36 / Abell 1682'
z         = 0.2259
width     = [3.,3.]
makeimagewrapper(prefixpath, 'PSZ2G114.99+70.36', center, width, z, name, rmsfactor=750., pexponent=0.3, plotcode ='A1682', lfont=24, tfont=22)




center    =  [198.77156,+51.81702]
name      = 'PSZ2 G114.31+64.89 / Abell 1703'
z         = 0.2836
width     = [2.3,2.3]
makeimagewrapper(prefixpath, 'PSZ2G114.31+64.89', center, width, z, name, rmsfactor=150., pexponent=0.35, plotcode='A1703', lfont=24, tfont=22)





center    = [168.8165798,53.33216888]
name      = 'PSZ2 G150.56+58.32 / MACS J1115.2+5320'
z         = 0.470
width     = [2.5,2.5]
makeimagewrapper(prefixpath, 'PSZ2G150.56+58.32', center, width, z, name, rmsfactor=150.,plotcode='MACSJ1115.2+5320', lfont=24, tfont=22, pexponent=0.35)




center    = [173.70601,+49.07690]
name      = 'Abell 1314'
z         = 0.0335
width     = [1.3,1.3]
makeimagewrapper(prefixpath, 'Abell1314', center, width, z, name, rmsfactor=150., plotcode='A1314', xs=0.25, lfont=24, tfont=22, pexponent=0.35)




#-------



center    = [202.62266667, +49.14669444] 
name      = 'PSZ2 G106.61+66.71'
z         =  0.331400
width     = [2.0,2.0]
makeimagewrapper('/net/lofar7/data1/botteon/DR2_sample/', 'PSZ2G106.61+66.71', center, width, z, name, rmsfactor=30., lfont=24, tfont=22)




center    = [228.5024535, +52.80381756]
name      = 'PSZ2 G086.93+53.18'
z         = 0.6752
width     = [2.0,2.0]
makeimagewrapper(prefixpath, 'PSZ2G086.93+53.18', center, width, z, name, rmsfactor=30.,  lfont=24, tfont=22, xs=0.25/6)



center    = [222.25533, +48.55664]
name      = 'PSZ2 G084.10+58.72'
z         = 0.731000
width     = [2.0,2.0]
makeimagewrapper(prefixpath, 'PSZ2G084.10+58.72', center, width, z, name, rmsfactor=30.,  lfont=24, tfont=22, xs=0.25/6)



center    = [231.63821,+54.15206]
name      = 'PSZ2 G087.39+50.92'
z         = 0.7480
width     = [3.0,3.0]
makeimagewrapper('/net/bovenrijn/data1/digennaro/HighRedshiftClusters/LOFAR/', 'PSZ2G087.39+50.92', center, width, z, name, rmsfactor=30.)







center    = [228.29890, +54.02127]
name      = 'WHL J151310.6+540116'
z         = 0.4366
width     = [2.5,2.5]
makeimagewrapper(prefixpath, 'WHLJ151310.6+540116', center, width, z, name, rmsfactor=30.)





center    = [227.53222,+53.95423]
name      = 'WHL J151001.6+535722'
z         = 0.4923
width     = [2.,2.]
makeimagewrapper(prefixpath, 'WHLJ151001.6+535722', center, width, z, name, rmsfactor=30.)

center    = [206.97006,+47.87677]
name      = 'WHL J134746.8+475214'
z         = 0.1695
width     = [2.0,2.0]
makeimagewrapper(prefixpath, 'WHLJ134746.8+475214', center, width, z, name, rmsfactor=30.)


center    = [201.56565,+48.87447]
name      = 'WHL J132615.8+485229'
z         = 0.2800
width     = [2.5,2.5]
makeimagewrapper(prefixpath, 'WHLJ132615.8+485229', center, width, z, name, rmsfactor=30.)




center    = [198.86497,+48.67361]
name      = 'WHL J131527.6+484025'
z         = 0.5158
width     = [2.5,2.5]
makeimagewrapper(prefixpath, 'WHLJ131527.6+484025', center, width, z, name, rmsfactor=30.)

center    = [194.44273,+48.91299]
name      = 'WHL J125746.3+485446'
z         = 0.2636
width     = [1.2,1.2]
makeimagewrapper(prefixpath, 'WHLJ125746.3+485446', center, width, z, name, rmsfactor=30.)

center    = [193.39217,+50.90403]
name      = 'WHL J125334.1+505414'
z         = 0.1226
width     = [1.2,1.2]
makeimagewrapper(prefixpath, 'WHLJ125334.1+505414', center, width, z, name, rmsfactor=30.)



center    = [192.13760,+55.80526]
name      = 'WHL J124833.0+554820'
z         = 0.27 # approx
width     = [1.75,1.75]
makeimagewrapper(prefixpath, 'WHLJ124833.0+554820', center, width, z, name, rmsfactor=30.)


center    = [187.06727,+49.83951]
name      = 'WHL J122816.1+495021'
z         = 0.2637
width     = [1.5,1.5]
makeimagewrapper(prefixpath, 'WHLJ122816.1+495021', center, width, z, name, rmsfactor=30.)




center    = [205.91115,+49.77400]
name      = 'GMBCG J205.93744+49.75911'
z         = 0.271
width     = [2.,2.]
makeimagewrapper(prefixpath, 'GMBCGJ205.93744+49.75911', center, width, z, name, rmsfactor=30.)

center    = [182.73972,+55.50167]
name      = 'SDSS J121057.53+553005.8'
z         = 0.3366
width     = [1.1,1.1]
makeimagewrapper(prefixpath, 'SDSSJ121057.53+553005.8', center, width, z, name, rmsfactor=30.)

center    = [177.54706,+53.72236]
name      = 'SDSS-C4-DR3 3106'
z         = 0.0603 
width     = [1.5,1.5]
makeimagewrapper(prefixpath, 'SDSS-C4-DR33106', center, width, z, name, rmsfactor=30.)


center    = [211.72963,+55.06747]
name      = 'GMBCG J211.77332+55.09968'
z         = 0.2506
width     = [2.5,2.5]
makeimagewrapper(prefixpath, 'GMBCGJ211.77332+55.09968', center, width, z, name, rmsfactor=30.)


center    = [181.90365,+52.91709]
name      = 'GMBCG J181.88181+52.89922'
z         = 0.2759
width     = [1.5,1.5]
makeimagewrapper(prefixpath, 'GMBCGJ181.88181+52.89922', center, width, z, name, rmsfactor=30.)



center    = [172.56418,+47.79317]
name      = 'ZwCl 1127.5+4804'
z         = 0.1267
width     = [1.5,1.5]
makeimagewrapper(prefixpath, 'ZwCl1127.5+4804', center, width, z, name, rmsfactor=30.)




center    = [189.71596,+55.64610]
name      = 'NSCS J123842+553825'
z         = 0.27834
width     = [2.1,2.1]
makeimagewrapper(prefixpath, 'NSCSJ123842+553825', center, width, z, name, rmsfactor=30.)



center    = [185.36489,+49.30730]
name      = 'MCXC J1221.4+4918'
z         = 0.700
width     = [1.5,1.5]
makeimagewrapper(prefixpath, 'MCXCJ1221.4+4918', center, width, z, name, rmsfactor=30.)



center    = [184.42961,+47.48647]
name      = 'MCXC J1217.7+4729'
z         = 0.2700
width     = [1.5,1.5]
makeimagewrapper(prefixpath, 'MCXCJ1217.7+4729', center, width, z, name, rmsfactor=30.)


center    = [163.43380,+49.49889]
name      = 'MCXC J1053.7+4929'
z         = 0.1400
width     = [1.5,1.5]
makeimagewrapper(prefixpath, 'MCXCJ1053.7+4929', center, width, z, name, rmsfactor=30.)


center    = [173.04772,+47.81041]
name      = 'MaxBCG J173.04772+47.81041'
z         = 0.2261
width     = [2.,2.]
makeimagewrapper(prefixpath, 'MaxBCGJ173.04772+47.81041', center, width, z, name, rmsfactor=30.)


center    = [182.49625, +49.89769]
name      = 'ClG J120958.9+495352'
z         = 0.9020
width     = [2.0,2.0]
makeimagewrapper(prefixpath, 'ClGJ120958.9+495352', center, width, z, name, rmsfactor=30.)


center    = [224.94504, +49.76864]
name      = 'Abell 2011'
z         = 0.1697
width     = [1.5,1.5]
makeimagewrapper(prefixpath, 'Abell2011', center, width, z, name, rmsfactor=30.)


center    = [223.59811,+54.32444]
name      = 'Abell 2000 / Abell 1999'
z         = 0.1012
width     = [2.0,2.0]
makeimagewrapper(prefixpath, 'Abell2000', center, width, z, name, rmsfactor=30.)

center    = [211.21535,+47.08447]
name      = 'Abell 1855'
z         = 0.2379
width     = [2.5,2.5]
makeimagewrapper(prefixpath, 'Abell1855', center, width, z, name, rmsfactor=30.)


center    = [207.25215, +49.31121]
name      = 'Abell 1804'
z         = 0.1665
width     = [1.6,1.6]
makeimagewrapper(prefixpath, 'Abell1804', center, width, z, name, rmsfactor=30.)


center    = [206.23556,+53.75058]
name      = 'Abell 1788'
z         = 0.1637
width     = [2.,2.]
makeimagewrapper(prefixpath, 'Abell1788', center, width, z, name, rmsfactor=30.)


center    = [201.69882,+53.82430]
name      = 'Abell 1745'
z         = 0.3665
width     = [2.,2.]
makeimagewrapper(prefixpath, 'Abell1745', center, width, z, name, rmsfactor=30.)


center    = [193.95482,+44.08842]
name      = 'Abell 1643'
z         = 0.2328
width     = [1.1,1.1]
makeimagewrapper(prefixpath, 'Abell1643', center, width, z, name, rmsfactor=30.)


center    = [191.92974,+48.86555]
name      = 'Abell 1615'
z         = 0.2106
width     = [1.5,1.5]
makeimagewrapper(prefixpath, 'Abell1615', center, width, z, name, rmsfactor=30.)

center    = [186.88430,+49.47882]
name      = 'Abell 1542'
z         = 0.1218
width     = [1.1,1.1]
makeimagewrapper(prefixpath, 'Abell1542', center, width, z, name, rmsfactor=30.)


center    = [176.83938,+55.72973]
name      = 'Abell 1377'
z         = 0.0514
width     = [1.5,1.5]
makeimagewrapper(prefixpath, 'Abell1377', center, width, z, name, rmsfactor=30.)


center    = [176.22393,+51.26713]
name      = 'Abell 1368'
z         = 0.1291
width     = [1.1,1.1]
makeimagewrapper(prefixpath, 'Abell1368', center, width, z, name, rmsfactor=30.)


center    = [173.09558, +55.96818]
name      = 'Abell 1291'
z         = 0.0527
width     = [1.1,1.1]
makeimagewrapper(prefixpath, 'Abell1291', center, width, z, name, rmsfactor=30.)


center    = [174.56139,+49.54541]
name      = 'Abell 1330'
z         = 0.2805
width     = [2.0,2.0]
makeimagewrapper(prefixpath, 'Abell1330', center, width, z, name, rmsfactor=30.)


center    = [171.83958,+48.29855]
name      = 'Abell 1261'
z         = 0.1562
width     = [1.1,1.1]
makeimagewrapper(prefixpath, 'Abell1261', center, width, z, name, rmsfactor=30.)


center    = [168.40705,+47.48700]
name      = 'Abell 1202'
z         = 0.1121
width     = [2.0,2.0]
makeimagewrapper(prefixpath, 'Abell1202', center, width, z, name, rmsfactor=30.)


center    = [203.15183,+54.31651]
name      = 'RXC J1332.6+5419'
z         = 0.10655
width     = [2.,2.]
makeimagewrapper(prefixpath, 'RXCJ1332.6+5419', center, width, z, name, rmsfactor=30.)


center    = [177.27167,+51.58205]
name      = 'PSZ2 G144.33+62.85 / Abell 1387'
z         = 0.1320
width     = [2.0,2.0]
makeimagewrapper(prefixpath, 'PSZ2G144.33+62.85', center, width, z, name, rmsfactor=30.)


center    = [216.852083,+55.750528]
name      = 'PSZ2 G098.44+56.59 / Abell 1920'
z         = 0.1318
width     = [2.6,2.6]
makeimagewrapper(prefixpath, 'PSZ2G098.44+56.59', center, width, z, name, rmsfactor=30.)


center    = [218.86564,+55.12877]
name      = 'PSZ2 G096.14+56.24 / Abell 1940'
z         = 0.139773
width     = [2.,2.]
makeimagewrapper(prefixpath, 'PSZ2G096.14+56.24', center, width, z, name, rmsfactor=30.)



center    = [224.7450589,+52.8173732]
name      = 'PSZ2 G088.98+55.07'
z         = 0.702346
width     = [2.0,2.0]
makeimagewrapper(prefixpath, 'PSZ2G088.98+55.07', center, width, z, name, rmsfactor=30.)


center    =  [163.71686, +55.35387]
name      = 'PSZ2 G151.62+54.78'
z         = 0.4864
width     = [2.0,2.0]
makeimagewrapper(prefixpath, 'PSZ2G151.62+54.78', center, width, z, name, rmsfactor=30.)




center    = [180.06851,56.26328]
name      = 'PSZ2 G136.92+59.46 / Abell 1436'
z         = 0.0650
width     = [1.5,1.5]
makeimagewrapper(prefixpath, 'PSZ2G136.92+59.46', center, width, z, name, rmsfactor=30.)


center    =  [203.1765617, 50.48]
name      = 'PSZ2 G107.10+65.32 / Abell 1758'
z         = 0.2799
width     = [4.3,4.3]
makeimagewrapper(prefixpath, 'PSZ2G107.10+65.32', center, width, z, name, rmsfactor=30., plotcode='A1758',  lfont=24, tfont=22)


center    = [213.69655,+54.78442]
name      = 'PSZ2 G099.86+58.45'
z         = 0.6305
width     = [3.,3.]
makeimagewrapper(prefixpath, 'PSZ2G099.86+58.45', center, width, z, name, rmsfactor=30.,  lfont=24, tfont=22)


center    = [207.92681, +46.36694]
name      = 'PSZ2 G095.22+67.41'
z         = 0.062500
width     = [2.0,2.0]
makeimagewrapper(prefixpath, 'PSZ2G095.22+67.41', center, width, z, name, rmsfactor=30.)


center    = [173.1767513,54.21999516]
name      = 'PSZ2 G145.65+59.30 / Abell 1294'
z         = 0.3475
width     = [2.0,2.0]
makeimagewrapper(prefixpath, 'PSZ2G145.65+59.30', center, width, z, name, rmsfactor=30., plotcode='A1294', lfont=24, tfont=22)









































































