import scipy
import os.path
import pickle
import scipy.integrate as si
from astropy.io import ascii
import matplotlib.pyplot as plt
import numpy as np
import bces.bces
import nmmn.stats
import astroquery
from astroquery.vizier import Vizier
from astropy.coordinates import Angle
from astroquery.ned import Ned
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import sys
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

def zmassplot():
  astroquery.vizier.Vizier.ROW_LIMIT = 100000
  astroquery.vizier.Vizier.columns= ['**']
  if not os.path.isfile('psz2.p'):
    catalog1 = Vizier.get_catalogs('J/A+A/594/A27')
    print(catalog1)
    psz2 = catalog1[0]
    f = open('psz2.p', 'wb')
    pickle.dump(psz2, f)
    f.close()
  else:
    print('Open pickle file from disk')
    f = open('psz2.p', 'rb')   # 'r' for reading; can be omitted
    psz2 = pickle.load(f)         # load file content as mydict
    f.close()   
    
  hetdexlist = ['PSZ2 G080.16+57.65', 'PSZ2 G084.10+58.72', 'PSZ2 G086.93+53.18', 'PSZ2 G087.39+50.92', \
                'PSZ2 G088.98+55.07', 'PSZ2 G089.52+62.34', 'PSZ2 G095.22+67.41', 'PSZ2 G096.14+56.24', \
                'PSZ2 G098.44+56.59', 'PSZ2 G099.86+58.45', 'PSZ2 G106.61+66.71', 'PSZ2 G107.10+65.32', \
                'PSZ2 G111.75+70.37', 'PSZ2 G114.31+64.89', 'PSZ2 G114.99+70.36', 'PSZ2 G118.34+68.79', \
                'PSZ2 G123.66+67.25', 'PSZ2 G133.60+69.04', 'PSZ2 G135.17+65.43', 'PSZ2 G136.92+59.46', \
                'PSZ2 G143.26+65.24', 'PSZ2 G144.33+62.85', 'PSZ2 G145.65+59.30', 'PSZ2 G150.56+58.32', \
                'PSZ2 G151.62+54.78', 'PSZ2 G156.26+59.64']

  hetdexcmass = []
  hetdexz    = []
  print(psz2.keys())
  print(psz2['Name'])
  plt.plot(  psz2['z'], psz2['MSZ']*1e14, '.', markersize=4, label='Full PSZ2 sample', c='cornflowerblue')
  plt.ylabel('M$_{500}$ [M$_\odot$]',fontsize=13)
  plt.xlabel('Redshift',fontsize=13)
  plt.xlim(0.,1.0)
  for cluster in hetdexlist:
     #idx = list(psz2['Name']).index(cluster)
     idx = np.where(psz2['Name'] == cluster)
     hetdexcmass.append(psz2['MSZ'][idx])
     hetdexz.append(psz2['z'][idx])
     
  plt.plot( np.array(hetdexz), np.array(hetdexcmass)*1e14, 'o', markersize=4, label='PSZ2 clusters this work', c='darkred')   
  plt.semilogy()
  plt.legend(loc='lower right')
  plt.show()

  return

# zmassplot()

linestyle_str = [
     ('solid', 'solid'),      # Same as (0, ()) or '-'
     ('dotted', 'dotted'),    # Same as (0, (1, 1)) or '.'
     ('dashed', 'dashed'),    # Same as '--'
     ('dashdot', 'dashdot')]  # Same as '-.'

linestyle_tuple = [
     ('loosely dotted',        (0, (1, 10))),
     ('dotted',                (0, (1, 1))),
     ('densely dotted',        (0, (1, 1))),

     ('loosely dashed',        (0, (5, 10))),
     ('dashed',                (0, (5, 5))),
     ('densely dashed',        (0, (5, 1))),

     ('loosely dashdotted',    (0, (3, 10, 1, 10))),
     ('dashdotted',            (0, (3, 5, 1, 5))),
     ('densely dashdotted',    (0, (3, 1, 1, 1))),

     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))]


def Y5R500_to_Y500(Y5R500):
# Arnaud+ 2010; Cassano+ 2013
    return Y5R500/1.79

def Y500_from_arcmin2_to_Mpc2(Y500, z):
  D_A = cosmo.angular_diameter_distance(z)
  factor = ((np.pi/10800.)**2)*(D_A**2)
  print(D_A)

  return Y500*factor



# convenient function for a line
def func(x): return x[1]*x[0]+x[2]

def bcesfit(x, xerr, y, yerr, i=3):
  # Values for i documentation
  #0 y|x Assumes x as the independent variable
  #1 x|y Assumes y as the independent variable
  #2 bissector Line that bisects the y|x and x|y. This approach is self-inconsistent.
  #3 orthogonal Orthogonal least squares: line that minimizes orthogonal distances. Should be used when it is not clear which variable should be treated as the independent one

  if i == 0:
      label = 'Y | X'
  if i == 1:
      label = 'X | Y'
  if i == 2:
      label = 'bisector'
  if i == 3:
      label = 'orthogonal'      

  cov = cov = np.zeros_like(x)

  
  # number of bootstrapping trials
  nboot=100000
  
  xlog = np.log10(x)
  ylog = np.log10(y)
  xerrlog = xerr/(x*np.log(10.))
  yerrlog = yerr/(y*np.log(10.))

  slope,offset,slopeerr,offseterr,covab=bces.bces.bcesp(xlog,xerrlog,ylog,yerrlog,cov, nboot)
  
  # array with best-fit parameters
  fitm=np.array([ slope[i],offset[i] ])
  # covariance matrix of parameter uncertainties
  covm=np.array([ (slopeerr[i]**2,covab[i]), (covab[i],offseterr[i]**2) ])
  # Gets lower and upper bounds on the confidence band 
  

  xvec=np.linspace(xlog.min()-1.,xlog.max()+1.)
  #lcb,ucb,xcb=nmmn.stats.confbandnl(xlog,ylog,func,fitm,covm,2,0.954,xvec)

  lcb,ucb,xcb=nmmn.stats.confbandnl(xlog,ylog,func,fitm,covm,2,0.997,xvec)  
  
  return slope[i], slopeerr[i], offset[i], offseterr[i], label, lcb, ucb, xcb

def radiopower150(flux, z, alpha, freq):
  '''
  compute 150 MHz radio power
    - flux: integrated flux density in units of Jy
    - z: redshift
    - alpha: spectra index (F \propto nu^alpha)
    - freq: freuqency in MHz
    - returns: radio power in astropy units of W/Hz
  '''
  flux2 = ((150./freq)**alpha) * flux
  fluxinJyunits = flux2 * u.Jansky
  D_l = cosmo.luminosity_distance(z) # pc
  p = 4.*np.pi*(D_l)**2 * fluxinJyunits * ((1.+z)**((-1.*alpha) - 1.))
  return p.to(u.Watt/u.Hz)


alpha = -1.2
candidates_plotcolor = 'lightsteelblue'
candidates = False
ftype=3

data  = ascii.read("750kpc-fluxmeasurementsRW_Jort.txt",\
                   fast_reader=False,  format='csv', comment='#',header_start=0)  


p150     = radiopower150(data['Flux'], data['z'], alpha, data['Freq'])
p150_upp = radiopower150(data['Flux']+data['Flux_err'],  data['z'], alpha, data['Freq'])
p150_low =  radiopower150(data['Flux']-data['Flux_err'], data['z'], alpha, data['Freq'])
p150err  = np.abs(0.5*(p150_upp - p150_low))

idx = np.where(data['haloorcandidate'] == 1)
print(idx)

p150_h     = radiopower150(data['Flux'][idx], data['z'][idx], alpha, data['Freq'][idx])
p150_upp_h = radiopower150(data['Flux'][idx]+data['Flux_err'][idx],  data['z'][idx], alpha, data['Freq'][idx])
p150_low_h =  radiopower150(data['Flux'][idx]-data['Flux_err'][idx], data['z'][idx], alpha, data['Freq'][idx])
p150err_h  = np.abs(0.5*(p150_upp_h - p150_low_h))


#plt.plot(data['Mass']*1e14,p150,'ro')


#plt.plot(literature['Mass']*1e14,power_lit,'bo')

#print 'Abell697', np.where(literature['Name'] == 'Abell697')
#plt.plot(literature['Mass'][np.where(literature['Name'] == 'Abell697')]*1e14,power_lit[np.where(literature['Name'] == 'Abell697')],'go')
#plt.plot(literature['Mass'][np.where(literature['Name'] == 'Abell521')]*1e14,power_lit[np.where(literature['Name'] == 'Abell521')],'go')
#plt.plot(literature['Mass'][np.where(literature['Name'] == 'Abell1132')]*1e14,power_lit[np.where(literature['Name'] == 'Abell1132')],'go')

plt.style.use('seaborn')
plt.rcParams.update({'axes.facecolor':'white'})

#line
#B, A = 3.70, 0.09 #BCES bisector Cassano et al. (2013); $\alpha$=-1.3
B, A = 4.51, 0.129 #5.05, 0.02 #BCES orthogonal Cassano et al. (2013); $\alpha$=-1.3
massvec = np.arange(0.1,100,0.1)*1e14
p150sc12  = ((150./1400.)**alpha) * (10**(24.5)) * 10**((B*np.log10(massvec/(10**(14.9)))) + A)

p150sc8  = ((150./1400.)**(-0.8)) * (10**(24.5)) * 10**((B*np.log10(massvec/(10**(14.9)))) + A)
p150sc10  = ((150./1400.)**(-1.0)) * (10**(24.5)) * 10**((B*np.log10(massvec/(10**(14.9)))) + A)
p150sc11  = ((150./1400.)**(-1.1)) * (10**(24.5)) * 10**((B*np.log10(massvec/(10**(14.9)))) + A)
p150sc13  = ((150./1400.)**(-1.3)) * (10**(24.5)) * 10**((B*np.log10(massvec/(10**(14.9)))) + A)
p150sc14  = ((150./1400.)**(-1.4)) * (10**(24.5)) * 10**((B*np.log10(massvec/(10**(14.9)))) + A)
p150sc15  = ((150./1400.)**(-1.5)) * (10**(24.5)) * 10**((B*np.log10(massvec/(10**(14.9)))) + A)
p150sc16  = ((150./1400.)**(-1.6)) * (10**(24.5)) * 10**((B*np.log10(massvec/(10**(14.9)))) + A)
p150sc18  = ((150./1400.)**(-1.8)) * (10**(24.5)) * 10**((B*np.log10(massvec/(10**(14.9)))) + A)
#p150sc  = ((150./1400.)**alpha)*p150sc

print (A, B)
# plt.plot(massvec, p150sc12, color='b', label='C2013; '+ r'$\alpha=-1.2$')

#plt.plot(massvec, p150sc8, color='k', label='C2013; '+ r'$\alpha=-0.8$',alpha=0.5,ls='--')
#plt.plot(massvec, p150sc10, color='k', label='C2013; '+ r'$\alpha=-1.0$',alpha=0.5,ls='dashdot')
#plt.plot(massvec, p150sc14, color='k', label='C2013; '+ r'$\alpha=-1.4$',alpha=0.5, ls=(0, (3, 1, 1, 1, 1, 1)))
#plt.plot(massvec, p150sc16, color='k', label='C2013; '+ r'$\alpha=-1.6$',alpha=0.5,ls='dotted' )


plt.errorbar(np.array(data['Mass'])*1e14, np.array(p150), xerr=np.array([data['Mass_err_down'],data['Mass_err_up']])*1e14, yerr=np.array(p150err), fmt='o',capsize=2,markersize=4,  label='van Weeren et al. 2021 (candidates)', c=candidates_plotcolor)
plt.errorbar(np.array(data['Mass'][idx])*1e14, np.array(p150_h), xerr=np.array([data['Mass_err_down'][idx],data['Mass_err_up'][idx]])*1e14, yerr=np.array(p150err_h), fmt='o',capsize=2,markersize=4,  label='van Weeren et al. 2021')

literature = ascii.read("literature150mhz.txt",fast_reader=False,  format='csv', comment='#',header_start=0)

p150lit     = radiopower150(1e-3*literature['Flux'], literature['z'], alpha, literature['Freq'])
p150lit_upp = radiopower150(1e-3*(literature['Flux']+literature['Flux_err']),  literature['z'], alpha, literature['Freq'])
p150lot_low =  radiopower150(1e-3*(literature['Flux']-literature['Flux_err']), literature['z'], alpha, literature['Freq'])
p150literr  = np.abs(0.5*(p150lit_upp - p150lot_low))

plt.errorbar(np.array(literature['Mass'])*1e14, np.array(p150lit), xerr=np.array([literature['Mass_err_down'],literature['Mass_err_up']])*1e14, yerr=np.array(p150literr), fmt='ok',capsize=2,markersize=4, label='literature')


# -- DO BCES fit --
# concatenate all arrays first
if candidates:
  massall =  np.concatenate((np.array(data['Mass']), np.array(literature['Mass'])), axis=0)
  p150all =  np.concatenate((np.array(p150), np.array(p150lit)), axis=0)
  p150allerr =  np.concatenate((np.array(p150err), np.array(p150literr)), axis=0)
  mass_avg_err    =  np.array(0.5*(data['Mass_err_up'] + data['Mass_err_down']))
  masslit_avg_err =  np.array(0.5*(literature['Mass_err_up'] + literature['Mass_err_down']))
  massallerr = np.concatenate((mass_avg_err,masslit_avg_err), axis=0)
else:
  massall =  np.concatenate((np.array(data['Mass'][idx]), np.array(literature['Mass'])), axis=0)
  p150all =  np.concatenate((np.array(p150_h), np.array(p150lit)), axis=0)
  p150allerr =  np.concatenate((np.array(p150err_h), np.array(p150literr)), axis=0)
  mass_avg_err    =  np.array(0.5*(data['Mass_err_up'][idx] + data['Mass_err_down'][idx]))
  masslit_avg_err =  np.array(0.5*(literature['Mass_err_up'] + literature['Mass_err_down']))
  massallerr = np.concatenate((mass_avg_err,masslit_avg_err), axis=0)


# DO fit
slope, slopeerr, offset, offseterr, plotlabel,  lc, uc, xc = bcesfit(massall*1e14/(10**(14.9)), massallerr*1e14/(10**(14.9)), p150all/(10**(24.5)), p150allerr/(10**(24.5)), i=ftype)

print ('M-P (slope, slope_err, offset, offset_err)',  slope, slopeerr, offset, offseterr)

#plot
#0 y|x Assumes x as the independent variable
#1 x|y Assumes y as the independent variable
#2 bissector Line that bisects the y|x and x|y. This approach is self-inconsistent.
#3 orthogonal Orthogonal least squares: line that minimizes orthogonal distances. Should be used when it is not clear which variable should be treated as the independent one

planck_A399, planck_A399_err = 5.2e14, 2.6e13
planck_A401, planck_A401_err = 6.7e14, 2.0e13
# xray_A399, xray_A399_err = 5.7e14, 0
# xray_A401, xray_A401_err = 9.3e14, 0
bridge_hincks, bridge_hincks_err = 3.3e14, 7e13
radio_A399, radio_A399_err = 1.63e+25, 5e+23
radio_A401, radio_A401_err = 1.28e+25, 4e+23
radio_bridge, radio_bridge_err = 8.3e24, 2e23

p150sc = (10**(24.5)) * 10**((slope*np.log10(massvec/(10**(14.9)))) + offset)
plt.grid(False)

plt.plot(massvec, p150sc, color='darkred', label='van Weeren et al. 2021 fit')

xc = (10**xc)*10**(14.9)
lc = (10**lc)*10**(24.5)
uc = (10**uc)*10**(24.5)
plt.fill_between(xc, lc, uc, alpha=0.2, facecolor='red')

plt.errorbar(np.array([planck_A399]), np.array([radio_A399*((144/150)**1.75)]), yerr=np.array([np.sqrt(radio_A399_err**2)]), xerr=np.array([planck_A399_err]), capsize=2,markersize=4, color='darkgreen', ecolor='darkgreen', label='A399')
plt.errorbar(np.array([planck_A401]), np.array([radio_A401*((144/150)**1.63)]), yerr=np.array([np.sqrt(radio_A401_err**2)]), xerr=np.array([planck_A401_err]), capsize=2,markersize=4, ecolor='black', color='black', label='A401')
# plt.errorbar(np.array([bridge_hincks]), np.array([radio_bridge*((144/150)**1.5)]), yerr=np.array([np.sqrt(radio_bridge_err**2+(0.1*radio_bridge)**2)]), xerr=np.array([bridge_hincks_err]), capsize=2, markersize=4, ecolor='darkblue', color='darkblue', label='Bridge')


plt.legend(loc='upper left')
plt.ylim(3e23,1e27)
plt.xlim(2e14, 2e15)

plt.xlabel('M$_{500}$ [M$_\odot$]',fontsize=13)
plt.ylabel('P$_{150 \: \mathrm{MHz}}$  [W Hz$^{-1}$]',fontsize=13) 
# plt.title('cluster mass $-$ 150 MHz radio halo power diagram',fontsize=13)
plt.semilogx()
plt.semilogy()
plt.savefig('masspower.png', bbox_inches='tight')


# Y500 = Y500_from_arcmin2_to_Mpc2(Y5R500_to_Y500(1e-3*data['Y5R500']), data['z'])
# Y500err = Y500_from_arcmin2_to_Mpc2(Y5R500_to_Y500(1e-3*data['Y5R500_err']), data['z'])
# Y500lit = Y500_from_arcmin2_to_Mpc2(Y5R500_to_Y500(1e-3*literature['Y5R500']), literature['z'])
# Y500literr = Y500_from_arcmin2_to_Mpc2(Y5R500_to_Y500(1e-3*literature['Y5R500_err']), literature['z'])
#
# #print literature['Y5R500_err'], literature['Y5R500']
# #print Y500lit
# #print Y500literr
# #sys.exit()
#
# plt.errorbar(np.array(Y500), np.array(p150), xerr=np.array(Y500err), yerr=np.array(p150err), fmt='o',capsize=2,markersize=4,  label='this work (candidates)',c=candidates_plotcolor)
#
# plt.errorbar(np.array(Y500[idx]), np.array(p150_h), xerr=np.array(Y500err[idx]), yerr=np.array(p150err[idx]), fmt='o',capsize=2,markersize=4,  label='this work')
#
#
#
# plt.errorbar(np.array(Y500lit), np.array(p150lit), xerr=np.array(Y500literr), yerr=np.array(p150literr), fmt='ok',capsize=2,markersize=4, label='literature')
#
#
# B, A = 2.28, -0.027 #2.48 , -0.167 #BCES orthogonal Cassano et al. (2013); $\alpha$=-1.3
# Yvec = 10**(np.arange(-6,-2,0.1))
# p150sc12  = ((150./1400.)**alpha)*(10**(24.5)) * 10**((B*np.log10(Yvec/1e-4)) + A)
#
#
# p150sc8   = ((150./1400.)**(-0.8)) * (10**(24.5)) * 10**((B*np.log10(Yvec/1e-4)) + A)
# p150sc10  = ((150./1400.)**(-1.0)) * (10**(24.5)) * 10**((B*np.log10(Yvec/1e-4)) + A)
# p150sc11  = ((150./1400.)**(-1.1)) * (10**(24.5)) * 10**((B*np.log10(Yvec/1e-4)) + A)
# p150sc13  = ((150./1400.)**(-1.3)) * (10**(24.5)) * 10**((B*np.log10(Yvec/1e-4)) + A)
# p150sc14  = ((150./1400.)**(-1.4)) * (10**(24.5)) * 10**((B*np.log10(Yvec/1e-4)) + A)
# p150sc15  = ((150./1400.)**(-1.5)) * (10**(24.5)) * 10**((B*np.log10(Yvec/1e-4)) + A)
# p150sc16  = ((150./1400.)**(-1.6)) * (10**(24.5)) * 10**((B*np.log10(Yvec/1e-4)) + A)
# p150sc18  = ((150./1400.)**(-1.8)) * (10**(24.5)) * 10**((B*np.log10(Yvec/1e-4)) + A)
#
#
#
# print (A, B)
# plt.plot(Yvec, p150sc12, color='b', label='C2013; '+ r'$\alpha=-1.2$')
#
# #plt.plot(Yvec, p150sc8, color='k', label='C2013; '+ r'$\alpha=-0.8$',alpha=0.5,ls='--')
# #plt.plot(Yvec, p150sc10, color='k', label='C2013; '+ r'$\alpha=-1.0$',alpha=0.5,ls='dashdot')
# #plt.plot(Yvec, p150sc14, color='k', label='C2013; '+ r'$\alpha=-1.4$',alpha=0.5, ls=(0, (3, 1, 1, 1, 1, 1)))
# #plt.plot(Yvec, p150sc16, color='k', label='C2013; '+ r'$\alpha=-1.6$',alpha=0.5,ls='dotted' )
#
# if candidates:
#   Y500all =  np.concatenate((np.array(Y500), np.array(Y500lit)), axis=0)
#   Y500allerr =  np.concatenate((np.array(Y500err), np.array(Y500literr)), axis=0)
# else:
#   Y500all =  np.concatenate((np.array(Y500[idx]), np.array(Y500lit)), axis=0)
#   Y500allerr =  np.concatenate((np.array(Y500err[idx]), np.array(Y500literr)), axis=0)
#
# # DO fit
# slope, slopeerr, offset, offseterr, plotlabel,  lc, uc, xc = bcesfit(Y500all/1e-4, Y500allerr/1e-4, p150all/(10**(24.5)), p150allerr/(10**(24.5)), i=ftype)
# p150sc = (10**(24.5)) * 10**((slope*np.log10(Yvec/1e-4)) + offset)
# plt.plot(Yvec, p150sc, color='k', label='BCES ' + plotlabel)
# xc = (10**xc)*1e-4
# lc = (10**lc)*10**(24.5)
# uc = (10**uc)*10**(24.5)
# plt.fill_between(xc, lc, uc, alpha=0.3, facecolor='grey')
#
# print('Y-P (slope, slope_err, offset, offset_err)', slope, slopeerr, offset, offseterr)
#
#
# plt.xlabel('Y$_{500}$ [Mpc$^{2}$]',fontsize=13)
# plt.ylabel('P$_{150 \: \mathrm{MHz}}$  [W Hz$^{-1}$]',fontsize=13)
# plt.title('cluster Y$_{500}$ $-$ 150 MHz radio halo power diagram',fontsize=13)
# plt.ylim(3e23,1e27)
# plt.xlim(1e-5, 6e-4)
# #plt.yscale('log')
# #plt.xscale('log')
# plt.semilogx()
# plt.semilogy()
# plt.legend(loc='lower right')
# plt.show()
#
# print(Y500lit)
# print(Y500literr)
