from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import warnings
from scipy.stats.stats import pearsonr, spearmanr
from past.utils import old_div
import astropy.units as u
from astropy.wcs import WCS


# os.system('python /home/jurjen/Documents/Python/lofar_helpers/analysis/ppt.py '
#           '-radio1 /home/jurjen/Documents/Python/lofar_helpers/fits/80all.fits '
#           '-limits 604 650 815 859 '
#           '-xsou /home/jurjen/Documents/Python/lofar_helpers/fits/mosaic_a399_a401.fits '
#           '-xbkg /home/jurjen/Documents/Python/lofar_helpers/fits/mosaic_a399_a401_bkg.fits '
#           '-xexp /home/jurjen/Documents/Python/lofar_helpers/fits/mosaic_a399_a401_exp.fits '
#           '-cellsize 15 '
#           '-y /home/jurjen/Documents/Python/lofar_helpers/fits/a401_curdecmaps_0.2_1.5s_sz.fits '
#           '-excluderegion /home/jurjen/Documents/Python/lofar_helpers/regions/excluderegions.reg')

warnings.filterwarnings('ignore')
plt.style.use('ggplot')

def findrms(mIn, maskSup=1e-7):
    """
    find the rms of an array, from Cycil Tasse/kMS
    """
    m = mIn[np.abs(mIn) > maskSup]
    rmsold = np.std(m)
    diff = 1e-1
    cut = 3.
    bins = np.arange(np.min(m), np.max(m), (np.max(m) - np.min(m)) / 30.)
    med = np.median(m)
    for i in range(10):
        ind = np.where(np.abs(m - med) < rmsold * cut)[0]
        rms = np.std(m[ind])
        if np.abs(old_div((rms - rmsold), rmsold)) < diff: break
        rmsold = rms
    return rms

def calc_beamarea(hdu):
    # Given a fitsfile this calculates the beamarea in pixels

    bmaj = hdu[0].header['BMAJ']
    bmin = hdu[0].header['BMIN']

    beammaj = bmaj / (2.0 * (2 * np.log(2)) ** 0.5)  # Convert to sigma
    beammin = bmin / (2.0 * (2 * np.log(2)) ** 0.5)  # Convert to sigma
    pixarea = abs(hdu[0].header['CDELT1'] * hdu[0].header['CDELT2'])

    beamarea = 2 * np.pi * 1.0 * beammaj * beammin  # Note that the volume of a two dimensional gaus$
    beamarea_pix = beamarea / pixarea

    return beamarea_pix

f1 = fits.open('../fits/20median.fits')
wcs =WCS(f1[0].header, naxis=2)
header = wcs.to_header()
rms = findrms(f1[0].data)/calc_beamarea(f1)/((header['CDELT2']*u.deg).to(u.arcsec)**2).value
f1.close()

f = fits.open('../ptp_dir/grid_75x75/grid_75x75_results.fits')
header = f[0].header
t = f[1].data

t = t[(t['xray_sb']>0)
      # & (t['radio1_sb']>rms*3)
      & (t['xray_sb_err']>0)
      & (t['radio1_sb_err']>0)]

# print(t['xray_sb_err']/t['xray_sb'])
# print(t['radio1_sb_err']/t['radio1_sb'])
# print(t['y_sb_err']/t['y_sb'])

def objective(x, a, b):
    return a * x + b

def fit(x, y):
    popt, _ = curve_fit(objective, x, y)
    a, b = popt
    x_line = np.arange(min(x)-1, max(x)+1, 0.01)
    y_line = objective(x_line, a, b)
    print('y = %.5f * x + %.5f' % (a, b))
    return x_line,y_line

fitlinex = fit(np.log10(t['xray_sb']), np.log10(t['radio1_sb']))
#
# plt.plot([np.min(np.log10(t['xray_sb'])-(0.434*t['xray_sb_err']/t['xray_sb'])),np.max(np.log10(t['xray_sb'])+(0.434*t['xray_sb_err']/t['xray_sb']))],
#          [np.log10(rms*3), np.log10(rms*3)], linestyle='--', color='red')
# plt.errorbar(np.log10(t['xray_sb']), np.log10(t['radio1_sb']), xerr=(0.234*t['xray_sb_err']/t['xray_sb']), yerr=0.234*t['radio1_sb_err']/t['radio1_sb'], fmt='.', ecolor='blue', color='blue')
# # plt.plot(fitline[0], fitline[1], color='darkslateblue', linestyle='--')
# # plt.fill_between(fitline[0], fitline[1]-t['radio1_sb_err'], fitline[1]+t['radio1_sb_err'])
# plt.ylabel('log $I_{R}$ [Jy/arcsec$^2$]')
# plt.xlabel('log $I_{X}$ [counts/s/arcsec$^2$]')
# plt.xlim(np.min(np.log10(t['xray_sb'])-(0.434*t['xray_sb_err']/t['xray_sb'])),
#          np.max(np.log10(t['xray_sb'])+(0.434*t['xray_sb_err']/t['xray_sb'])))
# plt.ylim(np.min(np.log10(t['radio1_sb'])-(0.434*t['radio1_sb_err']/t['radio1_sb'])),
#          np.max(np.log10(t['radio1_sb'])+(0.434*t['radio1_sb_err']/t['radio1_sb'])))
# plt.grid(False)
# plt.legend(['$3\sigma$'])
# plt.savefig('xraycorr.png')
# plt.close()


fitliney = fit(np.log10(t['y_sb']), np.log10(t['radio1_sb']))
# plt.plot([np.min(np.log10(t['y_sb'])-(0.434*t['y_sb_err']/t['y_sb'])),np.max(np.log10(t['y_sb'])+(0.434*t['y_sb_err']/t['y_sb']))],
#          [np.log10(rms*3), np.log10(rms*3)], linestyle='--', color='red')
# plt.errorbar(np.log10(t['y_sb']), np.log10(t['radio1_sb']), xerr=(t['y_sb_err']/t['y_sb']), yerr=0.434*t['radio1_sb_err']/t['radio1_sb'], fmt='.', ecolor='blue', color='blue')
# # plt.plot(fitline[0], fitline[1], color='darkslateblue', linestyle='--')
# plt.ylabel('log $I_{R}$ [Jy/arcsec$^2$]')
# plt.xlabel('log $y$ [Compton-y/arcsec$^2$]')
# plt.xlim(np.min(np.log10(t['y_sb'])-(t['y_sb_err']/t['y_sb'])),
#          np.max(np.log10(t['y_sb'])+(t['y_sb_err']/t['y_sb'])))
# plt.ylim(np.min(np.log10(t['radio1_sb'])-(0.434*t['radio1_sb_err']/t['radio1_sb'])),
#          np.max(np.log10(t['radio1_sb'])+(0.434*t['radio1_sb_err']/t['radio1_sb'])))
# plt.grid(False)
# plt.legend(['$3\sigma$'])
# plt.savefig('ymapcorr.png')
# plt.close()

fig, ax = plt.subplots(constrained_layout=True)
# ax.plot([np.min(np.log10(t['xray_sb'])-(t['xray_sb_err']/t['xray_sb'])),np.max(np.log10(t['xray_sb'])+(t['xray_sb_err']/t['xray_sb']))],
#          [np.log10(rms*3), np.log10(rms*3)], linestyle='--', color='green')
ax.errorbar(np.log10(t['xray_sb']), np.log10(t['radio1_sb']), xerr=(0.434*t['xray_sb_err']/t['xray_sb']), yerr=0.434*t['radio1_sb_err']/t['radio1_sb'], fmt='.', ecolor='red', elinewidth=0.4, color='red')
ax.set_ylim(np.min(np.log10(t['radio1_sb'])-(0.434*t['radio1_sb_err']/t['radio1_sb'])),
         np.max(np.log10(t['radio1_sb'])+(0.434*t['radio1_sb_err']/t['radio1_sb'])))
ax.set_xlim(np.min(np.log10(t['xray_sb'])-(t['xray_sb_err']/t['xray_sb'])),
         np.max(np.log10(t['xray_sb'])+(t['xray_sb_err']/t['xray_sb'])))
plt.grid(False)
ax2 = ax.twiny()
ax2.errorbar(np.log10(t['y_sb']), np.log10(t['radio1_sb']), xerr=(t['y_sb_err']/t['y_sb']), yerr=0.434*t['radio1_sb_err']/t['radio1_sb'], fmt='.', ecolor='blue', elinewidth=0.4, color='blue')
ax2.set_ylim(np.min(np.log10(t['radio1_sb'])-(0.434*t['radio1_sb_err']/t['radio1_sb'])),
         np.max(np.log10(t['radio1_sb'])+(0.434*t['radio1_sb_err']/t['radio1_sb'])))
ax2.set_xlim(np.min(np.log10(t['y_sb'])-(t['y_sb_err']/t['y_sb'])),
         np.max(np.log10(t['y_sb'])+(t['y_sb_err']/t['y_sb'])))
ax.set_ylabel('log $I_{R}$ [Jy/arcsec$^2$]')
ax.set_xlabel('log $I_{X}$ [counts/s/arcsec$^2$]')
ax2.set_xlabel('log $y$ [Compton-y/arcsec$^2$]')
fig.legend(['Radio vs. X-ray', 'Radio vs. SZ'], loc='upper left')
# ax.plot(fitlinex[0], fitlinex[1], color='red', linestyle='--')
# ax2.plot(fitliney[0], fitliney[1], color='blue', linestyle='--')
plt.tight_layout()
plt.grid(False)
plt.savefig('corrcombiplot.png', bbox_inches='tight')

print('Number of cells used: '+str(len(t)))

print('Pearson R (x-ray vs radio): '+str(pearsonr(np.log(t['xray_sb']), np.log(t['radio1_sb']))))
print('Pearson R (ymap vs radio): '+str(pearsonr(np.log(t['y_sb']), np.log(t['radio1_sb']))))

print('Spearman R (x-ray vs radio): '+str(spearmanr(np.log(t['xray_sb']), np.log(t['radio1_sb']))))
print('Spearman R (ymap vs radio): '+str(spearmanr(np.log(t['y_sb']), np.log(t['radio1_sb']))))