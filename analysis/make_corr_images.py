from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import warnings
from scipy.stats.stats import pearsonr

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

f = fits.open('../regions/ptp_dir/grid_13x13/grid_13x13_results.fits')
header = f[0].header
t = f[1].data
t = t[(t['xray_sb']>0)
      & (t['radio1_sb']>0)
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

fitline = fit(np.log10(t['xray_sb']), np.log10(t['radio1_sb']))

plt.errorbar(np.log10(t['xray_sb']), np.log10(t['radio1_sb']), xerr=(0.234*t['xray_sb_err']/t['xray_sb']).clip(0,0.4), yerr=0.234*t['radio1_sb_err']/t['radio1_sb'], fmt='.', ecolor='coral', color='coral')
plt.plot(fitline[0], fitline[1], color='darkslateblue', linestyle='--')
# plt.fill_between(fitline[0], fitline[1]-t['radio1_sb_err'], fitline[1]+t['radio1_sb_err'])
plt.ylabel('log $I_{R}$ [Jy/arcsec$^2$]')
plt.xlabel('log $I_{X}$ [counts/s/arcsec$^2$]')
plt.xlim(np.min(np.log10(t['xray_sb'])-(0.434*t['xray_sb_err']/t['xray_sb']).clip(0,0.4)),
         np.max(np.log10(t['xray_sb'])+(0.434*t['xray_sb_err']/t['xray_sb']).clip(0,0.4)))
plt.ylim(np.min(np.log10(t['radio1_sb'])-(0.434*t['radio1_sb_err']/t['radio1_sb']).clip(0,0.4)),
         np.max(np.log10(t['radio1_sb'])+(0.434*t['radio1_sb_err']/t['radio1_sb']).clip(0,0.4)))
plt.grid(False)
plt.savefig('xraycorr.png')


fitline = fit(np.log10(t['y_sb']), np.log10(t['radio1_sb']))
plt.errorbar(np.log10(t['y_sb']), np.log10(t['radio1_sb']), xerr=(0.234*t['y_sb_err']/t['y_sb']).clip(0,0.4), yerr=0.234*t['radio1_sb_err']/t['radio1_sb'], fmt='.', ecolor='coral', color='coral')
plt.plot(fitline[0], fitline[1], color='darkslateblue', linestyle='--')
plt.ylabel('log $I_{R}$ [Jy/arcsec$^2$]')
plt.xlabel('log $y$ [Compton-y/arcsec$^2$]')
plt.xlim(np.min(np.log10(t['y_sb'])-(0.434*t['y_sb_err']/t['y_sb']).clip(0,0.4)),
         np.max(np.log10(t['y_sb'])+(0.434*t['y_sb_err']/t['y_sb']).clip(0,0.4)))
plt.ylim(np.min(np.log10(t['radio1_sb'])-(0.434*t['radio1_sb_err']/t['radio1_sb']).clip(0,0.4)),
         np.max(np.log10(t['radio1_sb'])+(0.434*t['radio1_sb_err']/t['radio1_sb']).clip(0,0.4)))
plt.grid(False)
plt.savefig('ymapcorr.png')

print('Pearson R (x-ray vs radio): '+str(pearsonr(np.log(t['xray_sb']), np.log(t['radio1_sb']))))
print('Pearson R (ymap vs radio): '+str(pearsonr(np.log(t['y_sb']), np.log(t['radio1_sb']))))