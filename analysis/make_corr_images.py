from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import warnings
from scipy.stats.stats import pearsonr, spearmanr
from past.utils import old_div
import astropy.units as u
from astropy.wcs import WCS
import argparse

parser = argparse.ArgumentParser(
    description='Perform point-to-point analysis (radio/X or radio/radio/X) and save the results in a fits table. If (i) the X-ray counts < 0, (ii) the cell is not totally inside the X-ray FoV, and (ii) the radio flux density is < 0, it saves NaNs.')
parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
# REQUIRED
required.add_argument('-filein', type=str, required=True)
required.add_argument('-fileout', type=str, required=True)
args = parser.parse_args()




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

f1 = fits.open('fits/20median.fits')
wcs =WCS(f1[0].header, naxis=2)
header = wcs.to_header()
rms = findrms(f1[0].data)/calc_beamarea(f1)/((header['CDELT2']*u.deg).to(u.arcsec)**2).value

f1.close()

f = fits.open(args.filein)
header = f[0].header
t = f[1].data

t = t[(t['xray_sb']>0)
      & (t['radio1_sb']>rms*3)
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

# fig, ax = plt.subplots(constrained_layout=True)
# # ax.plot([np.min(np.log10(t['xray_sb'])-(t['xray_sb_err']/t['xray_sb'])),np.max(np.log10(t['xray_sb'])+(t['xray_sb_err']/t['xray_sb']))],
# #          [np.log10(rms*3), np.log10(rms*3)], linestyle='--', color='green')
# ax.errorbar(np.log10(t['xray_sb']), np.log10(t['radio1_sb']), xerr=(0.434*t['xray_sb_err']/t['xray_sb']), yerr=0.434*t['radio1_sb_err']/t['radio1_sb'], fmt='.', ecolor='red', elinewidth=0.4, color='red')
# ax.set_ylim(np.min(np.log10(t['radio1_sb'])-(0.434*t['radio1_sb_err']/t['radio1_sb'])),
#          np.max(np.log10(t['radio1_sb'])+(0.434*t['radio1_sb_err']/t['radio1_sb'])))
# ax.set_xlim(np.min(np.log10(t['xray_sb'])-(t['xray_sb_err']/t['xray_sb'])),
#          np.max(np.log10(t['xray_sb'])+(t['xray_sb_err']/t['xray_sb'])))
# plt.grid(False)
# ax2 = ax.twiny()
# ax2.errorbar(np.log10(t['y_sb']), np.log10(t['radio1_sb']), xerr=(t['y_sb_err']/t['y_sb']), yerr=0.434*t['radio1_sb_err']/t['radio1_sb'], fmt='.', ecolor='blue', elinewidth=0.4, color='blue')
# ax2.set_ylim(np.min(np.log10(t['radio1_sb'])-(0.434*t['radio1_sb_err']/t['radio1_sb'])),
#          np.max(np.log10(t['radio1_sb'])+(0.434*t['radio1_sb_err']/t['radio1_sb'])))
# ax2.set_xlim(np.min(np.log10(t['y_sb'])-(t['y_sb_err']/t['y_sb'])),
#          np.max(np.log10(t['y_sb'])+(t['y_sb_err']/t['y_sb'])))
# ax.set_ylabel('log $I_{R}$ [Jy/arcsec$^2$]')
# ax.set_xlabel('log $I_{X}$ [counts/s/arcsec$^2$]')
# ax2.set_xlabel('log $y$ [Compton-y/arcsec$^2$]')
# fig.legend(['Radio vs. X-ray', 'Radio vs. SZ'], loc='upper left')
# # ax.plot(fitlinex[0], fitlinex[1], color='red', linestyle='--')
# # ax2.plot(fitliney[0], fitliney[1], color='blue', linestyle='--')
# plt.tight_layout()
# plt.grid(False)
# plt.savefig('a399.png', bbox_inches='tight')

radio = t['radio1_sb']
xray = t['xray_sb']
y = t['y_sb']

radio_err = t['radio1_sb_err']*2
xray_err = t['xray_sb_err']
y_err = t['y_sb_err']*2

radio_err/=np.mean(radio)
xray_err/=np.mean(xray)
y_err/=np.mean(y)

radio/=np.mean(radio)
xray/=np.mean(xray)
y/=np.mean(y)

fig, ax = plt.subplots(constrained_layout=True)
ax.errorbar(np.log10(xray), np.log10(radio), xerr=(0.434*xray_err/xray), yerr=0.434*radio_err/radio, fmt='.', ecolor='red', elinewidth=0.4, color='darkred')
ax.set_ylim(np.min(np.log10(radio)-(0.434*radio_err/radio)),
         np.max(np.log10(radio)+(0.434*radio_err/radio)))
ax.set_xlim(np.min(np.log10(xray)-(xray_err/xray)),
         np.max(np.log10(xray)+(xray_err/xray)))
plt.grid(False)
ax.errorbar(np.log10(y), np.log10(radio), xerr=(y_err/y), yerr=0.434*radio_err/radio, fmt='.', ecolor='blue', elinewidth=0.4, color='darkblue')
ax.set_xlim(-1, 0.75)
ax.set_ylim(-0.7, 0.6)
ax.set_ylabel('log($I_{R}$) [SB/mean(SB)]')
# ax.set_xlabel('X-ray [SB/mean(SB)]')
ax.set_xlabel('log($I_{X}$) [SB/mean(SB)] and log(y) [SZ/mean(SZ)]')
fig.legend(['Radio vs. X-ray', 'Radio vs. SZ'], loc='upper left')
plt.tight_layout()
plt.grid(False)
plt.savefig(args.fileout, bbox_inches='tight')

print('Number of cells used: '+str(len(t)))

print('Pearson R (x-ray vs radio): '+str(pearsonr(np.log(xray), np.log(radio))))
print('Pearson R (ymap vs radio): '+str(pearsonr(np.log(y), np.log(radio))))

print('Spearman R (x-ray vs radio): '+str(spearmanr(np.log(xray), np.log(radio))))
print('Spearman R (ymap vs radio): '+str(spearmanr(np.log(y), np.log(radio))))