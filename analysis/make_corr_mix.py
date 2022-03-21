from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import warnings
from scipy.stats.stats import pearsonr, spearmanr, linregress
from past.utils import old_div
import astropy.units as u
from astropy.wcs import WCS
from scipy import stats
import argparse

parser = argparse.ArgumentParser(
    description='Perform point-to-point analysis (radio/X or radio/radio/X) and save the results in a fits table. If (i) the X-ray counts < 0, (ii) the cell is not totally inside the X-ray FoV, and (ii) the radio flux density is < 0, it saves NaNs.')
parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
# REQUIRED
required.add_argument('-obj', type=str, required=True)
args = parser.parse_args()

obj = args.obj

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

f1 = fits.open('fits/60rudnick.fits')
wcs =WCS(f1[0].header, naxis=2)
header = wcs.to_header()
rms = findrms(f1[0].data)/calc_beamarea(f1)/((header['CDELT2']*u.deg).to(u.arcsec)**2).value
f1.close()
f = fits.open(f'analysis/{obj}_results_rudnick.fits')
header = f[0].header
t1 = f[1].data
t1 = t1[(t1['radio1_sb']>2*rms)]

f1 = fits.open('fits/60cleanbridge_200kpc.fits')
wcs =WCS(f1[0].header, naxis=2)
header = wcs.to_header()
rms = findrms(f1[0].data)/calc_beamarea(f1)/((header['CDELT2']*u.deg).to(u.arcsec)**2).value
f1.close()
f = fits.open(f'analysis/{obj}_results_cb.fits')
header = f[0].header
t2 = f[1].data
t2 = t2[(t2['radio1_sb']>2*rms)]


def pearsonr_ci(x,y,alpha=0.05):
    ''' calculate Pearson correlation along with the confidence interval using scipy and numpy
    Parameters
    ----------
    x, y : iterable object such as a list or np.array
      Input for correlation calculation
    alpha : float
      Significance level. 0.05 by default
    Returns
    -------
    r : float
      Pearson's correlation coefficient
    pval : float
      The corresponding p value
    lo, hi : float
      The lower and upper bound of confidence intervals
    '''

    r, p = pearsonr(x,y)
    r_z = np.arctanh(r)
    se = 1/np.sqrt(x.size-3)
    z = stats.norm.ppf(1-alpha/2)
    lo_z, hi_z = r_z-z*se, r_z+z*se
    lo, hi = np.tanh((lo_z, hi_z))
    return r, p, lo, hi

def spearmanr_ci(x,y,alpha=0.05):
    ''' calculate Spearman correlation along with the confidence interval using scipy and numpy
    Parameters
    ----------
    x, y : iterable object such as a list or np.array
      Input for correlation calculation
    alpha : float
      Significance level. 0.05 by default
    Returns
    -------
    r : float
      Spearman's correlation coefficient
    pval : float
      The corresponding p value
    lo, hi : float
      The lower and upper bound of confidence intervals
    '''

    r, p = spearmanr(x,y)
    r_z = np.arctanh(r)
    se = 1/np.sqrt(x.size-3)
    z = stats.norm.ppf(1-alpha/2)
    lo_z, hi_z = r_z-z*se, r_z+z*se
    lo, hi = np.tanh((lo_z, hi_z))
    return r, p, lo, hi

# print(t['xray_sb_err']/t['xray_sb'])
# print(t['radio1_sb_err']/t['radio1_sb'])
# print(t['y_sb_err']/t['y_sb'])

def objective(x, a, b):
    return a * x + b


def fit(x, y):
    popt, _ = curve_fit(objective, x, y)
    a, b = popt
    x_line = np.arange(min(x) - 10, max(x) + 10, 0.01)
    y_line = objective(x_line, a, b)
    print('y = %.5f * x + %.5f' % (a, b))
    return x_line, y_line


def linreg(x, y):
    res = linregress(x, y)
    print(f'Slope is {res.slope} +- {res.stderr}')
    return res.slope, res.stderr

radio1 = t1['radio1_sb']
xray1 = t1['xray_sb']

radio_err1 = t1['radio1_sb_err']
xray_err1 = t1['xray_sb_err']

radio2 = t2['radio1_sb']
xray2 = t2['xray_sb']

radio_err2 = t2['radio1_sb_err']
xray_err2 = t2['xray_sb_err']

print('Number of cells used: '+str(len(t1)))
fitlinex1 = fit(np.log10(xray1), np.log10(radio1))
slopex1, errx1 = linreg(np.log10(xray1), np.log10(radio1))
pr = pearsonr_ci(np.log10(xray1), np.log10(radio1))
sr = spearmanr_ci(np.log10(xray1), np.log10(radio1))
print(pr)
print(f'Pearson R (x-ray vs radio): {pr[0]} +- {pr[-1]-pr[0]}')
print(sr)
print(f'Spearman R (x-ray vs radio): {sr[0]} +- {sr[-1]-sr[0]}')

print('Number of cells used: '+str(len(t2)))
fitlinex2 = fit(np.log10(xray2), np.log10(radio2))
slopex2, errx2 = linreg(np.log10(xray2), np.log10(radio2))
pr = pearsonr_ci(np.log10(xray2), np.log10(radio2))
sr = spearmanr_ci(np.log10(xray2), np.log10(radio2))
print(pr)
print(f'Pearson R (x-ray vs radio): {pr[0]} +- {pr[-1]-pr[0]}')
print(sr)
print(f'Spearman R (x-ray vs radio): {sr[0]} +- {sr[-1]-sr[0]}')


fig, ax = plt.subplots(constrained_layout=True)
ax.errorbar(np.log10(xray1), np.log10(radio1), xerr=(0.434*xray_err1/xray1), yerr=0.434*radio_err1/radio1, fmt='.', ecolor='red', elinewidth=0.4, color='darkred')
ax.errorbar(np.log10(xray2), np.log10(radio2), xerr=(0.434*xray_err2/xray2), yerr=0.434*radio_err2/radio2, fmt='.', ecolor='blue', elinewidth=0.4, color='darkblue')

ax.plot(fitlinex1[0], fitlinex1[1], color='darkred', linestyle='--')
ax.plot(fitlinex2[0], fitlinex2[1], color='darkblue', linestyle='--')

ax.set_ylim(np.min([np.min(np.log10(radio2) - (0.434*radio_err2 / radio2)),
                    np.min(np.log10(radio2) - (0.434*radio_err2 / radio2)),
                    np.min(np.log10(radio1) - (0.434*radio_err1 / radio1)),
                    np.min(np.log10(radio1) - (0.434*radio_err1 / radio1))]) - 0.1,
            np.max([np.max(np.log10(radio2) + (0.434*radio_err2 / radio2)),
                    np.max(np.log10(radio2) + (0.434*radio_err2 / radio2)),
                    np.max(np.log10(radio1) + (0.434 * radio_err1 / radio1)),
                    np.max(np.log10(radio1) + (0.434 * radio_err1 / radio1))
                    ]) + 0.1)
ax.set_xlim(np.min([np.min(np.log10(xray1) - (0.434 * xray_err1 / xray1)),
                    np.min(np.log10(xray1) - (0.434 * xray_err1 / xray1)),
                    np.min(np.log10(xray2) - (0.434 * xray_err2 / xray2)),
                    np.min(np.log10(xray2) - (0.434 * xray_err2 / xray2))]) - 0.1,
            np.max([np.max(np.log10(xray2) + (0.434 * xray_err2 / xray2)),
                    np.max(np.log10(xray2) + (0.434 * xray_err2/ xray2)),
                    np.max(np.log10(xray1) + (0.434 * xray_err1 / xray1)),
                    np.max(np.log10(xray1) + (0.434 * xray_err1 / xray1))
                    ]) + 0.1)
plt.grid(False)

ax.set_ylabel('log($I_{R}$) (Jy/arcsec$^{2}$)', fontsize=14)
# ax.set_xlabel('X-ray [SB/mean(SB)]')
ax.set_xlabel('log($I_{R}$) (counts/s/arcsec$^{2}$)', fontsize=14)
ax.legend(['Rudnick', 'UV-cut'], loc='upper left', fontsize=14)

plt.setp(ax.get_xticklabels(), fontsize=14)
plt.setp(ax.get_yticklabels(), fontsize=14)
plt.tight_layout()
plt.grid(False)
plt.savefig('analysis/combicorr'+obj+'.png', bbox_inches='tight')