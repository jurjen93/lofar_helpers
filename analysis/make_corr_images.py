from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import warnings
from scipy.stats.stats import pearsonr, spearmanr, linregress
from past.utils import old_div
import astropy.units as u
from astropy.wcs import WCS
import argparse
from scipy import stats
import linmix

parser = argparse.ArgumentParser(
    description='Perform point-to-point analysis (radio/X or radio/radio/X) and save the results in a fits table. If (i) the X-ray counts < 0, (ii) the cell is not totally inside the X-ray FoV, and (ii) the radio flux density is < 0, it saves NaNs.')
parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
# REQUIRED
required.add_argument('-filein', type=str, required=True)
required.add_argument('-fileout', type=str, required=True)
required.add_argument('-noisefits', type=str, required=True)
required.add_argument('-no_y', action='store_true')
args = parser.parse_args()


warnings.filterwarnings('ignore')
plt.style.use('seaborn-deep')

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

f1 = fits.open(args.noisefits)
wcs =WCS(f1[0].header, naxis=2)
header = wcs.to_header()
rms = findrms(f1[0].data)/calc_beamarea(f1)/((header['CDELT2']*u.deg).to(u.arcsec)**2).value

f1.close()

f = fits.open(args.filein)
header = f[0].header
t = f[1].data

t = t[(t['radio1_sb']>2*rms)]


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

def linear(x, a, b):
    return a * x + b

def linearfit(x, y):
    popt, _ = curve_fit(linear, x, y)
    a, b = popt
    x_line = np.arange(min(x) - 10, max(x) + 10, 0.01)
    y_line = linear(x_line, a, b)
    print('y = %.5f * x + %.5f' % (a, b))
    return x_line, y_line

def wlinear_fit(x,y,w):
    """
    Fit (x,y,w) to a linear function, using exact formulae for weighted linear
    regression. This code was translated from the GNU Scientific Library (GSL),
    it is an exact copy of the function gsl_fit_wlinear.
    """
    # compute the weighted means and weighted deviations from the means
    # wm denotes a "weighted mean", wm(f) = (sum_i w_i f_i) / (sum_i w_i)
    W = np.sum(w)
    wm_x = np.average(x,weights=w)
    wm_y = np.average(y,weights=w)
    dx = x-wm_x
    dy = y-wm_y
    wm_dx2 = np.average(dx**2,weights=w)
    wm_dxdy = np.average(dx*dy,weights=w)
    # In terms of y = a + b x
    b = wm_dxdy / wm_dx2
    a = wm_y - wm_x*b
    cov_00 = (1.0/W) * (1.0 + wm_x**2/wm_dx2)
    cov_11 = 1.0 / (W*wm_dx2)
    cov_01 = -wm_x / (W*wm_dx2)
    # Compute chi^2 = \sum w_i (y_i - (a + b * x_i))^2
    chi2 = np.sum(((y-(a+b*x))**2)/(a+b*x))
    return a,b,cov_00,cov_11,cov_01,chi2

# def linearfit_w(x, y, err):
#     popt = wlinear_fit(x, y, err)
#     a, b = popt[1], popt[0]
#     x_line = np.arange(min(x) - 10, max(x) + 10, 0.01)
#     y_line = linear(x_line, a, b)
#     print('y = %.5f * x + %.5f' % (a, b))
#     return x_line, y_line


def linreg(x, y):
    res = linregress(x, y)
    print(f'Linear regression slope is {res.slope} +- {res.stderr}')
    return res.slope, res.stderr

def linmix_f(x, y, xerr, yerr):
    lm = linmix.LinMix(np.log10(x), np.log10(y), 0.434*xerr/x, 0.434*yerr/y)
    lm.run_mcmc(silent=True)
    x_line = np.arange(min(x) - 10, max(x) + 10, 0.01)
    y_line = linear(x_line, lm.chain['beta'].mean(), lm.chain['alpha'].mean())
    # print("{}, {}".format(lm.chain['alpha'].mean(), lm.chain['alpha'].std()))
    # print("{}, {}".format(lm.chain['beta'].mean(), lm.chain['beta'].std()))
    # print("{}, {}".format(lm.chain['sigsqr'].mean(), lm.chain['sigsqr'].std()))
    print('y = %.5f * x + %.5f' % (lm.chain['beta'].mean(), lm.chain['alpha'].mean()))
    print(f"Linear regression slope is {lm.chain['beta'].mean()} +- {lm.chain['beta'].std()}")

    return lm.chain['beta'].mean(), lm.chain['beta'].std(), lm.chain['alpha'].mean(), lm.chain['alpha'].std(), x_line, y_line


radio = t['radio1_sb']
xray = t['xray_sb']

radio_err = t['radio1_sb_err']
xray_err = t['xray_sb_err']

if not args.no_y:
    y = np.power(t['y_sb'], 2)
    y_err = 2*np.sqrt(y)*t['y_sb_err']/35

    radio_err/=np.mean(radio)
    xray_err/=np.mean(xray)
    y_err/=np.mean(y)

    radio/=np.mean(radio)
    xray/=np.mean(xray)
    y/=np.mean(y)

    mask = np.log10(y)>-1.25
    radio=radio[mask]
    xray=xray[mask]
    y=y[mask]
    radio_err=radio_err[mask]
    xray_err=xray_err[mask]
    y_err=y_err[mask]

print('Number of cells used: '+str(len(t)))

print('XRAY VERSUS RADIO')
# fitlinex = linearfit(np.log10(xray), np.log10(radio))
slopex, errx, offset, offseterr, fitlinex, fitliney = linmix_f(xray, radio, xray_err, radio_err)
pr = pearsonr_ci(np.log10(xray), np.log10(radio))
sr = spearmanr_ci(np.log10(xray), np.log10(radio))
print(pr)
print(f'Pearson R (x-ray vs radio): {pr[0]} +- {pr[-1]-pr[0]}')
print(sr)
print(f'Spearman R (x-ray vs radio): {sr[0]} +- {sr[-1]-sr[0]}')
if not args.no_y:
    print('Y VERSUS RADIO')
    fitliney = linearfit(np.log10(y), np.log10(radio))
    slopex, errx = linreg(np.log10(y), np.log10(radio))
    pr = pearsonr_ci(np.log10(y), np.log10(radio))
    sr = spearmanr_ci(np.log10(y), np.log10(radio))
    print(pr)
    print(f'Pearson R (x-ray vs radio): {pr[0]} +- {pr[-1]-pr[0]}')
    print(sr)
    print(f'Spearman R (x-ray vs radio): {sr[0]} +- {sr[-1] - sr[0]}')

fig, ax = plt.subplots(constrained_layout=True)
ax.errorbar(np.log10(xray), np.log10(radio), xerr=(0.434*xray_err/xray), yerr=0.434*radio_err/radio, fmt='.', ecolor='red', elinewidth=0.4, color='darkred', capsize=2, markersize=4)
if not args.no_y:
    ax.errorbar(np.log10(y), np.log10(radio), xerr=(0.434*y_err/y), yerr=0.434*radio_err/radio, fmt='.', ecolor='blue', elinewidth=0.4, color='darkblue', capsize=2, markersize=4)
ax.plot(fitlinex, fitliney, color='darkred', linestyle='--')
if not args.no_y:
    ax.plot(fitliney[0], fitliney[1], color='darkblue', linestyle='--')
ax.set_ylim(np.min([np.min(np.log10(radio) - (0.434*radio_err / radio)),
                    np.min(np.log10(radio) - (0.434*radio_err / radio))]) - 0.1,
            np.max([np.max(np.log10(radio) + (0.434*radio_err / radio)),
                    np.max(np.log10(radio) + (0.434*radio_err / radio))]) + 0.1)
if not args.no_y:
    ax.set_xlim(np.min([np.min(np.log10(y) - (0.434*y_err / y)),
                        np.min(np.log10(xray) - (0.434*xray_err / xray))]) - 0.1,
                np.max([np.max(np.log10(y) + (0.434*y_err / y)),
                        np.max(np.log10(xray) + (0.434*xray_err / xray))]) + 0.1)
else:
    ax.set_xlim(np.min([np.min(np.log10(xray) - (0.434 * xray_err / xray)),
                        np.min(np.log10(xray) - (0.434 * xray_err / xray))]) - 0.1,
                np.max([np.max(np.log10(xray) + (0.434 * xray_err / xray)),
                        np.max(np.log10(xray) + (0.434 * xray_err / xray))]) + 0.1)
plt.grid(False)
if not args.no_y:

    ax.set_ylabel('log($I_{R}$) [SB/mean(SB)]', fontsize=14)
    # ax.set_xlabel('X-ray [SB/mean(SB)]')
    ax.set_xlabel('log($I_{X}$) [SB/mean(SB)] and log($I_{SZ}^{2}$) [SZ/mean(SZ)]', fontsize=14)
    ax.legend(['Radio vs. X-ray', 'Radio vs. SZ'], loc='upper left', fontsize=14)
else:
    ax.set_ylabel('log($I_{R}$) (Jy/arcsec$^{2}$)', fontsize=14)
    # ax.set_xlabel('X-ray [SB/mean(SB)]')
    ax.set_xlabel('log($I_{R}$) (counts/s/arcsec$^{2}$)', fontsize=14)
plt.setp(ax.get_xticklabels(), fontsize=14)
plt.setp(ax.get_yticklabels(), fontsize=14)
plt.tight_layout()
plt.grid(False)
plt.savefig(args.fileout, bbox_inches='tight')