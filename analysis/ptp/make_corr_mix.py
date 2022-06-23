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
from astropy.modeling import models, fitting
import nmmn.stats
import bces.bces
import linmix

np.random.seed(10)

parser = argparse.ArgumentParser(
    description='make corr plots')
parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
# REQUIRED
required.add_argument('-obj', type=str, required=True)
args = parser.parse_args()

obj = args.obj

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

def linear(x, a, b):
    return a * x + b

def linearfit(x, y):
    popt, _ = curve_fit(linear, x, y)
    a, b = popt
    x_line = np.arange(min(x) - 10, max(x) + 10, 0.01)
    y_line = linear(x_line, a, b)
    print('y = %.5f * x + %.5f' % (a, b))
    return x_line, y_line

def linearfit_w(x, y, err):
    popt = wlinear_fit(x, y, err)
    a, b = popt[1], popt[0]
    x_line = np.arange(min(x) - 10, max(x) + 10, 0.01)
    y_line = linear(x_line, a, b)
    print('y = %.5f * x + %.5f' % (a, b))
    return x_line, y_line


def linreg(x, y):
    res = linregress(x, y)
    print(f'Slope is {res.slope} +- {res.stderr}')
    return res.slope, res.stderr

def linreg_unc(x,y, y_err):
    fit = fitting.LinearLSQFitter()
    line_init = models.Linear1D()
    fitted_line = fit(line_init, x, y, weights=1/y_err)
    return fitted_line(x)

def func(x):
    return x[1]*x[0]+x[2]


def bcesfit(x, xerr, y, yerr, i=3):
    # Values for i documentation
    # 0 y|x Assumes x as the independent variable
    # 1 x|y Assumes y as the independent variable
    # 2 bissector Line that bisects the y|x and x|y. This approach is self-inconsistent.
    # 3 orthogonal Orthogonal least squares: line that minimizes orthogonal distances. Should be used when it is not clear which variable should be treated as the independent one

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
    nboot = 10000

    xlog = np.log10(x)
    ylog = np.log10(y)
    xerrlog = xerr / (x * np.log(10.))
    yerrlog = yerr / (y * np.log(10.))

    slope, offset, slopeerr, offseterr, covab = bces.bces.bcesp(xlog, xerrlog, ylog, yerrlog, cov, nboot)

    # array with best-fit parameters
    fitm = np.array([slope[i], offset[i]])
    # covariance matrix of parameter uncertainties
    covm = np.array([(slopeerr[i] ** 2, covab[i]), (covab[i], offseterr[i] ** 2)])
    # Gets lower and upper bounds on the confidence band

    xvec = np.linspace(xlog.min() - 1., xlog.max() + 1.)
    # lcb,ucb,xcb=nmmn.stats.confbandnl(xlog,ylog,func,fitm,covm,2,0.954,xvec)

    lcb, ucb, xcb = nmmn.stats.confbandnl(xlog, ylog, func, fitm, covm, 2, 0.95, xvec)

    return slope[i], slopeerr[i], offset[i], offseterr[i], lcb, ucb, xcb

def linmix_f(x, y, xerr, yerr):
    alpha = 1.0
    beta = 0.5
    lm = linmix.LinMix(np.log10(x), np.log10(y), 0.434*xerr/x, 0.434*yerr/y)
    lm.run_mcmc(silent=True)
    print("{}, {}".format(lm.chain['alpha'].mean(), lm.chain['alpha'].std()))
    print("{}, {}".format(lm.chain['beta'].mean(), lm.chain['beta'].std()))
    print("{}, {}".format(lm.chain['sigsqr'].mean(), lm.chain['sigsqr'].std()))
    # plt.scatter(x, y, alpha=0.5)
    # plt.errorbar(x, y, xerr=xerr, yerr=yerr, ls=' ', alpha=0.5)
    # for i in range(0, len(lm.chain), 25):
    #     xs = np.arange(-10, 11)
    #     ys = lm.chain[i]['alpha'] + xs * lm.chain[i]['beta']
    #     plt.plot(xs, ys, color='r', alpha=0.02)
    # ys = alpha + xs * beta
    # plt.plot(xs, ys, color='k')
    # plt.show
    return lm.chain['beta'].mean(), lm.chain['beta'].std(), lm.chain['alpha'].mean(), lm.chain['alpha'].std()

f1 = fits.open('fits/60rudnick.fits')
wcs =WCS(f1[0].header, naxis=2)
header = wcs.to_header()
rms = findrms(f1[0].data)/calc_beamarea(f1)/((header['CDELT2']*u.deg).to(u.arcsec)**2).value
f1.close()
f = fits.open(f'ptp_results_41/{obj}_results_rudnick_1.fits')
header = f[0].header
t1 = f[1].data
t1 = t1[(t1['radio1_sb']>2*rms)]

f1 = fits.open('fits/60cleanbridge_200kpc.fits')
wcs =WCS(f1[0].header, naxis=2)
header = wcs.to_header()
rms = findrms(f1[0].data)/calc_beamarea(f1)/((header['CDELT2']*u.deg).to(u.arcsec)**2).value
f1.close()
f = fits.open(f'ptp_results_41/{obj}_results_cb_1.fits')
header = f[0].header
t2 = f[1].data
t2 = t2[(t2['radio1_sb']>2*rms)]

radio1 = t1['radio1_sb']
xray1 = t1['xray_sb']

radio_err1 = t1['radio1_sb_err']
xray_err1 = t1['xray_sb_err']

radio2 = t2['radio1_sb']
xray2 = t2['xray_sb']

radio_err2 = t2['radio1_sb_err']
xray_err2 = t2['xray_sb_err']

# print('Number of cells used: '+str(len(t1)))
# fitlinex1 = linearfit_w(np.log10(xray1), np.log10(radio1), 0.434*radio_err1/radio1)
# slopex1, errx1 = linreg(np.log10(xray1), np.log10(radio1))
#
# pr = pearsonr_ci(np.log10(xray1), np.log10(radio1))
# sr = spearmanr_ci(np.log10(xray1), np.log10(radio1))
# print(f'Pearson R (x-ray vs radio): {pr[0]} +- {pr[-1]-pr[0]}')
# print(f'Spearman R (x-ray vs radio): {sr[0]} +- {sr[-1]-sr[0]}')

xray=np.append(xray1,xray2)
radio=np.append(radio1,radio2)
radio_err =np.append(radio_err1, radio_err2)
xray_err = np.append(xray_err1, xray_err2)

print('Number of cells used: '+str(len(t2)))
fitlinex2 = linearfit_w(np.log10(xray), np.log10(radio), 0.434*radio_err/radio)
slopex2, errx2 = linreg(np.log10(xray), np.log10(radio))
# print(wlinear_fit(np.log(xray1), np.log(radio1), 0.434*radio_err1/radio1))
print(wlinear_fit(np.log(xray), np.log(radio), 0.434*radio_err/radio))
pr = pearsonr_ci(np.log10(xray), np.log10(radio))
sr = spearmanr_ci(np.log10(xray), np.log10(radio))
print(f'Pearson R (x-ray vs radio): {pr[0]} +- {pr[-1]-pr[0]}')
print(f'Spearman R (x-ray vs radio): {sr[0]} +- {sr[-1]-sr[0]}')

slope, slopeerr, offset, offseterr,  lc, uc, xc = bcesfit(xray, xray_err, radio, radio_err, i=3)
# slope, slopeerr, offset, offseterr = linmix_f(np.log10(xray), np.log10(radio), 0.434*xray_err/xray, 0.434*radio_err/radio)


# print ('(slope, slope_err, offset, offset_err)',  slope, slopeerr, offset, offseterr)


fig, ax = plt.subplots(constrained_layout=True)
ax.errorbar(np.log10(xray1), np.log10(radio1), xerr=(0.434*xray_err1/xray1), yerr=0.434*radio_err1/radio1, fmt='.', ecolor='red', elinewidth=0.4, color='darkred', capsize=2, markersize=4, label='R02')
ax.errorbar(np.log10(xray2), np.log10(radio2), xerr=(0.434*xray_err2/xray2), yerr=0.434*radio_err2/radio2, fmt='.', ecolor='blue', elinewidth=0.4, color='darkblue', capsize=2, markersize=4, label='uv-subtract')

# ax.plot(fitlinex1[0], fitlinex1[1], color='darkred', linestyle='--')
# ax.plot(fitlinex2[0], fitlinex2[1], color='darkblue', linestyle='--', label='Best fit')
x_line = np.arange(- 10, 10, 0.01)
y_line = linear(x_line, slope, offset)
ax.plot(x_line, y_line, color='black', linestyle='--', label='Best fit')

ax.set_ylim(np.log10(rms), -5)
# ax.set_xlim(np.min([np.min(np.log10(xray1) - (0.434 * xray_err1 / xray1)),
#                     np.min(np.log10(xray1) - (0.434 * xray_err1 / xray1)),
#                     np.min(np.log10(xray2) - (0.434 * xray_err2 / xray2)),
#                     np.min(np.log10(xray2) - (0.434 * xray_err2 / xray2))]) - 0.1,
#             np.max([np.max(np.log10(xray2) + (0.434 * xray_err2 / xray2)),
#                     np.max(np.log10(xray2) + (0.434 * xray_err2/ xray2)),
#                     np.max(np.log10(xray1) + (0.434 * xray_err1 / xray1)),
#                     np.max(np.log10(xray1) + (0.434 * xray_err1 / xray1))
#                     ]) + 0.1)
ax.set_xlim(-6.8, -4.4)
plt.grid(False)

ax.plot((-10, -3), (np.log10(2*rms), np.log10(2*rms)), linestyle='--', color='orange', label='$2\sigma$')
# ax.fill_between(xc, lc, uc, alpha=0.2, facecolor='green')

ax.set_ylabel('log($I_{R}$) (Jy/arcsec$^{2}$)', fontsize=14)
# ax.set_xlabel('X-ray [SB/mean(SB)]')
ax.set_xlabel('log($I_{X}$) (counts/s/arcsec$^{2}$)', fontsize=14)
ax.legend(loc='upper left', fontsize=14)

plt.setp(ax.get_xticklabels(), fontsize=14)
plt.setp(ax.get_yticklabels(), fontsize=14)
plt.tight_layout()
plt.grid(False)
# plt.xticks([-6.7, -6.4, -6.1,-5.8, -5.5, -5.1, -4.8])
if 'trail' in args.obj:
    plt.title('A399 northwest extension')
elif 'bridge' in args.obj:
    plt.title('Radio bridge')
else:
    plt.title(args.obj.title())
plt.savefig('analysis/combicorr'+obj+'.png', bbox_inches='tight')