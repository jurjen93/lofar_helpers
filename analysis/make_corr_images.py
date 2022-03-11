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

f1 = fits.open('fits/60rudnick.fits')
wcs =WCS(f1[0].header, naxis=2)
header = wcs.to_header()
rms = findrms(f1[0].data)/calc_beamarea(f1)/((header['CDELT2']*u.deg).to(u.arcsec)**2).value

f1.close()

f = fits.open(args.filein)
header = f[0].header
t = f[1].data

t = t[(t['radio1_sb']>3*rms)]

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

radio = t['radio1_sb']
xray = t['xray_sb']
y = np.power(t['y_sb'], 2)

radio_err = t['radio1_sb_err']
xray_err = t['xray_sb_err']
y_err = 2*np.sqrt(y)*t['y_sb_err']

radio_err/=np.mean(radio)
xray_err/=np.mean(xray)/6
y_err/=np.mean(y)*31

radio/=np.mean(radio)
xray/=np.mean(xray)
y/=np.mean(y)

mask = np.log(y)>-1.25
radio=radio[mask]
xray=xray[mask]
y=y[mask]
radio_err=radio_err[mask]
xray_err=xray_err[mask]
y_err=y_err[mask]

print('Number of cells used: '+str(len(t)))

print('XRAY VERSUS RADIO')
fitlinex = fit(np.log10(xray), np.log10(radio))
slopex, errx = linreg(np.log10(xray), np.log10(radio))
print('Pearson R (x-ray vs radio): '+str(pearsonr(np.log(xray), np.log(radio))))
print('Spearman R (x-ray vs radio): '+str(spearmanr(np.log(xray), np.log(radio))))

print('Y VERSUS RADIO')
fitliney = fit(np.log10(y), np.log10(radio))
slopex, errx = linreg(np.log10(y), np.log10(radio))
print('Pearson R (ymap vs radio): '+str(pearsonr(np.log(y), np.log(radio))))
print('Spearman R (ymap vs radio): '+str(spearmanr(np.log(y), np.log(radio))))

fig, ax = plt.subplots(constrained_layout=True)
ax.errorbar(np.log10(xray), np.log10(radio), xerr=(0.434*xray_err/xray), yerr=0.434*radio_err/radio, fmt='.', ecolor='red', elinewidth=0.4, color='darkred')
ax.errorbar(np.log10(y), np.log10(radio), xerr=(0.434*y_err/y), yerr=0.434*radio_err/radio, fmt='.', ecolor='blue', elinewidth=0.4, color='darkblue')
ax.plot(fitlinex[0], fitlinex[1], color='darkred', linestyle='--')
ax.plot(fitliney[0], fitliney[1], color='darkblue', linestyle='--')
ax.set_ylim(np.min([np.min(np.log10(radio) - (0.434*radio_err / radio)),
                    np.min(np.log10(radio) - (0.434*radio_err / radio))]) - 0.1,
            np.max([np.max(np.log10(radio) + (0.434*radio_err / radio)),
                    np.max(np.log10(radio) + (0.434*radio_err / radio))]) + 0.1)
ax.set_xlim(np.min([np.min(np.log10(y) - (0.434*y_err / y)),
                    np.min(np.log10(xray) - (0.434*xray_err / xray))]) - 0.1,
            np.max([np.max(np.log10(y) + (0.434*y_err / y)),
                    np.max(np.log10(xray) + (0.434*xray_err / xray))]) + 0.1)
plt.grid(False)
ax.set_ylabel('log($I_{R}$) [SB/mean(SB)]', fontsize=14)
# ax.set_xlabel('X-ray [SB/mean(SB)]')
ax.set_xlabel('log($I_{X}$) [SB/mean(SB)] and log(y) [SZ/mean(SZ)]', fontsize=14)
ax.legend(['Radio vs. X-ray', 'Radio vs. SZ'], loc='upper left', fontsize=14)
plt.setp(ax.get_xticklabels(), fontsize=14)
plt.setp(ax.get_yticklabels(), fontsize=14)
plt.tight_layout()
plt.grid(False)
plt.savefig(args.fileout, bbox_inches='tight')