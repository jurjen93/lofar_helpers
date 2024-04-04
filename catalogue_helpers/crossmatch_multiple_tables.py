"""
Code to crossmatch catalogues

Note that this code is in some parts of the main function hardcoded as it was used for crossmatching
4 catalogues obtained after reducing data from ELAIS-N1 with LOFAR.
So, feel free to adapt and use it for your own purposes.
"""

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table, vstack, hstack
from astropy.coordinates import match_coordinates_sky
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import warnings
import scienceplots
import math
from scipy.special import erf
from scipy.optimize import curve_fit

plt.style.use('science')

# Disable all warnings temporarily
warnings.filterwarnings("ignore")

def merge_with_table(catalog1, catalog2, res_2=0.6):
    """Merge with other table"""

    sep = res_2/2

    coords1 = SkyCoord(ra=catalog1['RA'], dec=catalog1['DEC'], unit=(u.deg, u.deg))
    coords2 = SkyCoord(ra=catalog2['RA'], dec=catalog2['DEC'], unit=(u.deg, u.deg))
    idx_catalog2, separation, _ = match_coordinates_sky(coords1, coords2)

    match_idxs = np.where(separation < sep * u.arcsec)[0]

    catalog2 = catalog2[idx_catalog2][match_idxs]

    for col in ['Total_flux', 'Peak_flux', 'E_Total_flux', 'E_Peak_flux', 'Maj', 'Min', 'PA', 'E_PA', 'S_Code' , 'Isl_rms']:
        if col != 'S_Code':
            catalog1[col + f"_{res_2}"] = np.nan
        else:
            catalog1[col + f"_{res_2}"] = ''
        catalog1[col + f"_{res_2}"][match_idxs] = catalog2[col]

    # print(f"Number of cross-matches between two catalogs {len(match_idxs)} ({int(len(match_idxs)/len(catalog1)*100)}% matched)")
    if len(catalog1)>0:
        perc = len(match_idxs)/len(catalog1)
    else:
        perc = 0
    return catalog1[match_idxs], perc, len(catalog1)

def detectibility(cats, outputfolder=None):
    DR1 = '/home/jurjen/Documents/ELAIS/catalogues/en1_final_cross_match_catalogue-v1.0.fits'
    DR2 = '/home/jurjen/Documents/ELAIS/catalogues/pybdsf_sources_6asec.fits'
    lotss = Table.read(DR2, format='fits')
    lotss['dist'] = list(map(dist_pointing_center, lotss['RA', 'DEC']))
    for n in range(len(cats)):
        # cat = cat[(cat['Total_flux']*1000>10) & (cat['S_Code']=='S')]
        cats[n]['dist'] = list(map(dist_pointing_center, cats[n]['RA', 'DEC']))
        cats[n] = cats[n][(cats[n]['dist'] < 1.25) & (cats[n]['Peak_flux'] > 5*cats[n]['Isl_rms'])]
    lotss = lotss[(lotss['dist'] < 1.25) & (lotss['Peak_flux'] > 5*lotss['Isl_rms'])]

    resolutions = ['0.3"', '0.6"', '1.2"', '6"']
    colors = ['darkred', 'darkblue', 'darkgreen', 'orange']
    linestyle = ['-.', ':', '-', '--']

    plt.close()
    xp = np.logspace(-1, 3, 17)
    for n, cat in enumerate(cats):
        output = []
        e_output = []
        x_axis = []
        for m, x in enumerate(xp):
            try:

                _, out, N = merge_with_table(lotss[(lotss['Total_flux']*1000>x) & (lotss['Total_flux']*1000<xp[m+1])], cat, 12)
                err = np.sqrt((out*(1-out))/N)
                err = np.sqrt(1/N)
                if out==0:
                    break
                output.append(out)
                e_output.append(err)
                x_axis.append((x+xp[m+1])/2)
            except:
                pass
        plt.fill_between(x_axis, output, np.array(output)+np.array(e_output), alpha=0.1, color=colors[n])
        plt.fill_between(x_axis, output, np.array(output)-np.array(e_output), alpha=0.1, color=colors[n])
        plt.plot(x_axis, output, label=resolutions[n], color=colors[n], linestyle=linestyle[n])
        plt.xscale('log')
    plt.ylim(0, 1.2)
    plt.legend()
    plt.xlabel('$S_{6}$ [mJy]')
    plt.ylabel('Detectability')
    plt.savefig(f'{outputfolder}/detectability.png', dpi=150)
    plt.close()

    lotss = lotss[(lotss['S_Code']=='S')]

    resolutions = ['0.3"', '0.6"', '1.2"', '6"']
    colors = ['darkred', 'darkblue', 'darkgreen', 'orange']
    linestyle = ['-.', ':', '-', '--']

    plt.close()
    xp = np.logspace(-1, 1, 15)
    for n, cat in enumerate(cats):
        output = []
        e_output = []
        x_axis = []
        for m, x in enumerate(xp):
            try:

                _, out, N = merge_with_table(lotss[(lotss['Total_flux']*1000>x) & (lotss['Total_flux']*1000<xp[m+1])], cat, 12)
                err = np.sqrt((out*(1-out))/N)
                err = np.sqrt(1/N)
                if out==0:
                    break
                output.append(out)
                e_output.append(err)
                x_axis.append((x+xp[m+1])/2)
            except:
                pass
        plt.fill_between(x_axis, output, np.array(output)+np.array(e_output), alpha=0.1, color=colors[n])
        plt.fill_between(x_axis, output, np.array(output)-np.array(e_output), alpha=0.1, color=colors[n])
        plt.plot(x_axis, output, label=resolutions[n], color=colors[n], linestyle=linestyle[n])
        plt.xscale('log')
    plt.plot([0.1, 10], [1, 1], linestyle='--', color='black')
    plt.ylim(0, 1.2)
    plt.xlim(0.1, 10)
    plt.legend()
    plt.xlabel('$S_{6}$ [mJy]')
    plt.ylabel('Detectability')
    plt.savefig(f'{outputfolder}/detectability_compact.png', dpi=150)
    plt.close()


def get_compact(c):
    snr_factor = 15
    compact = c[(c['S_Code'] == 'S')
                & (c['S_Code_0.6'] == 'S')
                & (c['S_Code_1.2'] == 'S')
                & (c['Total_flux']>snr_factor*c['Isl_rms'])
                & (c['Total_flux_0.6']>snr_factor*c['Isl_rms_0.6'])
                & (c['Total_flux_1.2']>snr_factor*c['Isl_rms_1.2'])]
    return compact

def get_not_compact(c):
    snr_factor = 15
    compact = c[(c['S_Code'] == 'M')
                & (c['S_Code_0.6'] == 'M')
                & (c['S_Code_1.2'] == 'M')
                & (c['Total_flux']>snr_factor*c['Isl_rms'])
                & (c['Total_flux_0.6']>snr_factor*c['Isl_rms_0.6'])
                & (c['Total_flux_1.2']>snr_factor*c['Isl_rms_1.2'])]
    return compact

def fluxration_smearing(central_freq_hz, bandwidth_hz, integration_time_s, resolution,
                                         distance_from_phase_center_deg):
    """
    Calculate the expected peak flux over integrated flux ratio due to smearing,
    as a function of distance from the pointing center.
    """

    # Convert distance from degrees to radians
    distance_from_phase_center_rad = np.deg2rad(distance_from_phase_center_deg)

    # Calculate angular resolution (radians)
    angular_resolution_rad = resolution*4.8481*1e-6
    gamma = np.sqrt(np.log(2))*2
    beta = (bandwidth_hz / central_freq_hz * distance_from_phase_center_rad / angular_resolution_rad)

    # Bandwidth Smearing Factor (formula 18-24 from Bridle 1999)
    bw_smearing = erf(gamma*beta/2)*np.sqrt(np.pi)/(gamma*beta)
    # Time Smearing Factor (formula 18-43 from Bridle 1999)
    time_smearing = 1-1.2288710615597145e-09*(integration_time_s*(distance_from_phase_center_rad) /
                            angular_resolution_rad)**2

    total_smearing = bw_smearing*time_smearing

    return time_smearing, bw_smearing, total_smearing

def ratio_err(total_flux, e_total_flux, peak_flux, e_peak_flux):
    return peak_flux/total_flux * np.sqrt((e_total_flux / total_flux) ** 2
                        + (e_peak_flux / peak_flux) ** 2)

def dist_pointing_center(pos):
    # HARD CODED FOR ELAIS-N1
    pointing_center = SkyCoord(242.75 * u.degree, 54.95 * u.degree, frame='icrs')
    return pointing_center.separation(SkyCoord(pos['RA'] * u.deg, pos['DEC'] * u.deg, frame='icrs')).value



def make_plots_combine(subcat, outputfolder=None):
    """Make plots"""


    print(len(subcat))

    def model(x, m, a, b):
        return m * x ** 2 + a * x + b

    ############# Peak flux over Total flux ##############
    R_03 = subcat['Peak_flux'] / subcat['Total_flux']
    E_R03 = ratio_err(subcat['Total_flux'], subcat['E_Total_flux'], subcat['Peak_flux'], subcat['E_Peak_flux'])
    R_06 = subcat['Peak_flux_0.6'] / subcat['Total_flux_0.6']
    E_R06 = ratio_err(subcat['Total_flux_0.6'], subcat['E_Total_flux_0.6'], subcat['Peak_flux_0.6'], subcat['E_Peak_flux_0.6'])
    R_12 = subcat['Peak_flux_1.2'] / subcat['Total_flux_1.2']
    E_R12 = ratio_err(subcat['Total_flux_1.2'], subcat['E_Total_flux_1.2'], subcat['Peak_flux_1.2'], subcat['E_Peak_flux_1.2'])
    # R_6 = subcat['Peak_flux_6'] / subcat['Total_flux_6']

    subcat['dist'] = list(map(dist_pointing_center, subcat['RA', 'DEC']))
    R_03 *=1.1


    # plt.figure(figsize=(5,4))
    plt.scatter(subcat['RA'] % 360, subcat['DEC'] % 360, c=R_03, s=10, alpha=0.75, vmin=0, vmax=1)
    plt.xlabel("Right Ascension (degrees)")
    plt.ylabel("Declination (degrees)")
    plt.colorbar(label='Peak / integrated flux')
    plt.savefig(f'{outputfolder}/peak_total_total_im_03.png', dpi=150)
    plt.close()

    # plt.figure(figsize=(5,4))
    plt.scatter(subcat['RA'] % 360, subcat['DEC'] % 360, c=R_06, s=10, alpha=0.75, vmin=0, vmax=1)
    plt.xlabel("Right Ascension (degrees)")
    plt.ylabel("Declination (degrees)")
    plt.colorbar(label='Peak / integrated flux')
    plt.savefig(f'{outputfolder}/peak_total_total_im_06.png', dpi=150)
    plt.close()

    # plt.figure(figsize=(5,4))
    plt.scatter(subcat['RA'] % 360, subcat['DEC'] % 360, c=R_12, s=10, alpha=0.75, vmin=0, vmax=1)
    plt.xlabel("Right Ascension (degrees)")
    plt.ylabel("Declination (degrees)")
    plt.colorbar(label='Peak / integrated flux')
    plt.savefig(f'{outputfolder}/peak_total_total_im_12.png', dpi=150)
    plt.close()

    # plt.figure(figsize=(5,4))
    # plt.scatter(subcat['RA'] % 360, subcat['DEC'] % 360, c=R_6, s=20, alpha=0.75, vmin=0, vmax=1)
    # plt.xlabel("Right Ascension (degrees)")
    # plt.ylabel("Declination (degrees)")
    # plt.colorbar(label='Peak / integrated flux')
    # plt.savefig(f'{outputfolder}/peak_total_total_im_6.png', dpi=150)
    # plt.close()

    subcat['dist'] = list(map(dist_pointing_center, subcat['RA', 'DEC']))
    # degree = 2

    centralhz = 140232849.121094
    bandwidth = 12207
    inttime = 1
    xdist = np.linspace(subcat['dist'].min(), 1.25, 100)
    ts, bws, totals = fluxration_smearing(centralhz, bandwidth, inttime, 0.3, xdist)


    params_03, cov_03 = curve_fit(model, subcat['dist'], R_03, sigma=E_R03, absolute_sigma=True)
    m_fit_03, a_fit_03, b_fit_03 = params_03
    params_03_upper, cov_03_upper = curve_fit(model, subcat['dist'], R_03+E_R03, sigma=E_R03, absolute_sigma=True)
    m_fit_03_upper,a_fit_03_upper, b_fit_03_upper = params_03_upper
    params_03_lower, cov_03_lower = curve_fit(model, subcat['dist'], R_03-E_R03, sigma=E_R03, absolute_sigma=True)
    m_fit_03_lower,a_fit_03_lower, b_fit_03_lower = params_03_lower

    params_06, cov_06 = curve_fit(model, subcat['dist'], R_06, sigma=E_R06, absolute_sigma=True)
    m_fit_06, a_fit_06, b_fit_06 = params_06
    params_06_upper, cov_06_upper = curve_fit(model, subcat['dist'], R_06+E_R06, sigma=E_R06, absolute_sigma=True)
    m_fit_06_upper, a_fit_06_upper, b_fit_06_upper = params_06_upper
    params_06_lower, cov_06_lower = curve_fit(model, subcat['dist'], R_06-E_R06, sigma=E_R06, absolute_sigma=True)
    m_fit_06_lower,a_fit_06_lower, b_fit_06_lower = params_06_lower

    params_12, cov_12 = curve_fit(model, subcat['dist'], R_12, sigma=E_R12, absolute_sigma=True)
    m_fit_12, a_fit_12, b_fit_12 = params_12
    params_12_upper, cov_12_upper = curve_fit(model, subcat['dist'], R_12+E_R12, sigma=E_R12, absolute_sigma=True)
    m_fit_12_upper, a_fit_12_upper, b_fit_12_upper = params_12_upper
    params_12_lower, cov_12_lower = curve_fit(model, subcat['dist'], R_12-E_R12, sigma=E_R12, absolute_sigma=True)
    m_fit_12_lower, a_fit_12_lower, b_fit_12_lower = params_12_lower

    x_fit = np.linspace(subcat['dist'].min(), subcat['dist'].max(), len(R_03))
    y_fit_03 = model(x_fit, a_fit_03, m_fit_03, b_fit_03)
    y_fit_03_upper = model(x_fit,a_fit_03_upper, m_fit_03_upper, b_fit_03_upper)
    y_fit_03_lower = model(x_fit,a_fit_03_lower, m_fit_03_lower, b_fit_03_lower)

    y_fit_06 = model(x_fit, a_fit_06, m_fit_06, b_fit_06)
    y_fit_06_upper = model(x_fit, a_fit_06_upper, m_fit_06_upper, b_fit_06_upper)
    y_fit_06_lower = model(x_fit,a_fit_06_lower, m_fit_06_lower, b_fit_06_lower)

    y_fit_12 = model(x_fit, a_fit_12, m_fit_12, b_fit_12)
    y_fit_12_upper = model(x_fit, a_fit_12_lower, m_fit_12_upper, b_fit_12_upper)
    y_fit_12_lower = model(x_fit, a_fit_12_lower, m_fit_12_lower, b_fit_12_lower)

    # plt.figure(figsize=(5,4))
    # plt.plot(x_fit, totals, color='darkblue', label='Theoretical peak response', linestyle='-.')
    plt.plot(x_fit, y_fit_03, color='darkred', label='0.3"', linestyle='-.')
    plt.fill_between(x_fit, y_fit_03_lower, y_fit_03_upper,
                             color='red', alpha=0.2)

    plt.plot(x_fit, y_fit_06, color='darkblue', label='0.6"', linestyle=':')
    plt.fill_between(x_fit, y_fit_06_lower, y_fit_06_upper,
                             color='blue', alpha=0.2)

    plt.plot(x_fit, y_fit_12, color='darkgreen', label='1.2"', linestyle='-')
    plt.fill_between(x_fit, y_fit_12_lower, y_fit_12_upper,
                             color='green', alpha=0.2)

    plt.plot(xdist, totals, color='black', label='$\Delta$t=1s', linestyle='--')
    ts, bws, totals = fluxration_smearing(centralhz, bandwidth, inttime*2, 0.33, xdist)
    plt.plot(xdist, totals, color='black', label='$\Delta$t=2s', linestyle=':')

    # plt.plot(x_fit, y_fit_6, color='orange', label='6"', linestyle='--')

    # plt.scatter(subcat['dist'], R, s=6, c=np.clip(subcat['Peak_flux']/subcat['Isl_rms'], a_min=5,  a_max=80), alpha=0.75)
    plt.xlabel("Distance from pointing center (degrees)")
    plt.ylabel("Peak / total intensity")
    plt.ylim(0.2, 1)

    plt.xlim(0, 1.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{outputfolder}/peak_total_total_dist.png', dpi=150)
    plt.close()


    centralhz = 140232849.121094
    bandwidth = 12207
    inttime = 1
    xdist = np.linspace(0, 1.25, 100)
    ts, bws, totals = fluxration_smearing(centralhz, bandwidth, inttime, 0.2232305194805194, xdist)
    plt.plot(xdist, ts, color='black', label='$\Delta$t=1s (0.2")', linestyle='-')
    ts, bws, totals = fluxration_smearing(centralhz, bandwidth, inttime*2, 0.2232305194805194, xdist)
    plt.plot(xdist, ts, color='black', label='$\Delta$t=2s (0.2")', linestyle='--')

    ts, bws, totals = fluxration_smearing(centralhz, bandwidth, inttime, 0.3, xdist)
    plt.plot(xdist, ts, color='black', label='$\Delta$t=1s (0.3")', linestyle='-.')
    ts, bws, totals = fluxration_smearing(centralhz, bandwidth, inttime*2, 0.3, xdist)
    plt.plot(xdist, ts, color='black', label='$\Delta$t=2s (0.3")', linestyle=':')

    plt.ylim(0, 1)
    plt.xlim(0, 1.25)
    plt.xlabel("Distance from pointing center (degrees)", size=14)
    plt.ylabel("I / I$_{0}$", size=14)
    plt.xticks(size=14)
    plt.yticks(size=14)
    plt.legend(fontsize=14)
    plt.tight_layout()
    plt.savefig(f'{outputfolder}/theoretical_smearing.png', dpi=150)

def circ_size(dist):
    return dist**2*np.pi

def make_plots_compare(cats, outputfolder=''):

    resolutions = ['0.3"', '0.6"', '1.2"', '6"']
    colors = ['darkred', 'darkblue', 'darkgreen', 'orange']
    linestyle = ['-.', ':', '-', '--']

    for n in range(len(cats)):
        # cat = cat[(cat['Total_flux']*1000>10) & (cat['S_Code']=='S')]
        cats[n]['dist'] = list(map(dist_pointing_center, cats[n]['RA', 'DEC']))

    plt.close()
    for n, cat in enumerate(cats):
        kwargs = dict(histtype='step', alpha=1, bins=np.logspace(-1.2, 1, 30), color=colors[n], linestyle=linestyle[n], label=resolutions[n], linewidth=1.3)
        df = cat.to_pandas()
        df = df[df['dist']<1.25]
        plt.hist(df['Total_flux']*1000, **kwargs)
    # plt.plot([10**-0.323, 10**-0.323], [0, 850], linestyle='--', color='black')
    plt.legend()
    plt.xlabel('$S$ [mJy]')
    plt.ylabel('Counts')
    # plt.ylim(0, 850)
    plt.xscale('log')
    plt.savefig(f'{outputfolder}/source_count_flux.png', dpi=150)
    plt.close()

    # plt.figure(figsize=(5, 4))
    for n, cat in enumerate(cats):
        counts, bin_edges = np.histogram(cat['dist'], bins=100)
        cumsum = np.cumsum(counts)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        plt.plot(bin_centers, cumsum, linestyle=linestyle[n], label=resolutions[n], color=colors[n])
    plt.xlim(0, 1.75)
    plt.ylim(0, )
    # plt.plot([1.25, 1.25], [0, 9500], color='black', linestyle='--')
    plt.xlabel("Distance from pointing center (degrees)")
    plt.ylabel("Cumulative count")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{outputfolder}/cumulative_dist_count.png', dpi=150)
    plt.close()

    for n, cat in enumerate(cats):
        counts, bin_edges = np.histogram(cat['dist'], bins=np.arange(0, 1.75, 0.2))
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        density = [c/(circ_size(bin_edges[n+1])-circ_size(bin_edges[n])) for n, c in enumerate(counts)]
        plt.plot(bin_centers[1:], density[1:], linestyle=linestyle[n], label=resolutions[n], color=colors[n], linewidth=1)
    plt.xlim(0, 1.75)
    plt.ylim(0, )
    # plt.plot([1.25, 1.25], [0, 9500], color='black', linestyle='--')
    plt.xlabel("Distance from pointing center (degrees)")
    plt.ylabel("Source density (degree$^{-2})$")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{outputfolder}/density.png', dpi=150)
    plt.close()

    # plt.figure(figsize=(5, 4))
    for n, cat in enumerate(cats):
        counts, bin_edges = np.histogram(cat['dist'], bins=500)
        cumsum = np.cumsum(counts)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        plt.plot(bin_centers, cumsum/max(cumsum), linestyle=linestyle[n], label=resolutions[n], color=colors[n])
    plt.xlim(0, 1.75)
    plt.ylim(0, )
    # plt.plot([1.25, 1.25], [0, 9500], color='black', linestyle='--')
    plt.xlabel("Distance from pointing center (degrees)")
    plt.ylabel("Cumulative distribution")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{outputfolder}/cumulative_dist.png', dpi=150)
    plt.close()

    # plt.figure(figsize=(5, 4))
    for n, cat in enumerate(cats):
        kwargs = dict(histtype='step', alpha=1, bins=np.linspace(0, 1.75, 20), color=colors[n], linestyle=linestyle[n], label=resolutions[n], linewidth=1)
        plt.hist(cat.to_pandas()['dist'], **kwargs)
    plt.xlim(0, 1.75)
    plt.ylim(0, )
    # plt.plot([1.25, 1.25], [0, 9500], color='black', linestyle='--')
    plt.xlabel("Distance from pointing center (degrees)")
    plt.ylabel("Count")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{outputfolder}/hist_dist.png', dpi=150)
    plt.close()

    # plt.figure(figsize=(5, 4))
    for n, cat in enumerate(cats):
        _, bin_edges = np.histogram(cat['dist'], bins=30)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        rms = []
        for i in range(len(bin_centers)):
            if i==0:
                bi = 0
            else:
                bi = bin_edges[i-1]
            rms.append(np.median(cat[(cat['dist']>bi) & (cat['dist']<bin_edges[i])]['Isl_rms'])*1000)
        params = np.polyfit(bin_centers[1:], rms[1:], 4)
        polynomial_function = np.poly1d(params)
        # plt.plot(bin_edges[1:], rms, linestyle=linestyle[n], label=resolutions[n], color=colors[n])
        xp = np.linspace(bin_centers[1], 1.75, 400)
        data = polynomial_function(xp)
        # data = [data[0]]+list(data[100:])
        # data = [d if d>data[m-1] else  for m, d in enumerate(data)][1:]
        plt.plot(xp, data, linestyle=linestyle[n], label=resolutions[n], color=colors[n], linewidth=1)
    plt.legend()
    plt.xlabel("Distance from pointing center (degree)")
    plt.ylabel("RMS (mJy/beam)")
    plt.xlim(0, 1.5)
    plt.ylim(0, )
    plt.tight_layout()
    plt.savefig(f'{outputfolder}/rms_dist.png', dpi=150)
    plt.close()

    # plt.figure(figsize=(5, 4))
    xp = np.linspace(bin_centers[1], 1.75, 400)
    for n, cat in enumerate(cats):
        _, bin_edges = np.histogram(cat['dist'], bins=30)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        rms = []
        for i in range(len(bin_centers)):
            if i==0:
                bi = 0
            else:
                bi = bin_edges[i-1]
            rms.append(np.median(cat[(cat['dist']>bi) & (cat['dist']<bin_edges[i])]['Isl_rms'])*1000)
        params = np.polyfit(bin_centers[1:], rms[1:], 4)
        polynomial_function = np.poly1d(params)
        # plt.plot(bin_edges[1:], rms, linestyle=linestyle[n], label=resolutions[n], color=colors[n])
        data = polynomial_function(xp)
        # data = [data[0]]+list(data[100:])
        # data = [d if d>data[m-1] else  for m, d in enumerate(data)][1:]
        plt.plot(xp, data/min(data), linestyle=linestyle[n], label=resolutions[n], color=colors[n], linewidth=1)
    plt.plot(xp, primary_beam_intensity_140(xp, 1.02), linestyle='--', label='$\\alpha=1.02$', color='black')
    plt.plot(xp, primary_beam_intensity_140(xp, 1.3), linestyle=':', label='$\\alpha=1.3$', color='black')
    plt.legend()
    plt.xlabel("Distance from pointing center (degree)")
    plt.ylabel("Relative RMS")
    plt.xlim(0, 1.5)
    plt.ylim(0.9, 3)
    plt.tight_layout()
    plt.savefig(f'{outputfolder}/rms_dist_relative.png', dpi=150)

def flux_ratios(T):
    T = T[T['Total_flux'] > 20 * T['Isl_rms']]
    R = T['Total_flux_6'] / T['Total_flux']
    print(f'median ratio: '+str(np.median(R[~np.isnan(R)])))
    print(f'std ratio: '+str(np.std(R[~np.isnan(R)])))
    print('Percentage above 1')
    print(len(R[R>1])/len(R))
    print('Percentage below 1')
    print(len(R[R<1])/len(R))

def primary_beam_intensity_140(theta, alpha):
    """
    Calculate the primary beam intensity as a function of angular distance.
    Parameters:
    - theta: Angular distance from the pointing center (in degrees).
    - fwhm: Full Width at Half Maximum of the beam (in degrees).
    Returns:
    - Intensity at the given angular distance, relative to the peak intensity.
    """
    # fwhm = 3/3600
    fwhm = alpha*2.173
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    intensity = np.exp(-0.5 * (theta / sigma) ** 2)
    return 1/intensity


def primary_beam_intensity_168(theta, alpha):
    """
    Calculate the primary beam intensity as a function of angular distance.
    Parameters:
    - theta: Angular distance from the pointing center (in degrees).
    - fwhm: Full Width at Half Maximum of the beam (in degrees).
    Returns:
    - Intensity at the given angular distance, relative to the peak intensity.
    """
    # fwhm = 3/3600
    fwhm = alpha*1.810998735777497
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    intensity = np.exp(-0.5 * (theta / sigma) ** 2)
    return 1/intensity


def parse_args():
    """Parse input arguments"""

    parser = argparse.ArgumentParser(description='Catalogue matching')
    parser.add_argument('--cat1', type=str, help='Catalogue 0.3', default=None)
    parser.add_argument('--cat2', type=str, help='Catalogue 0.6', default=None)
    parser.add_argument('--cat3', type=str, help='Catalogue 1.2', default=None)

    return parser.parse_args()



def main():
    """Main"""

    args = parse_args()
    catalog1 = Table.read(args.cat1)
    catalog2 = Table.read(args.cat2)
    catalog3 = Table.read(args.cat3)
    # flux_ratios(catalog1)
    # flux_ratios(catalog2)
    # flux_ratios(catalog3)
    # #
    # total = merge_with_table(catalog1, catalog2, res_2=0.6)[0]
    # total = merge_with_table(total, catalog3, res_2=1.2)[0]
    # compact = get_compact(total)
    # detectibility([catalog1, catalog2, catalog3], outputfolder='/home/jurjen/Documents/ELAIS/paperplots/')
    make_plots_combine(compact, outputfolder='/home/jurjen/Documents/ELAIS/paperplots/')
    make_plots_compare([catalog1, catalog2, catalog3], outputfolder='/home/jurjen/Documents/ELAIS/paperplots/')




if __name__ == '__main__':
    main()

# python catalogue_helpers/crossmatch_multiple_tables.py --cat1 final_merged_03.fits --cat2 final_merged_06.fits --cat3 final_merged_12.fits
