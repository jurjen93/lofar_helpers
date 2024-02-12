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

    print(f"Number of cross-matches between two catalogs {len(match_idxs)} ({int(len(match_idxs)/len(catalog1)*100)}% matched)")

    return catalog1[match_idxs]

def get_compact(c):
    snr_factor = 15
    compact = c[(c['S_Code'] == 'S')
                & (c['S_Code_0.6'] == 'S')
                & (c['S_Code_1.2'] == 'S')
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


def make_plots(subcat, outputfolder=None):
    """Make plots"""

    def dist_pointing_center(pos):
        # HARD CODED FOR ELAIS-N1
        pointing_center = SkyCoord(242.75 * u.degree, 54.95 * u.degree, frame='icrs')
        return pointing_center.separation(SkyCoord(pos['RA'] * u.deg, pos['DEC'] * u.deg, frame='icrs')).value

    print(len(subcat))

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


    plt.figure(figsize=(5,4))
    plt.scatter(subcat['RA'] % 360, subcat['DEC'] % 360, c=R_03, s=20, alpha=0.75, vmin=0, vmax=1)
    plt.xlabel("Right Ascension (degrees)")
    plt.ylabel("Declination (degrees)")
    plt.colorbar(label='Peak / integrated flux')
    plt.savefig(f'{outputfolder}/peak_total_total_im_03.png', dpi=150)
    plt.close()

    plt.figure(figsize=(5,4))
    plt.scatter(subcat['RA'] % 360, subcat['DEC'] % 360, c=R_06, s=20, alpha=0.75, vmin=0, vmax=1)
    plt.xlabel("Right Ascension (degrees)")
    plt.ylabel("Declination (degrees)")
    plt.colorbar(label='Peak / integrated flux')
    plt.savefig(f'{outputfolder}/peak_total_total_im_06.png', dpi=150)
    plt.close()

    plt.figure(figsize=(5,4))
    plt.scatter(subcat['RA'] % 360, subcat['DEC'] % 360, c=R_12, s=20, alpha=0.75, vmin=0, vmax=1)
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

    def model(x, m, a, b):
        return m * x**2 + a*x + b

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

    centralhz = 140232849.121094
    bandwidth = 12207
    inttime = 1
    ts, bws, totals = fluxration_smearing(centralhz, bandwidth, inttime, 0.3, np.linspace(subcat['dist'].min(), subcat['dist'].max(), 100))

    plt.figure(figsize=(5,4))
    # plt.plot(x_fit, totals, color='darkblue', label='Theoretical peak response', linestyle='-.')
    plt.plot(x_fit, y_fit_03, color='darkred', label='0.3"', linestyle='-.')
    plt.fill_between(x_fit, y_fit_03_lower, y_fit_03_upper,
                             color='red', alpha=0.2)

    plt.plot(x_fit, y_fit_06, color='darkblue', label='0.6"', linestyle=':')
    plt.fill_between(x_fit, y_fit_06_lower, y_fit_06_upper,
                             color='blue', alpha=0.2)

    plt.plot(x_fit, y_fit_12, color='darkgreen', label='1.2"', linestyle='--')
    plt.fill_between(x_fit, y_fit_12_lower, y_fit_12_upper,
                             color='green', alpha=0.2)
    # plt.plot(x_fit, y_fit_6, color='orange', label='6"', linestyle='--')

    # plt.scatter(subcat['dist'], R, s=6, c=np.clip(subcat['Peak_flux']/subcat['Isl_rms'], a_min=5,  a_max=80), alpha=0.75)
    plt.xlabel("Distance from pointing center (degrees)")
    plt.ylabel("Peak / integrated flux")
    plt.ylim(0.2, 1)

    plt.xlim(0, 1.25)
    plt.legend()
    plt.savefig(f'{outputfolder}/peak_total_total_dist.png', dpi=150)
    plt.close()

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
    total = merge_with_table(catalog1, catalog2, res_2=0.6)
    total = merge_with_table(total, catalog3, res_2=1.2)
    total = get_compact(total)
    make_plots(total, outputfolder='/home/jurjen/Documents/ELAIS/paperplots/')




if __name__ == '__main__':
    main()

# python catalogue_helpers/crossmatch.py --cat1 /home/jurjen/Documents/ELAIS/catalogues/finalcat03/*.fits --cat2 /home/jurjen/Documents/ELAIS/catalogues/pybdsf_sources_6asec.fits --out_table final_merged_03.fits
# python catalogue_helpers/crossmatch.py --cat1 /home/jurjen/Documents/ELAIS/catalogues/finalcat06/*.fits --cat2 /home/jurjen/Documents/ELAIS/catalogues/pybdsf_sources_6asec.fits --out_table final_merged_06.fits --resolution 0.6
# python catalogue_helpers/crossmatch.py --cat1 /home/jurjen/Documents/ELAIS/catalogues/finalcat12/*.fits --cat2 /home/jurjen/Documents/ELAIS/catalogues/pybdsf_sources_6asec.fits --out_table final_merged_12.fits --resolution 1.2
