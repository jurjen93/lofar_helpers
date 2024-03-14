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

plt.style.use('science')

# Disable all warnings temporarily
warnings.filterwarnings("ignore")

def find_matches(cat1, cat2, separation_asec):
    """
    Find crossmatches with two catalogues
    """

    catalog1 = Table.read(cat1, format='fits')
    catalog2 = Table.read(cat2, format='fits')

    # Define the celestial coordinates for each catalog
    coords1 = SkyCoord(ra=catalog1['RA'], dec=catalog1['DEC'], unit=(u.deg, u.deg))
    coords2 = SkyCoord(ra=catalog2['RA'], dec=catalog2['DEC'], unit=(u.deg, u.deg))

    # Perform the crossmatch using Astropy's match_coordinates_sky function
    idx_catalog2, separation, _ = match_coordinates_sky(coords1, coords2)

    # Define a maximum separation threshold (adjust as needed)
    max_sep_threshold = separation_asec * u.arcsec
    matched_sources_mask = separation < max_sep_threshold

    # You can use the mask to filter the matched sources in catalog1
    matched_sources_catalog1 = catalog1[matched_sources_mask]
    matched_sources_catalog2 = catalog2[idx_catalog2[matched_sources_mask]]

    return matched_sources_catalog1, matched_sources_catalog2


def separation_match(cat1, cat2, separation_asec):
    """
    Find non-crossmatches between two catalogues and remove those from catalogue 1 based on distance in arcsec

    :param cat1: catalogue 1
    :param cat2: catalogue 2
    :param separation_asec: max separation between catalogue matches in arcsec

    :return: corrected catalogue 1, removed sources from catalogue 1
    """

    catalog1 = Table.read(cat1, format='fits')
    catalog2 = Table.read(cat2, format='fits')

    # Define the celestial coordinates for each catalog
    coords1 = SkyCoord(ra=catalog1['RA'], dec=catalog1['DEC'], unit=(u.deg, u.deg))
    coords2 = SkyCoord(ra=catalog2['RA'], dec=catalog2['DEC'], unit=(u.deg, u.deg))

    # Perform the crossmatch using Astropy's match_coordinates_sky function
    idx_catalog2, separation, _ = match_coordinates_sky(coords1, coords2)
    catalog1['separation_cat2'] = separation

    # Define a maximum separation threshold (adjustplots as needed)
    max_sep_threshold_large = separation_asec * u.arcsec

    # Sources below RMS threshold and above separation threshold
    non_matches_large = catalog1[(catalog1['separation_cat2'] > max_sep_threshold_large)]

    non_match_mask_large = np.sum([non_matches_large['Source_id'] == i for i in catalog1['Source_id']], axis=1).astype(bool)

    non_match_mask = non_match_mask_large

    catalog1_corrected = catalog1[~non_match_mask]
    catalog1_removed = catalog1[non_match_mask]

    print(f'Source count after cross-match deep field DR2: {len(catalog1_corrected)} '
          f'({int(round(1 - len(catalog1_corrected) / len(catalog1),2)*100)}% removed)')

    return catalog1_corrected, catalog1_removed


def remove_snr(catalog, snr=5, cat=None):
    """Remove sources based on Signal-To-Noise"""

    print(f"Source count before SNR cut {len(catalog)}")
    newcat = catalog[catalog['Peak_flux'] > snr * catalog['Isl_rms']]
    print(f'Source count after SNR cut {len(newcat)} ({int(abs(1-len(catalog)/len(newcat))*100)}% removed)')
    if cat is not None:
        for name, dens in density.items():
            if name in cat:
                print(name, cat)
                print(f'Density: {len(newcat)/dens} sources/arcsec**2')

    return newcat


def crossmatch_itself(catalog, min_sep=0.15):
    """Crossmatch with itself to find nearest neighbour"""

    coords = SkyCoord(ra=catalog['RA'], dec=catalog['DEC'], unit=(u.deg, u.deg))
    idx_catalog, separation, _ = match_coordinates_sky(coords, coords, nthneighbor=2)
    nearest_neighbour = separation<min_sep*u.arcsec
    print(catalog[nearest_neighbour]['Cat_id'])


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


def make_plots(cat, res=0.3, outputfolder=None):
    """Make plots"""

    def dist_pointing_center(pos):
        # HARD CODED FOR ELAIS-N1
        pointing_center = SkyCoord(242.75 * u.degree, 54.95 * u.degree, frame='icrs')
        return pointing_center.separation(SkyCoord(pos['RA'] * u.deg, pos['DEC'] * u.deg, frame='icrs')).value

    # make peak flux plot
    plt.hist(np.log10(cat['Peak_flux'] * 1000), bins=30)
    plt.xlabel('Peak flux (mJy/beam)')
    plt.savefig('peak_flux.png')
    plt.close()

    # make total flux plot
    plt.hist(np.log10(cat['Total_flux'] * 1000), bins=30)
    plt.xlabel('Total flux (mJy)')
    plt.savefig('total_flux.png', dpi=150)
    plt.close()

    ############### 2D HISTOGRAM ###############
    # 6arcsec offset
    subcat = cat[(cat[f'dRA_{res}'] == cat[f'dRA_{res}'])
                 & (cat['S_Code'] == 'S')
                 & (cat['S_Code_6'] == 'S')
                 & (cat['Total_flux_6']*1000 > 0.5)
                 & (cat['Total_flux']*1000 > 0.5)]
    print(f"Number of sources for hist2d plot: {len(subcat)}")

    plt.hexbin(subcat[f'dRA_{res}'] * 3600,
               subcat[f'dDEC_{res}'] * 3600,
               gridsize=30, cmap='viridis', norm=LogNorm())
    plt.xlabel("dRA (arcsec)")
    plt.ylabel("dDEC (arcsec)")
    plt.colorbar(label='Source count')

    medianedRA = round(np.mean(subcat[f'E_dRA_{res}']) * 3600, 4)
    medianedDEC = round(np.mean(subcat[f'E_dDEC_{res}']) * 3600, 4)
    print(f"Mean E_dRA: {medianedRA} "
              f"Mean E_dDEC: {medianedDEC}")
    mediandRA = round(np.median(subcat[f'dRA_{res}']) * 3600, 4)
    mediandDEC = round(np.median(subcat[f'dDEC_{res}']) * 3600, 4)
    print(f"Mean dRA: {mediandRA} "
              f"Mean dDEC: {mediandDEC}")

    axs = np.arange(-16, 17)
    plt.plot(axs, [mediandDEC]*len(axs), color='black', linestyle='--')
    plt.plot([mediandRA]*len(axs), axs, color='black', linestyle='--')
    plt.xlim(-4, 4)
    plt.ylim(-4, 4)
    plt.savefig(f'{outputfolder}/dRA_dDEC_{res}.png', dpi=150)
    plt.close()

    ############# Flux ratio 6" ##############
    subcat = cat[(cat['S_Code'] == 'S')
                 & (cat['S_Code_6'] == 'S')
                 & (cat['Peak_flux'] * 1000 > 0.5)
                 & (cat['Peak_flux_6'] * 1000 > 0.5)]
    subcat = subcat[(subcat[f'dDEC_{res}']*3600 < res/2) & (subcat[f'dRA_{res}']*3600 < res/2)]
    # random_index = np.random.choice(len(subcat), size=50)
    # subcat = subcat[random_index]
    print(f"Number of sources for flux ratio 6'': {len(subcat)}")

    R = subcat['Total_flux_6'] / subcat['Total_flux']
    R_err = R * np.sqrt((subcat['E_Total_flux_6'] / subcat['Total_flux_6']) ** 2
                        + (subcat[f'E_Total_flux'] / subcat[f'Total_flux']) ** 2)
    plt.errorbar(subcat['Total_flux_6'], R, xerr=subcat['E_Total_flux_6'], yerr=R_err, color='darkred',
                 ecolor='darkred', markersize=1, fmt='o', capsize=2)
    plt.xscale('log')
    plt.xlabel("Total flux")
    plt.ylabel("Flux ratio")
    plt.plot([subcat['Total_flux_6'].min(), subcat['Total_flux_6'].max()], [1, 1], color='black', linestyle='--')
    plt.ylim(0.5, 1.5)
    plt.savefig(f'{outputfolder}/lotssdeep_ratio_{res}.png', dpi=150)
    print(f'Median error on the ratio Total_flux_6/Total_flux: {round(np.mean(R_err), 2)}')
    print(f'Median ratio Total_flux_6/Total_flux: {round(np.median(R[np.isfinite(R)]), 2)}')
    print(f'Percentage of sources with ratio above 1 {round(len(R[R>1])/len(R), 3)} and below 1 {round(len(R[R<1])/len(R), 3)}')
    plt.close()

    ############# Peak flux over Total flux ##############
    subcat = cat[(cat['S_Code'] == 'S')
                 & (cat['Peak_flux'] > cat['Isl_rms'] * 30)]
    print(f"Number of sources for peak flux over total flux: {len(subcat)}")
    R = subcat['Peak_flux'] / subcat['Total_flux']
    plt.scatter(subcat['Total_flux'], R, s=5)
    plt.xscale('log')
    plt.xlabel("Total flux")
    plt.ylabel("Peakflux/Totalflux")
    plt.savefig(f'{outputfolder}/peak_total_{res}.png', dpi=150)
    print(f'Peak_flux/Total_flux: {round(np.median(R[np.isfinite(R)]), 2)}')
    plt.close()

    # subcat = cat[(cat['S_Code'] == 'S')
    #              & (cat['Peak_flux'] > cat['Isl_rms'] * 30)]
    R = subcat['Peak_flux'] / subcat['Total_flux']
    subcat['dist'] = list(map(dist_pointing_center, subcat['RA', 'DEC']))
    plt.figure(figsize=(5,4))

    plt.scatter(subcat['RA'] % 360, subcat['DEC'] % 360, c=R, s=20, alpha=0.75)
    plt.xlabel("Right Ascension (degrees)")
    plt.ylabel("Declination (degrees)")
    plt.colorbar(label='Peak / integrated flux')
    plt.savefig(f'{outputfolder}/peak_total_{res}_im.png', dpi=150)
    plt.close()

    # subcat = cat[(cat['S_Code'] == 'S')
    #              & (cat['Peak_flux'] > cat['Isl_rms'] * 30)]
    R = subcat['Peak_flux'] / subcat['Total_flux']
    subcat['dist'] = list(map(dist_pointing_center, subcat['RA', 'DEC']))
    degree = 2
    coeffs = np.polyfit(subcat['dist'], R, degree)
    poly_func = np.poly1d(coeffs)
    x_fit = np.linspace(subcat['dist'].min(), subcat['dist'].max(), 100)
    y_fit = poly_func(x_fit)

    centralhz = 140232849.121094
    bandwidth = 12207*res/0.3
    inttime = 1*res/0.3
    ts, bws, totals = fluxration_smearing(centralhz, bandwidth, inttime, res, np.linspace(subcat['dist'].min(), subcat['dist'].max(), 100))

    plt.figure(figsize=(5,4))
    plt.plot(x_fit, totals, color='darkblue', label='Theoretical peak response', linestyle='-.')
    plt.plot(x_fit, y_fit, color='black', label='Peak/integrated flux fit', linestyle='--')
    plt.scatter(subcat['dist'], R, s=6, c=np.clip(subcat['Peak_flux']/subcat['Isl_rms'], a_min=5,  a_max=80), alpha=0.75)
    plt.xlabel("Distance from pointing center (degrees)")
    plt.ylabel("Peak / integrated flux")
    plt.colorbar(label='Peak/RMS')
    plt.ylim(0.4, 1)
    plt.xlim(0, x_fit.max())
    plt.legend()
    plt.savefig(f'{outputfolder}/peak_total_{res}_dist.png', dpi=150)
    plt.close()

def merge_with_table(catalog1, catalog2, sep=6, res=0.3):
    """Merge with other table"""

    coords1 = SkyCoord(ra=catalog1['RA'], dec=catalog1['DEC'], unit=(u.deg, u.deg))
    coords2 = SkyCoord(ra=catalog2['RA'], dec=catalog2['DEC'], unit=(u.deg, u.deg))
    idx_catalog2, separation, _ = match_coordinates_sky(coords1, coords2)

    match_idxs = np.where(separation < sep * u.arcsec)[0]

    catalog2 = catalog2[idx_catalog2][match_idxs]

    for col in ['Total_flux', 'Peak_flux', 'E_Total_flux', 'E_Peak_flux', 'Maj', 'Min', 'PA', 'E_PA', 'S_Code', 'Isl_rms']:
        if col != 'S_Code':
            catalog1[col + "_6"] = np.nan
        else:
            catalog1[col + "_6"] = ''
        catalog1[col + "_6"][match_idxs] = catalog2[col]

    catalog1[f'dRA_{res}'] = np.nan
    catalog1[f'dDEC_{res}'] = np.nan
    catalog1[f'E_dRA_{res}'] = np.nan
    catalog1[f'E_dDEC_{res}'] = np.nan

    catalog1[f'dRA_{res}'][match_idxs] = catalog1[match_idxs] ['RA'] % 360 - catalog2['RA'] % 360
    catalog1[f'dDEC_{res}'][match_idxs] = catalog1[match_idxs] ['DEC'] % 360 - catalog2['DEC'] % 360
    catalog1[f'E_dRA_{res}'][match_idxs] = np.sqrt(catalog1[match_idxs] ['E_RA'] ** 2 + catalog2['E_RA'] ** 2)
    catalog1[f'E_dDEC_{res}'][match_idxs] = np.sqrt(catalog1[match_idxs] ['E_DEC'] ** 2 + catalog2['E_DEC'] ** 2)

    print(f"Number of cross-matches between two catalogs {len(match_idxs)} ({int(len(match_idxs)/len(catalog1)*100)}% matched)")

    return catalog1


def parse_args():
    """Parse input arguments"""

    parser = argparse.ArgumentParser(description='Catalogue matching')
    parser.add_argument('--cat1', nargs='+', help='Catalogue 1 (can be multiple)', default=None)
    parser.add_argument('--cat2', type=str, help='Catalogue 2', default=
                        "/project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/"
                          "imaging/split_facets2/catalogue/"
                          "6asec/pybdsf_sources_6asec.fits")
    parser.add_argument('--separation_asec', type=float, default=6, help=
                        'minimal separation between catalogue 1 and catalogue 2')
    parser.add_argument('--source_id_prefix', type=str)
    parser.add_argument('--out_table', type=str, default='final_merged.fits')
    parser.add_argument('--resolution', type=float, default=0.3)
    parser.add_argument('--only_plot', action='store_true', help='make only plot')

    return parser.parse_args()


def main():
    """Main"""


    outcols = ['Cat_id', 'RA', 'E_RA', 'DEC', 'E_DEC', 'Total_flux', 'E_Total_flux', 'Peak_flux', 'E_Peak_flux',
                     'Maj', 'E_Maj', 'Min', 'E_Min', 'PA', 'E_PA', 'S_Code', 'Isl_rms']

    args = parse_args()

    if args.only_plot:
        totalcat = Table.read(args.out_table, format='fits')

    else:
        for n, cat in enumerate(args.cat1):
            print(cat)
            catalog1_new, _ = separation_match(cat, args.cat2, args.separation_asec)
            remove_snr(catalog1_new, snr=5, cat=cat)
            if args.source_id_prefix is not None:
                catalog1_new['Cat_id'] = [f'{args.source_id_prefix}_{id}' for id in list(catalog1_new['Source_id'])]
            else:
                catalog1_new['Cat_id'] = [f'{cat.split("/")[-1].split("_")[1]}_{id}' for id in list(catalog1_new['Source_id'])]
            if n==0:
                totalcat = catalog1_new
            else:
                totalcat = vstack([totalcat, catalog1_new], join_type='exact')

        totalcat = totalcat[outcols]

        # signal to noise cut
        totalcat = remove_snr(totalcat, snr=5)

        # add 6asec columns
        DR1 = '/home/jurjen/Documents/ELAIS/catalogues/en1_final_cross_match_catalogue-v1.0.fits'
        DR2 = '/home/jurjen/Documents/ELAIS/catalogues/pybdsf_sources_6asec.fits'
        totalcat = merge_with_table(totalcat, Table.read(DR1, format='fits'), sep=args.separation_asec, res=args.resolution)

        totalcat.write(args.out_table, format='fits', overwrite=True)
        totalcat[outcols].write('publication_'+args.out_table, format='fits', overwrite=True)

    print(len(totalcat))
    make_plots(totalcat, res=args.resolution, outputfolder='/home/jurjen/Documents/ELAIS/paperplots/')


if __name__ == '__main__':
    density = {'facet_0': 0.24143183950617234,
               'facet_1': 0.13690175308641947,
               'facet_10': 0.20427424691357982,
               'facet_11': 0.2855608888888883,
               'facet_12': 0.22253392592592544,
               'facet_13': 0.19539255555555515,
               'facet_14': 0.16753998765432063,
               'facet_15': 0.10965179012345656,
               'facet_16': 0.12828291358024665,
               'facet_17': 0.34949097530864126,
               'facet_18': 0.26597167901234514,
               'facet_19': 0.09615233333333313,
               'facet_2': 0.12172407407407382,
               'facet_20': 0.11137129629629607,
               'facet_21': 0.14146774074074045,
               'facet_22': 0.20576508641975266,
               'facet_23': 0.4636706543209867,
               'facet_24': 0.07083577777777762,
               'facet_25': 0.23814687654320937,
               'facet_26': 0.47480614814814714,
               'facet_27': 0.219969469135802,
               'facet_28': 0.22957644444444394,
               'facet_29': 0.23002959259259212,
               'facet_3': 0.23208735802469085,
               'facet_4': 0.19441218518518477,
               'facet_5': 0.22606349382716,
               'facet_6': 0.12530912345678985,
               'facet_7': 0.1945322962962959,
               'facet_8': 0.15238839506172808,
               'facet_9': 0.21153599999999956}
    main()

# python catalogue_helpers/crossmatch.py --cat1 /home/jurjen/Documents/ELAIS/catalogues/finalcat03/*.fits --cat2 /home/jurjen/Documents/ELAIS/catalogues/pybdsf_sources_6asec.fits --out_table final_merged_03.fits
# python catalogue_helpers/crossmatch.py --cat1 /home/jurjen/Documents/ELAIS/catalogues/finalcat06/*.fits --cat2 /home/jurjen/Documents/ELAIS/catalogues/pybdsf_sources_6asec.fits --out_table final_merged_06.fits --resolution 0.6
# python catalogue_helpers/crossmatch.py --cat1 /home/jurjen/Documents/ELAIS/catalogues/finalcat12/*.fits --cat2 /home/jurjen/Documents/ELAIS/catalogues/pybdsf_sources_6asec.fits --out_table final_merged_12.fits --resolution 1.2