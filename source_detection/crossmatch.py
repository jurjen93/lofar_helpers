"""Code to crossmatch catalogues"""

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

    # idx_catalog2 contains the indices of matched sources in catalog2
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

    return catalog1_corrected, catalog1_removed


def remove_snr(catalog, snr=5):
    """Remove sources based on Signal-To-Noise"""

    print(f"Source count before SNR cut {len(catalog)}")
    newcat = catalog[catalog['Peak_flux'] > snr * catalog['Isl_rms']]
    print(f'Source count after SNR cut {len(newcat)} ({int(abs(1-len(catalog)/len(newcat))*100)}% removed)')

    return newcat


def crossmatch_itself(catalog, min_sep=0.15):
    """Crossmatch with itself to find nearest neighbour"""

    coords = SkyCoord(ra=catalog['RA'], dec=catalog['DEC'], unit=(u.deg, u.deg))
    idx_catalog, separation, _ = match_coordinates_sky(coords, coords, nthneighbor=2)
    nearest_neighbour = separation<min_sep*u.arcsec
    print(catalog[nearest_neighbour]['Cat_id'])


def make_plots(cat, res=0.3):
    """Make plots"""

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
    subcat = cat[(cat[f'dRA_{res}'] == cat[f'dRA_{res}']) & (cat['S_Code'] == 'S') & (cat['Total_flux_6']*1000 > 0.5)]
    print(f"Number of sources for hist2d plot: {len(subcat)}")

    plt.hist2d(subcat[f'dRA_{res}'] * 3600,
               subcat[f'dDEC_{res}'] * 3600,
               bins=50, cmap=plt.cm.jet, norm=LogNorm())
    plt.xlabel("dRA (arcsec)")
    plt.ylabel("dDEC (arcsec)")
    plt.colorbar(label='Source count')

    mediandRA = round(np.median(subcat[f'dRA_{res}']) * 3600, 4)
    mediandDEC = round(np.median(subcat[f'dDEC_{res}']) * 3600, 4)
    plt.title(f"Mean dRA: {mediandRA} "
              f"Mean dDEC: {mediandDEC}")

    axs = np.arange(-16, 17)
    plt.plot(axs, [mediandDEC]*len(axs), color='black', linestyle='--')
    plt.plot([mediandRA]*len(axs), axs, color='black', linestyle='--')
    plt.xlim(-6, 6)
    plt.ylim(-6, 6)
    plt.savefig('dRA_dDEC.png', dpi=150)
    plt.close()


    ############# Flux ratio 6" ##############
    subcat = cat[(cat['S_Code'] == 'S') & (cat['Total_flux'] * 1000 > 0.1)]
    subcat = subcat[(subcat['dDEC_0.3']*3600 < 0.15) & (subcat['dRA_0.3']*3600 < 0.15)]
    R = subcat['Total_flux_6'] / subcat['Total_flux']
    plt.scatter(subcat['Total_flux_6'], R, color='darkred')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel("Total flux")
    plt.ylabel("Flux ratio")
    plt.title(f'Median ratio: {round(np.median(R[np.isfinite(R)]), 2)}')
    plt.savefig('lotssdeep_ratio.png', dpi=150)

    ############# Peak flux over Total flux ##############
    subcat = cat[(cat['S_Code'] == 'S') & (cat['Total_flux'] * 1000 > 1)]
    R = subcat['Peak_flux'] / subcat['Total_flux']
    plt.scatter(subcat['Total_flux'], R)
    plt.xscale('log')
    plt.xlabel("Total flux")
    plt.ylabel("Peakflux/Totalflux")
    plt.title(f'Median ratio: {round(np.median(R[np.isfinite(R)]), 2)}')
    plt.savefig('peak_total.png', dpi=150)


def merge_with_table(catalog1, catalog2, sep=6, res=0.3):
    """Merge with other table"""

    coords1 = SkyCoord(ra=catalog1['RA'], dec=catalog1['DEC'], unit=(u.deg, u.deg))
    coords2 = SkyCoord(ra=catalog2['RA'], dec=catalog2['DEC'], unit=(u.deg, u.deg))
    idx_catalog2, separation, _ = match_coordinates_sky(coords1, coords2)

    match_idxs = np.where(separation < sep * u.arcsec)[0]
    nonmatch_idxs = np.where(separation >= sep * u.arcsec)[0]

    for col in ['Total_flux', 'Peak_flux', 'E_Total_flux', 'E_Peak_flux', 'Maj', 'Min', 'PA', 'E_PA', 'S_Code']:
        if col != 'S_Code':
            catalog1[col + "_6"] = np.nan
        else:
            catalog1[col + "_6"] = ''
        for idx in match_idxs:
            catalog1[idx][col + "_6"] = catalog2[idx][col]

    catalog1[f'dRA_{res}'] = np.nan
    catalog1[f'dDEC_{res}'] = np.nan
    catalog1[f'E_dRA_{res}'] = np.nan
    catalog1[f'E_dDEC_{res}'] = np.nan

    catalog1[f'dRA_{res}'] = catalog1['RA'] % 360 - catalog2[idx_catalog2]['RA'] % 360
    catalog1[f'dDEC_{res}'] = catalog1['DEC'] % 360 - catalog2[idx_catalog2]['DEC'] % 360
    catalog1[f'E_dRA_{res}'] = np.sqrt(catalog1['E_RA'] ** 2 + catalog2[idx_catalog2]['E_RA'] ** 2)
    catalog1[f'E_dDEC_{res}'] = np.sqrt(catalog1['E_DEC'] ** 2 + catalog2[idx_catalog2]['E_DEC'] ** 2)

    for idx in nonmatch_idxs:
        catalog1[idx][f'dRA_{res}'] = np.nan
        catalog1[idx][f'dDEC_{res}'] = np.nan
        catalog1[idx][f'E_dRA_{res}'] = np.nan
        catalog1[idx][f'E_dDEC_{res}'] = np.nan

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

    outcols = ['Cat_id', 'Isl_id', 'RA', 'E_RA', 'DEC','E_DEC', 'Total_flux', 'E_Total_flux', 'Peak_flux', 'E_Peak_flux',
               'Maj', 'E_Maj', 'Min', 'E_Min', 'PA', 'E_PA', 'S_Code', 'Isl_rms']

    args = parse_args()

    if args.only_plot:
        totalcat = Table.read(args.out_table, format='fits')

    else:
        for n, cat in enumerate(args.cat1):
            print(cat)
            catalog1_new, _ = separation_match(cat, args.cat2, args.separation_asec)
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

    print(len(totalcat))
    make_plots(totalcat, res=args.resolution)


if __name__ == '__main__':
    main()

# python source_detection/crossmatch.py --cat1 /home/jurjen/Documents/ELAIS/catalogues/finalcat03/*.fits --cat2 /home/jurjen/Documents/ELAIS/catalogues/pybdsf_sources_6asec.fits --out_table final_merged_03.fits
# python source_detection/crossmatch.py --cat1 /home/jurjen/Documents/ELAIS/catalogues/finalcat06/*.fits --cat2 /home/jurjen/Documents/ELAIS/catalogues/pybdsf_sources_6asec.fits --out_table final_merged_06.fits
