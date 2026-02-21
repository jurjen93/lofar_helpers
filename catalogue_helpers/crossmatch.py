"""
Code to crossmatch catalogues
"""

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
from astropy.coordinates import match_coordinates_sky

__author__ = "Jurjen de Jong"


def find_matches(cat1, cat2, separation_asec):
    """
    Find crossmatches with two catalogues
    """

    catalog1 = Table.read(cat1, format='fits')
    catalog2 = Table.read(cat2, format='fits')

    coords1 = SkyCoord(ra=catalog1['RA'], dec=catalog1['DEC'], unit=(u.deg, u.deg))
    coords2 = SkyCoord(ra=catalog2['RA'], dec=catalog2['DEC'], unit=(u.deg, u.deg))

    idx_catalog2, separation, _ = match_coordinates_sky(coords1, coords2)

    max_sep_threshold = separation_asec * u.arcsec
    matched_sources_mask = separation < max_sep_threshold

    matched_sources_catalog1 = catalog1[matched_sources_mask]
    matched_sources_catalog2 = catalog2[idx_catalog2[matched_sources_mask]]

    return matched_sources_catalog1, matched_sources_catalog2


def crossmatch_itself(catalog, min_sep=0.15):
    """Crossmatch with itself to find nearest neighbour"""

    coords = SkyCoord(ra=catalog['RA'], dec=catalog['DEC'], unit=(u.deg, u.deg))
    idx_catalog, separation, _ = match_coordinates_sky(coords, coords, nthneighbor=2)
    nearest_neighbour = separation<min_sep*u.arcsec

    return catalog[nearest_neighbour]