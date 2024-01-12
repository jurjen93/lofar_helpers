"""
This script runs pybdsf on a fits file image to extract sources and components.
From the source table, it makes cut out images of all resolved sources (fits and png images).

TODO: Make contours and see which sources are contained in same contour
"""

import bdsf
import argparse
from astropy.nddata import Cutout2D
import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.wcs.utils import skycoord_to_pixel
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm, PowerNorm
from astropy.visualization.wcsaxes import WCSAxes
import matplotlib.path as mpath
import matplotlib.patches as patches
from matplotlib import ticker
import os
from astropy.table import Table
import pyregion
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial import distance


def get_rms(image_data):
    """
    from Cyril Tasse/kMS

    :param image_data: image data array
    :return: rms (noise measure)
    """
    from past.utils import old_div

    maskSup = 1e-7
    m = image_data[np.abs(image_data) > maskSup]
    rmsold = np.std(m)
    diff = 1e-1
    cut = 3.
    med = np.median(m)
    for _ in range(10):
        ind = np.where(np.abs(m - med) < rmsold * cut)[0]
        rms = np.std(m[ind])
        if np.abs(old_div((rms - rmsold), rmsold)) < diff: break
        rmsold = rms
    print(f'Noise : {str(round(rms * 1000, 4))} {u.mJy / u.beam}')
    return rms


def get_beamarea(hdu):
    """
    Get beam area in pixels
    """

    bmaj = hdu[0].header['BMAJ']
    bmin = hdu[0].header['BMIN']

    beammaj = bmaj / (2.0 * (2 * np.log(2)) ** 0.5)  # Convert to sigma
    beammin = bmin / (2.0 * (2 * np.log(2)) ** 0.5)  # Convert to sigma
    pixarea = abs(hdu[0].header['CDELT1'] * hdu[0].header['CDELT2'])

    beamarea = 2 * np.pi * 1.0 * beammaj * beammin  # Note that the volume of a two dimensional gaus$
    beamarea_pix = beamarea / pixarea

    return beamarea_pix


def make_cutout(fitsfile=None, pos: tuple = None, size: tuple = (1000, 1000), savefits=None):
    """
    Make cutout from your image with this method.
    pos (tuple) -> position in pixels
    size (tuple) -> size of your image in pixel size, default=(1000,1000)
    """
    fts = fits.open(fitsfile)
    image_data = fts[0].data
    wcs = WCS(fts[0].header, naxis=2)

    out = Cutout2D(
        data=image_data,
        position=pos,
        size=size,
        wcs=wcs,
        mode='partial'
    )

    wcs = out.wcs
    header = wcs.to_header()

    image_data = out.data
    # rms = get_rms(image_data)
    # hdu = [fits.PrimaryHDU(image_data, header=header)]

    if savefits:
        image_data = out.data
        image_data = np.expand_dims(np.expand_dims(image_data, axis=0), axis=0)
        fits.writeto(savefits, image_data, header, overwrite=True)


def make_image(fitsfile=None, cmap: str = 'RdBu_r', components: str = None):
    """
    Image your data with this method.
    fits file -> fits file
    cmap -> choose your preferred cmap
    """

    hdu = fits.open(fitsfile)
    image_data = hdu[0].data
    while image_data.ndim > 2:
        image_data = image_data[0]
    header = hdu[0].header

    rms = get_rms(image_data)
    vmin = rms/10
    vmax = rms * 10

    if hdu is None:
        wcs = WCS(header, naxis=2)
    else:
        wcs = WCS(hdu[0].header, naxis=2)

    fig = plt.figure(figsize=(7, 10), dpi=200)
    plt.subplot(projection=wcs)
    WCSAxes(fig, [0.1, 0.1, 0.8, 0.8], wcs=wcs)
    im = plt.imshow(image_data, origin='lower', cmap=cmap)
    im.set_norm(PowerNorm(gamma=0.5, vmin=vmin, vmax=vmax))
    plt.xlabel('Right Ascension (J2000)', size=14)
    plt.ylabel('Declination (J2000)', size=14)
    plt.tick_params(axis='both', which='major', labelsize=12)

    def fixed_color(shape, saved_attrs):
        from pyregion.mpl_helper import properties_func_default
        attr_list, attr_dict = saved_attrs
        attr_dict["color"] = 'pink'
        kwargs = properties_func_default(shape, (attr_list, attr_dict))

        return kwargs

    if components is not None:
        r = pyregion.open(components).as_imagecoord(header=hdu[0].header)
        patch_list, artist_list = r.get_mpl_patches_texts(fixed_color)

        # fig.add_axes(ax)
        for patch in patch_list:
            plt.gcf().gca().add_patch(patch)
        for artist in artist_list:
            plt.gca().add_artist(artist)


    # p0 = axes[0].get_position().get_points().flatten()
    # p2 = axes[2].get_position().get_points().flatten()

    orientation = 'horizontal'
    ax_cbar1 = fig.add_axes([0.22, 0.15, 0.73, 0.02])
    cb = plt.colorbar(im, cax=ax_cbar1, orientation=orientation)
    cb.set_label('Surface brightness [mJy/beam]', size=16)
    cb.ax.tick_params(labelsize=16)

    cb.outline.set_visible(False)

    # Extend colorbar
    bot = -0.05
    top = 1.05

    # Upper bound
    xy = np.array([[0, 1], [0, top], [1, top], [1, 1]])
    if orientation == "horizontal":
        xy = xy[:, ::-1]

    Path = mpath.Path

    # Make Bezier curve
    curve = [
        Path.MOVETO,
        Path.CURVE4,
        Path.CURVE4,
        Path.CURVE4,
    ]

    color = cb.cmap(cb.norm(cb._values[-1]))
    patch = patches.PathPatch(
        mpath.Path(xy, curve),
        facecolor=color,
        linewidth=0,
        antialiased=False,
        transform=cb.ax.transAxes,
        clip_on=False,
    )
    cb.ax.add_patch(patch)

    # Lower bound
    xy = np.array([[0, 0], [0, bot], [1, bot], [1, 0]])
    if orientation == "horizontal":
        xy = xy[:, ::-1]

    color = cb.cmap(cb.norm(cb._values[0]))
    patch = patches.PathPatch(
        mpath.Path(xy, curve),
        facecolor=color,
        linewidth=0,
        antialiased=False,
        transform=cb.ax.transAxes,
        clip_on=False,
    )
    cb.ax.add_patch(patch)

    # Outline
    xy = np.array(
        [[0, 0], [0, bot], [1, bot], [1, 0], [1, 1], [1, top], [0, top], [0, 1], [0, 0]]
    )
    if orientation == "horizontal":
        xy = xy[:, ::-1]

    Path = mpath.Path

    curve = [
        Path.MOVETO,
        Path.CURVE4,
        Path.CURVE4,
        Path.CURVE4,
        Path.LINETO,
        Path.CURVE4,
        Path.CURVE4,
        Path.CURVE4,
        Path.LINETO,
    ]
    path = mpath.Path(xy, curve, closed=True)

    patch = patches.PathPatch(
        path, facecolor="None", lw=1, transform=cb.ax.transAxes, clip_on=False
    )
    cb.ax.add_patch(patch)
    tick_locator = ticker.MaxNLocator(nbins=3)
    cb.locator = tick_locator
    cb.update_ticks()

    fig.tight_layout(pad=1.0)
    plt.grid(False)
    plt.grid('off')
    plt.savefig(fitsfile.replace('.fits', '.png'), dpi=250, bbox_inches='tight')
    plt.close()

    hdu.close()


def run_pybdsf(fitsfile, rmsbox):
    """
    Run pybdsf

    :param fitsfile: fits file
    :param rmsbox: rms box first parameter

    :return: source catalogue
    """

    prefix = fitsfile.replace('.fits', '')
    img = bdsf.process_image(fitsfile,
                             thresh_isl=3.5,
                             thresh_pix=5.5,
                             atrous_do=True,
                             rms_box=(int(rmsbox), int(rmsbox // 4)),
                             rms_box_bright=(int(rmsbox//2), int(rmsbox//8)),
                             adaptive_rms_box=True,
                             group_tol=10.0)  # , rms_map=True, rms_box = (160,40))

    img.write_catalog(clobber=True, outfile=prefix + '_source_catalog.fits', format='fits', catalog_type='srl')
    img.write_catalog(clobber=True, outfile=prefix + '_gaussian_catalog.fits', format='fits', catalog_type='gaul')
    for type in ['island_mask', 'gaus_model', 'gaus_resid', 'mean', 'rms']:
        img.export_image(clobber=True, img_type=type, outfile=prefix + f'_{type}.fits')
    return prefix + '_source_catalog.fits'


def get_pix_coord(table):
    """
    Get pixel coordiantes from table RA/DEC
    :param table: fits table
    :return: pixel coordinates
    """

    f = fits.open(table)
    t = f[1].data
    # res = t[t['S_Code'] != 'S']
    pos = list(zip(t['RA'], t['DEC'], t['Source_id']))
    pos = [[SkyCoord(f'{c[0]}deg', f'{c[1]}deg', frame='icrs'), c[2]] for c in pos]
    fts = fits.open(table.replace("_source_catalog", ""))
    w = WCS(fts[0].header, naxis=2)
    pix_coord = [([int(c) for c in skycoord_to_pixel(sky[0], w, 0, 'all')], sky[1]) for sky in pos]
    return pix_coord


def parse_args():
    """
    Parse input arguments
    """

    parser = argparse.ArgumentParser(description='Source detection')
    parser.add_argument('--rmsbox', type=int, help='rms box pybdsf', default=160)
    parser.add_argument('--no_pybdsf', action='store_true', help='Skip pybdsf')
    parser.add_argument('fitsf', nargs='+', help='fits files')
    # parser.add_argument('--ref_catalogue', help='fits table')
    return parser.parse_args()


def make_point_file(t):
    """
    Make ds9 file with ID in it
    """

    header = """# Region file format: DS9 version 4.1
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5
"""
    t = Table.read(t, format='fits')
    file = open('components.reg', 'w')
    file.write(header)
    for n, c in enumerate(zip(t['RA'], t['DEC'])):
        file.write(f'\n# text({c[0]},{c[1]}) text=' + '{' + f'{t["Source_id"][n]}' + '}')
    return


def get_table_index(t, source_id):
    return int(np.argwhere(t['Source_id'] == source_id).squeeze())


def get_clusters(t):
    """
    Get clusters of sources based on euclidean distance
    """
    deg_dist = 0.01
    ra_dec = np.stack((list(t['RA']),list(t['DEC'])),axis=1)
    Z = linkage(ra_dec, method='complete', metric='euclidean')
    return fcluster(Z, deg_dist, criterion='distance')


def cluster_idx(clusters, idx):
    return np.argwhere(clusters==clusters[idx]).squeeze(axis=1)


def max_dist(coordinates):
    """
    Get the longest distance between coordinates

    :param coordinates: indices from table
    """
    return np.max(distance.cdist(coordinates, coordinates, 'euclidean'))



def main():
    """
    Main function
    """

    args = parse_args()
    for m, fts in enumerate(args.fitsf):
        if not args.no_pybdsf:
            tbl = run_pybdsf(fts, args.rmsbox)
        else:
            tbl = fts.replace('.fits', '') + '_source_catalog_final.fits'
        # make ds9 region file with sources in it
        make_point_file(tbl)
        # loop through resolved sources and make images
        coord = get_pix_coord(tbl)

        T = Table.read(tbl, format='fits')
        T['Peak_flux_min'] = T['Peak_flux'] - T['E_Peak_flux']
        T['Total_flux_min'] = T['Total_flux'] - T['E_Total_flux']

        clusters = get_clusters(T)

        f = fits.open(fts)
        beamarea = get_beamarea(f)
        imdat = f[0].data
        rms = get_rms(imdat)
        f.close()

        os.system('mkdir -p bright_sources')
        os.system('mkdir -p weak_sources')
        os.system('mkdir -p cluster_sources')

        to_delete = []
        to_ignore = []

        for c, n in coord:
            table_idx = get_table_index(T, n)
            if table_idx in to_ignore:
                continue

            cluster_indices = cluster_idx(clusters, table_idx)
            if T[T['Source_id'] == n]['Peak_flux'][0] < rms/3:
                to_delete.append(table_idx)
                print("Delete Source_id: "+str(n))

            elif len(cluster_indices) > 1:
                pix_coord = np.array([p[0] for p in coord])[cluster_indices]
                imsize = int(max_dist(pix_coord)*2*1.15)
                idxs = '-'.join([str(p) for p in cluster_indices])
                make_cutout(fitsfile=fts, pos=tuple(c), size=(imsize, imsize), savefits=f'cluster_sources/source_{m}_{idxs}.fits')
                make_image(f'cluster_sources/source_{m}_{idxs}.fits', 'RdBu_r', 'components.reg')
                for i in cluster_indices:
                    to_ignore.append(i)

            elif (T[T['Source_id'] == n]['Peak_flux_min'][0] < rms or
                  T[T['Source_id'] == n]['Peak_flux'][0] < 2*rms or
                  T[T['Source_id'] == n]['Total_flux_min'][0]/beamarea < rms/10):
                make_cutout(fitsfile=fts, pos=tuple(c), size=(300, 300), savefits=f'weak_sources/source_{m}_{n}.fits')
                make_image(f'weak_sources/source_{m}_{n}.fits', 'RdBu_r', 'components.reg')

            else:
                make_cutout(fitsfile=fts, pos=tuple(c), size=(300, 300), savefits=f'bright_sources/source_{m}_{n}.fits')
                make_image(f'bright_sources/source_{m}_{n}.fits', 'RdBu_r', 'components.reg')

        for i in sorted(to_delete)[::-1]:
            del T[i]

if __name__ == '__main__':
    main()