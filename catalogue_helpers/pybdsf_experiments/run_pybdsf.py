"""
This script runs pybdsf on a fits file image to extract sources and components.

This has been used for catalogue reduction of the ELAIS-N1 field.
Feel free to adapt to your own needs. (Watch out for hardcoded parameters or paths)

"""

import bdsf
import argparse

def run_pybdsf(fitsfile, rmsbox, frequency):
    """
    Run pybdsf

    :param fitsfile: fits file
    :param rmsbox: rms box first parameter

    :return: source catalogue
    """

    prefix = fitsfile.replace('.fits', '')
    img = bdsf.process_image(fitsfile,
                             thresh_isl=3,
                             thresh_pix=5,
                             atrous_do=True,
                             rms_box=(int(rmsbox), int(rmsbox // 8)),
                             rms_box_bright=(int(rmsbox//3), int(rmsbox//12)),
                             adaptive_rms_box=True,
                             group_tol=10.0,
                             frequency=frequency)  # , rms_map=True, rms_box = (160,40))

    img.write_catalog(clobber=True, outfile=prefix + '_source_catalog.fits', format='fits', catalog_type='srl')
    img.write_catalog(clobber=True, outfile=prefix + '_gaussian_catalog.fits', format='fits', catalog_type='gaul')
    for type in ['island_mask', 'gaus_model', 'gaus_resid', 'mean', 'rms']:
        img.export_image(clobber=True, img_type=type, outfile=prefix + f'_{type}.fits')
    return prefix + '_source_catalog.fits'


def parse_args():
    """
    Parse input arguments
    """

    parser = argparse.ArgumentParser(description='Source detection')
    parser.add_argument('--rmsbox', type=int, help='rms box pybdsf', default=120)
    parser.add_argument('--frequency', help='frequency')
    parser.add_argument('fits', help='fits files')
    return parser.parse_args()


def main():
    """
    Main function
    """

    args = parse_args()
    run_pybdsf(args.fits, args.rmsbox, args.frequency)


if __name__ == '__main__':
    main()
