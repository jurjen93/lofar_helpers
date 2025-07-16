"""
This script runs pybdsf on a fits file image to extract sources and components.
"""

import bdsf
import argparse
from glob import glob

__author__ = "Jurjen de Jong"


def run_pybdsf(fitsfile, rmsbox, detection_image):
    """
    Run pybdsf

    :param fitsfile: fits file
    :param rmsbox: rms box first parameter

    :return: source catalogue
    """

    prefix = fitsfile.replace('.fits', '')
    img = bdsf.process_image(fitsfile,
                             thresh_isl=3.0,
                             thresh_pix=5.0,
                             atrous_do=True,
                             rms_box=(int(rmsbox), int(rmsbox // 8)),
                             rms_box_bright=(int(rmsbox//3), int(rmsbox//12)),
                             adaptive_rms_box=True,
                             group_tol=10.0,
                             advanced_opts=True,
                             detection_image=detection_image,
                             group_by_isl=True)  # , rms_map=True, rms_box = (160,40))

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
    parser.add_argument('--rmsbox', type=int, help='RMS box pybdsf', default=120)
    parser.add_argument('--detection_image', type=int, help='Use alternative (for instance non-pb) image for detection', default=None)
    parser.add_argument('fits', help='FITS image')
    return parser.parse_args()


def main():
    """
    Main function
    """

    args = parse_args()
    if args.detection_image is not None:
        detection_image = args.detection_images
    else:
        try:
            detection_image = glob(args.fits.replace("-pb",""))[0]
        except:
            detection_image = None

    print(f"Science image: {args.fits}")
    print(f"Detection image: {detection_image}")

    run_pybdsf(args.fits, args.rmsbox, detection_image)


if __name__ == '__main__':
    main()
