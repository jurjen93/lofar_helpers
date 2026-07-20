"""Make a cutout with a ds9 region file"""

import numpy as np
from argparse import ArgumentParser
from glob import glob
from astropy.io import fits
import pyregion
import numpy as np



def parse_args():
    """
    Command line argument parser
    :return: parsed arguments
    """
    parser = ArgumentParser()
    parser.add_argument('dir', help='dir', type=str)
    return parser.parse_args()


def main():
    """ Main function"""
    args = parse_args()

    ifits, qfits, ufits, vfits = sorted(glob(args.dir+"/0.3arcsec-MFS-?-image.fits"))

    regionfile = "/net/rijn9/data2/jurjendejong/ELAIS/ILTJ160607.62+552135.4_DI_polarisation/ds9.reg"

    with fits.open(ifits) as hdu:
        ifitsdata = hdu[0].data
        header = hdu[0].header
    with fits.open(qfits) as hdu:
        qfitsdata = hdu[0].data
    with fits.open(ufits) as hdu:
        ufitsdata = hdu[0].data
    with fits.open(vfits) as hdu:
        vfitsdata = hdu[0].data


    r = pyregion.open(regionfile).as_imagecoord(header=header)
    mask = r.get_mask(hdu=hdu[0], shape=(header["NAXIS2"], header["NAXIS1"])).astype(int)
    ifitsdata *= mask
    qfitsdata *= mask
    ufitsdata *= mask
    vfitsdata *= mask

    print(f"Leakage linear: {np.sqrt(np.sum(qfitsdata)**2+np.sum(ufitsdata)**2)/np.sum(ifitsdata)*100}%")
    print(f"Leakage V: {np.sum(vfitsdata)/np.sum(ifitsdata)*100}%")



if __name__ == '__main__':
    main()
