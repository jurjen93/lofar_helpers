#!/usr/bin/env python
import argparse
import os

from astropy.io import fits # type: ignore
import numpy as np

def blank_image(imname: str, beamlvl: float, MFS: bool) -> str:
    """ Blanks images created by WSClean beyond a specified level of the primary beam.

    Args:
        imname : str
            Name of the images to blank as passed to WSClean's -name option.
        beamlvl : float
            Level of the primary beam to blank below.
        MFS : bool
            Set to True if the image was made in MFS mode (i.e. channels-out > 1).

    Returns:
        imname_blanked : str
            Base name of the blanked images.
    """
    if MFS:
        im_app = imname + '-MFS-image.fits'
        im_pb = imname + '-MFS-image-pb.fits'
        im_beam = imname + '-MFS-beam-0.fits'
    else:
        im_app = imname + '-image.fits'
        im_pb = imname + '-image-pb.fits'
        im_beam = imname + '-beam-0.fits'

    if not os.path.isfile(im_app):
        raise FileNotFoundError(f"Image {im_app} does not exist.")
    if not os.path.isfile(im_pb):
        raise FileNotFoundError(f"Image {im_pb} does not exist.")
    if not os.path.isfile(im_beam):
        raise FileNotFoundError(f"Image {im_beam} does not exist.")
    beam = fits.getdata(im_beam)
    mask = beam > beamlvl

    data_app = fits.getdata(im_app)
    data_pb = fits.getdata(im_pb)

    im_app_masked = np.where(mask, data_app, np.nan)
    im_pb_masked = np.where(mask, data_pb, np.nan)

    print(f'Writing primary-beam-blanked images for {imname}...')
    fits.writeto(f"{imname}.pbblanked-MFS-image.fits", data=im_app_masked, header=fits.getheader(im_app))
    fits.writeto(f"{imname}.pbblanked-MFS-image-pb.fits", data=im_pb_masked, header=fits.getheader(im_pb))
    print('Done')
    return imname + ".pbblanked"


if __name__ == '__main__':
    parser = argparse.ArgumentParser("Blank a WSClean image to a certain value of the primary beam.")
    parser.add_argument("--imagename", type=str, help="Image name as passed to WSClean.")
    parser.add_argument("--beam_cut", type=float, help="Fractional value of the primary beam to blank to relative to the centre of the image.")
    parser.add_argument("--no-MFS", action="store_true", required=False, help="Set if image was made using a single output channel (i.e. -channels-out=1 or not given).")
    args = parser.parse_args()

    blank_image(args.imagename, args.beam_cut, not args.no_MFS)
