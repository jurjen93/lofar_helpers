from argparse import ArgumentParser
from os.path import basename
from astropy.io import fits
import numpy as np
from scipy.ndimage import median_filter, binary_dilation, label

def remove_stripes(image, median_size=5, min_size=30, dilation=3):
    """
    Detect elongated stripe artefacts and replace those pixels with their
    local median, while protecting compact real sources.

    Parameters
    ----------
    image : 2D ndarray
        Input image data.
    median_size : int, optional
        Size (in px) of the local median filter used both to estimate the
        smooth background (for residual computation). Default 5.
    min_size : int, optional
        Minimum connected-component size (in px) for a candidate region to
        be kept; smaller specks are treated as isolated noise and dropped
        before shape analysis. Default 30.
    dilation : int, optional
        Number of iterations to grow the final stripe mask by, to catch
        faint edge pixels just below the detection threshold. Also acts as
        a switch: only when dilation > 0 is the stripe mask additionally
        restricted to source_mask (pixels below 5x the image RMS), which
        excludes bright compact sources from being overwritten even if
        they were misclassified as "thin". Default 3.

    Returns
    -------
    cleaned : ndarray
        Copy of `image` with detected stripe pixels replaced by their
        local median value.
    """

    # Make candidate stripe mask
    med = median_filter(image, size=median_size)
    residual = image - med
    mad = np.median(np.abs(residual))
    sigma_est = 1.4826 * mad
    candidate_mask = np.abs(residual) > sigma_est

    # Drop tiny specks
    labels, n = label(candidate_mask)
    sizes = np.bincount(labels.ravel())
    candidate_mask = candidate_mask & (sizes[labels] >= min_size)

    if dilation > 0:
        stripe_mask = binary_dilation(candidate_mask, iterations=dilation)

    # Add source mask
    rms = get_rms(image)
    stripe_mask = stripe_mask & (image < rms*5)

    cleaned = image.copy()
    cleaned[stripe_mask] = med[stripe_mask]
    return cleaned


def get_rms(image_data):
    """
    from Cyril Tasse/kMS

    :param image_data: image data array
    :return: rms (noise measure)
    """

    maskSup = 1e-7
    m = image_data[np.abs(image_data) > maskSup]
    rmsold = np.std(m)
    diff = 1e-1
    cut = 3.
    med = np.median(m)
    for _ in range(10):
        ind = np.where(np.abs(m - med) < rmsold * cut)[0]
        rms = np.std(m[ind])
        if np.abs(np.divide((rms - rmsold), rmsold)) < diff: break
        rmsold = rms
    return rms


def parse_args():
    """
    Command line argument parser
    :return: parsed arguments
    """
    parser = ArgumentParser(description='Remove stripes through image from bad subtraction')
    parser.add_argument('fits_input', help='fits input file', type=str)
    parser.add_argument('--output_name', help='fits output file', type=str)
    return parser.parse_args()


def main():
    """ Main function"""
    args = parse_args()

    hdu = fits.open(args.fits_input)
    header = hdu[0].header
    data = hdu[0].data

    destriped_data = remove_stripes(data)

    if args.output_name is None:
        output_name = basename(args.fits_input)
    else:
        output_name = args.output_name

    hdu = fits.PrimaryHDU(header=header, data=destriped_data)
    hdu.writeto(output_name, overwrite=True)


if __name__ == '__main__':
    main()
