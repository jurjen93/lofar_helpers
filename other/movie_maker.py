import os
from glob import glob
import imghdr
import sys
from argparse import ArgumentParser
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm, PowerNorm, LogNorm


def findrms(mIn,maskSup=1e-7):
    """
    find the rms of an array, from Cycil Tasse/kMS
    """
    m=mIn[np.abs(mIn)>maskSup]
    rmsold=np.std(m)
    diff=1e-1
    cut=3.
    bins=np.arange(np.min(m),np.max(m),(np.max(m)-np.min(m))/30.)
    med=np.median(m)
    for i in range(10):
        ind=np.where(np.abs(m-med)<rmsold*cut)[0]
        rms=np.std(m[ind])
        if np.abs((rms-rmsold)/rmsold)<diff: break
        rmsold=rms
    return rms


def from_png(frames, output='movie.mp4'):
    """
    Make video or gif from png images

    :param frames: list of frame names
    :param output: output name + type (possibility .mp4 or .gif)
    """
    print('\n'.join(sorted(frames)))
    os.system('mkdir -p framestmp')
    for n, frame in enumerate(sorted(frames)):
        if imghdr.what(frame)!='png' and imghdr.what(frame)!='jpg':
            sys.exit('ERROR: '+frame+' is not jpg or png file')
        os.system(f'cp {frame} framestmp/frame_{str(n).rjust(5, "0")}.png')
    os.system(f'ffmpeg -y -f image2 -r 2 -lavfi "fps=10,scale=720:-1:flags=lanczos" -i framestmp/frame_%05d.png '
              f'{output} && rm -rf framestmp')
    return

def png_from_fits(fitsfiles, cmap, norm):
    """
    Make png image from fits image

    :param fitsfiles: fits file names
    :return: list of png frames
    """
    os.system('if [ -d framesfits ]; then rm -rf framesfits; fi && mkdir -p framesfits')
    for m, fitsfile in enumerate(fitsfiles):
        print(f'Make framesfits/frame_{str(m).rjust(5, "0")}.png')
        hdu = fits.open(fitsfile)
        imdata = hdu[0].data
        while imdata.ndim>2:
            imdata = imdata[0]
        # hdr = hdu[0].header
        fig, axes = plt.subplots(figsize=(12, 10), nrows=1, ncols=1)
        rms = findrms(imdata)
        if m==0:
            vmin, vmax = rms, 30 * rms
        if norm=='SquareRoot':
            nrm = PowerNorm(gamma=1/2, vmin=vmin, vmax=vmax)
        elif norm=='SymLogNorm':
            nrm = SymLogNorm(linthresh=abs(vmin), vmin=0, vmax=vmax)
        elif norm=='LogNorm':
            nrm = SymLogNorm(linthresh=abs(vmin), vmin=abs(vmin), vmax=vmax)
        else:
            sys.exit("ERROR: unexpected norm (expecting SquareRoot, SymLogNorm, or LogNorm")
        axes.imshow(imdata, origin='lower', cmap=cmap,
                            norm=nrm)
        axes.set_yticks([])
        axes.set_xticks([])
        fig.tight_layout(pad=1.0)
        plt.savefig(f'framesfits/frame_{str(m).rjust(5, "0")}.png', dpi=300)

    return glob('framesfits/frame_*.png')


def parse_arg():
    """
    Command line argument parser
    """

    parser = ArgumentParser(description='make videos or gifs from png, jpg or fits files')
    parser.add_argument('--input_frames', help='png, jpg, or fits', required=True, nargs='+')
    parser.add_argument('--output_name', type=str, help='output name (f.e.: movie.mp4)', default='movie.mp4')
    parser.add_argument('--cmap', type=str, help='Color map', default='RdBu_r')
    parser.add_argument('--norm', type=str, help='Color scale norm (SquareRoot or SymLogNorm or LogNorm)', default='SquareRoot')

    return parser.parse_args()

def main():
    """
    Main function
    """

    args = parse_arg()

    if imghdr.what(args.input_frames[0])=='png' or imghdr.what(args.input_frames[0])=='jpg':
        from_png(args.input_frames, args.output_name)
    elif args.input_frames[0].endswith('.fits'):
        fits_frames = png_from_fits(args.input_frames, args.cmap, args.norm)
        from_png(fits_frames, args.output_name)
        os.system("if [ -d framesfits ]; then rm -rf framesfits; fi")
    else:
        sys.exit("ERROR: unexpected input files (expecting png, jpg, or fits files)")


if __name__ == '__main__':
    main()
