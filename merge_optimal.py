import tables
from astropy.wcs import WCS
from astropy.io import fits
from argparse import ArgumentParser
from math import pi, cos, sin, acos

def degree_to_radian(inp):
    """degree to radian"""
    return float(inp) / 360 * pi * 2

def radian_to_degree(inp):
    """degree to radian"""
    return float(inp) * 360 / (pi * 2)

def angular_distance(p1, p2):
    """angular distance for points in ra and dec"""
    return acos(sin(p1[1])*sin(p2[1])+cos(p1[1])*cos(p2[1])*cos(p1[0]-p2[0]))


parser = ArgumentParser()
parser.add_argument('-f', '--file', type=str, help='fitsfile name')
parser.add_argument('-ac', '--angular_cutoff', type=float, default=None, help='angular distances higher than this value from the center will be excluded from the box selection')
args = parser.parse_args()



hdu = fits.open(args.file)[0]
header = WCS(hdu.header, naxis=2).to_header()

center= (header['CRVAL1'], header['CRVAL2'])

directions = []
H = tables.open_file('lotss.h5')
for dir in H.root.sol000.source[:]:
    position = [radian_to_degree(i) for i in dir[1]]
    if angular_distance(center, position)>args.angular_cutoff:
        print(dir)
H.close()
