"""
Split a ds9 region file with multiple facets into separate facet region files.
You need this if you want to make images of individual facets instead of full facet-imaging in wsclean.
The script also returns a polygon_info.csv containing the center of the polygon and the calibrator source direction with a polygon area
and an estimate for how many times you can average based on the measurement sets from a 2.5x2.5 degree wide-field (for long-baselines)

Example:
    python split_polygon_facets.py --reg all_facets.reg --h5 merged.h5 --extra_boundary 0.1

    This example takes a region file 'all_facets.reg' with X-number of polygons and splits into X separate polygon files.
    It loops over the directions from the corresponding h5 solution file with all solution per facet.

IMPORTANT: The script assumes that the polygons and h5 solutions are sorted in the same order and the coordinates
are in degrees!
"""


from shapely import geometry
import numpy as np
import tables
from glob import glob
import csv
from argparse import ArgumentParser
from astropy.coordinates import SkyCoord


def make_utf8(inp):
    """
    Convert input to utf8 instead of bytes

    :param inp: string input
    """
    try:
        inp = inp.decode('utf8')
        return inp
    except (UnicodeDecodeError, AttributeError):
        return inp


def split_polygons_ds9(regionfile, extra_boundary=0.):
    """
    Split polygons

    :param regionfile: region file
    :param extra_boundary: adding extra boundary layer
    :return:
    """

    regionfile = open(regionfile, 'r')
    lines = regionfile.readlines()
    header = lines[0:4]
    polygons = lines[4:]
    for n, poly in enumerate(polygons):
        if extra_boundary>0:
            poly_file = open('poly_' + str(n)+'.reg', 'w')
        else:
            poly_file = open('poly_' + str(n) + '.reg', 'w')
        poly_file.writelines(header)
        polyp = [float(p) for p in poly.replace('polygon(', '').replace(')', '').replace('\n', '').split(',')]
        poly_geo = geometry.Polygon(tuple(zip(polyp[0::2], polyp[1::2])))
        if extra_boundary!=0.:
            poly_geo = poly_geo.buffer(extra_boundary, resolution=len(polyp[0::2]), join_style=2)
        poly = 'polygon'+str(tuple(item for sublist in poly_geo.exterior.coords[:] for item in sublist))
        poly_file.writelines(poly)
    regionfile.close()


def distance(c1, c2):
    """
    Get distance based on coordiantes in degrees

    :param c1: coordinate 1
    :param c2: coordinate 2
    :return:
    """
    c1 = SkyCoord(f'{c1[0]}deg', f'{c1[1]}deg', frame='icrs')
    c2 = SkyCoord(f'{c2[0]}deg', f'{c2[1]}deg', frame='icrs')
    return c1.separation(c2).value


def point_in_polygon(point, poly_reg):
    """
    Is point in polygon region file

    :param point: list or tuple with 2D coordinate
    :param poly_reg: polygon region file
    :return: point in geo, polygon area
    """
    polyregion = open(poly_reg, 'r')
    lines = polyregion.readlines()
    poly = lines[4]
    polyregion.close()
    polyp = [float(p) for p in poly.replace('polygon(', '').replace(')', '').replace('\n', '').split(',')]
    poly_geo = geometry.Polygon(tuple(zip(polyp[0::2], polyp[1::2])))
    point_geo = geometry.Point(point)
    if poly_geo.contains(point_geo):
        c_x, c_y = poly_geo.centroid.x, poly_geo.centroid.y

        list_x = poly_geo.boundary.coords.xy[0]
        list_y = poly_geo.boundary.coords.xy[1]

        # max distance from center
        max_dist = 0
        for i in range(len(list_x)):
            dist = distance([c_x, c_y], [list_x[i], list_y[i]])
            if dist > max_dist:
                max_dist = dist
                max_point = [list_x[i], list_y[i]]

        # calculate averaging factor based on 2.5 by 2.5 degrees field size
        avg = max(int(2.5//(2*max(distance([c_x, c_y], [max_point[0], c_y]), distance([c_x, c_y], [c_x, max_point[1]])))), 1) + 1

        print(c_x, c_y, max_point, avg)
        return poly_geo.contains(point_geo), poly_geo.area, [c_x, c_y], avg

    else:
        return False, None, None, None


def parse_args():
    """Argument parser"""

    parser = ArgumentParser(description='Split multi-facet region file with polygon regions out into multiple region files')
    parser.add_argument('--reg', help='region file', type=str, required=True)
    parser.add_argument('--h5', help='h5 file to write directions from', type=str, required=True)
    parser.add_argument('--extra_boundary', help='make polygons with extra boundaries', type=float, default=0.)
    return parser.parse_args()


def main():
    """Main function"""

    args = parse_args()

    reg = args.reg
    solutionfile = args.h5

    split_polygons_ds9(regionfile=reg, extra_boundary=0)

    H = tables.open_file(solutionfile)
    dirs = H.root.sol000.source[:]['dir']
    dirname = H.root.sol000.source[:]['name']
    H.close()

    #TODO: dangerous! Convert to degrees
    if np.all(np.abs(dirs) < 2 * np.pi):
        # Converting radians to degrees
        dirs = np.degrees(dirs)
        dirs[:,0]=np.mod(dirs[:,0], 360)

    f = open('polygon_info.csv', 'w')
    writer = csv.writer(f)
    writer.writerow(['idx', 'dir_name', 'polygon_file', 'dir', 'poly_center', 'area', 'avg'])
    for n, dir in enumerate(dirs):
        for polygonregion_file in glob('poly_*.reg'):
            point_in_poly, poly_area, poly_center, avg = point_in_polygon(dir, polygonregion_file)
            if point_in_poly:
                print(n, make_utf8(dirname[n]), polygonregion_file)
                writer.writerow([n,
                                 make_utf8(dirname[n]),
                                 polygonregion_file,
                                 '['+str(dir[0])+'deg'+','+str(dir[1])+'deg'+']',
                                 '['+str(round(poly_center[0], 5))+'deg'+','+str(round(poly_center[1], 5))+'deg'+']',
                                 poly_area,
                                 avg])

    f.close()

    if args.extra_boundary > 0:
        split_polygons_ds9(regionfile=reg, extra_boundary=args.extra_boundary)

if __name__ == "__main__":
    main()
