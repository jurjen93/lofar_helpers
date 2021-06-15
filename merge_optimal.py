import tables
from astropy.wcs import WCS
from astropy.io import fits
from argparse import ArgumentParser, ArgumentTypeError
from math import pi, cos, sin, acos
from losoto.h5parm import h5parm
from numpy import ones, zeros

def str2bool(v):
    v = str(v)
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise ArgumentTypeError('Boolean value expected.')

def degree_to_radian(inp):
    """degree to radian"""
    return float(inp) / 360 * pi * 2

def radian_to_degree(inp):
    """radion to degree"""
    return float(inp) * 360 / (pi * 2)

def angular_distance(p1, p2):
    """angular distance for points in ra and dec"""
    return acos(sin(p1[1])*sin(p2[1])+cos(p1[1])*cos(p2[1])*cos(p1[0]-p2[0]))

def create_new_soltab(h5_in, h5_out, solset, soltab, directions, sources):
    """
    Create a new dataset in the h5 table
    :param filename: name of ourput file
    :param solset: solution set name
    :param soltab: solution table name
    :param dirs: directions to include
    """

    print('Filter {solset}/{soltab} from {h5_in} into {h5_out}'.format(solset=solset, soltab=soltab, h5_in=h5_in, h5_out=h5_out))

    h5_in = h5parm(h5_out, readonly=True)
    h5_out = h5parm(h5_out, readonly=False)

    solutiontable = h5_in.getSolset(solset).getSoltab(soltab)
    axes = solutiontable.getValues()[1]
    values_in = solutiontable.getValues()[0]
    indexes = [list(axes['dir']).index(dir.decode('UTF-8')) for dir in directions]
    axes['dir']=directions
    dir_index = solutiontable.getAxesNames().index('dir')
    h5_in.close()
    shape = list(values_in.shape)
    shape[dir_index]=len(directions)
    if 'amplitude' in soltab:
        values_new = zeros(shape)
    elif 'phase' in soltab:
        values_new = ones(shape)
    else:
        print('Skip {soltab}'.format(soltab=soltab))
        return

    for idx_new, idx_old in enumerate(indexes):
        if dir_index == 0:
            if 'amplitude' in soltab:
                values_new[idx_new,...] *= values_in[idx_old, ...]
            elif 'phase' in soltab:
                values_new[idx_new,...] += values_in[idx_old, ...]
        elif dir_index == 1:
            if 'amplitude' in soltab:
                values_new[:, idx_new,...] *= values_in[:, idx_old, ...]
            elif 'phase' in soltab:
                values_new[:, idx_new,...] += values_in[:, idx_old, ...]
        elif dir_index == 2:
            if 'amplitude' in soltab:
                values_new[:, :, idx_new,...] *= values_in[:, :, idx_old, ...]
            elif 'phase' in soltab:
                values_new[:, :, idx_new,...] += values_in[:, :, idx_old, ...]
        elif dir_index == 3:
            if 'amplitude' in soltab:
                values_new[:, :, :, idx_new,...] *= values_in[:, :, :, idx_old, ...]
            elif 'phase' in soltab:
                values_new[:, :, :, idx_new,...] += values_in[:, :, :, idx_old, ...]
        elif dir_index == 4:
            if 'amplitude' in soltab:
                values_new[:, :, :, :, idx_new,...] *= values_in[:, :, :, :, idx_old,...]
            elif 'phase' in soltab:
                values_new[:, :, :, :, idx_new,...] += values_in[:, :, :, :, idx_old,...]

    print('New number of sources {num}'.format(num=len(sources)))
    print('Filtered output shape {shape}'.format(shape=values_new.shape))

    solsetout = h5_out.getSolset(solset)
    current_sources = [source[0].decode('UTF-8') for source in solsetout.obj.source[:]]
    new_sources = [source for source in sources if source[0] not in current_sources]
    if len(new_sources) > 0:
        solsetout.obj.source.append(new_sources)


    weights = ones(values_new.shape)
    solsetout.makeSoltab(soltab, axesNames=list(axes.keys()), axesVals=list(axes.values()), vals=values_new,
                         weights=weights)

    h5_out.close()

parser = ArgumentParser()
parser.add_argument('-f', '--fits', type=str, help='fitsfile name')
parser.add_argument('-h5o', '--output_h5', type=str, help='name of output h5')
parser.add_argument('-ac', '--angular_cutoff', type=float, default=None, help='angular distances higher than this value from the center will be excluded from the box selection')
parser.add_argument('-in', '--inside', type=str2bool, default=False, help='keep directions inside the angular cutoff')
parser.add_argument('-h5out', '--h5_file_out', type=str, help='h5 files with directions outside of the angular cutoff')
parser.add_argument('-h5in', '--h5_file_in', type=str, help='h5 files with directions inside the angular cutoff')

args = parser.parse_args()

hdu = fits.open(args.fits)[0]
header = WCS(hdu.header, naxis=2).to_header()

center= (header['CRVAL1'], header['CRVAL2'])

H = tables.open_file(args.h5_file_out)
sources = []
directions = []
for dir in H.root.sol000.source[:]:
    position = [radian_to_degree(i) for i in dir[1]]
    print(angular_distance(center, position))
    if args.inside and angular_distance(center, position)<args.angular_cutoff:
        print('Keep {dir}'.format(dir=dir))
        directions.append(dir[0])
        sources.append(dir)
    elif not args.inside and angular_distance(center, position)>=args.angular_cutoff:
        print('Keep {dir}'.format(dir=dir))
        directions.append(dir[0])
        sources.append(dir)
    else:
        print('Remove {dir}'.format(dir=dir))
H.close()

h5 = h5parm(args.h5_file_out)
for ss in h5.getSolsetNames():
    for st in h5.getSolset(ss).getSoltabNames():
        print(ss, st)
        create_new_soltab(h5_in, h5_out, ss, st, directions, sources)
h5.close()