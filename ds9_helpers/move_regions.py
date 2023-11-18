"""
Move ds9 region files from existing ds9.
Save them under a new name in $PWD/boxes such that this script will recognize you made changes.
When closing, the new region files will be generated.

Note that this script was originally written to move calibrator boxes, but it will probably also work for other type of
regions (polygons, lines, etc..)
"""

from glob import glob
import os
from argparse import ArgumentParser


def move_regions(file, folder):
    """Move ds9 region files"""

    try:
        print('Opening ds9...'
              '\nIf you wish to make changes, please save the new full regions file under a new name in ' + folder + '/boxes.')
        current_boxes = glob(folder + '/boxes/*')
        os.system("ds9 {FILE} -regions load all '{DATALOC}/boxes/*.reg'".format(FILE=file, DATALOC=folder))
        new_box = [b for b in glob(folder + '/boxes/*') if b not in current_boxes]
        if len(new_box) == 1:
            os.system('mkdir ' + folder + '/boxestemp')
            for nb in new_box:
                with open(nb, 'r') as f:
                    for line in f:
                        if '{box' in line:
                            g = open(folder + '/boxestemp/' + line.strip().split()[-1].
                                     replace('text={', '').replace('}', '') + '.reg', "a")
                            g.write("fk5\n"+line)
                            g.close()
            os.system('rm -rf ' + folder + '/boxes && mv ' + folder + '/boxestemp ' + folder + '/boxes')

        print('Closed ds9.')
    except:
        print("Failing to open ds9..."
              "\nds9 {FILE} -regions load all '{DATALOC}/boxes/*.reg'".format(FILE=file, DATALOC=folder))


def parse_args():
    """Command line parser"""

    parser = ArgumentParser()
    parser.add_argument('-p', '--folder_path', type=str, help='path to folder', default='.')
    parser.add_argument('-f', '--file', type=str, help='fits file name',
                        default='image_full_ampphase_di_m.NS.app.restored.fits')
    args = parser.parse_args()

    folder = args.folder_path
    if folder[-1] == '/':
        folder = folder[0:-1]

    return args, folder


def main():
    """Main function"""

    args, folder = parse_args()

    move_regions(args.file, folder)


if __name__ == '__main__':
    main()
