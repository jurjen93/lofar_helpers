"""
With this script you can move ds9 region files.
Save them under a new name in $PWD/boxes such that this script will recognize you made changes.
When closing, the new region files will be generated.
"""

from glob import glob
import os

def move_boxes(file, folder):
    try:
        print('Opening ds9 to verify box selections and make manual changes if needed.'
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
                            g = open(folder + '/boxestemp/' + line.strip().split()[-1].replace('text={', '').replace('}', '') + '.reg', "a")
                            g.write("fk5\n"+line)
                            g.close()
            os.system('rm -rf ' + folder + '/boxes && mv ' + folder + '/boxestemp ' + folder + '/boxes')

        print('Closed ds9.')
    except:
        print("Failing to open ds9 to verify box selection, check if installed and try to run on the commandline"
              "\nds9 {FILE} -regions load all '{DATALOC}/boxes/*.reg'".format(FILE=file, DATALOC=folder))

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('-l', '--location', type=str, help='data location folder name', default='.')
    parser.add_argument('-f', '--file', type=str, help='fitsfile name',
                        default='image_full_ampphase_di_m.NS.app.restored.fits')
    args = parser.parse_args()

    folder = args.location
    if folder[-1] == '/':
        folder = folder[0:-1]

    move_boxes(args.file, folder)
