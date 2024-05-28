from pathlib import Path

import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt


def main(root_folder):
    data = []
    for fitsfile in Path(root_folder).rglob('*.fits'):
        with fits.open(fitsfile) as fts:
            data.append(fts[0].data)

    data = np.concatenate(data).flatten()
    print(data.min(), data.max())

    plt.hist(data, density=True, bins=25)
    plt.yscale('log')
    plt.title('marginal distribution of observation values')
    plt.savefig('marginal.png')

if __name__ == '__main__':
    main('/dev/shm/scratch-shared/CORTEX/calibrator_selection_robertjan/autosettings/')
