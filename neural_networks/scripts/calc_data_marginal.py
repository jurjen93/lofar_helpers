import itertools
import math
import os
from pathlib import Path

import numpy as np
import torch
import tqdm
from astropy.io import fits
from matplotlib import pyplot as plt

from pre_processing_for_ml import FitsDataset


def main_fits(root_folder):
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

def main_dataloader(root_folder):
    num_workers = 6
    prefetch_factor, persistent_workers = (
        (2, True) if num_workers > 0 else
        (None, False)
    )

    train_dataloader, val_dataloader = (
        torch.utils.data.DataLoader(
            dataset=FitsDataset(root_folder, mode=mode),
            batch_size=1,
            num_workers=num_workers,
            prefetch_factor=prefetch_factor,
            persistent_workers=persistent_workers,
            pin_memory=False,
            drop_last=False,  # needed for torch.compile,
            shuffle=False,  # needed so that val AP is non nan
        )
        for mode in ('train', 'validation')
    )

    arr_size = 1, 3, 2048, 2048
    train_arr = np.empty([len(train_dataloader), *arr_size])
    val_arr = np.empty([len(val_dataloader), *arr_size])

    for i, (data, label) in tqdm.tqdm(enumerate(train_dataloader), total=len(train_dataloader)):
        train_arr[i] = data.numpy()


    for i, (data, label) in tqdm.tqdm(enumerate(val_dataloader), total=len(val_dataloader)):
        val_arr[i] = data.numpy()

    print(train_arr.mean(axis=(0, 1, 3, 4)), train_arr.std(axis=(0, 1, 3, 4)))
    print(val_arr.mean(axis=(0, 1, 3, 4)), val_arr.std(axis=(0, 1, 3, 4)))

    exit(0)
    flattened_val = val_arr.flatten()

    plt.hist([flattened_train, flattened_val], density=True, bins=25, label=['train', 'val'])
    plt.legend()
    plt.title('marginal distribution of observation values')
    plt.savefig('marginal_train_val_dataloader_exped.png')



if __name__ == '__main__':
    # main_fits('/dev/shm/scratch-shared/CORTEX/calibrator_selection_robertjan/autosettings/')

    main_dataloader('/dev/shm/scratch-shared/CORTEX/public.spider.surfsara.nl/project/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/cnn_data/')
