import random
from pathlib import Path

import torch
from torch.utils.data import Dataset
import os
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from PIL import Image
import torchvision.transforms as T


def get_rms(data: np.ndarray, maskSup=1e-7):
    """
    find the rms of an array, from Cycil Tasse/kMS

    :param data: numpy array
    :param maskSup: mask threshold

    :return: rms --> rms of image
    """

    mIn = np.ndarray.flatten(data)
    m = mIn[np.abs(mIn) > maskSup]
    rmsold = np.std(m)
    diff = 1e-1
    cut = 3.
    med = np.median(m)

    for i in range(10):
        ind = np.where(np.abs(m - med) < rmsold * cut)[0]
        rms = np.std(m[ind])
        if np.abs((rms - rmsold) / rmsold) < diff:
            break
        rmsold = rms

    return rms  # jy/beam


def make_image(imdat, plotname):
    """
    Make an image of data
    """
    plt.imshow(imdat)
    # plt.colorbar()
    # plt.clim(0, 1)
    plt.savefig(plotname)
    plt.close()


def crop(data):
    """
    Crop image
    """
    width, height = list(np.array(data.shape) // 2)

    half_width = width // 2
    half_height = height // 2
    x, y = half_width * 2, half_height * 2

    x_min = max(x - half_width, 0)
    x_max = min(x + half_width, data.shape[1])
    y_min = max(y - half_height, 0)
    y_max = min(y + half_height, data.shape[0])

    cutout = data[y_min:y_max, x_min:x_max]
    return cutout


class FitsDataset(Dataset):
    def __init__(self, root_dir, mode='train', train_split=0.85):
        """
        Args:
            root_dir (string): Directory with good/bad folders in it.
        """

        assert mode in ('train', 'validation')
        assert 0 <= train_split <= 1

        random.seed(42)

        classes = {'stop': 0, 'continue': 1}
        # Yes this code is way overengineered. Yes I also derive pleasure from writing it :) - RJS
        #
        # Actual documentation:
        # You want all 'self.x' variables to be non-python objects such as numpy arrays,
        # otherwise you get memory leaks in the PyTorch dataloader
        data_paths, labels = map(np.asarray, list(zip(*(
            (str(file), val)
            for cls, val in classes.items()
            for file in (Path(root_dir, cls).rglob('*.fits'))
        ))))

        # Filter found data_paths and labels on train/val splits
        # TODO: make this a bit cleaner

        # First randomly permute so that the train/val splits are randomly distributed
        rng = np.random.default_rng(42)
        p = rng.permutation(len(labels))
        data_paths, labels = data_paths[p], labels[p]

        filtered_paths, filtered_labels = [], []
        for val in classes.values():
            indices = np.nonzero(labels == val)[0]
            split = train_split if mode == 'train' else 1 - train_split
            num_examples = np.int32((np.floor if mode == 'train' else np.ceil)(split * len(indices)))

            if num_examples == 0:
                print(f"Warning: num remaining examples is 0 for class {val} and split {split}")

            filtered_idx = (indices[:num_examples] if mode == 'train' else indices[-num_examples:],)
            filtered_paths.append(data_paths[filtered_idx])
            filtered_labels.append(labels[filtered_idx])

        self.data_paths, self.labels = map(np.concatenate, (filtered_paths, filtered_labels))

        sources = ", ".join(sorted([str(elem).split('/')[-1].strip('.fits') for elem in self.data_paths]))
        print(f'{mode}: using the following sources: {sources}')


    @staticmethod
    def transform_data(image_data):
        """
        Transform data for preprocessing
        """
        while image_data.ndim > 2:
            image_data = image_data[0]

        # crop data (half data size)
        image_data = crop(image_data)

        # re-normalize data (such that values are between 0 and 1)
        rms = get_rms(image_data)
        norm = SymLogNorm(linthresh=rms * 2, linscale=2, vmin=-rms, vmax=rms * 50000, base=10)

        image_data = norm(image_data)
        image_data = np.clip(image_data - image_data.min(), a_min=0, a_max=1)

        # make RGB image
        cmap = plt.get_cmap('RdBu_r')
        image_data = cmap(image_data)
        image_data = np.delete(image_data, 3, 2)

        image_data = -image_data + 1  # make the peak exist at zero
        # image_data = np.e ** image_data
        #
        # image_data = image_data - np.array([1.34517796, 1.2179354 , 1.15546612])
        # image_data = image_data / np.array([0.37540418, 0.18891336, 0.13838379])


        # if random.randint(0, 1):
        #     image_data = -image_data

        image_data = torch.from_numpy(image_data).to(dtype=torch.float32)
        image_data = torch.movedim(image_data, -1, 0)
        # image_data = T.Resize((256, 256))(image_data)
        # Just some regularization
        image_data = image_data + 0.1*torch.randn_like(image_data)
        image_data = image_data

        return image_data

    def __len__(self):
        return len(self.data_paths)

    def __getitem__(self, idx):
        fits_path = self.data_paths[idx]
        with fits.open(fits_path) as hdul:
            image_data = hdul[0].data

        # Pre-processing
        image_data = self.transform_data(image_data.astype(np.float32))

        label = self.labels[idx]

        return image_data, label


if __name__ == '__main__':
    root = '/project/lofarvwf/Public/jdejong/CORTEX/calibrator_selection_robertjan/cnn_data'

    Idat = FitsDataset(root)
    for n in range(Idat.__len__()):
        imdat, label = Idat.__getitem__(n)
        make_image(imdat, f'{label}_{Idat.data[n].split("/")[-1].replace(".fits", "")}.png')
