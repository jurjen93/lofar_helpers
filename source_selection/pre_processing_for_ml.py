import random

import torch
from torch.utils.data import Dataset
import os
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from PIL import Image
import torchvision.transforms as T

def get_rms(data: np.ndarray, maskSup = 1e-7):
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
    width, height = list(np.array(data.shape)//2)

    half_width = width // 2
    half_height = height // 2
    x, y = half_width*2, half_height*2

    x_min = max(x - half_width, 0)
    x_max = min(x + half_width, data.shape[1])
    y_min = max(y - half_height, 0)
    y_max = min(y + half_height, data.shape[0])

    cutout = data[y_min:y_max, x_min:x_max]
    return cutout


class FitsDataset(Dataset):
    def __init__(self, root_dir):
        """
        Args:
            root_dir (string): Directory with good/bad folders in it.
        """
        self.root_dir = root_dir
        self.data = []
        self.labels = []
        self.classes = os.listdir(root_dir)

        # Load dataset
        for cls in self.classes:
            cls_path = os.path.join(self.root_dir, cls)
            if not os.path.isdir(cls_path):
                continue

            for file in os.listdir(cls_path):
                if file.endswith('.fits'):
                    self.data.append(os.path.join(cls_path, file))
                    label = 'stop' in cls
                    self.labels.append(int(label))

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
        norm = SymLogNorm(linthresh=rms * 2, linscale=2, vmin=-rms, vmax=rms*50000, base=10)
        image_data = norm(image_data)
        image_data = np.clip(image_data - image_data.min(), a_min=0, a_max=1)

        # make RGB image
        cmap = plt.get_cmap('RdBu_r')
        image_data = cmap(image_data)
        image_data = np.delete(image_data, 3, 2)

        image_data = -image_data + 1  # make the peak exist at zero

        if random.randint(0, 1):
            image_data = -image_data

        image_data = torch.from_numpy(image_data)
        image_data = torch.movedim(image_data, -1, 0)
        image_data = T.Resize((256, 256))(image_data)

        return image_data

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        fits_path = self.data[idx]
        with fits.open(fits_path) as hdul:
            image_data = hdul[0].data

        # Pre-processing
        image_data = self.transform_data(image_data.astype(np.float64))

        label = self.labels[idx]

        return image_data, label


if __name__ == '__main__':
    root = '/project/lofarvwf/Public/jdejong/CORTEX/calibrator_selection_robertjan/cnn_data'

    Idat = FitsDataset(root)
    for n in range(Idat.__len__()):
        imdat, label = Idat.__getitem__(n)
        make_image(imdat, f'{label}_{Idat.data[n].split("/")[-1].replace(".fits","")}.png')
