import itertools
import os
from functools import lru_cache
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import torch
from astropy.io import fits
from joblib import Parallel, delayed, Memory
import joblib
from matplotlib.colors import SymLogNorm
from torch.utils.data import Dataset


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
    cut = 3.0
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


def normalize_fits(image_data: np.ndarray):
    image_data = image_data.squeeze()

    # Pre-processing
    rms = get_rms(image_data)
    norm_f = SymLogNorm(
        linthresh=rms * 2, linscale=2, vmin=-rms, vmax=rms * 50000, base=10
    )

    image_data = norm_f(image_data)

    image_data = image_data - image_data.min()
    image_data = np.clip(image_data, a_min=0, a_max=1)

    return image_data[..., None]


def transform_data(root_dir, classes=("continue", "stop"), modes=("", "_val")):

    def process_fits(fits_path):
        with fits.open(fits_path) as hdul:
            image_data = hdul[0].data
        if image_data.shape[2] == 2048:
            transformed = normalize_fits(image_data)

            np.savez_compressed(
                fits_path.with_suffix(".npz"), transformed.astype(np.float32)
            )
        else:
            print(f"Skipping {fits_path}. Improper image size {image_data.shape}")
        # print(fits_path)

    root_dir = Path(root_dir)
    assert root_dir.exists()

    Parallel(n_jobs=len(os.sched_getaffinity(0)))(
        delayed(process_fits)(fits_path)
        for cls, mode in itertools.product(classes, modes)
        for fits_path in (root_dir / (cls + mode)).glob("*.fits")
    )


class FitsDataset(Dataset):
    def __init__(self, root_dir, mode="train"):
        """
        Args:
            root_dir (string): Directory with good/bad folders in it.
        """

        modes = ("train", "val")
        assert mode in modes

        classes = {"stop": 0, "continue": 1}

        root_dir = Path(root_dir)
        assert root_dir.exists(), f"'{root_dir}' doesn't exist!"

        ext = ".npz"
        glob_ext = "*" + ext

        self.root_dir = root_dir

        for folder in (
            root_dir / (cls + ("" if mode == "train" else "_val")) for cls in classes
        ):
            assert (
                folder.exists()
            ), f"root folder doesn't exist, got: '{str(folder.resolve())}'"
            assert (
                len(list(folder.glob(glob_ext))) > 0
            ), f"no '{ext}' files were found in '{str(folder.resolve())}'"

        # Yes this code is way overengineered. Yes I also derive pleasure from writing it :) - RJS
        #
        # Actual documentation:
        # You want all 'self.x' variables to be non-python objects such as numpy arrays,
        # otherwise you get memory leaks in the PyTorch dataloader
        self.data_paths, self.labels = map(
            np.asarray,
            list(
                zip(
                    *(
                        (str(file), val)
                        for cls, val in classes.items()
                        for file in (
                            root_dir / (cls + ("" if mode == "train" else "_val"))
                        ).glob(glob_ext)
                    )
                )
            ),
        )

        assert len(self.data_paths) > 0
        self.sources = ", ".join(
            sorted([str(elem).split("/")[-1].strip(ext) for elem in self.data_paths])
        )
        self.mode = mode
        _, counts = np.unique(self.labels, return_counts=True)
        self.label_ratio = counts[0] / counts[1]
        # print(f'{mode}: using the following sources: {sources}')

    def compute_statistics(self, normalize):
        cache = Memory(location=self.root_dir / "_cache")
        cached_compute = cache.cache(FitsDataset._compute_statistics)
        self.mean, self.std = cached_compute(self, normalize)
        return self.mean, self.std

    # For caching
    def __reduce__(self):
        return (
            self.__class__,
            (self.labels, self.mode, self.label_ratio, self.sources),
        )

    @staticmethod
    def _compute_statistics(loader, normalize, verbose=True):
        if verbose:
            print("Computing dataset statistics")
        if not normalize:
            return torch.asarray([0]), torch.asarray([1])

        means = []
        sums_of_squares = []
        f = (lambda x: torch.log(x + 1e-10)) if normalize == 2 else lambda x: x
        for i, (imgs, _) in enumerate(loader):
            imgs = f(imgs)
            means.append(torch.mean(imgs, dim=(1, 2)))
            sums_of_squares.append((imgs**2).sum(dim=(1, 2)))

        mean = torch.stack(means).mean(0)
        sums_of_squares = torch.stack(sums_of_squares).sum(0)
        variance = (sums_of_squares / (len(loader) * imgs.shape[1] * imgs.shape[2])) - (
            mean**2
        )
        if verbose:
            print("Finished Computing dataset statistics")
        return mean, torch.sqrt(variance)

    @staticmethod
    def transform_data(image_data):
        """
        Transform data for preprocessing
        """

        # FIXME: this should really be a parameter
        image_data = torch.from_numpy(image_data).to(torch.bfloat16)
        image_data = torch.movedim(image_data, -1, 0)

        return image_data

    @lru_cache(maxsize=1)
    def __len__(self):
        return len(self.data_paths)

    def __getitem__(self, idx):

        npy_path = self.data_paths[idx]
        label = self.labels[idx]
        image_data = np.load(npy_path)["arr_0"]  # there is always only one array

        # Pre-processing
        image_data = self.transform_data(image_data)

        return image_data, label


def make_histogram(root_dir):
    root_dir = Path(root_dir)
    assert root_dir.exists(), f"'{root_dir}' doesn't exist!"

    ext = ".npz"
    glob_ext = "*" + ext

    classes = {"stop": 0, "continue": 1}
    file_paths = []
    for mode in ["train", "val"]:
        for folder in (
            root_dir / (cls + ("" if mode == "train" else "_val")) for cls in classes
        ):
            assert (
                folder.exists()
            ), f"root folder doesn't exist, got: '{str(folder.resolve())}'"
            assert (
                len(list(folder.glob(glob_ext))) > 0
            ), f"no '{ext}' files were found in '{str(folder.resolve())}'"

            # Yes this code is way overengineered. Yes I also derive pleasure from writing it :) - RJS
            #
            # Actual documentation:
            # You want all 'self.x' variables to be non-python objects such as numpy arrays,
            # otherwise you get memory leaks in the PyTorch dataloader
        file_paths += [
            str(file)
            for cls, _ in classes.items()
            for file in (root_dir / (cls + ("" if mode == "train" else "_val"))).glob(
                glob_ext
            )
        ]

    result = np.asarray(
        list(
            map(
                lambda file_path: file_path.split("/")[-1].split("_")[0].split("-")[0],
                file_paths,
            )
        )
    )
    unique, inverse, counts = np.unique(result, return_inverse=True, return_counts=True)
    print(f"number of unique stations: {len(unique)}")
    plt.ylabel("Count")
    plt.xlabel("Number of images per station")
    plt.hist(counts, bins=len(np.unique(counts)))
    plt.savefig("image_count_histogram.png")

    plt.hist(inverse, bins=len(unique))
    plt.ylabel("Number of images")
    plt.xlabel("Station")
    plt.savefig("images_per_station.png")


if __name__ == "__main__":
    root = f"/scratch-shared/CORTEX/public.spider.surfsara.nl/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/cnn_data"
    transform_data(root)

    # make_histogram(root)
    # dataset.compute_statistics(normalize=1)

    # dataset = FitsDataset(root, mode="train")
    # sources = dataset.sources
    # hash(dataset)
    # print(sources)
    # imgs, label = dataset[0]
    # from PIL import Image
    # plt.imshow(imgs.permute(1, 2, 0).to(torch.float32).numpy())
    # plt.savefig('test.png')

    # for img, label in dataset:
    #     if img.shape[2] != 2048:
    #         exit()

    # images = np.concatenate([image.flatten() for image, label in Idat])
    # print("creating hist")
    # plt.hist(images)
    # plt.savefig('preselect_fig.png')
    # make_image(imdat, f'{label}_{Idat.data[n].split("/")[-1].replace(".fits", "")}.png')
