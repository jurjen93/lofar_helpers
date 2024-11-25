from cortexchange.architecture import get_architecture, Architecture
from pathlib import Path
import sys
import os

SCRIPT_DIR = Path(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(SCRIPT_DIR))
from pre_processing_for_ml import normalize_fits
import matplotlib.pyplot as plt
import numpy as np
import torch
from functools import lru_cache
from torch.utils.data import Dataset, DataLoader
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
from astropy.io import fits


class RawFitsDataset(Dataset):
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

        ext = ".fits"
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

        fits_path = self.data_paths[idx]
        label = self.labels[idx]

        image_data = process_fits(fits_path)
        # there is always only one array

        # Pre-processing
        image_data = self.transform_data(image_data)

        return image_data, label


def load_model(architecture_name, model_name, device="cpu"):
    StopPredictor: type(Architecture) = get_architecture(architecture_name)
    predictor = StopPredictor(device=device, model_name=model_name)
    return predictor


@torch.no_grad()
def get_confusion_matrix(predictor, dataloader, mean, std, thresholds):
    conf_matrix_dict = {thres: np.zeros((2, 2)) for thres in thresholds}
    thresholds_tensor = torch.tensor(thresholds)
    for i, (img, label) in enumerate(dataloader):
        data = predictor.prepare_batch(img, mean=mean, std=std)
        pred = torch.sigmoid(predictor.model(data)).to("cpu")
        preds_thres = pred >= thresholds_tensor
        for i, thres in enumerate(thresholds):
            conf_matrix_dict[thres] += confusion_matrix(
                label, preds_thres[:, i], labels=[0, 1]
            )

    return conf_matrix_dict


def plot_conf_matrices(savedir, conf_matrix_dict, thresholds):
    os.makedirs(savedir, exist_ok=True)
    for thres, conf_matrix in conf_matrix_dict.items():

        disp = ConfusionMatrixDisplay(
            # Normalization
            conf_matrix / np.sum(conf_matrix, axis=1, keepdims=True),
            display_labels=["stop", "continue"],
        )
        # print(conf_matrix)
        disp.plot()

        plt.savefig(f"{savedir}/confusion_thres_{thres:.3f}.png")


def process_fits(fits_path):
    with fits.open(fits_path) as hdul:
        image_data = hdul[0].data

    return normalize_fits(image_data)


def get_dataloader(data_root, mode="val", batch_size=32):
    dataset = RawFitsDataset(data_root, mode="val")
    num_workers = min(12, len(os.sched_getaffinity(0)))

    prefetch_factor, persistent_workers = (
        (2, True) if num_workers > 0 else (None, False)
    )
    dataloader = DataLoader(
        dataset,
        batch_size=32,
        shuffle=True,
        num_workers=num_workers,
        persistent_workers=persistent_workers,
        prefetch_factor=prefetch_factor,
        drop_last=False,
    )

    return dataloader


if __name__ == "__main__":
    # Latest model
    model_name = "surf/dinov2_09739_rotations"
    TESTING = True
    architecture_name = "surf/TransferLearning"
    # Set Device here
    DEVICE = "cuda"
    # Thresholds to consider for classification
    thresholds = [0.2, 0.3, 0.4, 0.5]
    # Change to directory of files. Should have subfolders 'continue_val' and 'stop_val'
    data_root = "/scratch-shared/CORTEX/public.spider.surfsara.nl/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/cnn_data"
    # Uses cached confusion matrix for testing the plotting functionalities
    savedir = model_name.split("/")[-1]
    dict_fname = f"{savedir}/conf_matrix.npydict"
    print(Path(dict_fname).exists(), dict_fname)
    if Path(dict_fname).exists() and TESTING:
        conf_matrix_dict = np.load(dict_fname, allow_pickle=True)[()]
    else:

        dataloader = get_dataloader(data_root, mode="val")

        predictor = load_model(architecture_name, model_name, device=DEVICE)

        mean, std = predictor.args["dataset_mean"], predictor.args["dataset_std"]

        conf_matrix_dict = get_confusion_matrix(
            predictor, dataloader, mean, std, thresholds
        )
        with open(dict_fname, "wb") as f:
            np.save(f, conf_matrix_dict)

    print(conf_matrix_dict)

    plot_conf_matrices(savedir, conf_matrix_dict, thresholds)
