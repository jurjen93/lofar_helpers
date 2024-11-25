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
from sklearn.metrics import precision_recall_curve
from astropy.io import fits
import torcheval.metrics.functional as tef


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
def get_dropout_output(predictor, dataloader, mean, std, vi_iters_list):
    labels = []
    vi_dict = {
        vi_iters: {"std": [], "pred": [], "labels": []} for vi_iters in vi_iters_list
    }
    for i, vi_iters in enumerate(vi_iters_list):
        for img, label in dataloader:
            data = predictor.prepare_batch(img, mean=mean, std=std)
            if not i:
                labels += label.numpy().tolist()

            predictor.variational_dropout = vi_iters
            pred, stds = predictor.predict(data.clone())
            vi_dict[vi_iters]["pred"] += pred.cpu().to(torch.float32).numpy().tolist()
            vi_dict[vi_iters]["std"] += stds.cpu().to(torch.float32).numpy().tolist()
        vi_dict[vi_iters]["pred"] = torch.asarray(vi_dict[vi_iters]["pred"])
        vi_dict[vi_iters]["std"] = torch.asarray(vi_dict[vi_iters]["std"])
        vi_dict[vi_iters]["labels"] = torch.asarray(labels)
    return vi_dict


def plot_pr_curves(savedir, vi_dict, vi_iters_list):
    os.makedirs(savedir, exist_ok=True)
    # for vi_iter, pred_dict in vi_dict.items():
    for vi_iter in sorted(vi_iters_list):
        pred_dict = vi_dict[vi_iter]
        preds, labels = pred_dict["pred"], pred_dict["labels"]
        # Reverse labels to compute pr curve for predicting "stop"
        precision, recall, thresholds = precision_recall_curve(
            np.asarray(labels), np.asarray(preds)
        )
        auprc = tef.binary_auprc(torch.asarray(preds), torch.asarray(labels))
        print(f"auprc vi_iters {vi_iter}", auprc)

        plt.plot(
            recall, precision, label=f"VI Iters: {vi_iter} auprc: {auprc.item():.3f}"
        )
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title("PR Curves for varying Variational Inference iteration counts")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.7)

    plt.tight_layout()
    plt.ylim(0, 1)
    plt.savefig(f"{savedir}/precision_recall_curves.png", dpi=300)
    plt.clf()


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
        shuffle=False,
        num_workers=num_workers,
        persistent_workers=persistent_workers,
        prefetch_factor=prefetch_factor,
        drop_last=False,
    )

    return dataloader


if __name__ == "__main__":
    # Latest model
    model_name = "surf/dinov2_09739_rotations"
    # model_name = "surf/dinov2_097_rotations_dropout_01"
    TESTING = True
    architecture_name = "surf/TransferLearning"
    # Set Device here
    DEVICE = "cuda"
    # Thresholds to consider for classification
    vi_iters_list = [0, 1, 2, 4, 8, 16, 32, 64, 128]
    # vi_iters_list = [0]
    # Change to directory of files. Should have subfolders 'continue_val' and 'stop_val'
    data_root = "/scratch-shared/CORTEX/public.spider.surfsara.nl/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/cnn_data"
    # Uses cached results for testing the plotting functionalities
    savedir = model_name.split("/")[-1]
    dict_fname = f"{savedir}/pr_curve_dict.npydict"
    VI_DICT, LABELS = None, None
    # Load cached results
    if Path(dict_fname).exists() and TESTING:
        VI_DICT = np.load(dict_fname, allow_pickle=True)[()]
        existing_vi_iters = list(VI_DICT.keys())
        vi_iters_list = list(set(vi_iters_list) - set(existing_vi_iters))

    if vi_iters_list != []:

        dataloader = get_dataloader(data_root, mode="val")

        predictor = load_model(architecture_name, model_name, device=DEVICE)

        mean, std = predictor.args["dataset_mean"], predictor.args["dataset_std"]

        vi_dict = get_dropout_output(predictor, dataloader, mean, std, vi_iters_list)
        # Add new results to cached results and save
        if VI_DICT is not None:
            VI_DICT = vi_dict | VI_DICT
        else:
            VI_DICT = vi_dict

        os.makedirs(savedir, exist_ok=True)
        with open(dict_fname, "wb") as f:
            np.save(f, VI_DICT)

    plot_pr_curves(savedir, VI_DICT, vi_iters_list)
