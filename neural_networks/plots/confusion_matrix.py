from cortexchange.architecture import get_architecture, Architecture
from pathlib import Path
import sys
import os

SCRIPT_DIR = Path(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(SCRIPT_DIR))
from train_nn import MultiEpochsDataLoader
from pre_processing_for_ml import FitsDataset
import matplotlib.pyplot as plt
import numpy as np
import torch
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay


def load_model(architecture_name, model_name):
    StopPredictor: type(Architecture) = get_architecture(architecture_name)
    predictor = StopPredictor(device="cuda", model_name=model_name)
    return predictor


def get_dataloader(data_root, mode, batch_size):
    num_workers = min(12, len(os.sched_getaffinity(0)))

    prefetch_factor, persistent_workers = (
        (2, True) if num_workers > 0 else (None, False)
    )

    return MultiEpochsDataLoader(
        dataset=FitsDataset(
            data_root,
            mode=mode,
        ),
        batch_size=batch_size,
        num_workers=num_workers,
        prefetch_factor=prefetch_factor,
        persistent_workers=persistent_workers,
        pin_memory=True,
        shuffle=False,
        drop_last=False,
    )


def get_statistics(data_root, mode):
    return FitsDataset(
        data_root,
        mode=mode,
    ).compute_statistics(1)


@torch.no_grad()
def get_confusion_matrix(
    predictor, dataloader, mean, std, thresholds=[0.2, 0.3, 0.4, 0.5]
):
    confusion_matrices = np.zeros((len(thresholds), 2, 2))
    thresholds = torch.tensor(thresholds)
    for img, label in dataloader:
        data = predictor.prepare_batch(img, mean=mean, std=std)
        pred = torch.sigmoid(predictor.model(data)).to("cpu")
        preds_thres = pred >= thresholds
        for i, _ in enumerate(thresholds):
            confusion_matrices[i] += confusion_matrix(
                label, preds_thres[:, i], labels=[0, 1]
            )

    for i, conf_matrix in enumerate(confusion_matrices):

        disp = ConfusionMatrixDisplay(
            # Normalization
            conf_matrix / np.sum(conf_matrix, axis=1, keepdims=True),
            display_labels=["continue", "stop"],
        )
        disp.plot()

        plt.savefig(f"confusion_thres_{thresholds[i]:.3f}.png")


if __name__ == "__main__":
    model_name = "surf/dinov2_09814"
    architecture_name = "surf/TransferLearning"
    predictor = load_model(architecture_name, model_name)
    if hasattr(predictor, "args") and "dataset_mean" in predictor.args:
        mean, std = predictor.args["dataset_mean"], predictor.args["dataset_std"]
    else:
        mean, std = get_statistics(
            "/scratch-shared/CORTEX/public.spider.surfsara.nl/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/cnn_data/",
            mode="train",
        )

    dataloader = get_dataloader(
        "/scratch-shared/CORTEX/public.spider.surfsara.nl/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/cnn_data/",
        mode="val",
        batch_size=32,
    )
    get_confusion_matrix(predictor, dataloader, mean, std)
