from collections import defaultdict
from pickle import Pickler
import sys
import os

import tensorboard
import numpy as np
from sklearn.metrics import precision_recall_curve
import matplotlib.pyplot as plt

sys.path.append(
    ".."
)  # yes this is mega ugly, but otherwise I need to restructure the whole project...

from inference import dataset_inference_vi, is_file
from train_nn import ImagenetTransferLearning
from pre_processing_for_ml import FitsDataset

DATASET_ROOT = "/dev/shm/scratch-shared/CORTEX/public.spider.surfsara.nl/project/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/cnn_data/"
CHECKPOINT = "../grid_search/version_7944362_0__model_efficientnet_v2_l__lr_0.0001__normalize_0__dropout_p_0.25__use_compile_True/ckpt_step=1643.pth"
PICKLE_DICT_FNAME = "precision_recall_curves.npy"


if os.path.isfile(PICKLE_DICT_FNAME):
    pickle_dict = np.load(PICKLE_DICT_FNAME, allow_pickle=True)[
        ()
    ]  # Don't ask me why this is the syntax
else:
    # Yes I'm aware of the existence of collections.defaultdict, but that doesn't work with np.save...
    pickle_dict = {}

variational_iter_vals = (0, 1, 2, 4, 16)
modes = ("val", "train")

# variational_iter_vals = (0, 1,)
# modes = ('val',)

new_vals = False


for mode in modes:
    if mode not in pickle_dict:
        pickle_dict[mode] = {}

    for variational_iters in variational_iter_vals:
        variational_iters_str = f"variational_iters_{variational_iters}"

        if variational_iters_str in pickle_dict[mode]:
            print(f"{mode}/{variational_iters_str} already in saved pickle; skipping.")

        else:
            print(f"{mode}/{variational_iters_str} not found; calculating")
            new_vals = True

            pickle_dict[mode][variational_iters_str] = {}

            dataset = FitsDataset(DATASET_ROOT, mode=mode)
            preds, stds = dataset_inference_vi(
                dataset, CHECKPOINT, variational_iters=variational_iters
            )

            precision, recall, thresholds = precision_recall_curve(
                dataset.labels, preds
            )

            for name, value in (
                ("precision", precision),
                ("recall", recall),
                ("thresholds", thresholds),
                ("preds", preds),
                ("stds", stds),
                ("labels", dataset.labels),
                ("sources", dataset.data_paths),
            ):
                pickle_dict[mode][variational_iters_str][name] = value

        plt.plot(
            pickle_dict[mode][variational_iters_str]["recall"],
            pickle_dict[mode][variational_iters_str]["precision"],
            label=f"VI Iters: {variational_iters}",
        )

    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title("PR Curves for varying Variational Inference iteration counts")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.7)

    plt.tight_layout()
    plt.savefig(f"precision_recall_curves_{mode}.png", dpi=300)
    plt.clf()

if new_vals:
    print("saving new/updated pickle_dict")

    # pylance or w/e can complain, but this is valid syntax.
    np.save(PICKLE_DICT_FNAME, pickle_dict)
