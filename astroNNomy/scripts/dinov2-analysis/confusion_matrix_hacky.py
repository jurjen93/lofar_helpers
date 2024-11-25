from pathlib import Path
import sys
import os

SCRIPT_DIR = Path(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.dirname(SCRIPT_DIR))
import matplotlib.pyplot as plt
import numpy as np
import torch
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay



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


def plot_conf_matrices(savedir, conf_matrix_dict):
    os.makedirs(savedir, exist_ok=True)
    for thres, conf_matrix in conf_matrix_dict.items():

        disp = ConfusionMatrixDisplay(
            conf_matrix / np.sum(conf_matrix, axis=1, keepdims=True),  # Normalization
            display_labels=["stop", "continue"]
        )

        fig, ax = plt.subplots()
        disp.plot(cmap='Blues', ax=ax)
        im = ax.images[0]
        im.set_clim(0, 1)

        ax.set_xlabel("Predicted Label", fontsize=16)
        ax.set_ylabel("True Label", fontsize=16)

        ax.tick_params(axis='both', which='major', labelsize=16)

        for text in disp.text_.ravel():
            text.set_fontsize(16)

        colorbar = im.colorbar
        colorbar.set_label("Accuracy", fontsize=16)
        colorbar.ax.tick_params(labelsize=16)

        fig.tight_layout(pad=2.0)

        plt.savefig(f"{savedir}/confusion_thres_{thres:.3f}.png", dpi=150)
        plt.close(fig)



if __name__ == "__main__":
    # Latest model
    model_name = "surf/dinov2_09739_rotations"
    TESTING = True
    architecture_name = "surf/TransferLearning"
    DEVICE = "cuda"
    thresholds = [0.2, 0.3, 0.4, 0.5]
    dict_fname = f"plots/dinov2_09739_rotations/conf_matrix.npydict"
    print(Path(dict_fname).exists(), dict_fname)
    if TESTING:
        conf_matrix_dict = np.load(dict_fname, allow_pickle=True)[()]

    print(conf_matrix_dict)

    plot_conf_matrices('/home/jurjen/Documents/ELAIS/paperplots', conf_matrix_dict, thresholds)
