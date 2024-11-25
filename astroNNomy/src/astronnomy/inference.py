import argparse
import os
from functools import partial
from typing import Callable, Any, Sequence

import matplotlib.image
import torch
from torch.utils.data import DataLoader
from tqdm import tqdm

from train_nn import ImagenetTransferLearning, load_checkpoint
from pre_processing_for_ml import FitsDataset


def variational_dropout(model, batch, variational_iters: int):

    model.feature_extractor.eval()

    if variational_iters == 0:
        # Disable dropout
        model.classifier.eval()
        variational_iters = 1
    else:
        # Enable dropout
        model.classifier.train()

    model.cuda()

    with torch.amp.autocast(device_type="cuda", dtype=torch.bfloat16):
        sample = batch[0].cuda(non_blocking=True).expand(-1, 3, -1, -1)

        means, stds = variational_dropout_step(model, sample, variational_iters)

    return means, stds


@torch.no_grad()
def variational_dropout_step(model: torch.nn.Module, sample, variational_iters: int):
    preds = torch.sigmoid(
        torch.concat([model(sample).clone() for _ in range(variational_iters)], dim=1)
    )

    means = preds.mean(dim=1)

    if preds.shape[1] == 1:
        stds = torch.zeros_like(means)
    else:
        stds = preds.std(dim=1)

    return means, stds


def save_images(dataset, out_path, preds, stds):
    os.makedirs(out_path, exist_ok=True)

    for elem, in_path, pred, std in tqdm(
        zip(dataset, dataset.data_paths, preds, stds), total=len(dataset)
    ):
        batch, label = elem
        name = in_path.strip(".npz").split("/")[-1]
        matplotlib.image.imsave(
            fname=f"{out_path}/{std:.3f}_{pred:.3f}_{label}_{name}.png",
            arr=batch.to(torch.float).movedim(0, -1).numpy(),
        )


def gen_output_from_dataset(
    dataset: torch.utils.data.Dataset, inference_f: Any, **dataloader_kwargs
):

    num_workers = dataloader_kwargs.pop(
        "num_workers", min(18, len(os.sched_getaffinity(0)))
    )
    batch_size = dataloader_kwargs.pop("batch_size", 64)

    dataloader = DataLoader(
        dataset=dataset,
        batch_size=batch_size,
        num_workers=num_workers,
        pin_memory=True,
        shuffle=False,
        drop_last=False,
        **dataloader_kwargs,
    )

    outputs = [inference_f(batch=batch) for batch in dataloader]

    return tuple(
        map(
            lambda x: torch.concat(x)
            .to(torch.device("cpu"), dtype=torch.float)
            .numpy(),
            zip(*outputs),
        )
    )


def dataset_inference_vi(
    dataset: torch.utils.data.Dataset,
    checkpoint_path: str,
    variational_iters: int,
):
    torch.set_float32_matmul_precision("high")

    model = load_checkpoint(checkpoint_path)["model"]

    inference_f = partial(
        variational_dropout, model=model, variational_iters=variational_iters
    )

    preds, stds = gen_output_from_dataset(dataset=dataset, inference_f=inference_f)

    return preds, stds


def is_dir(path):
    if not os.path.isdir(path):
        raise argparse.ArgumentTypeError(f"'{path}' is not a valid directory.")
    return path


def is_file(path):
    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError(f"'{path}' is not a valid file.")
    return path


def positive_int(value):
    def err(value):
        raise argparse.ArgumentTypeError(f"'{value}' is not a valid positive integer.")

    try:
        ivalue = int(value)
    except:
        err(value)

    if ivalue < 1:
        err(value)

    return ivalue


def get_args():
    parser = argparse.ArgumentParser(
        description="Argument parser for dataset and variational iterations"
    )

    parser.add_argument(
        "--dataset_root",
        type=is_dir,
        required=True,
        help="Path to the dataset root directory",
    )
    parser.add_argument(
        "--checkpoint_path",
        type=is_file,
        required=True,
        help="Path to the checkpoint root directory",
    )
    parser.add_argument(
        "--save_images_path",
        type=is_dir,
        default=None,
        help="Path to save images (optional)",
    )
    parser.add_argument(
        "--variational_iters",
        type=positive_int,
        default=5,
        help="Number of variational iterations (must be >= 1)",
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()

    dataset = FitsDataset(dataset_root, mode=mode)
    dataset_inference(dataset, args.variational_iters, args.save_images_path)

    if save_images_path is not None:
        out_path = save_images_path + f"/preds_{mode}"
        save_images(dataset, out_path, preds, stds)

    # root = f'{os.environ["TMPDIR"]}/public.spider.surfsara.nl/project/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/cnn_data/'
    # root = f'/dev/shm/scratch-shared/CORTEX/public.spider.surfsara.nl/project/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/cnn_data/'
    # checkpoint_path = './gridsearch_efficientnet/version_6665871_5__model_efficientnet_v2_l__lr_0.0001__normalize_0__dropout_p_0.25/ckpt_step=7999.pth'
    # out_path = f'/scratch-shared/CORTEX/public.spider.surfsara.nl/project/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/preds_{mode}/'
