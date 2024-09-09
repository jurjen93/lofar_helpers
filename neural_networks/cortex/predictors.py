import argparse
import functools
import logging
import os
import tarfile
from typing import Optional

import torch
from tqdm import tqdm
from webdav3.client import Client

from cortex.pre_processing_for_ml import process_fits
from cortex.train_nn import load_checkpoint
from cortex.train_nn import ImagenetTransferLearning  # noqa


def download_model(cache, model):
    tarred_file = f'{model}.tar.gz'
    full_path_tar = os.path.join(cache, tarred_file)

    options = {
        'webdav_hostname': "https://surfdrive.surf.nl/files/public.php/webdav/",
        'webdav_login': "5lnKaoagQi92y0j",
        'webdav_password': "1234"
    }

    client = Client(options)
    files = client.list()
    if tarred_file not in files:
        logging.error(f"Available files: {files}")
        raise ValueError("No file exists remotely with this name.")

    bar = None

    def progress(current, total):
        nonlocal bar
        if bar is None:
            bar = tqdm(total=total, unit_scale=True, unit="B", unit_divisor=1024)
        bar.n = current
        bar.refresh()

    os.makedirs(cache, exist_ok=True)
    client.download(tarred_file, full_path_tar, progress=progress)

    # Extract tar and remove original file.
    tar = tarfile.TarFile(full_path_tar)
    tar.extractall(cache)
    tar.close()
    os.remove(full_path_tar)


class StopPredictor:
    def __init__(self, cache, model, device: str, variational_dropout: int = 0):
        self.dtype = torch.float32
        self.device = device
        model_path = os.path.join(cache, model)
        if not os.path.exists(model_path):
            download_model(cache, model)

        checkpoint = load_checkpoint(model_path, device)
        self.model = checkpoint.get("model").to(self.dtype)
        self.model.eval()

        assert variational_dropout >= 0
        self.variational_dropout = variational_dropout

    @functools.lru_cache(maxsize=1)
    def _prepare_data(self, input_path):
        input_data: torch.Tensor = torch.from_numpy(process_fits(input_path))
        input_data = input_data.to(self.dtype)
        input_data = input_data.swapdims(0, 2).unsqueeze(0)
        return input_data

    @torch.no_grad()
    def predict(self, input_path):
        data = self._prepare_data(input_path)

        with torch.autocast(dtype=self.dtype, device_type=self.device):
            if self.variational_dropout > 0:
                self.model.feature_extractor.eval()
                self.model.classifier.train()

            predictions = torch.concat(
                [torch.sigmoid(self.model(data)).clone() for _ in range(self.variational_dropout)],
                dim=1
            )

            mean = predictions.mean()
            std = predictions.std()

        print(mean, std)
        return mean, std


def process_args():
    parser = argparse.ArgumentParser(
        description="Arguments for dolly inference/test code"
    )

    parser.add_argument("--cache", type=str, default=".cache", help="Where to store the downloaded model weights.")
    parser.add_argument(
        "--model",
        type=str,
        default="version_7743995_4__model_resnext101_64x4d__lr_0.001__normalize_0__dropout_p_0.25__use_compile_1",
        help="Name of the model."
    )
    parser.add_argument("--input", type=str, default=None, help="Path to the fits file.")
    parser.add_argument("--device", type=str, default="cpu", help="Device for inference, default=cpu.")
    parser.add_argument(
        "--variational_dropout",
        type=int,
        default=None,
        help="Optional: Amount of times to run the model to obtain a variational estimate of the stdev"
    )
    return parser.parse_args()


def main(args):
    predictor = StopPredictor(
        cache=args.cache,
        device=args.device,
        model=args.model,
        variational_dropout=args.variational_dropout
    )
    print("Initialized models")
    predictor.predict(input_path=args.input)


if __name__ == "__main__":
    main(process_args())
