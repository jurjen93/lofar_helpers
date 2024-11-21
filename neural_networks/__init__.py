import argparse
import functools

import torch
from torch.nn.functional import interpolate
import os

from cortexchange.architecture import Architecture
import __main__
from astropy.io import fits

from .train_nn import (
    ImagenetTransferLearning,
    load_checkpoint,
    normalize_inputs,
)  # noqa
from .pre_processing_for_ml import normalize_fits

setattr(__main__, "TransferLearning", ImagenetTransferLearning)


def process_fits(fits_path):
    with fits.open(fits_path) as hdul:
        image_data = hdul[0].data

    return normalize_fits(image_data)


class TransferLearning(Architecture):
    def __init__(
        self,
        model_name: str = None,
        device: str = None,
        variational_dropout: int = 0,
        **kwargs,
    ):
        super().__init__(model_name, device)

        self.dtype = torch.bfloat16

        self.model = self.model.to(self.dtype)
        self.model.eval()

        assert variational_dropout >= 0
        self.variational_dropout = variational_dropout

    def load_checkpoint(self, path) -> torch.nn.Module:
        # To avoid errors on CPU
        if "gpu" not in self.device and self.device != "cuda":
            os.environ["XFORMERS_DISABLED"] = "1"
        (
            model,
            _,
            self.args,
        ) = load_checkpoint(path, self.device).values()
        self.resize = self.args["resize"]
        self.lift = self.args["lift"]
        return model

    @functools.lru_cache(maxsize=1)
    def prepare_data(self, input_path: str) -> torch.Tensor:
        input_data: torch.Tensor = torch.from_numpy(process_fits(input_path))
        input_data = input_data.to(self.dtype)
        input_data = input_data.swapdims(0, 2).unsqueeze(0)
        return self.prepare_batch(input_data)

    def prepare_batch(self, batch: torch.Tensor, mean=None, std=None) -> torch.Tensor:
        batch = batch.to(self.dtype).to(self.device)
        if self.resize != 0:
            batch = interpolate(
                batch, size=self.resize, mode="bilinear", align_corners=False
            )
        if mean is None:
            mean = self.args["dataset_mean"]
        if std is None:
            std = self.args["dataset_std"]
        batch = normalize_inputs(batch, mean, std, normalize=1)
        return batch

    @torch.no_grad()
    def predict(self, data: torch.Tensor):
        with torch.autocast(dtype=self.dtype, device_type=self.device):
            if self.variational_dropout > 0:
                self.model.train()
            else:
                self.model.eval()

            predictions = torch.concat(
                [
                    torch.sigmoid(self.model(data)).clone()
                    for _ in range(max(self.variational_dropout, 1))
                ],
                dim=1,
            )

            mean = predictions.mean(dim=1)
            std = predictions.std(dim=1)

        print(mean, std)
        return mean, std

    @staticmethod
    def add_argparse_args(parser: argparse.ArgumentParser) -> None:
        parser.add_argument(
            "--variational_dropout",
            type=int,
            default=0,
            help="Optional: Amount of times to run the model to obtain a variational estimate of the stdev",
        )
