import os
from typing import Any

import numpy as np
import torch
import torchmetrics
from lightning import LightningModule, Trainer
from matplotlib import pyplot as plt
from torch import nn
from torchvision import models
from torchmetrics.functional import accuracy

from pre_processing_for_ml import FitsDataset


class ImagenetTransferLearning(LightningModule):
    def __init__(self):
        super().__init__()

        # init a pretrained resnet
        backbone = models.resnet50(weights="DEFAULT")

        layers = list(backbone.children())[:-1]
        self.feature_extractor = nn.Sequential(*layers)
        self.feature_extractor.eval()

        num_target_classes = 1
        num_filters = backbone.fc.in_features

        self.classifier = nn.Sequential(
            nn.Linear(num_filters, num_filters),
            nn.ReLU(),
            nn.Linear(num_filters, num_target_classes),
        )

    def forward(self, x):
        with torch.no_grad():
            representations = self.feature_extractor(x).flatten(1)
        x = self.classifier(representations)
        return x

    def training_step(self, batch):
        inputs, target = batch

        output = self(inputs)

        target = target[None].T.to(dtype=output.dtype)  # massaging the targets into the correct dtype
        loss = torch.nn.functional.binary_cross_entropy(torch.sigmoid(output), target)
        train_accuracy = accuracy(output, target, task="binary")

        self.log("train_bce_loss", loss, on_step=True, on_epoch=False, prog_bar=True, logger=True)
        self.log("train_accuracy", train_accuracy, on_step=True, on_epoch=False, prog_bar=True, logger=True)
        return loss

    def validation_step(self, batch):
        inputs, target = batch

        output = self(inputs)

        target = target[None].T.to(dtype=output.dtype)  # massaging the targets into the correct dtype
        val_loss = torch.nn.functional.binary_cross_entropy(torch.sigmoid(output), target)
        val_accuracy = accuracy(output, target, task="binary")

        self.log("val_bce_loss", val_loss, on_step=False, on_epoch=True, prog_bar=True, logger=True)
        self.log("val_accuracy", val_accuracy, on_step=False, on_epoch=True, prog_bar=True, logger=True)

        return val_loss

    def configure_optimizers(self):
        return torch.optim.AdamW(self.classifier.parameters(), lr=1e-3)


def main(root: str):
    torch.set_float32_matmul_precision('high')

    trainer = Trainer()

    model = ImagenetTransferLearning()

    num_workers = len(os.sched_getaffinity(0))
    # num_workers = 0
    prefetch_factor, persistent_workers = (
        (2, True) if num_workers > 0 else
        (None, False)
    )

    train_dataloader, val_dataloader = (
        torch.utils.data.DataLoader(
            dataset=FitsDataset(root, mode=mode),
            batch_size=32,
            num_workers=num_workers,
            prefetch_factor=prefetch_factor,
            persistent_workers=persistent_workers,
            pin_memory=False,
            drop_last=(True if mode == 'train' else False),  # needed for torch.compile,
            shuffle=(True if mode == 'train' else False),
        )
        for mode in ('train', 'validation')
    )

    trainer.fit(model, train_dataloaders=train_dataloader, val_dataloaders=val_dataloader)


def plot_marginal(root):
    Idat = FitsDataset(root)
    data = np.concatenate([elem[0].flatten() for elem in Idat])
    data = data / data.std()

    plt.hist(data, density=True, bins=50)
    # plt.yscale('log')
    plt.title('marginal distribution of observation values')
    plt.savefig('marginal_flipped_doubled.png')


if __name__ == '__main__':
    root = '/dev/shm/public.spider.surfsara.nl/project/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/cnn_data/'

    main(root)
