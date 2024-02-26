import os
from typing import Any

import numpy as np
import torch
import torchmetrics
from lightning import LightningModule, Trainer
from matplotlib import pyplot as plt
from torch import nn
from torchmetrics import AveragePrecision
from torchmetrics.classification import BinaryAveragePrecision
from torchvision import models
from torchmetrics.functional import accuracy

from pre_processing_for_ml import FitsDataset




class ImagenetTransferLearning(LightningModule):
    def __init__(self):
        super().__init__()

        # init a pretrained resnet
        backbone = models.resnet50(weights="DEFAULT")
        # backbone = models.resnet152(weights="DEFAULT")
        # backbone = models.resnext101_64x4d(weights="DEFAULT")
        # backbone = models.efficientnet_v2_l

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

        self.average_precision = nn.ModuleDict({
            'train_ap': BinaryAveragePrecision(),
            'val_ap': BinaryAveragePrecision()
        })

    def forward(self, x):
        with torch.no_grad():
            representations = self.feature_extractor(x).flatten(1)
        x = self.classifier(representations)
        return x

    def training_step(self, batch):
        loss = self.shared_step(batch, mode='train')
        return loss

    def validation_step(self, batch):
        loss = self.shared_step(batch, 'val')
        return loss

    def shared_step(self, batch, mode):
        assert mode in ('train', 'val', 'test')

        inputs, target = batch

        if torch.isnan(inputs).any():
            breakpoint()

        output = torch.sigmoid(self(inputs))

        target = target[None].T.to(dtype=output.dtype)  # massaging the targets into the correct dtype
        loss = torch.nn.functional.binary_cross_entropy(output, target)
        acc = accuracy(output, target, task="binary")
        ap = self.average_precision[f"{mode}_ap"](output, target.int())

        self.log_dict(
            {
                f"{mode}_{name}": val for name, val in [
                    ("bce_loss", loss),
                    ("accuracy", acc),
                    ("average_precision", ap)
                ]
            },
            on_step=mode == 'train',
            on_epoch=mode != 'train',
            prog_bar=True,
            logger=True,
        )

        return loss

    def configure_optimizers(self):
        return torch.optim.AdamW(self.classifier.parameters(), lr=1e-4)


def main(root: str):
    torch.set_float32_matmul_precision('high')

    trainer = Trainer()

    model = ImagenetTransferLearning()

    num_workers = len(os.sched_getaffinity(0))
    # num_workers = 0
    # num_workers = 5
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
            pin_memory=True,
            drop_last=(True if mode == 'train' else False),  # needed for torch.compile,
            shuffle=True,  # needed so that val AP is non nan
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
    root = '/scratch-shared/CORTEX/public.spider.surfsara.nl/project/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/cnn_data/'

    main(root)
