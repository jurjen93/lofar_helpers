from typing import Any

import numpy as np
import torch
from lightning import LightningModule, Trainer
from matplotlib import pyplot as plt
from torch import nn
from torchvision import models

from pre_processing_for_ml import FitsDataset


class ImagenetTransferLearning(LightningModule):
    def __init__(self):
        super().__init__()

        # init a pretrained resnet
        backbone = models.resnet50(weights="DEFAULT")
        num_filters = backbone.fc.in_features
        layers = list(backbone.children())[:-1]
        self.feature_extractor = nn.Sequential(*layers)
        self.feature_extractor.eval()

        # use the pretrained model to classify cifar-10 (10 image classes)
        num_target_classes = 1
        self.classifier = nn.Linear(num_filters, num_target_classes)

    def forward(self, x):
        with torch.no_grad():
            representations = self.feature_extractor(x).flatten(1)
        x = self.classifier(representations)
        return x

    def training_step(self, batch):
        inputs, target = batch

        output = self(inputs.to(torch.float32))
        loss = torch.nn.functional.binary_cross_entropy(torch.sigmoid(output.reshape(-1)), target.float())

        self.log("bce_loss", loss, on_step=True, on_epoch=True, prog_bar=True, logger=True)
        return loss

    def configure_optimizers(self):
        return torch.optim.AdamW(self.classifier.parameters(), lr=1e-4)





def main(root: str):
    trainer = Trainer()

    model = ImagenetTransferLearning()

    dataloader = torch.utils.data.DataLoader(
        dataset=FitsDataset(root),
        num_workers=0,
        batch_size=32,
        pin_memory=True,
        drop_last=True,
        persistent_workers=False
    )

    trainer.fit(model, train_dataloaders=dataloader)


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