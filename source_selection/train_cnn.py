import os
from functools import partial, lru_cache

import numpy as np
import torch
import torcheval.metrics.functional as tef
from lightning import LightningModule, Trainer
from matplotlib import pyplot as plt
from torch import nn
from torchvision import models
from torchvision.transforms import v2

from pre_processing_for_ml import FitsDataset


def setup_transferred_model(name: str, num_target_classes: int):
    # use partial to prevent loading all models at once
    model_map = {
        'resnet50': partial(models.resnet50, weights="DEFAULT"),
        'resnet152': partial(models.resnet152, weights="DEFAULT"),
        'resnext50_32x4d': partial(models.resnext50_32x4d, weights="DEFAULT"),
        'resnext101_64x4d': partial(models.resnext101_64x4d, weights="DEFAULT"),
        'efficientnet_v2_l': partial(models.efficientnet_v2_l, weights="DEFAULT"),
    }

    backbone = model_map[name]()

    feature_extractor = nn.Sequential(*list(backbone.children())[:-1])
    num_out_features = (
        backbone.fc if name in ('resnet50', 'resnet152', 'resnext101_64x4d')
        else backbone.classifier[-1]  # efficientnet
    ).in_features

    classifier = nn.Sequential(
        nn.Flatten(),
        nn.Dropout1d(p=0.25),
        nn.Linear(num_out_features, num_out_features),
        nn.ReLU(),
        nn.Linear(num_out_features, num_target_classes),
    )

    return feature_extractor, classifier


class ImagenetTransferLearning(LightningModule):
    def __init__(self, model: str):
        super().__init__()
        self.save_hyperparameters()

        self.feature_extractor, self.classifier = setup_transferred_model(name=model, num_target_classes=1)
        self.feature_extractor.eval()

        # self.average_precision = nn.ModuleDict({
        #     'train_ap': BinaryAveragePrecision(),
        #     'val_ap': BinaryAveragePrecision()
        # })

    @partial(torch.compile, mode='reduce-overhead')
    def forward(self, x):
        with torch.no_grad():
            representations = self.feature_extractor(x)
        x = self.classifier(representations)
        return x

    @torch.no_grad()
    def augmentation(self, inputs):
        inputs = get_transforms()(inputs)
        inputs = inputs + 0.01 * torch.randn_like(inputs)

        return inputs

    @torch.no_grad()
    def normalize(self, inputs):
        inputs = (inputs - inputs.mean(dim=(0, 2, 3), keepdim=True)) / inputs.std(dim=(0, 2, 3), keepdim=True)
        return inputs

    def training_step(self, batch):
        loss = self.shared_step(batch, mode='train')
        return loss

    def validation_step(self, batch):
        loss = self.shared_step(batch, mode='val')
        return loss

    def predict_step(self, batch):
        pass

    def shared_step(self, batch, mode):
        assert mode in ('train', 'val', 'test')

        inputs, target = batch

        # inputs = self.normalize(inputs)

        if mode == 'train':
            inputs = self.augmentation(inputs)

        # TODO: do this from the global mean

        output = self(inputs).flatten()

        target = target.to(dtype=output.dtype)  # massaging the targets into the correct dtype
        loss = torch.nn.functional.binary_cross_entropy_with_logits(output, target, pos_weight=torch.as_tensor(2.698717948717949))

        with torch.no_grad():
            pseudo_probs = torch.sigmoid(output)

            ap = tef.binary_auprc(pseudo_probs, target)
            acc = tef.binary_accuracy(pseudo_probs, target)

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
        return torch.optim.AdamW(self.classifier.parameters(), lr=3e-5)


def main(root: str):
    torch.set_float32_matmul_precision('high')

    trainer = Trainer(benchmark=True, precision='bf16-mixed', num_sanity_val_steps=0)

    model = ImagenetTransferLearning(model='resnet50')

    # num_workers = len(os.sched_getaffinity(0))
    # num_workers = 0
    num_workers = 12
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


@lru_cache(maxsize=1)
def get_transforms():
    return v2.Compose([
        # v2.Resize(size=1024),
        v2.ColorJitter(brightness=.5, hue=.3, saturation=0.1, contrast=0.1),
        v2.RandomInvert(),
        v2.RandomEqualize(),
        v2.RandomVerticalFlip(p=0.5),
        v2.RandomHorizontalFlip(p=0.5),
    ])


if __name__ == '__main__':
    # root = f'{os.environ["TMPDIR"]}/public.spider.surfsara.nl/project/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/cnn_data/'
    root = f'/dev/shm/public.spider.surfsara.nl/project/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/cnn_data/'

    main(root)
