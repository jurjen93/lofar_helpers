import argparse
import os
from functools import partial, lru_cache
from pathlib import Path

import torch
import torcheval.metrics.functional as tef
from lightning import LightningModule, Trainer
from lightning.pytorch.profilers import PyTorchProfiler
from pytorch_lightning.loggers import TensorBoardLogger
from torch import nn
from torchvision import models
from torchvision.transforms import v2

from pre_processing_for_ml import FitsDataset


PROFILE = False

def get_feature_extractor(name: str):
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
        backbone.fc if name in ('resnet50', 'resnet152', 'resnext50_32x4d', 'resnext101_64x4d')
        else backbone.classifier[-1]  # efficientnet
    ).in_features

    return feature_extractor, num_out_features


def get_classifier(dropout_p: float, n_features: int, num_target_classes: int):
    assert 0 <= dropout_p <= 1

    classifier = nn.Sequential(
        nn.Flatten(),
        nn.Dropout1d(p=dropout_p),
        nn.Linear(n_features, n_features),
        nn.ReLU(),
        nn.Linear(n_features, num_target_classes),
    )

    return classifier


class ImagenetTransferLearning(LightningModule):
    def __init__(
            self,
            model: str = 'resnet50',
            lr: float = 1e-4,
            normalize: int = 0,
            dropout_p: float = 0.25
    ):
        super().__init__()
        self.save_hyperparameters()

        self.feature_extractor, num_features = get_feature_extractor(name=model)
        self.feature_extractor.eval()

        self.classifier = get_classifier(dropout_p, n_features=num_features, num_target_classes=1)

        self.lr = lr
        self.normalize = normalize

        # self.average_precision = nn.ModuleDict({
        #     'train_ap': BinaryAveragePrecision(),
        #     'val_ap': BinaryAveragePrecision()
        # })

    @torch.no_grad()
    def normalize_inputs(self, inputs):
        assert self.normalize in range(3)

        if self.normalize == 0:
            # Actual op instead of simple return because of debugging reasons
            means, stds = [0, 0, 0], [1, 1, 1]
        elif self.normalize == 1:
            # Inputs are lognormal -> log to make normal
            means, stds = [-1.55642344, -1.75137082, -2.13795913], [1.25626133, 0.79308821, 0.7116124]
            inputs = inputs.log()
        else:
            # Inputs are lognormal
            means, stds = [0.35941373, 0.23197646, 0.15068751], [0.28145176, 0.17234328, 0.10769559]

        # Resulting shape of means and stds: [1, 3, 1, 1]
        means, stds = map(
            lambda x: torch.tensor(x, device=self.device).unsqueeze(0).unsqueeze(-1).unsqueeze(-1),
            (means, stds)
        )

        return (inputs - means) / stds

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

        inputs = self.normalize_inputs(inputs)

        if mode == 'train':
            inputs = self.augmentation(inputs)

        # TODO: do this from the global mean

        output = self(inputs).flatten()

        target = target.to(dtype=output.dtype)  # massaging the targets into the correct dtype
        loss = torch.nn.functional.binary_cross_entropy_with_logits(output, target,
                                                                    pos_weight=torch.as_tensor(2.698717948717949))

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

        if PROFILE:
            profiler.step()

        return loss

    def configure_optimizers(self):
        return torch.optim.AdamW(self.classifier.parameters(), lr=self.lr)


def main(dataset_root: str, model: str, lr: float, normalize=False, dropout_p=0.1):
    torch.set_float32_matmul_precision('high')

    logging_root = 'lightning_logs/'

    version = int(os.getenv('SLURM_ARRAY_JOB_ID', os.getenv('SLURM_JOB_ID', 0)))
    version_appendix = int(os.getenv('SLURM_ARRAY_TASK_ID', 0))
    while True:
        version_dir = f"version_{version}_{version_appendix}__model_{model}__lr_{lr}__normalize_{normalize}___dropoutP_{dropout_p}"
        logging_dir = Path.cwd() / logging_root / version_dir
        if not logging_dir.exists():
            break

        version_appendix += 1

    logger = TensorBoardLogger(
        save_dir=os.getcwd(),
        name=logging_root,
        version=version_dir
    )

    profiler_kwargs = {}
    if PROFILE:
        profiling_dir = str(logging_dir)

        global profiler
        profiler = torch.profiler.profile(
            schedule=torch.profiler.schedule(wait=0, warmup=1, active=10, repeat=0),
            on_trace_ready=torch.profiler.tensorboard_trace_handler(profiling_dir),
            record_shapes=True,
            with_stack=True,
        )

        profiler_kwargs = {
            'limit_train_batches': 5,
            'limit_val_batches': 5,
            'max_epochs': 1
        }

        profiler.start()


    trainer = Trainer(
        **profiler_kwargs,

        benchmark=True,
        precision='bf16-mixed',
        num_sanity_val_steps=0,
        logger=logger

    )

    model = ImagenetTransferLearning(model=model, lr=lr, normalize=normalize, dropout_p=dropout_p)

    # num_workers = len(os.sched_getaffinity(0))
    # num_workers = 0
    num_workers = 12
    prefetch_factor, persistent_workers = (
        (2, True) if num_workers > 0 else
        (None, False)
    )

    train_dataloader, val_dataloader = (
        torch.utils.data.DataLoader(
            dataset=FitsDataset(dataset_root, mode=mode),
            batch_size=24,
            num_workers=num_workers,
            prefetch_factor=prefetch_factor,
            persistent_workers=persistent_workers,
            pin_memory=False,
            drop_last=(True if mode == 'train' else False),  # needed for torch.compile,
            shuffle=(True if mode == 'train' else False),  # needed so that val AP is non nan
        )
        for mode in ('train', 'validation')
    )

    trainer.fit(model, train_dataloaders=train_dataloader, val_dataloaders=val_dataloader)

    if PROFILE:
        profiler.stop()

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


def get_argparser():
    """
    Create and return an argument parser for hyperparameter tuning.
    """
    # Create the parser
    parser = argparse.ArgumentParser(description="Hyperparameter tuning for a machine learning model.")

    # Add arguments
    parser.add_argument('dataset_root', type=Path)
    parser.add_argument('--lr', type=float, help='Learning rate for the model.', default=1e-4)
    parser.add_argument('--model', type=str, help='The model to use.', default='resnet50')
    parser.add_argument('--normalize', type=int, help='Whether to do normalization', default=0, choices=[0, 1, 2])
    parser.add_argument('--dropout_p', type=float, help='Dropout probability', default=0.1)
    parser.add_argument('--profile', action='store_true')

    return parser.parse_args()


if __name__ == '__main__':
    kwargs = vars(get_argparser())
    print(kwargs)

    if kwargs['profile']:
        PROFILE = True
    del kwargs['profile']

    main(**kwargs)
    # '/dev/shm/scratch-shared/CORTEX/public.spider.surfsara.nl/project/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/cnn_data/'

    # # root = f'{os.environ["TMPDIR"]}/public.spider.surfsara.nl/project/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/cnn_data/'
    #
    # main(root)
