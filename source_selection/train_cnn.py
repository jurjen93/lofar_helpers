import argparse
import itertools
import os
from functools import partial, lru_cache
from pathlib import Path

import torch
import torcheval.metrics.functional as tef
from torch import nn, binary_cross_entropy_with_logits
from torch.utils.data import SequentialSampler, Sampler, RandomSampler
from torch.utils.tensorboard import SummaryWriter
from torchvision import models
from torchvision.transforms import v2
from tqdm import tqdm

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


@torch.no_grad()
def normalize_inputs(inputs, mode=0):
    assert mode in range(3)

    if mode == 0:
        # Actual op instead of simple return because of debugging reasons
        means, stds = [0, 0, 0], [1, 1, 1]
    elif mode == 1:
        # Inputs are lognormal -> log to make normal
        means, stds = [-1.55642344, -1.75137082, -2.13795913], [1.25626133, 0.79308821, 0.7116124]
        inputs = inputs.log()
    else:
        # Inputs are lognormal
        means, stds = [0.35941373, 0.23197646, 0.15068751], [0.28145176, 0.17234328, 0.10769559]

    # Resulting shape of means and stds: [1, 3, 1, 1]
    means, stds = map(
        lambda x: torch.tensor(x, device=inputs.device).reshape(1, 3, 1, 1),
        (means, stds)
    )

    return (inputs - means) / stds


@torch.no_grad()
def augmentation(inputs):
    inputs = get_transforms()(inputs)
    inputs = inputs + 0.01 * torch.randn_like(inputs)

    return inputs



class ImagenetTransferLearning(nn.Module):
    def __init__(
            self,
            model_name: str = 'resnet50',
            dropout_p: float = 0.25
    ):
        super().__init__()

        # For saving in the state dict
        self.kwargs = {'model': model_name, 'dropout_p': dropout_p}

        self.feature_extractor, num_features = get_feature_extractor(name=model_name)
        self.feature_extractor.eval()

        self.classifier = get_classifier(dropout_p, n_features=num_features, num_target_classes=1)

    @partial(torch.compile, mode='reduce-overhead')
    def forward(self, x):
        with torch.no_grad():
            representations = self.feature_extractor(x)
        x = self.classifier(representations)
        return x

    def step(self, inputs, targets):
        logits = self(inputs).flatten()

        loss = binary_cross_entropy_with_logits(
            logits,
            targets,
            pos_weight=torch.as_tensor(2.698717948717949)
        )

        if PROFILE:
            global profiler
            profiler.step()

        return logits, loss


def get_dataloaders(dataset_root, batch_size):
    num_workers = min(12, len(os.sched_getaffinity(0)))
    prefetch_factor, persistent_workers = (
        (2, True) if num_workers > 0 else
        (None, False)
    )

    loaders = tuple(
        MultiEpochsDataLoader(
            dataset=FitsDataset(dataset_root, mode=mode),
            batch_size=batch_size,
            num_workers=num_workers,
            prefetch_factor=prefetch_factor,
            persistent_workers=persistent_workers,
            pin_memory=True,
            shuffle=True if mode == 'train' else False,
            drop_last=True if mode == 'train' else False,
        )
        for mode in ('train', 'val')
    )

    return loaders


def get_logging_dir(logging_root: str, /, **kwargs):
    # As version string, prefer $SLURM_ARRAY_JOB_ID, then $SLURM_JOB_ID, then 0.
    version = int(os.getenv('SLURM_ARRAY_JOB_ID', os.getenv('SLURM_JOB_ID', 0)))
    version_appendix = int(os.getenv('SLURM_ARRAY_TASK_ID', 0))
    while True:
        version_dir = "__".join((
            f"version_{version}_{version_appendix}",
            *(f"{k}_{v}" for k, v in kwargs.items())
        ))

        logging_dir = Path(logging_root) / version_dir

        if not logging_dir.exists():
            break
        version_appendix += 1

    return str(logging_dir.resolve())

def get_tensorboard_logger(logging_dir):
    writer = SummaryWriter(log_dir=str(logging_dir))

    # writer.add_hparams()

    return writer


def write_metrics(writer, metrics: dict, global_step: int):
    for metric_name, value in metrics.items():
        if isinstance(value, tuple):
            probs, labels = value
            writer_fn = partial(
                writer.add_pr_curve, labels=labels, predictions=probs,
            )
        else:
            writer_fn = partial(
                writer.add_scalar, scalar_value=value
            )

        writer_fn(tag=metric_name, global_step=global_step)


def get_optimizer(parameters: list[torch.Tensor], lr: float):
    return torch.optim.AdamW(parameters, lr=lr)


def merge_metrics(suffix, **kwargs):
    return {
        f"{k}/{suffix}": v for k, v in kwargs.items()
    }


@torch.no_grad()
def prepare_data(data, labels, normalize, device):
    data, labels = (
        data.to(device, non_blocking=True, memory_format=torch.channels_last),
        labels.to(device, non_blocking=True, dtype=data.dtype)
    )

    data = normalize_inputs(data, mode=normalize)

    return data, labels

def main(dataset_root: str, model_name: str, lr: float, normalize: int, dropout_p: float, batch_size: int):
    torch.set_float32_matmul_precision('high')

    profiler_kwargs = {}

    if PROFILE:
        pass
        # profiling_dir = str(logging_dir)
        #
        # global profiler
        # profiler = torch.profiler.profile(
        #     schedule=torch.profiler.schedule(wait=0, warmup=1, active=9, repeat=0),
        #     on_trace_ready=torch.profiler.tensorboard_trace_handler(profiling_dir),
        # )
        #
        # profiler.start()

    logging_dir = get_logging_dir(
        str(Path.cwd() / 'runs'),
        # kwargs
        model=model_name,
        lr=lr,
        normalize=normalize,
        dropout_p=dropout_p
    )

    writer = get_tensorboard_logger(logging_dir)

    device = torch.device('cuda')

    model: nn.Module = ImagenetTransferLearning(model_name=model_name, dropout_p=dropout_p)
    model.feature_extractor.eval()

    # noinspection PyArgumentList
    model.to(device=device, memory_format=torch.channels_last)

    optimizer = get_optimizer(model.classifier.parameters(), lr=lr)

    train_dataloader, val_dataloader = get_dataloaders(dataset_root, batch_size)

    logging_interval = 10

    train_step_f, val_step_f = (
        partial(
            step_f,
            model=model,
            prepare_data_f=partial(prepare_data, device=device, normalize=normalize),
            metrics_logger=partial(log_metrics, write_metrics_f=partial(write_metrics, writer=writer))
        )
        for step_f in (train_step, val_step)
    )

    train_step_f = partial(train_step_f, train_dataloader=train_dataloader, optimizer=optimizer, logging_interval=logging_interval)
    val_step_f = partial(val_step_f, val_dataloader=val_dataloader)

    checkpoint_saver = partial(save_checkpoint, logging_dir=logging_dir, model=model, optimizer=optimizer, normalize=normalize, batch_size=batch_size)

    best_val_loss = torch.inf
    global_step = 0  # make it a tensor so we can do in-place edits

    with torch.cuda.amp.autocast(dtype=torch.bfloat16):

        n_epochs = 250
        for epoch in range(n_epochs):

            global_step = train_step_f(global_step=global_step)
            val_loss = val_step_f(global_step=global_step)

            if val_loss < best_val_loss:
                checkpoint_saver(global_step=global_step)
                best_val_loss = val_loss

    if PROFILE:
        profiler.stop()

    writer.flush()
    writer.close()


@torch.no_grad()
def log_metrics(loss, logits, targets, log_suffix, global_step, write_metrics_f):
    probs = torch.sigmoid(logits)
    ap = tef.binary_auprc(probs, targets)

    metrics = merge_metrics(
        bce_loss=loss,
        au_pr_curve=ap,
        pr_curve=(probs, targets),
        suffix=log_suffix
    )

    write_metrics_f(metrics=metrics, global_step=global_step)

@torch.no_grad()
def val_step(model, val_dataloader, global_step, metrics_logger, prepare_data_f):
    val_losses, val_logits, val_targets = [], [], []

    model.classifier.eval()
    for i, (data, labels) in tqdm(enumerate(val_dataloader), desc='Validation', total=len(val_dataloader)):
        # print("validation start")

        data, labels = prepare_data_f(data, labels)

        logits, loss = model.step(data, labels)

        val_losses.append(loss)
        val_logits.append(logits.clone())
        val_targets.append(labels)

    losses, logits, targets = map(torch.concatenate, (val_losses, val_logits, val_targets))

    mean_loss = losses.mean()
    metrics_logger(loss=mean_loss, logits=logits, targets=targets, global_step=global_step, log_suffix='Validation')

    return mean_loss



def train_step(model, optimizer, train_dataloader, prepare_data_f, global_step, logging_interval, metrics_logger):
    # print("training")
    model.classifier.train()

    for i, (data, labels) in tqdm(enumerate(train_dataloader), desc='Training', total=len(train_dataloader)):
        global_step += 1

        data, labels = prepare_data_f(data, labels)

        data = augmentation(data)

        optimizer.zero_grad(set_to_none=True)

        logits, loss = model.step(data, labels)
        mean_loss = loss.mean()

        mean_loss.backward()
        optimizer.step()


        if i % logging_interval == 0:
            with torch.no_grad():
                metrics_logger(loss=mean_loss.detach(), logits=logits.detach(), targets=labels, global_step=global_step, log_suffix='training')

    return global_step


class MultiEpochsDataLoader(torch.utils.data.DataLoader):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._DataLoader__initialized = False
        if self.batch_sampler is None:
            self.sampler = _RepeatSampler(self.sampler)
        else:
            self.batch_sampler = _RepeatSampler(self.batch_sampler)
        self._DataLoader__initialized = True
        self.iterator = super().__iter__()

    def __len__(self):
        return len(self.sampler) if self.batch_sampler is None else len(self.batch_sampler.sampler)

    def __iter__(self):
        for i in range(len(self)):
            yield next(self.iterator)


class _RepeatSampler(object):
    """ Sampler that repeats forever.

    Args:
        sampler (Sampler)
    """

    def __init__(self, sampler):
        self.sampler = sampler

    def __iter__(self):
        while True:
            yield from iter(self.sampler)


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


def save_checkpoint(logging_dir, model, optimizer, global_step, **kwargs):
    old_checkpoints = Path(logging_dir).glob('*.pth')
    for old_checkpoint in old_checkpoints:
        Path.unlink(old_checkpoint)

    torch.save(
        {
            'model': type(model),
            'model_state_dict': model.state_dict(),
            'optimizer': type(optimizer),
            'optimizer_state_dict': optimizer.state_dict(),
            'step': global_step,
            **kwargs
        },
        f=(logging_dir + f'/ckpt_step={global_step}.pth')
    )

def load_checkpoint(ckpt_path):

    ckpt_dict = torch.load(ckpt_path)

    # ugh, this is so ugly, something something hindsight something something 20-20
    # FIXME: probably should do a pattern match, but this works for now
    kwargs = str(Path(ckpt_path).parent).split('/')[-1].split('__')

    # strip 'model_' from the name
    model_name = kwargs[1][6:]
    lr = float(kwargs[2].split('_')[-1])
    normalize = int(kwargs[3].split('_')[-1])
    dropout_p = float(kwargs[4].split('_')[-1])

    model = ckpt_dict['model'](model=model_name, dropout_p=dropout_p)
    model.load_state_dict(ckpt_dict['model_state_dict'])

    # FIXME: add optim class and args to state dict
    optim = ckpt_dict.get('optimizer', torch.optim.AdamW)(
        lr=lr,
        params=model.classifier.parameters()
    ).load_state_dict(ckpt_dict['optimizer_state_dict'])

    return {'model': model, 'optim': optim, 'normalize': normalize}

def get_argparser():
    """
    Create and return an argument parser for hyperparameter tuning.
    """
    # Create the parser
    parser = argparse.ArgumentParser(description="Hyperparameter tuning for a machine learning model.")

    # Add arguments
    parser.add_argument('dataset_root', type=Path)
    parser.add_argument('--lr', type=float, help='Learning rate for the model.', default=1e-4)
    parser.add_argument('--batch_size', type=int, help='Batch size', default=64)
    parser.add_argument('--model_name', type=str, help='The model to use.', default='resnet50')
    parser.add_argument('--normalize', type=int, help='Whether to do normalization', default=0, choices=[0, 1, 2])
    parser.add_argument('--dropout_p', type=float, help='Dropout probability', default=0.25)
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
