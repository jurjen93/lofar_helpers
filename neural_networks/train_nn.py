import argparse
import os
from functools import partial, lru_cache
from pathlib import Path

import torch
import torcheval.metrics.functional as tef
from torch import nn, binary_cross_entropy_with_logits
from torch.nn.functional import interpolate
from torch.utils.data import SequentialSampler, Sampler
from torch.utils.tensorboard import SummaryWriter
from torchvision import models
from torchvision.transforms import v2
from tqdm import tqdm
import joblib

import numpy as np
import random

from pre_processing_for_ml import FitsDataset

PROFILE = False
SEED = None

cache = joblib.Memory(location='_cache', verbose=0)

def init_vit(model_name):
    assert model_name == 'vit_l_16'

    backbone = models.vit_l_16(weights='ViT_L_16_Weights.IMAGENET1K_SWAG_E2E_V1')
    for param in backbone.parameters():
        param.requires_grad_(False)

    # backbone.class_token[:] = 0
    backbone.class_token.requires_grad_(True)

    hidden_dim = backbone.heads[0].in_features

    del backbone.heads

    return backbone, hidden_dim

def init_cnn(name: str):
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
def normalize_inputs(inputs, means, stds, normalize=1):
    if normalize == 2:
        inputs = torch.log(inputs)

    return (
        (inputs - means[None, :, None, None].to(inputs.device))
        / stds[None, :, None, None].to(inputs.device)
    )



@torch.no_grad()
def augmentation(inputs):
    inputs = get_transforms()(inputs)
    inputs = inputs + 0.01 * torch.randn_like(inputs)

    return inputs



class ImagenetTransferLearning(nn.Module):
    def __init__(
            self,
            model_name: str = 'resnet50',
            dropout_p: float = 0.25,
            use_compile: bool = True
    ):
        super().__init__()

        get_classifier_f = partial(get_classifier, dropout_p=dropout_p, num_target_classes=1)

        # For saving in the state dict
        self.kwargs = {'model_name': model_name, 'dropout_p': dropout_p}

        if model_name == 'vit_l_16':
            self.vit, num_features = init_vit(model_name)
            self.vit.eval()

            classifier = get_classifier_f(n_features=num_features)
            self.vit.heads = classifier

            self.forward = self.vit_forward

        else:
            self.feature_extractor, num_features = init_cnn(name=model_name)
            self.feature_extractor.eval()

            self.classifier = get_classifier_f(n_features=num_features)

            self.forward = self.cnn_forward

        if use_compile:
            self.forward = torch.compile(model=self.forward, mode='reduce-overhead')


    # @partial(torch.compile, mode='reduce-overhead')
    def cnn_forward(self, x):
        with torch.no_grad():
            representations = self.feature_extractor(x)
        x = self.classifier(representations)
        return x

    # @partial(torch.compile, mode='reduce-overhead')
    def vit_forward(self, x):

        x = self.vit.forward(x)

        return x

    def step(self, inputs, targets, ratio=1):
        logits = self(inputs).flatten()

        loss = binary_cross_entropy_with_logits(
            logits,
            targets,
            pos_weight=torch.as_tensor(ratio)
        )

        if PROFILE:
            global profiler
            profiler.step()

        return logits, loss

    def eval(self):
        if self.kwargs['model_name'] == 'vit_l_16':
            self.vit.heads.eval()
        else:
            self.classifier.eval()

    def train(self):
        if self.kwargs['model_name'] == 'vit_l_16':
            self.vit.heads.train()
        else:
            self.classifier.train()

def get_dataloaders(dataset_root, batch_size, normalize):
    num_workers = min(12, len(os.sched_getaffinity(0)))
    #num_workers = 0
    prefetch_factor, persistent_workers = (
        (2, True) if num_workers > 0 else
        (None, False)
    )
    generators = {}
    for mode in ('val', 'train'):
        generators[mode] = torch.Generator()
        if SEED is not None:
            generators[mode].manual_seed(SEED)


    loaders = tuple(
        MultiEpochsDataLoader(
            dataset=FitsDataset(dataset_root, mode=mode, normalize=normalize),
            batch_size=batch_size,
            num_workers=num_workers,
            prefetch_factor=prefetch_factor,
            persistent_workers=persistent_workers,
            worker_init_fn=seed_worker,
            generator=generators[mode],
            pin_memory=True,
            shuffle=True if mode == 'train' else False,
            drop_last=True if mode == 'train' else False,
        )
        for mode in ('train', 'val')
    )

    return loaders

@cache.cache(ignore=['loader'])
def _compute_statistics(loader, normalize, _):
    if not normalize:
        return torch.asarray([0, 0, 0]), torch.asarray([1, 1, 1])
    means = []
    sums_of_squares = []
    f = torch.log if normalize == 2 else lambda x: x
    for i, (imgs, _) in enumerate(loader):
        imgs = f(imgs)
        means.append(torch.mean(imgs, dim=(0, 2, 3)))
        sums_of_squares.append((imgs**2).sum(dim=(0, 2, 3)))

    mean = torch.stack(means).mean(0)
    sums_of_squares = torch.stack(sums_of_squares).sum(0)
    variance = (sums_of_squares / (len(loader) * imgs.shape[0] * imgs.shape[2] * imgs.shape[3])) - (mean ** 2)
    return mean, torch.sqrt(variance)


def compute_statistics(loader, normalize):
    return _compute_statistics(loader, normalize, loader.dataset.sources)



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
def prepare_data(data: torch.Tensor, labels: torch.Tensor, resize: int, normalize: int, device: torch.device, mean:torch.Tensor, std: torch.Tensor):

    data, labels = (
        data.to(device, non_blocking=True, memory_format=torch.channels_last),
        labels.to(device, non_blocking=True, dtype=data.dtype)
    )

    if resize:
      data = interpolate(data, size=resize, mode='bilinear', align_corners=False)

    data = normalize_inputs(data, mean, std, normalize)
    # data = data.expand(-1, 3)

    return data, labels


def main(dataset_root: str, model_name: str, lr: float, resize: int, normalize: int, dropout_p: float, batch_size: int, use_compile: bool):
    torch.set_float32_matmul_precision('high')
    torch.backends.cudnn.benchmark = True

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
        str(Path.cwd() / 'grid_search'),
        # kwargs
        model=model_name,
        lr=lr,
        normalize=normalize,
        dropout_p=dropout_p, 
        use_compile=use_compile
    )

    writer = get_tensorboard_logger(logging_dir)

    device = torch.device('cuda')

    model: nn.Module = ImagenetTransferLearning(model_name=model_name, dropout_p=dropout_p, use_compile=use_compile)

    # noinspection PyArgumentList
    model.to(device=device, memory_format=torch.channels_last)

    optimizer = get_optimizer([param for param in model.parameters() if param.requires_grad], lr=lr)

    train_dataloader, val_dataloader = get_dataloaders(dataset_root, batch_size, normalize)

    mean, std = compute_statistics(train_dataloader, normalize)

    logging_interval = 10

    train_step_f, val_step_f = (
        partial(
            step_f,
            prepare_data_f=partial(prepare_data, resize=resize, normalize=normalize, device=device, mean=mean, std=std),
            metrics_logger=partial(log_metrics, write_metrics_f=partial(write_metrics, writer=writer))
        )
        for step_f in (train_step, val_step)
    )

    train_step_f = partial(train_step_f, train_dataloader=train_dataloader, optimizer=optimizer, logging_interval=logging_interval)
    val_step_f = partial(val_step_f, val_dataloader=val_dataloader)

    checkpoint_saver = partial(save_checkpoint, logging_dir=logging_dir, model=model, optimizer=optimizer, normalize=normalize, batch_size=batch_size)

    best_val_loss = torch.inf
    global_step = 0  # make it a tensor so we can do in-place edits

    best_results = {}

    n_epochs = 250
    for epoch in range(n_epochs):

        global_step = train_step_f(global_step=global_step, model=model)
        val_loss, logits, targets = val_step_f(global_step=global_step, model=model)
        if val_loss < best_val_loss:
            best_results['logits'] = logits.clone()
            best_results['targets'] = targets.clone()
            checkpoint_saver(global_step=global_step)
            best_val_loss = val_loss

        with torch.no_grad():
            log_metrics(loss=best_val_loss, logits=best_results['logits'], targets=best_results['targets'], global_step=global_step, log_suffix='validation_best', write_metrics_f=partial(write_metrics, writer=writer))

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
        pr_curve=(probs.to(torch.float32), targets.to(torch.float32)),
        suffix=log_suffix
    )

    write_metrics_f(metrics=metrics, global_step=global_step)

@torch.no_grad()
def val_step(model, val_dataloader, global_step, metrics_logger, prepare_data_f):
    val_losses, val_logits, val_targets = [], [], []

    model.eval()
    for i, (data, labels) in tqdm(enumerate(val_dataloader), desc='Validation', total=len(val_dataloader)):
        # print("validation start")

        data, labels = prepare_data_f(data, labels)
        with torch.autocast('cuda', dtype=torch.bfloat16):
            logits = model(data).flatten()
            loss = binary_cross_entropy_with_logits(
            logits,
            labels,
            pos_weight=torch.as_tensor(val_dataloader.dataset.label_ratio)
            )
        val_losses.append(loss)
        val_logits.append(logits.clone())
        val_targets.append(labels)

    losses, logits, targets = map(torch.concatenate, (val_losses, val_logits, val_targets))

    mean_loss = losses.mean()
    metrics_logger(loss=mean_loss, logits=logits, targets=targets, global_step=global_step, log_suffix='validation')

    return mean_loss, logits, targets


def train_step(model, optimizer, train_dataloader, prepare_data_f, global_step, logging_interval, metrics_logger):
    # print("training")
    model.train()

    for i, (data, labels) in tqdm(enumerate(train_dataloader), desc='Training', total=len(train_dataloader)):
        global_step += 1

        data, labels = prepare_data_f(data, labels)

        data = augmentation(data)

        optimizer.zero_grad(set_to_none=True)
        with torch.autocast('cuda', dtype=torch.bfloat16):
            logits, loss = model.step(data, labels, ratio=train_dataloader.dataset.label_ratio)
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
    
    def __hash__(self):
        return hash(self.dataset) + 10000


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

    ckpt_dict = torch.load(ckpt_path, weights_only=False)

    # ugh, this is so ugly, something something hindsight something something 20-20
    # FIXME: probably should do a pattern match, but this works for now
    kwargs = str(Path(ckpt_path).parent).split('/')[-1].split('__')

    # strip 'model_' from the name
    model_name = kwargs[1][6:]
    lr = float(kwargs[2].split('_')[-1])
    normalize = int(kwargs[3].split('_')[-1])
    dropout_p = float(kwargs[4].split('_')[-1])

    model = ckpt_dict['model'](model_name=model_name, dropout_p=dropout_p)
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
    parser.add_argument('--batch_size', type=int, help='Batch size', default=32)
    parser.add_argument('--model_name', type=str, help='The model to use.', default='resnet50',
                        choices=['resnet50', 'resnet152', 'resnext50_32x4d', 'resnext101_64x4d', 'efficientnet_v2_l', 'vit_l_16'])
    parser.add_argument('--normalize', type=int, help='Whether to do normalization', default=0, choices=[0, 1, 2])
    parser.add_argument('--dropout_p', type=float, help='Dropout probability', default=0.25)
    parser.add_argument('--resize', type=int, default=0, help="size to resize to. Will be set to 512 for ViT.")
    parser.add_argument('--use_compile', action='store_true')
    parser.add_argument('--profile', action='store_true', help="[DISABLED] profile the training and validation loop")
    parser.add_argument('-d', '--deterministic', action='store_true', help="use deterministic training", default=False)

    return parser.parse_args()


def sanity_check_args(parsed_args):
    assert parsed_args.lr >= 0
    assert parsed_args.batch_size >= 0
    assert 0 <= parsed_args.dropout_p <= 1
    assert parsed_args.resize >= 0
    # ViT always needs the input size to be 512x512
    if parsed_args.model_name == 'vit_l_16' and parsed_args.resize != 512:
        print("Setting resize to 512 since vit_16_l is being used")
        parsed_args.resize = 512

    return parsed_args

def set_seed(seed):
    np.random.seed(seed)
    torch.manual_seed(seed)
    random.seed(seed)
    # torch.use_deterministic_algorithms(True)
    # torch.utils.deterministic.fill_uninitialized_memory = False

def seed_worker(worker_id):
    worker_seed = torch.initial_seed() % 2**32
    np.random.seed(worker_id)
    random.seed(worker_seed)

if __name__ == '__main__':
    args = get_argparser()
    parsed_args = sanity_check_args(args)
    kwargs = vars(args)
    print(kwargs)

    if kwargs['profile']:
        PROFILE = True
    if kwargs['deterministic']:
        SEED = 42
        set_seed(SEED)
    del kwargs['deterministic']
    del kwargs['profile']

    main(**kwargs)
