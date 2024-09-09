import argparse
import os
import subprocess
# os.environ['CUDA_LAUNCH_BLOCKING']="1"
from functools import partial, lru_cache
from pathlib import Path

import torch
import torcheval.metrics.functional as tef
import torch.distributed as dist
from torch.nn.parallel import DistributedDataParallel as DDP
import torch.multiprocessing as mp
from torch import nn, binary_cross_entropy_with_logits
from torch.nn.functional import interpolate
from torch.utils.data import SequentialSampler, Sampler, DistributedSampler
from torch.utils.tensorboard import SummaryWriter
from torchvision import models
from torchvision.transforms import v2
from tqdm import tqdm

from cortex.pre_processing_for_ml import FitsDataset

PROFILE = False

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
    for param in feature_extractor.parameters():
        param.requires_grad = False
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

    # def compile_forward(self):
    #     self.forward  = torch.compile(model=self.forward, mode='reduce-overhead')


    # @partial(torch.compile, mode='reduce-overhead')
    def cnn_forward(self, x):

        with torch.no_grad():
            representations = self.feature_extractor(x)
        classified = self.classifier(representations)
        return classified

    # @partial(torch.compile, mode='reduce-overhead')
    def vit_forward(self, x):

        x = self.vit.forward(x)

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

    def eval(self):
        if self.kwargs['model_name'] == 'vit_l_16':
            self.vit.heads.eval()
        else:
            self.classifier.eval()

    def train(self, mode=True):
        if self.kwargs['model_name'] == 'vit_l_16':
            self.vit.heads.train(mode)
        else:
            self.classifier.train(mode)

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
def prepare_data(data: torch.Tensor, labels: torch.Tensor, resize: int, normalize: int, device: torch.device):

    data, labels = (
        data.to(device, non_blocking=True, memory_format=torch.channels_last),
        labels.to(device, non_blocking=True, dtype=data.dtype)
    )

    if resize:
      data = interpolate(data, size=resize, mode='bilinear', align_corners=False)

    data = normalize_inputs(data, mode=normalize)

    return data, labels

def main(rank: int, local_rank: int, dataset_root: str, model_name: str, lr: float, resize: int, normalize: int, dropout_p: float, batch_size: int, use_compile: bool, world_size: int):
    torch.set_float32_matmul_precision('high')
    torch.cuda.set_device(local_rank)
    dist.init_process_group("NCCL", rank=rank, world_size=world_size)
    
    

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
        dropout_p=dropout_p, 
        use_compile=use_compile,
        gpus=world_size
    )
    writer = get_tensorboard_logger(logging_dir) if not rank else None

    # device = torch.device('cuda')
    model: nn.Module = ImagenetTransferLearning(model_name=model_name, dropout_p=dropout_p, use_compile=False)
    
    # noinspection PyArgumentList
    model.to(device=local_rank, memory_format=torch.channels_last)
    model = DDP(model, device_ids=[local_rank])


    if use_compile:
        model.forward = torch.compile(model=model.forward, mode='reduce-overhead')

    if model_name == 'vit_l_16':
        params = [param for param in model.parameters() if param.requires_grad]
        optimizer = get_optimizer(params, lr=lr)
    else:
        params = [param for param in model.parameters() if param.requires_grad]
        optimizer = get_optimizer(params, lr=lr)

    train_dataloader, val_dataloader = get_dataloaders(dataset_root, batch_size)

    logging_interval = 10

    train_step_f, val_step_f = (
        partial(
            step_f,
            model=model,
            prepare_data_f=partial(prepare_data, resize=resize, normalize=normalize, device=local_rank),
            metrics_logger=partial(log_metrics, write_metrics_f=partial(write_metrics, writer=writer) if writer is not None else None)
        )
        for step_f in (train_step, val_step)
    )

    train_step_f = partial(train_step_f, train_dataloader=train_dataloader, optimizer=optimizer, logging_interval=logging_interval, rank=rank, local_rank=local_rank)
    val_step_f = partial(val_step_f, val_dataloader=val_dataloader, rank=rank)
    if not rank:
        checkpoint_saver = partial(save_checkpoint, logging_dir=logging_dir, model=model, optimizer=optimizer, normalize=normalize, batch_size=batch_size)

    best_val_loss = torch.inf
    global_step = 0  # make it a tensor so we can do in-place edits
    
    n_epochs = 250
    for epoch in range(n_epochs):
        dist.barrier()
        global_step = train_step_f(global_step=global_step)
        if not rank:
            val_loss = val_step_f(global_step=global_step)
            if val_loss < best_val_loss:
                checkpoint_saver(global_step=global_step)
                best_val_loss = val_loss

    if PROFILE:
        profiler.stop()

    writer.flush()
    writer.close()
    dist.destroy_process_group()

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
def val_step(model, val_dataloader, global_step, metrics_logger, prepare_data_f, rank):
    val_losses, val_logits, val_targets = [], [], []

    model.eval()
    for i, (data, labels) in tqdm(enumerate(val_dataloader), desc='Validation', total=len(val_dataloader), disable=rank):
        # print("validation start")

        data, labels = prepare_data_f(data, labels)
        with torch.autocast('cuda', dtype=torch.bfloat16):
            logits = model(data).flatten()
            loss = binary_cross_entropy_with_logits(
            logits,
            labels,
            pos_weight=torch.as_tensor(2.698717948717949)
            )

        val_losses.append(loss)
        val_logits.append(logits.clone())
        val_targets.append(labels)

    losses, logits, targets = map(torch.concatenate, (val_losses, val_logits, val_targets))

    mean_loss = losses.mean()
    metrics_logger(loss=mean_loss, logits=logits, targets=targets, global_step=global_step, log_suffix='validation')

    return mean_loss



def train_step(model, optimizer, train_dataloader, prepare_data_f, global_step, logging_interval, metrics_logger, local_rank, rank):
    model.train()
    train_losses, train_logits, train_targets = [], [], []
    for i, (data, labels) in tqdm(enumerate(train_dataloader), desc='Training', total=len(train_dataloader), disable=rank):
        global_step += 1

        data, labels = prepare_data_f(data, labels, device=local_rank)

        data = augmentation(data)

        optimizer.zero_grad(set_to_none=True)
        with torch.autocast('cuda', dtype=torch.bfloat16, enabled=True):
            logits = model(data).flatten()
            loss = binary_cross_entropy_with_logits(
            logits,
            labels,
            pos_weight=torch.as_tensor(2.698717948717949)
            )
            mean_loss = loss.mean()
        
        mean_loss.backward()
        dist.barrier()
        # param_grad = next(model.module.classifier[2].parameters()).grad
        # grads = [torch.zeros_like(param_grad) for _ in range(dist.get_world_size())]
        # dist.all_gather(grads, param_grad)

        # Check if all gathered gradients are the same
        # synced = all(torch.equal(grads[0], g) for g in grads)
        # if synced:
        #     print(f"Process {rank}: Gradients are synchronized across processes.")
        # else:
        #     print(f"Process {rank}: Gradients are NOT synchronized!")
        optimizer.step()
        train_losses.append(loss.detach().clone())
        train_logits.append(logits.clone())
        train_targets.append(labels)
        # dist.barrier()
        # for grad in grads:
        #     print(grad)
        # dist.barrier()
        # dist.destroy_process_group()
        # exit()

        if i % logging_interval == 0 and not rank:
            with torch.no_grad():
                losses, logits, targets = map(torch.concatenate, (train_losses, train_logits, train_targets))
                metrics_logger(loss=losses.mean(), logits=logits, targets=targets, global_step=global_step, log_suffix='training')
                train_losses, train_logits, train_targets = [], [], []


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
    parser.add_argument('--dataset_root', type=Path, default=f'/scratch-shared/CORTEX/public.spider.surfsara.nl/project/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/cnn_data')
    parser.add_argument('--lr', type=float, help='Learning rate for the model.', default=1e-4)
    parser.add_argument('--batch_size', type=int, help='Batch size', default=32)
    parser.add_argument('--model_name', type=str, help='The model to use.', default='resnet50',
                        choices=['resnet50', 'resnet152', 'resnext50_32x4d', 'resnext101_64x4d', 'efficientnet_v2_l', 'vit_l_16'])
    parser.add_argument('--normalize', type=int, help='Whether to do normalization', default=0, choices=[0, 1, 2])
    parser.add_argument('--dropout_p', type=float, help='Dropout probability', default=0.25)
    parser.add_argument('--resize', type=int, default=0, help="size to resize to. Will be set to 512 for ViT.")
    parser.add_argument('--use_compile', action='store_true', default=True)
    parser.add_argument('--no_compile', dest='use_compile', action='store_false')

    parser.add_argument('--profile', action='store_true', help="[DISABLED] profile the training and validation loop")

    # parser.add_argument('-g', '--gpus', type=int, default=1, help='number of gpus')
    # parser.add_argument('-n', '--nodes', type=int, default=1, help='number of nodes')

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


if __name__ == '__main__':

    args = get_argparser()
    parsed_args = sanity_check_args(args)
    nodes = int(os.environ.get('SLURM_JOB_NUM_NODES', 1))
    parsed_args.world_size = int(os.environ.get('SLURM_NTASKS', 1))
    rank = int(os.environ.get('SLURM_PROCID', 0))
    gpus_per_node = torch.cuda.device_count()
    local_rank = rank - gpus_per_node * (rank //gpus_per_node)
    parsed_args.rank = rank
    parsed_args.local_rank = local_rank
    kwargs = vars(args)
    print(kwargs)

    if kwargs['profile']:
        PROFILE = True
    del kwargs['profile']

    main_args = (
        kwargs['dataset_root'],
        kwargs['model_name'],
        kwargs['lr'],
        kwargs['resize'],
        kwargs['normalize'],
        kwargs['dropout_p'],
        kwargs['batch_size'],
        kwargs['use_compile'],
        kwargs['world_size']
    )
    if "MASTER_ADDR" not in os.environ:
        hostnode = subprocess.getoutput("scontrol show hostnames $SLURM_NODELIST | head -n 1")
        if 'error' in hostnode:
            hostnode = 'localhost'
        os.environ['MASTER_ADDR'] = hostnode
    if "MASTER_PORT" not in os.environ:
        os.environ["MASTER_PORT"] = "29400"
    

    main(**kwargs)
    # mp.spawn(main, args=main_args, nprocs=args.world_size, join=True)
