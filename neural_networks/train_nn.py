import argparse
import os
from functools import partial, lru_cache
from pathlib import Path
import warnings

import torch
import torcheval.metrics.functional as tef
from torch import nn, binary_cross_entropy_with_logits
from torch.nn.functional import interpolate
from torch.utils.data import SequentialSampler, Sampler
from torch.utils.tensorboard import SummaryWriter
from torchvision import models
from torchvision.transforms import v2
from tqdm import tqdm
import torchvision.transforms.functional as TF
from torchvision.transforms.functional import InterpolationMode

import numpy as np
import random

from pre_processing_for_ml import FitsDataset
from dino_model import DINOV2FeatureExtractor

PROFILE = False
SEED = None


def init_vit(model_name):
    assert model_name == "vit_l_16"

    backbone = models.vit_l_16(weights="ViT_L_16_Weights.IMAGENET1K_SWAG_E2E_V1")
    for param in backbone.parameters():
        param.requires_grad_(False)

    # backbone.class_token[:] = 0
    backbone.class_token.requires_grad_(True)

    hidden_dim = backbone.heads[0].in_features

    del backbone.heads
    return backbone, hidden_dim


def init_dino(model_name, get_classifier_f, use_lora, rank=16, alpha=16):

    backbone = torch.hub.load("facebookresearch/dinov2", model_name)
    hidden_dim = backbone.cls_token.shape[-1]
    classifier = get_classifier_f(n_features=hidden_dim)

    dino_lora = DINOV2FeatureExtractor(
        encoder=backbone,
        decoder=classifier,
        r=rank,
        use_lora=use_lora,
        alpha=alpha,
    )

    return dino_lora, hidden_dim


def init_cnn(name: str, lift="stack"):
    # use partial to prevent loading all models at once
    model_map = {
        "resnet50": partial(models.resnet50, weights="DEFAULT"),
        "resnet152": partial(models.resnet152, weights="DEFAULT"),
        "resnext50_32x4d": partial(models.resnext50_32x4d, weights="DEFAULT"),
        "resnext101_64x4d": partial(models.resnext101_64x4d, weights="DEFAULT"),
        "efficientnet_v2_l": partial(models.efficientnet_v2_l, weights="DEFAULT"),
    }

    backbone = model_map[name]()

    feature_extractor = nn.Sequential(*list(backbone.children())[:-1])
    for param in feature_extractor.parameters():
        param.requires_grad_(False)

    if lift == "reinit_first":
        if name in ("resnet50", "resnet152", "resnext50_32x4d", "resnext101_64x4d"):
            conv = feature_extractor[0]
            new_conv = init_first_conv(conv)
            feature_extractor[0] = new_conv
            del conv
        elif name == "efficientnet_v2_l":
            conv = feature_extractor[0][0][0]
            new_conv = init_first_conv(conv)
            feature_extractor[0][0][0] = new_conv
            del conv

    num_out_features = (
        backbone.fc
        if name in ("resnet50", "resnet152", "resnext50_32x4d", "resnext101_64x4d")
        else backbone.classifier[-1]  # efficientnet
    ).in_features
    return feature_extractor, num_out_features


def init_first_conv(conv):
    kernel_size = conv.kernel_size
    stride = conv.stride
    padding = conv.padding
    bias = conv.bias
    out_channels = conv.out_channels
    return nn.Conv2d(
        in_channels=1,
        out_channels=out_channels,
        kernel_size=kernel_size,
        stride=stride,
        padding=padding,
        bias=bias,
    )


def get_classifier(dropout_p: float, n_features: int, num_target_classes: int):
    assert 0 <= dropout_p <= 1
    print(n_features)

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
        inputs = torch.log(inputs + 1e-10)
    return (inputs - means[None, :, None, None].to(inputs.device)) / stds[
        None, :, None, None
    ].to(inputs.device)


@torch.no_grad()
def augmentation(inputs):
    inputs = get_transforms()(inputs)
    inputs = inputs + 0.01 * torch.randn_like(inputs)

    return inputs


class ImagenetTransferLearning(nn.Module):
    def __init__(
        self,
        model_name: str = "resnet50",
        dropout_p: float = 0.25,
        use_compile: bool = True,
        lift: str = "stack",
        use_lora: bool = False,
        alpha: float = 16.0,
        rank: int = 16,
        pos_embed: bool = False,
    ):
        super().__init__()

        get_classifier_f = partial(
            get_classifier, dropout_p=dropout_p, num_target_classes=1
        )

        # For saving in the state dict
        self.kwargs = {
            "model_name": model_name,
            "dropout_p": dropout_p,
            "use_compile": use_compile,
            "lift": lift,
            "use_lora": use_lora,
            "rank": rank,
            "alpha": alpha,
            "pos_embed": pos_embed,
        }

        if lift == "stack":
            self.lift = partial(torch.repeat_interleave, repeats=3, dim=1)
        elif lift == "conv":
            self.lift = nn.Conv2d(1, 3, 1)
        elif lift == "reinit_first":
            self.lift = nn.Identity()

        if model_name == "vit_l_16":
            self.vit, num_features = init_vit(model_name)
            self.vit.eval()

            classifier = get_classifier_f(n_features=num_features)
            self.vit.heads = classifier

            self.forward = self.vit_forward

        elif "dinov2" in model_name:
            self.dino, num_features = init_dino(
                model_name, get_classifier_f, use_lora=use_lora, alpha=alpha, rank=rank
            )
            # self.classifier = get_classifier_f(n_features=num_features)
            if "zeros" in self.kwargs["pos_embed"]:
                self.dino.encoder.pos_embed[:, 1:, :] = torch.zeros_like(
                    self.dino.encoder.pos_embed[:, 1:, :]
                )
            self.dino.encoder.pos_embed.requires_grad = (
                True if "fine-tune" in self.kwargs["pos_embed"] else False
            )
            self.forward = self.dino_forward

        else:
            self.feature_extractor, num_features = init_cnn(name=model_name, lift=lift)
            self.feature_extractor.eval()

            self.classifier = get_classifier_f(n_features=num_features)

            self.forward = self.cnn_forward

        if use_compile:
            self.forward = torch.compile(model=self.forward, mode="reduce-overhead")

    # @partial(torch.compile, mode='reduce-overhead')
    def cnn_forward(self, x):
        x = self.lift(x)
        representations = self.feature_extractor(x)
        x = self.classifier(representations)
        return x

    # @partial(torch.compile, mode='reduce-overhead')
    def vit_forward(self, x):
        x = self.lift(x)
        x = self.vit.forward(x)

        return x

    def dino_forward(self, x):
        x = self.lift(x)
        return self.dino(x)

    def step(self, inputs, targets, ratio=1):
        logits = self(inputs).flatten()

        loss = binary_cross_entropy_with_logits(
            logits, targets, pos_weight=torch.as_tensor(ratio)
        )

        if PROFILE:
            global profiler
            profiler.step()

        return logits, loss

    def eval(self):
        if self.kwargs["model_name"] == "vit_l_16":
            self.vit.heads.eval()
        elif "dinov2" in self.kwargs["model_name"]:
            if self.kwargs["use_lora"]:
                self.dino.eval()
            else:
                self.dino.decoder.eval()
        else:
            self.classifier.eval()

    def train(self):
        if self.kwargs["model_name"] == "vit_l_16":
            self.vit.heads.train()
        elif "dinov2" in self.kwargs["model_name"]:
            if self.kwargs["use_lora"]:
                self.dino.train()
            else:
                self.dino.decoder.train()
            # Finetune learnable pos_embedding

        else:
            self.classifier.train()


def get_dataloaders(dataset_root, batch_size):
    num_workers = min(18, len(os.sched_getaffinity(0)))

    prefetch_factor, persistent_workers = (
        (2, True) if num_workers > 0 else (None, False)
    )
    generators = {}
    for mode in ("val", "train"):
        generators[mode] = torch.Generator()
        if SEED is not None:
            generators[mode].manual_seed(SEED)

    loaders = tuple(
        MultiEpochsDataLoader(
            dataset=FitsDataset(dataset_root, mode=mode),
            batch_size=batch_size,
            num_workers=num_workers,
            prefetch_factor=prefetch_factor,
            persistent_workers=persistent_workers,
            worker_init_fn=seed_worker,
            generator=generators[mode],
            pin_memory=True,
            shuffle=True if mode == "train" else False,
            drop_last=True if mode == "train" else False,
        )
        for mode in ("train", "val")
    )

    return loaders


def get_logging_dir(logging_root: str, /, **kwargs):
    # As version string, prefer $SLURM_ARRAY_JOB_ID, then $SLURM_JOB_ID, then 0.
    version = int(os.getenv("SLURM_ARRAY_JOB_ID", os.getenv("SLURM_JOB_ID", 0)))
    version_appendix = int(os.getenv("SLURM_ARRAY_TASK_ID", 0))
    while True:
        version_dir = "__".join(
            (
                f"version_{version}_{version_appendix}",
                *(f"{k}_{v}" for k, v in kwargs.items()),
            )
        )

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
                writer.add_pr_curve,
                labels=labels,
                predictions=probs,
            )
        else:
            writer_fn = partial(writer.add_scalar, scalar_value=value)

        writer_fn(tag=metric_name, global_step=global_step)


def get_optimizer(parameters: list[torch.Tensor], lr: float):
    return torch.optim.AdamW(parameters, lr=lr)


def merge_metrics(suffix, **kwargs):
    return {f"{k}/{suffix}": v for k, v in kwargs.items()}


@torch.no_grad()
def prepare_data(
    data: torch.Tensor,
    labels: torch.Tensor,
    resize: int,
    normalize: int,
    device: torch.device,
    mean: torch.Tensor,
    std: torch.Tensor,
):

    data, labels = (
        data.to(device, non_blocking=True, memory_format=torch.channels_last),
        labels.to(device, non_blocking=True, dtype=data.dtype),
    )

    if resize:
        data = interpolate(data, size=resize, mode="bilinear", align_corners=False)

    data = normalize_inputs(data, mean, std, normalize)

    return data, labels


def main(
    dataset_root: str,
    model_name: str,
    lr: float,
    resize: int,
    normalize: int,
    dropout_p: float,
    batch_size: int,
    use_compile: bool,
    label_smoothing: float,
    stochastic_smoothing: bool,
    lift: str,
    use_lora: bool,
    rank: int = 16,
    alpha: float = 16,
    log_path: Path = "runs",
    epochs: int = 120,
    pos_embed: str = "pre-trained",
):
    torch.set_float32_matmul_precision("high")
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
        str(Path.cwd() / log_path),
        # kwargs
        model=model_name,
        lr=lr,
        normalize=normalize,
        dropout_p=dropout_p,
        use_compile=use_compile,
        label_smoothing=label_smoothing,
        stochastic_smoothing=stochastic_smoothing,
        use_lora=use_lora,
        resize=resize,
        rank=rank,
        alpha=alpha,
        lift=lift,
        pos_embed=pos_embed,
    )

    writer = get_tensorboard_logger(logging_dir)

    device = torch.device("cuda")

    model: nn.Module = ImagenetTransferLearning(
        model_name=model_name,
        dropout_p=dropout_p,
        use_compile=use_compile,
        lift=lift,
        use_lora=use_lora,
        alpha=alpha,
        rank=rank,
        pos_embed=pos_embed,
    )

    # noinspection PyArgumentList
    model.to(device=device, memory_format=torch.channels_last)

    optimizer = get_optimizer(
        [param for param in model.parameters() if param.requires_grad], lr=lr
    )

    train_dataloader, val_dataloader = get_dataloaders(dataset_root, batch_size)

    mean, std = train_dataloader.dataset.compute_statistics(normalize)

    logging_interval = 10

    train_step_f, val_step_f = (
        partial(
            step_f,
            prepare_data_f=partial(
                prepare_data,
                resize=resize,
                normalize=normalize,
                device=device,
                mean=mean,
                std=std,
            ),
            metrics_logger=partial(
                log_metrics, write_metrics_f=partial(write_metrics, writer=writer)
            ),
        )
        for step_f in (train_step, val_step)
    )

    train_step_f = partial(
        train_step_f,
        train_dataloader=train_dataloader,
        optimizer=optimizer,
        logging_interval=logging_interval,
        smoothing_fn=partial(
            label_smoother,
            stochastic=stochastic_smoothing,
            smoothing_factor=label_smoothing,
        ),
    )
    val_step_f = partial(val_step_f, val_dataloader=val_dataloader)

    checkpoint_saver = partial(
        save_checkpoint,
        logging_dir=logging_dir,
        model=model,
        optimizer=optimizer,
        args={
            "normalize": normalize,
            "batch_size": batch_size,
            "use_compile": use_compile,
            "label_smoothing": label_smoothing,
            "stochastic_smoothing": stochastic_smoothing,
            "lift": lift,
            "use_lora": use_lora,
            "rank": rank,
            "alpha": alpha,
            "resize": resize,
            "lr": lr,
            "dropout_p": dropout_p,
            "model_name": model_name,
            "dataset_mean": mean,
            "dataset_std": std,
        },
        model_args={
            "model_name": model_name,
            "use_compile": use_compile,
            "lift": lift,
            "use_lora": use_lora,
            "rank": rank,
            "alpha": alpha,
            "dropout_p": dropout_p,
            "pos_embed": pos_embed,
        },
    )

    best_val_loss = torch.inf
    global_step = 0  # make it a tensor so we can do in-place edits

    best_results = {}

    n_epochs = epochs
    for epoch in range(n_epochs):

        global_step = train_step_f(global_step=global_step, model=model)
        val_loss, logits, targets = val_step_f(global_step=global_step, model=model)
        if val_loss < best_val_loss:
            best_results["logits"] = logits.clone()
            best_results["targets"] = targets.clone()
            checkpoint_saver(global_step=global_step)
            best_val_loss = val_loss

        with torch.no_grad():
            log_metrics(
                loss=best_val_loss,
                logits=best_results["logits"],
                targets=best_results["targets"],
                global_step=global_step,
                log_suffix="validation_best",
                write_metrics_f=partial(write_metrics, writer=writer),
            )

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
        suffix=log_suffix,
    )

    write_metrics_f(metrics=metrics, global_step=global_step)


@torch.no_grad()
def val_step(model, val_dataloader, global_step, metrics_logger, prepare_data_f):
    val_losses, val_logits, val_targets = [], [], []

    model.eval()
    for i, (data, labels) in tqdm(
        enumerate(val_dataloader), desc="Validation", total=len(val_dataloader)
    ):
        # print("validation start")

        data, labels = prepare_data_f(data, labels)
        with torch.autocast("cuda", dtype=torch.bfloat16):
            logits = model(data).flatten()
            loss = binary_cross_entropy_with_logits(
                logits,
                labels,
                pos_weight=torch.as_tensor(val_dataloader.dataset.label_ratio),
            )
        val_losses.append(loss)
        val_logits.append(logits.clone())
        val_targets.append(labels)

    losses, logits, targets = map(
        torch.concatenate, (val_losses, val_logits, val_targets)
    )

    mean_loss = losses.mean()
    metrics_logger(
        loss=mean_loss,
        logits=logits,
        targets=targets,
        global_step=global_step,
        log_suffix="validation",
    )

    return mean_loss, logits, targets


def label_smoother(
    labels: torch.tensor, smoothing_factor: float = 0.1, stochastic: bool = True
):
    smoothing_factor = smoothing_factor - (
        torch.rand_like(labels) * smoothing_factor * stochastic
    )
    smoothed_label = (1 - smoothing_factor) * labels + 0.5 * smoothing_factor
    return smoothed_label


def train_step(
    model,
    optimizer,
    train_dataloader,
    prepare_data_f,
    global_step,
    logging_interval,
    metrics_logger,
    smoothing_fn,
):
    # print("training")
    model.train()

    for i, (data, labels) in tqdm(
        enumerate(train_dataloader), desc="Training", total=len(train_dataloader)
    ):
        global_step += 1
        data, labels = prepare_data_f(data, labels)
        smoothed_label = smoothing_fn(labels)
        data = augmentation(data)

        optimizer.zero_grad(set_to_none=True)
        with torch.autocast("cuda", dtype=torch.bfloat16):
            logits, loss = model.step(
                data, smoothed_label, ratio=train_dataloader.dataset.label_ratio
            )
            mean_loss = loss.mean()

        mean_loss.backward()
        optimizer.step()
        if i % logging_interval == 0:
            with torch.no_grad():
                metrics_logger(
                    loss=mean_loss.detach(),
                    logits=logits.detach(),
                    targets=labels,
                    global_step=global_step,
                    log_suffix="training",
                )

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
        return (
            len(self.sampler)
            if self.batch_sampler is None
            else len(self.batch_sampler.sampler)
        )

    def __iter__(self):
        for i in range(len(self)):
            yield next(self.iterator)

    def __hash__(self):
        return hash(self.dataset) + 10000


class _RepeatSampler(object):
    """Sampler that repeats forever.

    Args:
        sampler (Sampler)
    """

    def __init__(self, sampler):
        self.sampler = sampler

    def __iter__(self):
        while True:
            yield from iter(self.sampler)


class Rotate90Transform:
    def __init__(self, angles=[0, 90, 180, 270]):
        self.angles = angles

    def __call__(self, x):
        angle = np.random.choice(self.angles)
        return v2.functional.rotate(x, int(angle), InterpolationMode.BILINEAR)


@lru_cache(maxsize=1)
def get_transforms():

    return v2.Compose(
        [
            (Rotate90Transform()),
            v2.RandomHorizontalFlip(p=0.5),
        ]
    )


def save_checkpoint(logging_dir, model, optimizer, global_step, **kwargs):
    old_checkpoints = Path(logging_dir).glob("*.pth")
    for old_checkpoint in old_checkpoints:
        Path.unlink(old_checkpoint)

    torch.save(
        {
            "model": type(model),
            "model_state_dict": model.state_dict(),
            "optimizer": type(optimizer),
            "optimizer_state_dict": optimizer.state_dict(),
            "step": global_step,
            **kwargs,
        },
        f=(logging_dir + f"/ckpt_step={global_step}.pth"),
    )


def load_checkpoint(ckpt_path, device="cuda"):
    if os.path.isfile(ckpt_path):
        ckpt_dict = torch.load(ckpt_path, weights_only=False)
    else:
        files = os.listdir(ckpt_path)
        possible_checkpoints = list(filter(lambda x: x.endswith(".pth"), files))
        if len(possible_checkpoints) != 1:
            raise ValueError(
                f"Too many checkpoint files in the given checkpoint directory. Please specify the model you want to load directly."
            )
        ckpt_path = f"{ckpt_path}/{possible_checkpoints[0]}"
        ckpt_dict = torch.load(ckpt_path, weights_only=False)

    # strip 'model_' from the name
    model_name = ckpt_dict["args"]["model_name"]
    if "model_args" in ckpt_dict["args"]:
        model = ckpt_dict["model"](**ckpt_dict["model_args"]).to(device)
    else:
        dropout_p = ckpt_dict["args"]["dropout_p"]
        use_lora = ckpt_dict["args"]["use_lora"]
        rank = ckpt_dict["args"]["rank"]
        alpha = ckpt_dict["args"]["alpha"]
        lift = ckpt_dict["args"]["lift"]
        model_name = ckpt_dict["args"]["model_name"]

        model = ckpt_dict["model"](
            model_name=model_name,
            dropout_p=dropout_p,
            use_lora=use_lora,
            lift=lift,
            alpha=alpha,
            rank=rank,
        ).to(device)
    model.load_state_dict(ckpt_dict["model_state_dict"])
    lr = ckpt_dict["args"]["lr"]
    try:
        # FIXME: add optim class and args to state dict
        optim = ckpt_dict.get("optimizer", torch.optim.AdamW)(
            lr=lr, params=model.classifier.parameters()
        ).load_state_dict(ckpt_dict["optimizer_state_dict"])
    except Exception as e:
        print(f"Could not load optim due to {e}; skipping.")
        optim = None

    return {"model": model, "optim": optim, "args": ckpt_dict["args"]}


def get_argparser():
    """
    Create and return an argument parser for hyperparameter tuning.
    """
    # Create the parser
    parser = argparse.ArgumentParser(
        description="Hyperparameter tuning for a machine learning model."
    )

    # Add arguments
    parser.add_argument("dataset_root", type=Path)
    parser.add_argument(
        "--lr", type=float, help="Learning rate for the model.", default=1e-4
    )
    parser.add_argument("--batch_size", type=int, help="Batch size", default=32)
    parser.add_argument(
        "--model_name",
        type=str,
        help="The model to use.",
        default="resnet50",
        choices=[
            "resnet50",
            "resnet152",
            "resnext50_32x4d",
            "resnext101_64x4d",
            "efficientnet_v2_l",
            "vit_l_16",
            "dinov2_vits14",
            "dinov2_vitb14",
            "dinov2_vitl14",
            "dinov2_vitg14",
            "dinov2_vits14_reg",
            "dinov2_vitb14_reg",
            "dinov2_vitl14_reg",
            "dinov2_vitg14_reg",
        ],
    )
    parser.add_argument(
        "--normalize",
        type=int,
        help="Whether to do normalization",
        default=1,
        choices=[0, 1, 2],
    )
    parser.add_argument(
        "--dropout_p", type=float, help="Dropout probability", default=0.25
    )
    parser.add_argument(
        "--resize",
        type=int,
        default=0,
        help="size to resize to. Will be set to 512 for ViT.",
    )
    parser.add_argument(
        "--epochs",
        type=int,
        default=120,
        help="number of epochs",
    )
    parser.add_argument("--use_compile", action="store_true")
    parser.add_argument(
        "--profile",
        action="store_true",
        help="[DISABLED] profile the training and validation loop",
    )
    parser.add_argument(
        "-d",
        "--deterministic",
        action="store_true",
        help="use deterministic training",
        default=False,
    )

    parser.add_argument(
        "--label_smoothing",
        type=float,
        default=0,
        help="Label smoothing factor",
    )

    parser.add_argument(
        "--stochastic_smoothing",
        action="store_true",
        help="use stochastic smoothing",
    )

    parser.add_argument(
        "--lift",
        type=str,
        default="stack",
        choices=["stack", "conv", "reinit_first"],
        help="How to lift single channel to 3 channels. Stacking stacks the single channel thrice. Conv adds a 1x1 convolution layer before the model. reinit_first re-initialises the first layer if applicable.",
    )

    parser.add_argument(
        "--use_lora",
        action="store_true",
        help="Whether to use LoRA if applicable.",
    )

    parser.add_argument(
        "--pos_embed",
        type=str,
        default="pre-trained",
        choices=["pre-trained", "fine-tune", "zeros", "zeros-fine-tune"],
        help="How to handle positional embeddings",
    )

    parser.add_argument("--rank", type=int, default=16, help="rank of LoRA")

    parser.add_argument(
        "--alpha",
        type=float,
        default=None,
        help="LoRA alpha scaling. Defaults to rank value if not set",
    )

    parser.add_argument("--log_path", type=str, default="runs")

    return parser.parse_args()


def sanity_check_args(parsed_args):
    assert parsed_args.lr >= 0
    assert parsed_args.batch_size >= 0
    assert 0 <= parsed_args.dropout_p <= 1
    assert parsed_args.resize >= 0
    assert 0 <= parsed_args.label_smoothing <= 1
    assert not (
        (parsed_args.model_name == "vit_l_16" or parsed_args.model_name == "dino_v2")
        and parsed_args.lift == "reinit_first"
    )
    # ViT always needs the input size to be 512x512
    if parsed_args.model_name == "vit_l_16" and parsed_args.resize != 512:
        print("Setting resize to 512 since vit_16_l is being used")
        parsed_args.resize = 512
    if "dinov2" in parsed_args.model_name and parsed_args.resize == 0:
        resize = 560
        print(f"\n#######\nSetting resize to {resize} \n######\n")
        parsed_args.resize = resize

    if parsed_args.use_lora and not "dinov2" in parsed_args.model_name:
        warnings.warn(
            "Warning: LoRA is only supported for Dino V2 models. Ignoring setting....\n",
            UserWarning,
        )

    if parsed_args.alpha is None:
        parsed_args.alpha = parsed_args.rank

    assert parsed_args.resize % 14 == 0 or parsed_args.model_name != "dino_v2"

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


if __name__ == "__main__":
    args = get_argparser()
    parsed_args = sanity_check_args(args)
    kwargs = vars(args)
    print(kwargs)

    if kwargs["profile"]:
        PROFILE = True
    if kwargs["deterministic"]:
        SEED = 42
        set_seed(SEED)
    del kwargs["deterministic"]
    del kwargs["profile"]

    main(**kwargs)
