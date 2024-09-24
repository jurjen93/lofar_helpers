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

from pre_processing_for_ml import FitsDataset

from train_nn import *

PROFILE = False


def get_dataloaders(dataset_root, batch_size):
    num_workers = min(12, len(os.sched_getaffinity(0)))
    prefetch_factor, persistent_workers = (
        (2, True) if num_workers > 0 else (None, False)
    )

    loaders = tuple(
        MultiEpochsDataLoader(
            dataset=FitsDataset(dataset_root, mode=mode),
            batch_size=batch_size,
            num_workers=num_workers,
            prefetch_factor=prefetch_factor,
            persistent_workers=persistent_workers,
            pin_memory=True,
            shuffle=True if mode == "train" else False,
            drop_last=True if mode == "train" else False,
        )
        for mode in ("train", "val")
    )

    return loaders


def main(
    rank: int,
    local_rank: int,
    dataset_root: str,
    model_name: str,
    lr: float,
    resize: int,
    normalize: int,
    dropout_p: float,
    batch_size: int,
    use_compile: bool,
    label_smoothing: float,
    world_size: int,
):
    torch.set_float32_matmul_precision("high")
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
        str(Path.cwd() / "runs"),
        # kwargs
        model=model_name,
        lr=lr,
        normalize=normalize,
        dropout_p=dropout_p,
        use_compile=use_compile,
        gpus=world_size,
    )
    writer = get_tensorboard_logger(logging_dir) if not rank else None

    # device = torch.device('cuda')
    model: nn.Module = ImagenetTransferLearning(
        model_name=model_name, dropout_p=dropout_p, use_compile=False
    )

    # noinspection PyArgumentList
    model.to(device=local_rank, memory_format=torch.channels_last)
    model = DDP(model, device_ids=[local_rank])

    if use_compile:
        model.forward = torch.compile(model=model.forward, mode="reduce-overhead")

    if model_name == "vit_l_16":
        params = [param for param in model.parameters() if param.requires_grad]
        optimizer = get_optimizer(params, lr=lr)
    else:
        params = [param for param in model.parameters() if param.requires_grad]
        optimizer = get_optimizer(params, lr=lr)

    train_dataloader, val_dataloader = get_dataloaders(dataset_root, batch_size)

    mean, std = train_dataloader.dataset.compute_statistics(normalize)

    logging_interval = 10

    train_step_f, val_step_f = (
        partial(
            step_f,
            model=model,
            prepare_data_f=partial(
                prepare_data,
                resize=resize,
                normalize=normalize,
                device=local_rank,
                mean=mean,
                mean=std,
            ),
            metrics_logger=partial(
                log_metrics,
                write_metrics_f=(
                    partial(write_metrics, writer=writer)
                    if writer is not None
                    else None
                ),
            ),
        )
        for step_f in (train_step, val_step)
    )

    train_step_f = partial(
        train_step_f,
        train_dataloader=train_dataloader,
        optimizer=optimizer,
        logging_interval=logging_interval,
        label_smoothing=label_smoothing,
        rank=rank,
        local_rank=local_rank,
    )
    val_step_f = partial(val_step_f, val_dataloader=val_dataloader, rank=rank)
    if not rank:
        checkpoint_saver = partial(
            save_checkpoint,
            logging_dir=logging_dir,
            model=model,
            optimizer=optimizer,
            normalize=normalize,
            batch_size=batch_size,
        )

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
def val_step(model, val_dataloader, global_step, metrics_logger, prepare_data_f, rank):
    val_losses, val_logits, val_targets = [], [], []

    model.eval()
    for i, (data, labels) in tqdm(
        enumerate(val_dataloader),
        desc="Validation",
        total=len(val_dataloader),
        disable=rank,
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

    return mean_loss


def train_step(
    model,
    optimizer,
    train_dataloader,
    prepare_data_f,
    global_step,
    logging_interval,
    metrics_logger,
    local_rank,
    rank,
    label_smoothing=0,
):
    model.train()
    train_losses, train_logits, train_targets = [], [], []
    for i, (data, labels) in tqdm(
        enumerate(train_dataloader),
        desc="Training",
        total=len(train_dataloader),
        disable=rank,
    ):
        global_step += 1

        data, labels = prepare_data_f(data, labels, device=local_rank)
        smoothed_label = (1 - label_smoothing) * labels + 0.5 * label_smoothing

        data = augmentation(data)

        optimizer.zero_grad(set_to_none=True)
        with torch.autocast("cuda", dtype=torch.bfloat16):
            logits = model(data).flatten()
            loss = binary_cross_entropy_with_logits(
                logits,
                labels,
                pos_weight=torch.as_tensor(train_dataloader.dataset.label_ratio),
            )
            mean_loss = loss.mean()

        mean_loss.backward()
        dist.barrier()
        optimizer.step()
        train_losses.append(loss.detach().clone())
        train_logits.append(logits.clone())
        train_targets.append(labels)

        if i % logging_interval == 0 and not rank:
            with torch.no_grad():
                losses, logits, targets = map(
                    torch.concatenate, (train_losses, train_logits, train_targets)
                )
                metrics_logger(
                    loss=losses.mean(),
                    logits=logits,
                    targets=targets,
                    global_step=global_step,
                    log_suffix="training",
                )
                train_losses, train_logits, train_targets = [], [], []

    return global_step


def get_argparser():
    """
    Create and return an argument parser for hyperparameter tuning.
    """
    # Create the parser
    parser = argparse.ArgumentParser(
        description="Hyperparameter tuning for a machine learning model."
    )

    # Add arguments
    parser.add_argument(
        "--dataset_root",
        type=Path,
        default=f"/scratch-shared/CORTEX/public.spider.surfsara.nl/project/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/cnn_data",
    )
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
        ],
    )
    parser.add_argument(
        "--normalize",
        type=int,
        help="Whether to do normalization",
        default=0,
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
    parser.add_argument("--use_compile", action="store_true")

    parser.add_argument(
        "--profile",
        action="store_true",
        help="[DISABLED] profile the training and validation loop",
    )

    parser.add_argument(
        "--label_smoothing",
        type=float,
        default=0,
        help="Label smoothing factor",
    )

    # parser.add_argument('-g', '--gpus', type=int, default=1, help='number of gpus')
    # parser.add_argument('-n', '--nodes', type=int, default=1, help='number of nodes')

    return parser.parse_args()


if __name__ == "__main__":

    args = get_argparser()
    parsed_args = sanity_check_args(args)
    nodes = int(os.environ.get("SLURM_JOB_NUM_NODES", 1))
    parsed_args.world_size = int(os.environ.get("SLURM_NTASKS", 1))
    rank = int(os.environ.get("SLURM_PROCID", 0))
    gpus_per_node = torch.cuda.device_count()
    local_rank = rank - gpus_per_node * (rank // gpus_per_node)
    parsed_args.rank = rank
    parsed_args.local_rank = local_rank
    kwargs = vars(args)
    print(kwargs)

    if kwargs["profile"]:
        PROFILE = True
    del kwargs["profile"]

    main_args = (
        kwargs["dataset_root"],
        kwargs["model_name"],
        kwargs["lr"],
        kwargs["resize"],
        kwargs["normalize"],
        kwargs["dropout_p"],
        kwargs["batch_size"],
        kwargs["use_compile"],
        kwargs["world_size"],
    )
    if "MASTER_ADDR" not in os.environ:
        hostnode = subprocess.getoutput(
            "scontrol show hostnames $SLURM_NODELIST | head -n 1"
        )
        if "error" in hostnode:
            hostnode = "localhost"
        os.environ["MASTER_ADDR"] = hostnode
    if "MASTER_PORT" not in os.environ:
        os.environ["MASTER_PORT"] = "29400"

    main(**kwargs)
    # mp.spawn(main, args=main_args, nprocs=args.world_size, join=True)
