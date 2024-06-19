import matplotlib.image
import torch
from pytorch_lightning import Trainer
from tqdm import tqdm

from train_cnn import ImagenetTransferLearning, load_checkpoint
from pre_processing_for_ml import FitsDataset


@torch.no_grad()
def variational_dropout(model, dataloader, variational_iters=25):
    model.feature_extractor.eval()
    model.classifier.train()

    model.cuda()

    for sample in tqdm(dataloader):
        with torch.cuda.amp.autocast(dtype=torch.bfloat16):
            batch, labels = sample[0].cuda(), sample[1].cuda()
            preds = torch.sigmoid(torch.concat([model(batch).clone() for _ in range(variational_iters)], dim=1))

            means = preds.mean(dim=1)
            stds = preds.std(dim=1)

        yield (
            means.cpu(),
            stds.cpu()
        )


def save_images(dataset, out_paths, preds, stds, mode):
    train_len, val_len = dataset.train_len, dataset.val_len
    if mode == 'val':
        idx = torch.arange(train_len, train_len+val_len)
        data_paths = dataset.val_paths
        dataset_len = val_len
    else:
        idx = torch.arange(train_len)
        data_paths = dataset.train_paths
        dataset_len = train_len

    assert dataset_len == len(preds)

    for index, in_path, pred, std in tqdm(zip(idx, data_paths, preds, stds), total=dataset_len):
        batch, label = dataset[index]
        name = in_path.strip('.npz').split('/')[-1]
        matplotlib.image.imsave(
            fname=f'{out_paths}/{std:.3f}_{pred:.3f}_{label}_{name}.png',
            arr=batch.to(torch.float).movedim(0, -1).numpy()
        )

def main(dataroot, checkpoint_path):
    torch.set_float32_matmul_precision('high')

    ckpt_dict = load_checkpoint(checkpoint_path)
    model = ckpt_dict['model']

    # num_workers = len(os.sched_getaffinity(0))
    # num_workers = 0
    num_workers = 12
    prefetch_factor, persistent_workers = (
        (2, True) if num_workers > 0 else
        (None, False)
    )

    sample_loader = torch.utils.data.DataLoader(
        dataset=FitsDataset(root),
        batch_size=1,
        num_workers=num_workers,
        prefetch_factor=prefetch_factor,
        persistent_workers=persistent_workers,
        pin_memory=False,
        drop_last=False,  # needed for torch.compile,
        shuffle=False,  # needed so that val AP is non nan
    )

    train_dataloader = get_subdataloader(sample_loader, mode='train')
    val_dataloader = get_subdataloader(sample_loader, mode='val')

    out_path = '/scratch-shared/CORTEX/public.spider.surfsara.nl/project/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/preds_val3/'

    val_pred, val_std = map(torch.concat, zip(*[elem for elem in variational_dropout(model, val_dataloader)]))
    save_images(FitsDataset(root), out_path, val_pred, val_std, mode='val')
    # train_pred, train_std = map(torch.concat, zip(*[elem for elem in variational_dropout(model, train_dataloader)]))
    # save_images(FitsDataset(root, mode='train'), out_path, train_pred, train_std)


def get_subdataloader(dataloader, mode):
    train_len = dataloader.dataset.train_len // dataloader.batch_size
    for i, sample in enumerate(dataloader):
        print(i)
        if (mode == 'train' and i < train_len) or (mode == 'val' and i >= train_len):
            yield sample


if __name__ == '__main__':
    # root = f'{os.environ["TMPDIR"]}/public.spider.surfsara.nl/project/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/cnn_data/'
    root = f'/dev/shm/scratch-shared/CORTEX/public.spider.surfsara.nl/project/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/cnn_data/'
    checkpoint_path = './gridsearch_efficientnet/version_6665871_5__model_efficientnet_v2_l__lr_0.0001__normalize_0__dropout_p_0.25/ckpt_step=7999.pth'

    main(root, checkpoint_path)
