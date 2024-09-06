import os

import matplotlib.image
import torch
from torch.utils.data import DataLoader
from tqdm import tqdm

from train_nn import ImagenetTransferLearning, load_checkpoint
from pre_processing_for_ml import FitsDataset


@torch.no_grad()
def variational_dropout(model, dataloader, variational_iters=25):
    model.feature_extractor.eval()
    model.classifier.train()

    model.cuda()

    for sample in tqdm(dataloader):
        with torch.cuda.amp.autocast(dtype=torch.bfloat16):
            batch, labels = sample[0].cuda(non_blocking=True), sample[1].cuda(non_blocking=True)
            preds = torch.sigmoid(torch.concat([model(batch).clone() for _ in range(variational_iters)], dim=1))

            means = preds.mean(dim=1)
            stds = preds.std(dim=1)

        yield (
            means.cpu(),
            stds.cpu()
        )


def save_images(dataset, out_paths, preds, stds, mode):
    for elem, in_path, pred, std in tqdm(zip(dataset, dataset.data_paths, preds, stds), total=len(dataset)):
        batch, label = elem
        name = in_path.strip('.npz').split('/')[-1]
        matplotlib.image.imsave(
            fname=f'{out_paths}/{std:.3f}_{pred:.3f}_{label}_{name}.png',
            arr=batch.to(torch.float).movedim(0, -1).numpy()
        )

def main(dataset_root, checkpoint_path):
    torch.set_float32_matmul_precision('high')

    ckpt_dict = load_checkpoint(checkpoint_path)
    model = ckpt_dict['model']
    breakpoint()

    num_workers = min(18, len(os.sched_getaffinity(0)))
    prefetch_factor, persistent_workers = (
        (2, True) if num_workers > 0 else
        (None, False)
    )
    batch_size = 64

    def gen_and_save(mode):
        dataset = FitsDataset(dataset_root, mode=mode)

        dataloader = DataLoader(
            dataset=dataset,
            batch_size=batch_size,
            num_workers=num_workers,
            prefetch_factor=prefetch_factor,
            persistent_workers=False,
            pin_memory=True,
            shuffle=False,
            drop_last=False,
        )

        out_path = f'/scratch-shared/CORTEX/public.spider.surfsara.nl/project/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/preds_{mode}/'

        preds, stds = map(torch.concat, zip(*[elem for elem in variational_dropout(model, dataloader)]))
        save_images(dataset, out_path, preds, stds, mode=mode)

    for mode in ('train', 'val'):
        gen_and_save(mode)



if __name__ == '__main__':
    # root = f'{os.environ["TMPDIR"]}/public.spider.surfsara.nl/project/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/cnn_data/'
    root = f'/dev/shm/scratch-shared/CORTEX/public.spider.surfsara.nl/project/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/cnn_data/'
    checkpoint_path = './gridsearch_efficientnet/version_6665871_5__model_efficientnet_v2_l__lr_0.0001__normalize_0__dropout_p_0.25/ckpt_step=7999.pth'

    main(root, checkpoint_path)
