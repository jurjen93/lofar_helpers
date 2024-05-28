import matplotlib.image
import torch
from pytorch_lightning import Trainer
from tqdm import tqdm

from train_cnn import ImagenetTransferLearning
from pre_processing_for_ml import FitsDataset


@torch.no_grad()
def variational_dropout(model, dataloader, variational_iters=25):
    model.feature_extractor.eval()
    model.classifier.train()

    model.cuda()

    for sample in tqdm(dataloader):
        batch, labels = sample[0].cuda(), sample[1].cuda()
        logits = torch.concat([model(batch).clone() for _ in range(variational_iters)], dim=1)
        stds = logits.std(dim=1)
        preds = torch.sigmoid(logits).mean(dim=1)

        yield (
            preds.cpu(),
            stds.cpu()
        )


def save_images(dataset, out_paths, preds, stds):
    assert len(dataset) == len(preds)

    for elem, in_path, pred, std in tqdm(zip(dataset, dataset.data_paths, preds, stds), total=len(dataset)):
        batch, label = elem
        name = in_path.strip('.npz').split('/')[-1]
        matplotlib.image.imsave(
            fname=f'{out_paths}/{std:.3f}_{pred:.3f}_{label}_{name}.png',
            arr=batch.movedim(0, -1).numpy()
        )

def main(dataroot, checkpoint_path):
    torch.set_float32_matmul_precision('high')

    model = ImagenetTransferLearning.load_from_checkpoint(checkpoint_path)

    # num_workers = len(os.sched_getaffinity(0))
    # num_workers = 0
    num_workers = 6
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
            drop_last=False,  # needed for torch.compile,
            shuffle=False,  # needed so that val AP is non nan
        )
        for mode in ('train', 'validation')
    )
    out_path = '/scratch-shared/CORTEX/public.spider.surfsara.nl/project/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/preds_train/'

    # val_pred, val_std = map(torch.concat, zip(*[elem for elem in variational_dropout(model, val_dataloader)]))
    #
    # save_images(FitsDataset(root, mode='validation'), out_path, val_pred, val_std)
    train_pred, train_std = map(torch.concat, zip(*[elem for elem in variational_dropout(model, train_dataloader)]))
    save_images(FitsDataset(root, mode='train'), out_path, train_pred, train_std)



if __name__ == '__main__':
    # root = f'{os.environ["TMPDIR"]}/public.spider.surfsara.nl/project/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/cnn_data/'
    root = f'/dev/shm/scratch-shared/CORTEX/public.spider.surfsara.nl/project/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/cnn_data/'
    checkpoint_path = '/gpfs/home5/robertsc/CORTEX/lofar_helpers/source_selection/lightning_logs/version_5471545/checkpoints/epoch=89-step=4680.ckpt'

    main(root, checkpoint_path)
