from ..train_nn import *
import matplotlib.pyplot as plt
import torchvision.transforms.functional as F

@lru_cache(maxsize=1)
def get_transforms():
    return v2.Compose([
        v2.ColorJitter(brightness=.5, hue=.3, saturation=0.1, contrast=0.1),
        v2.RandomInvert(),
        v2.RandomEqualize(),
        v2.RandomVerticalFlip(p=0.5),
        v2.RandomHorizontalFlip(p=0.5),
    ])

def compute_statistics(loader, normalize: int):
    if not normalize:
        return torch.asarray([0, 0, 0]), torch.asarray([1, 1, 1])
    means = []
    sums_of_squares = []
    f = torch.log if normalize==2 else lambda x: x
    for i, (imgs, _) in enumerate(loader):
        print(i, len(loader))
        imgs = imgs.to('cuda')
        imgs = f(imgs)
        means.append(torch.mean(imgs, dim=(0, 2, 3)))
        sums_of_squares.append((imgs**2).sum(dim=(0, 2, 3)))
    mean = torch.stack(means).mean(0)
    sums_of_squares = torch.stack(sums_of_squares).sum(0)
    variance = (sums_of_squares / (len(loader) * imgs.shape[0] * imgs.shape[2] * imgs.shape[3])) - (mean ** 2)
    return mean, torch.sqrt(variance)

def plot_image(img, fname):
    img = img[0].cpu().permute(1, 2, 0).to(torch.float32).numpy()
    plt.imshow((img - np.min(img))/(np.max(img)-np.min(img)))
    plt.savefig(fname)

if __name__=='__main__':
    dataset_root = 'public.spider.surfsara.nl/project/lofarvwf/jdejong/CORTEX/calibrator_selection_robertjan/cnn_data'
    batch_size = 32
    normalize = 2
    train_dataloader, val_loader = get_dataloaders(dataset_root, batch_size)
    print('computing_statistics')
    mean, std = compute_statistics(train_dataloader, normalize=normalize)
    print('done_with_statistics')

    prepare_data_f = partial(prepare_data, resize=0, device='cuda', mean=mean, std=std, normalize=normalize)

    output_folder = 'image_samples/'
    os.makedirs(output_folder, exist_ok=True)
    data, labels = next(iter(val_loader))

    data_norm, labels = prepare_data_f(data, labels)
    plot_image(data_norm, f"{output_folder}/normalize_{normalize}")
    transforms = {'brightness': [v2.ColorJitter(brightness=0.5)],
                  'hue': [v2.ColorJitter(hue=0.3)],
                  'saturation': [v2.ColorJitter(saturation=0.1)],
                  'contrast': [v2.ColorJitter(contrast=0.1)],
                  'colorjitter': [v2.ColorJitter(brightness=0.5, hue=0.3, saturation=0.1, contrast=0.1)], 
                  'invert': [v2.RandomInvert(1)],
                  'equalize': [v2.RandomEqualize(1)],
                  'all': [v2.ColorJitter(brightness=0.5, hue=0.3, saturation=0.1, contrast=0.1), v2.RandomInvert(1),v2.RandomEqualize(1)]}
    for t_name, transformations in transforms.items():
        transform_f = v2.Compose(transformations)
        data_transformed = transform_f(data_norm)
        plot_image(data_transformed, f"{output_folder}/normalize_{normalize}_{t_name}")
    
    print('done')
    exit() 



