import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from argparse import ArgumentParser
import os
from glob import glob

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-p', '--path', type=str, help='h5 file name for output')
    args = parser.parse_args()
    if args.path:
        path = args.path
    else:
        path = os.path.dirname(os.path.realpath(__name__))

    results = os.system('mkdir {path}/results'.format(path=path))
    boxes = glob('{path}/box_*'.format(path=path))

    for box in boxes:
        images = glob('{box_path}/image_*.png'.format(box_path=box))
        first_image, last_image = images[0], images[-1]
        fig, axs = plt.subplots(1, 2)
        im_0 = mpimg.imread(first_image)
        axs[0].imshow(im_0)
        axs[0].title.set_text(first_image)
        im_1 = mpimg.imread(last_image)
        axs[1].imshow(im_1)
        axs[1].title.set_text(last_image)
        fig.suptitle(box.split('/')[-1])
        [axi.set_axis_off() for axi in axs.ravel()]
        plt.savefig('{path}/results/{box}.png'.format(path=path, box=box.split('/')[-1]))
        plt.close()