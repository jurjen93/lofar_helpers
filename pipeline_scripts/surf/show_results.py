__author__ = "Jurjen de Jong (jurjendejong@strw.leidenuniv.nl)"

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from argparse import ArgumentParser
import os
from glob import glob

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-p', '--path', type=str, help='output folder')
    parser.add_argument('-in', '--include_boxes', help='Include only the following boxes (numbers only)')
    args = parser.parse_args()
    if args.path:
        path = args.path
    else:
        path = os.path.dirname(os.path.realpath(__name__))

    results = os.system('mkdir {path}/results'.format(path=path))
    boxes = glob('{path}/box_*'.format(path=path))

    if args.include_boxes:
        included_boxes = ['box_' + n for n in args.include_boxes.split(',')]
        boxes = [b for b in boxes if b.split('/')[-1] in included_boxes]

    for box in boxes:
        images = sorted(glob('{box_path}/image_*.png'.format(box_path=box)))
        if len(images)>0:
            first_image, last_image = images[0], images[-1]
            fig, axs = plt.subplots(1, 2, figsize=(18,6))
            im_0 = mpimg.imread(first_image)
            axs[0].imshow(im_0)
            axs[0].title.set_text(first_image.split('/')[-1])
            im_1 = mpimg.imread(last_image)
            axs[1].imshow(im_1)
            axs[1].title.set_text(last_image.split('/')[-1])
            fig.suptitle(box.split('/')[-1])
            [axi.set_axis_off() for axi in axs.ravel()]
            plt.savefig('{path}/results/{box}.png'.format(path=path, box=box.split('/')[-1]))
            scalarcomplexgainplots = glob(box+'/plotlosoto*/scalarcomplexgain*007*.png')
            for plot in scalarcomplexgainplots:
                for el in plot.split('.'):
                    if 'archive' in el:
                        archive=el
                        break
                os.system('cp '+plot+' '+path+'/results/'+box.split('/')[-1]+'_'+archive+'_'+plot.split('/')[-1])
            plt.close()