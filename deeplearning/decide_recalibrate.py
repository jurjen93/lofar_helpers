import warnings
warnings.filterwarnings("ignore")
import numpy as np
from tensorflow.keras.models import load_model
from skimage.transform import resize
from matplotlib.colors import SymLogNorm

def predict_recalibration(image, model_name):
    """
    Preprocessing image
    --> only points bigger than 5 times the noise are considered
    --> resize image
    --> divided by the max value of the image
    --> reshape image for the neural network
    :param image: input image
    :param model_name: name of the model
    """
    model = load_model(model_name)
    image = np.where(image > 5 * np.std(image), 0, image)
    image = resize(image, (150, 150))
    image /= np.max(image)
    image = image.reshape(list(image.shape) + [1])
    predict_score = model.predict(np.array([image]))
    return round(predict_score[0][0],3), int(round(predict_score[0][0],0))


if __name__ == '__main__':
    import pandas as pd
    from astropy.io import fits
    import matplotlib.pyplot as plt
    from glob import glob
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('--model_name', type=str, help='Name of the output model (h5 file)', default='best_model.h5')
    args = parser.parse_args()

    df_labels = pd.read_csv('DL_data/labels.csv').set_index('Name')
    labels = []
    images = []
    for im in glob('DL_data/numpy/*.npy') + glob('DL_data/fits/*'):
        if '.npy' in im:
            image = np.load(im)
        else:
            hdul = fits.open(im)[0]
            image = hdul.data
            while len(image.shape) != 2:
                image = image[0]
        images.append(image)
        labels.append(df_labels.loc[im.split('/')[-1].split('.npy')[0], 'Recalibrate'])

    for n, im in enumerate(images):
        score, label = predict_recalibration(im, model_name=args.model_name)

        if label!=labels[n]:
            plt.imshow(im, norm=SymLogNorm(linthresh=np.nanstd(im)/50, vmin=np.nanstd(im)/1000, vmax=np.nanstd(im)*20),
                       cmap='CMRmap', origin='lower')
            plt.title('Model score: '+str(score)+', Real: '+str(labels[n]))
            plt.show()
