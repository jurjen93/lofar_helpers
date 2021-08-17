"""
Make Deep Learning model
It dumps a pickle which should be used to decide if recalibration is needed
"""

__author__ = "Jurjen de Jong (jurjendejong@strw.leidenuniv.nl)"

import warnings
warnings.filterwarnings("ignore")
from tensorflow.keras import layers
from tensorflow.keras.models import Sequential
from tensorflow.keras import callbacks
from glob import glob
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc
from astropy.io import fits
from skimage.transform import resize
import matplotlib.pyplot as plt


def flip_image(image):
    """Flip images upside-down"""
    if len(image.shape) == 3:
        return image[::-1, :, :]
    elif len(image.shape) == 2:
        return image[::-1, :]
    else:
        print('Image shape is: ' + str(image.shape) + ' but should have 2 or 3 dimensions')
        return image


def mirror_image(image):
    """Mirror image left-to-right"""
    if len(image.shape) == 3:
        return image[:, ::-1, :]
    elif len(image.shape) == 2:
        return image[:, ::-1]
    else:
        print('Image shape is: ' + image.shape + ' but should have 2 or 3 dimensions')
        return image

class RecalRecognition:
    def __init__(self, images, labels, image_size=224, model_name='best_model.h5'):
        """
        :param images: input images
        :param labels: input labels
        :param image_size: image size (note that images are always squares)
        :param model_name: name of the output model
        """
        self.images = images
        self.labels = labels
        self.im_size = image_size
        self.model_name = model_name

    def _preprocessing(self, image):
        """
        Preprocessing image
        --> only points bigger than 5 times the noise are considered
        --> resize image
        --> divided by the max value of the image
        --> reshape image for the neural network
        :param image: input image
        """
        image = np.where(image > 5*np.std(image), 0, image)
        image = resize(image, (self.im_size, self.im_size))
        image /= np.max(image)
        image = image.reshape(list(image.shape) + [1])
        return image

    def _preprocessing_en(self, image):
        """
        Preprocessing image
        --> only points bigger than 5 times the noise are considered
        --> resize image
        --> divided by the max value of the image
        --> reshape image for the neural network
        :param image: input image
        """
        image = np.where(image > 5*np.std(image), 0, image)
        image = resize(image, (self.im_size, self.im_size))
        image /= np.max(image)
        image_new = np.zeros(list(image.shape) + [3])
        for i in range(3):
            image_new[:, :, i] += image
        return image_new

    def preprocess_data(self, EN=False):
        """
        Preprocessing all input images
        :param EN: Using efficient net or not (because they need a (224, 224, 3) input)
        """
        if EN:
            self.images = [self._preprocessing_en(i) for i in self.images]
        else:
            self.images = [self._preprocessing(i) for i in self.images]
        self.images = np.array(self.images)
        print("Shape of images: " + str(self.images.shape))
        self.labels = np.array(self.labels)
        return self

    def split_data(self):
        """Split data in train and test"""
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(self.images, self.labels, test_size=0.25, random_state=42)
        return self

    def make_model(self):
        """
        Set model layers
        """
        self.model = Sequential()
        self.model.add(layers.Conv2D(16, (2, 2), padding='same', activation='relu',
                                     input_shape=(self.im_size, self.im_size, 1)))
        self.model.add(layers.MaxPooling2D(2, 2)),
        self.model.add(layers.Flatten())
        self.model.add(layers.Dense(64, activation='relu'))
        self.model.add(layers.Dropout(.3))
        self.model.add(layers.Dense(32, activation='relu'))
        self.model.add(layers.Dropout(.3))
        self.model.add(layers.Dense(16, activation='relu'))
        self.model.add(layers.Dense(1, activation="sigmoid"))

        self.model.compile(optimizer='adam',
                      loss='binary_crossentropy',
                      metrics=['accuracy'])
        self.model.summary()

        return self

    def efficientnet_model(self):
        """
        Efficientnet is an optimized automatic neural network [But doesn't work yet for our case ...]
        """
        import efficientnet.keras as efn
        self.model = efn.EfficientNetB0(input_shape=(self.im_size, self.im_size, 3), weights='imagenet')
        return self


    def fit_model(self):
        self.history = self.model.fit(self.X_train,
                                     self.y_train,
                                     batch_size=len(self.X_train) // 3,
                                     epochs=20,
                                     verbose=2,
                                     validation_data=(self.X_test, self.y_test),
                                     callbacks=[callbacks.EarlyStopping(monitor='val_loss', patience=3),
                                                callbacks.ModelCheckpoint(self.model_name, save_best_only=True)])
        return self

    def evaluate_model(self):
        """
        Evaluate deep learning model
        --> Area under the curve
        --> Accuracy over epoch
        --> Loss over epoch
        """

        test_loss, accuracy = self.model.evaluate(self.X_test, self.y_test, verbose=0)


        preds = self.model.predict(self.X_test, verbose=1)
        fpr, tpr, _ = roc_curve(self.y_test, preds)
        roc_auc = auc(fpr, tpr)

        print('Label 1 in training data: '+str(round(np.sum(self.y_train)/len(self.y_train)*100, 2))+'%')
        print('Label 1 in test data: '+str(round(np.sum(self.y_test)/len(self.y_test)*100, 2))+'%')
        print('Label 1 in predictions: '+str(round(np.sum(preds)/len(preds)*100, 2))+'%')
        print('Accuracy: %.2f' % (accuracy * 100))
        print('Test_loss: %.2f' % (test_loss * 100))
        print('AUC: %.2f' %(round(roc_auc*100,2)))

        """
        True positive rate (tpr) is TP/(TP+FN)
        False positive rate (fpr) is FP/(FP+TN)
        The area under the curve (AUC) is the area under this tpr vs. fpr plot.
        For useful for imbalanced data sets (not equally 50/50)
        """
        f, ax = plt.subplots(1, 3)
        ax[0].plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve (area = %0.2f)'%roc_auc)
        ax[0].plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
        ax[0].set_xlim([0.0, 1.0])
        ax[0].set_ylim([0.0, 1.05])
        ax[0].set_title('ROC curve')
        ax[0].set_xlabel('False positive rate')
        ax[0].set_ylabel('True positive rate')
        ax[0].legend(loc='lower right')

        ax[1].plot(self.history.history['accuracy'], color='blue')
        ax[1].plot(self.history.history['val_accuracy'], color='red', linestyle='--')
        ax[1].set_title('Model Accuracy')
        ax[1].set_ylabel('accuracy')
        ax[1].set_xlabel('epoch')
        ax[1].legend(['train', 'test'], loc='upper left')

        ax[2].plot(self.history.history['loss'], color='blue')
        ax[2].plot(self.history.history['val_loss'], color='red', linestyle='--')
        ax[2].set_title('Model Loss')
        ax[2].set_ylabel('loss')
        ax[2].set_xlabel('epoch')
        ax[2].legend(['train', 'test'], loc='upper left')
        plt.savefig('evaluation_model.png')
        plt.show()

        return self


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('--model_name', type=str, help='Name of the output model (h5 file)', default='best_model.h5')
    args = parser.parse_args()

    # ----- LOAD DATA -----
    df_labels = pd.read_csv('DL_data/labels.csv').set_index('Name')
    labels = []
    images = []
    for im in glob('DL_data/numpy/*.npy')+glob('DL_data/fits/*'):
        if '.npy' in im:
            image = np.load(im)
        else:
            hdul = fits.open(im)[0]
            image = hdul.data
            while len(image.shape)!=2:
                image = image[0]
        images.append(image)
        images.append(flip_image(image))
        images.append(mirror_image(image))
        images.append(flip_image(mirror_image(image)))
        labels.extend([df_labels.loc[im.split('/')[-1].split('.npy')[0], 'Recalibrate']]*4)

    # ----- TRAIN AND TEST MODEL -----
    Model = RecalRecognition(images, labels, image_size=150, model_name=args.model_name)
    Model.preprocess_data().\
        split_data().\
        make_model().\
        fit_model().\
        evaluate_model()