"""
Make Deep Learning model
It dumps a pickle which should be used to decide if recalibration is needed
"""

__author__ = "Jurjen de Jong (jurjendejong@strw.leidenuniv.nl)"

import warnings
warnings.filterwarnings("ignore")
import tensorflow as tf
from tensorflow.keras import layers
from tensorflow.keras.models import Sequential
from glob import glob
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
import joblib

def flip_image(image):
    """Flip images upside-down"""
    return image[::-1, :, :]

def mirror_image(image):
    """Mirror image left-to-right"""
    return image[:, ::-1, :]

def train_model(X_train, y_train, save=False):
    """
    Train the model with Keras layers
    :param X_train: Training data (Numpy array)
    :param y_train: Training labels (Numpy array)
    """
    model = Sequential()
    # model.add(layers.Conv2D(16, (2, 2), padding='same', activation='relu', input_shape=(img_height, img_width, 1)))
    # model.add(layers.MaxPooling2D()),
    # model.add(layers.Conv2D(32, (2, 2), padding='same', activation='relu'))
    # model.add(layers.MaxPooling2D())
    # model.add(layers.Conv2D(64, (2, 2), padding='same', activation='relu'))
    # model.add(layers.MaxPooling2D())
    model.add(layers.Flatten(input_shape=(X_train.shape[1], X_train.shape[2], X_train.shape[3])))
    model.add(layers.Dense(128, activation='relu'))
    model.add(layers.Dropout(.3))
    model.add(layers.Dense(64, activation='relu'))
    model.add(layers.Dropout(.3))
    model.add(layers.Dense(32, activation='relu'))
    model.add(layers.Dense(2, activation="softmax"))

    model.compile(optimizer='adam',
                  loss=tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True),
                  metrics=['accuracy'])
    model.summary()

    model.fit(X_train, y_train, epochs=5, verbose=2)
    if save:
        joblib.dump(model, 'recal_dec_model.sav')
    return model

def evaluate_model(X_test, y_test, model):
    """
    Evaluate deep learning model
    :param X_test: Test data (Numpy array)
    :param y_test: Test labels (Numpy array)
    :param model: keras/tensorflow model
    """
    test_loss, accuracy = model.evaluate(X_test, y_test, verbose=0)
    print('Accuracy: %.2f' % (accuracy * 100))
    print('Test_loss: %.2f' % (test_loss * 100))

    return


if __name__ == '__main__':

    df_labels = pd.read_csv('DL_data/DL_data.csv').set_index('Name')
    labels = []
    images = []
    for im in glob('DL_data/numpy/*.npy'):
        image = np.load(im)
        image /= np.max(image)
        image = image.reshape(list(image.shape) + [1])
        images.append(image)
        # Image augmentation by flipping and mirroring
        images.append(flip_image(image))
        images.append(mirror_image(image))
        images.append(flip_image(mirror_image(image)))
        # Adding the label four times because of the augmentation
        labels.extend([df_labels.loc[im.split('/')[-1].split('.npy')[0], 'Recalibrate']]*4)

    images = np.array(images)
    labels = np.array(labels)
    images = images

    X_train, X_test, y_train, y_test = train_test_split(images, labels, test_size=0.33, random_state=42)

    print(df_labels.Recalibrate.value_counts() / len(df_labels))
    model = train_model(X_train, y_train, save=False)
    evaluate_model(X_test, y_test, model)