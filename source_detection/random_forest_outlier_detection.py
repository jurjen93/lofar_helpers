from astropy.table import Table
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score, accuracy_score
from sklearn.model_selection import RandomizedSearchCV as RSCV
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import tree
import os
from glob import glob


# Table columns
cols = ['Source_id', 'E_Total_flux', 'E_Peak_flux',
        'S_Code', 'E_pos', 'Total_flux/peak_flux',
        'Peak_flux/RMS', 'Total_flux/RMS', 'E_MajMin',
        'MajMin', 'PA', 'E_PA']

sep = '\n##########################################################################\n'


def gridsearch_test(df):
    """

    Grid search test

    :param df: pandas data frame
    :return: best model
    """

    df = df.dropna()

    param_grid = {'n_estimators': np.arange(50, 200, 15),
                  'max_features': np.arange(0.1, 1, 0.1),
                  'max_depth': [3, 4, 5],
                  'max_samples': [0.3, 0.4, 0.5, 0.6, 0.7, 0.8]}

    X_train, X_test, y_train, y_test = train_test_split(df[[c for c in cols if c not in ['Source_id', 'label']]],
                                                        df.label, test_size=0.3, random_state=42)
    model = RSCV(RandomForestClassifier(), param_grid, n_iter=15).fit(X_train,
                                                                                          y_train.astype('int'))
    model = model.best_estimator_
    print(sep + 'BEST MODEL\n')
    print(model)
    print(sep)
    predictions = model.predict(X_test)
    mse = mean_squared_error(y_test, predictions)
    print('MSE: ' + str(mse))
    r2 = r2_score(y_test, predictions.round(0))
    print('R2: ' + str(r2))
    indices = np.argwhere((y_test - predictions.round(0)).to_numpy() != 0)
    print('Differently classified sources:' + str(y_test.index[indices].squeeze()))
    print('Accuracy: '+str(accuracy_score(y_test.astype(int), np.array(predictions).astype(int))))
    print(sep)
    return model


def gridsearch(df):
    """

    Grid search

    :param df: pandas data frame
    :return: best model
    """
    df = df.dropna()

    param_grid = {'n_estimators': np.arange(30, 200, 5),
                  'max_features': np.arange(0.1, 1, 0.1),
                  'max_depth': [3, 4, 5],
                  'max_samples': [0.3, 0.4, 0.5, 0.6, 0.7, 0.8]}

    model = RSCV(RandomForestClassifier(), param_grid, n_iter=30). \
        fit(df[[c for c in cols if c not in ['Source_id', 'label']]], df.label.astype('int'))
    model = model.best_estimator_
    print(sep + 'BEST MODEL\n')
    print(model)
    print(sep)
    return model


def sort_feature_importance(feat_importance, labels):
    """
    Sort feature importance from model

    :param feat_importance: unsorted feature importance model
    :param labels: unsorted feature importance labels

    :return: sorted feature importance
    """
    important_features_dict = {}
    for idx, val in enumerate(feat_importance):
        important_features_dict[idx] = val

    important_features_list = sorted(important_features_dict,
                                     key=important_features_dict.get,
                                     reverse=True)
    return feat_importance[important_features_list], np.array(labels)[important_features_list]


def get_table_data(tbl):
    """
    Get table data for training/testing/predicting

    :param tbl: table name

    :return: pandas data frame

    """
    T = Table.read(tbl)
    df = T.to_pandas()

    le = LabelEncoder()
    le.fit([b'S', b'C', b'M'])
    df['S_Code'] = le.transform(df.S_Code)
    df['Total_flux/peak_flux'] = df['Total_flux'] / df['Peak_flux']
    df['Peak_flux/RMS'] = df['Peak_flux'] / df['Isl_rms']
    df['Total_flux/RMS'] = df['Total_flux'] / df['Isl_rms']
    df['E_pos'] = np.sqrt(df['E_RA']**2+df['E_DEC']**2)
    df["E_MajMin"] = np.sqrt(df["E_Maj"]**2 + df["E_Min"]**2)
    df["MajMin"] = df["Maj"] * df["Min"]

    df = df[cols]

    return df


def plot_feature_importance(feat_importance, labels):
    """
    Plot feature importance

    :param feat_importance: feature importance
    :param labels: labels
    """
    plt.style.use("bmh")

    hfont = {'fontname': 'monospace'}
    forest_importances = pd.Series(feat_importance,
                                   index=[l for l in labels])
    fig, ax = plt.subplots()
    forest_importances.plot.bar(ax=ax, color='darkred')
    # ax.set_title()
    ax.set_ylabel("Feature importance")
    plt.xticks(rotation=75, **hfont)
    fig.tight_layout()
    plt.grid(False)
    for pos in ['right', 'top', 'bottom', 'left']:
        plt.gca().spines[pos].set_visible(False)
    plt.savefig('feat_importance.png', dpi=300)
    # plt.show()


def predict_falsepositives(tbl, model):
    """
    Predict which sources are false positives with model

    :param tbl: table name
    :param model: RF model

    :return sources to delete
    """
    df1 = get_table_data(tbl).set_index("Source_id")[[c for c in cols if c not in ['Source_id', 'label']]]
    predicted = model.predict(df1)
    df1['prediction'] = predicted
    return list(df1[df1.prediction == 0].index), df1.columns


def plot_tree(model, tree_idx, columns):
    """
    Plot a tree
    """
    plt.figure(figsize=(200, 200))
    _ = tree.plot_tree(model.estimators_[tree_idx], feature_names=columns, filled=True)
    plt.show()

def get_table_index(t, source_id):
    return int(np.argwhere(t['Source_id'] == source_id).squeeze())

def clean_table(tbl, model):
    """
    Make clean table

    :param tbl: fits table name
    """
    to_remove, columns = predict_falsepositives(tbl, model)
    T = Table.read(tbl, format='fits')
    print(tbl)
    for i in sorted(to_remove)[::-1]:
        n = get_table_index(T, i)
        if T[n]['Peak_flux'] < 5.5*T[n]['Isl_rms']:
            print('Deleted source ID: ' + str(i))
            del T[n]

    T.write(tbl.replace('.fits', '_clean.fits'), format='fits', overwrite=True)


def main():
    """
    Main function
    """

    train_table = "finalcat03/facet_19_source_catalog_final.fits"
    df = get_table_data(train_table)

    # false positive indices from test data by eye
    falsepositives = [246, 248, 249, 252, 253, 254, 258, 330, 263, 264, 282, 288, 289, 290, 291, 292, 293, 294,
                      295, 296, 297, 298, 299, 300, 301, 303, 305, 308, 347, 310, 313, 315, 316, 318, 319, 320,
                      322, 323, 327, 330, 331, 333, 334, 335, 341, 255, 259, 262, 266, 267, 268, 279, 284, 367,
                      362, 358, 357, 356, 355, 351, 349, 344, 343, 342]

    falsepositives = sorted(list(set(falsepositives)))

    df['label'] = None
    df.loc[df.Source_id.isin(falsepositives), 'label'] = False
    df.loc[~df.Source_id.isin(falsepositives), 'label'] = True

    gridsearch_test(df)
    model = gridsearch(df)
    importances, labels = sort_feature_importance(model.feature_importances_,
                                                  [c for c in cols if c not in ['Source_id', 'label']])
    plot_feature_importance(importances, labels)

    for cat in glob("finalcat03/facet_*_source_catalog.fits"):
        clean_table(cat, model)

if __name__ == '__main__':
    main()
