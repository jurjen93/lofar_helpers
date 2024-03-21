import numpy as np
import sys


def replace_nan_with_avg(img, kernel_size=5):
    """
    Replace nan values with average in kernel of X by X size

    :param img: image (grayscale)
    :param kernel_size: kernel size (X by X)

    :return: smoothed image
    """

    if kernel_size % 2 == 0:
        sys.exit("Kernel needs to be uneven")
    else:
        k = (kernel_size-1)//2

    rows, cols = img.shape
    result_img = np.copy(img)

    coor_nan = np.argwhere(img != img)

    for c in coor_nan:
        row, col = c
        row_min = max(row - k, 0)
        row_max = min(row + k, rows - 1)
        col_min = max(col - k, 0)
        col_max = min(col + k, cols - 1)

        neighborhood = img[row_min:row_max + 1, col_min:col_max + 1]
        avg_val = np.nanmean(neighborhood)

        result_img[row, col] = avg_val

    return result_img