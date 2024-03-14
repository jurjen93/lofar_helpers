import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('TkAgg')

SPD_INDEX = 0
SOLINT_INDEX = 1
RA_IDX = 2
DEC_IDX = 3

MU = 0
SIG = 0.3


def gaussian(x, mu, sig):
    return 1.0 / (np.sqrt(2.0 * np.pi) * sig) * np.exp(-np.power((x - mu) / sig, 2.0) / 2)


def apply_gauss(work_arr, source_point, mu=MU, sig=SIG):
    temp_arr = np.copy(work_arr)
    ra_diff = temp_arr[:, RA_IDX] - source_point[RA_IDX]
    dec_diff = temp_arr[:, DEC_IDX] - source_point[DEC_IDX]

    distance = np.sqrt(ra_diff ** 2 + dec_diff ** 2)

    scale = gaussian(0, mu, sig)

    work_arr[:, SPD_INDEX] *= 1 - (gaussian(distance, mu, sig) / scale)


def plot_data(arr, selection=None, text=None):
    graph, plot = plt.subplots(1, 1)
    plt.scatter(arr[:, RA_IDX], arr[:, DEC_IDX], arr[:, SPD_INDEX], alpha=0.5)
    if selection is not None:
        plt.scatter(selection[:, RA_IDX], selection[:, DEC_IDX], selection[:, SPD_INDEX])

    if text is not None:
        for text, x, y in zip(text, selection[:, RA_IDX], selection[:, DEC_IDX]):
            plt.annotate(text, (x, y))
    plot.invert_xaxis()
    plt.xlabel("DEC")
    plt.ylabel("RA")
    plt.show()


def compare_selected_names(selection_names, output=False):
    selection = sorted(selection_names)
    good_selection = np.loadtxt("good-selection.txt", dtype=str)

    s_selection = set(selection)
    s_good = set(good_selection)

    if output:
        print(f"Selected: {len(s_selection)}. True: {len(s_good)}")
        print(f"Same in {len(s_selection.intersection(s_good))}/{len(s_selection.union(s_good))} cases.")

        print("name:   S T")
        for name in sorted(list(s_selection.union(s_good))):
            print(f"{name} {name in selection_names} {name in good_selection}")

    return len(s_selection.intersection(s_good)) / len(s_selection.union(s_good))


def main():
    arr, names = load_data()

    selection_indices = select_sources(arr, threshold=0.1, sig=0.12)
    selection_names = names[selection_indices]
    compare_selected_names(selection_names, output=True)

    plot_data(arr, arr[selection_indices])

    # Mask out all non-selected values
    mask = np.zeros(len(arr), bool)
    mask[selection_indices] = 1
    arr[:, SPD_INDEX] = mask * 10

    plot_data(arr, arr[selection_indices], text=selection_names)


def select_sources(arr, threshold=0.1, sig=SIG):
    selection_indices = []
    work_arr = np.copy(arr)
    while np.max(work_arr[:, SPD_INDEX]) > 0:
        idx = np.argmax(work_arr[:, SPD_INDEX])
        point = np.copy(work_arr[idx])
        if point[SPD_INDEX] < threshold:
            break

        selection_indices.append(idx)

        apply_gauss(work_arr, point, sig=sig)

        # Uncomment to find the order of selection.
    return selection_indices


def load_data():
    arr = np.loadtxt("phasediff_output.csv", delimiter=",", usecols=[1, 2, 3, 4], skiprows=1)
    names = np.loadtxt("phasediff_output.csv", delimiter=",", usecols=[0], skiprows=1, dtype=str)
    names = np.array([name.split("/", 1)[0] for name in names])
    arr[:, SPD_INDEX] = 1 / (arr[:, SPD_INDEX])
    return arr, names


if __name__ == "__main__":
    main()
