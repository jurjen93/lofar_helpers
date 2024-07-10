from casacore.tables import table
import sys
import matplotlib.pyplot as plt
from glob import glob
import os

# Set the MPLCONFIGDIR environment variable
os.system('mkdir -p ~/matplotlib_cache')
os.environ['MPLCONFIGDIR'] = os.path.expanduser('~/matplotlib_cache')


def get_station_id(ms):
    """
    Get station with corresponding id number

    :param:
        - ms: measurement set

    :return:
        - antenna names, IDs
    """

    t = table(ms+'::ANTENNA', ack=False)
    ants = t.getcol("NAME")
    t.close()

    t = table(ms+'::FEED', ack=False)
    ids = t.getcol("ANTENNA_ID")
    t.close()

    return ants, ids


def plot_baseline_track(t_final_name: str = None, t_input_names: list = None, baseline='0-1', UV=True, saveas=None):
    """
    Plot baseline track

    :param:
        - t_final_name: table with final name
        - t_input_names: tables to compare with
        - mappingfiles: baseline mapping files
    """

    if len(t_input_names) > 4:
        sys.exit("ERROR: Can just plot 4 inputs")

    colors = ['red', 'green', 'yellow', 'black']

    if not UV:
        print("MAKE UW PLOT")

    ant1, ant2 = baseline.split('-')

    for n, t_input_name in enumerate(t_input_names):

        ref_stats, ref_ids = get_station_id(t_final_name)
        new_stats, new_ids = get_station_id(t_input_name)

        id_map = dict(zip([ref_stats.index(a) for a in new_stats], new_ids))

        with table(t_final_name, ack=False) as f:
            fsub = f.query(f'ANTENNA1={ant1} AND ANTENNA2={ant2} AND NOT ALL(WEIGHT_SPECTRUM == 0)', columns='UVW')
            uvw1 = fsub.getcol("UVW")

        with table(t_input_name, ack=False) as f:
            fsub = f.query(f'ANTENNA1={id_map[int(ant1)]} AND ANTENNA2={id_map[int(ant2)]} AND NOT ALL(WEIGHT_SPECTRUM == 0)', columns='UVW')
            uvw2 = fsub.getcol("UVW")

        # Scatter plot for uvw1
        if n == 0:
            lbl = 'Final MS'
        else:
            lbl = None

        plt.scatter(uvw1[:, 0], uvw1[:, 2] if UV else uvw1[:, 3], label=lbl, color='blue', edgecolor='black', alpha=0.2, s=100, marker='o')

        # Scatter plot for uvw2
        plt.scatter(uvw2[:, 0], uvw2[:, 2] if UV else uvw2[:, 3], label=f'Dataset {n}', color=colors[n], edgecolor='black', alpha=0.7, s=40, marker='*')


    # Adding labels and title
    plt.xlabel("U (m)", fontsize=14)
    plt.ylabel("V (m)" if UV else "W (m)", fontsize=14)

    # Adding grid
    plt.grid(True, linestyle='--', alpha=0.6)

    # Adding legend
    plt.legend(fontsize=12)

    plt.xlim(561260, 563520)
    plt.ylim(192782, 194622)

    plt.tight_layout()

    if saveas is None:
        plt.show()
    else:
        plt.savefig(saveas, dpi=150)
        plt.close()

plot_baseline_track('test_1.ms', sorted(glob('a*.ms')), '0-71', saveas='test_1.png')
plot_baseline_track('test_2.ms', sorted(glob('a*.ms')), '0-71', saveas='test_2.png')
plot_baseline_track('test_3.ms', sorted(glob('a*.ms')), '0-71', saveas='test_3.png')
plot_baseline_track('test_4.ms', sorted(glob('a*.ms')), '0-71', saveas='test_4.png')
plot_baseline_track('test_6.ms', sorted(glob('a*.ms')), '0-71', saveas='test_6.png')
plot_baseline_track('test_8.ms', sorted(glob('a*.ms')), '0-71', saveas='test_8.png')