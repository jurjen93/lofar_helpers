from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(
    description='Make montecarlo results')
parser._action_groups.pop()
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
# REQUIRED
required.add_argument('-cellsize', type=str, required=True)
args = parser.parse_args()

cellsize = args.cellsize

def uncertainty(d_err):
    unc = 1/np.power(d_err, 2)
    return np.sqrt(len(unc)/np.sum(unc))

def weighted_mean(values, err):
    return np.sum(np.divide(values, np.power(err, 2))) / np.sum(1 / np.power(err, 2))

def bootstrap(values, err):
    res = []
    for i in zip(values, err):
        for n in range(len(values)):
            res.append(np.random.normal(i[0], i[1], 1))

    return np.mean(res), np.std(res)

a399trail_rudnick = glob(f'ptp_results_{cellsize}/A399trailresults_rudnick_*.txt')
a399_rudnick = glob(f'ptp_results_{cellsize}/a399results_rudnick_*.txt')
a401_rudnick = glob(f'ptp_results_{cellsize}/a401results_rudnick_*.txt')
bridge_rudnick = glob(f'ptp_results_{cellsize}/bridgeresults_rudnick_*.txt')
a399trail_cb = glob(f'ptp_results_{cellsize}/A399trailresults_cb_*.txt')
a399_cb = glob(f'ptp_results_{cellsize}/a399results_cb_*.txt')
a401_cb = glob(f'ptp_results_{cellsize}/a401results_cb_*.txt')
bridge_cb = glob(f'ptp_results_{cellsize}/bridgeresults_cb_*.txt')

a399_results_rudnick, a399_rudnick_err, a399_pearson, a399_pearson_err, a399_spearman, a399_spearman_err = [], [], [], [], [], []
a399trail_results_rudnick, a399trail_rudnick_err, a399trail_pearson, a399trail_pearson_err, a399trail_spearman, a399trail_spearman_err = [], [], [], [], [], []
a401_results_rudnick, a401_rudnick_err, a401_pearson, a401_pearson_err, a401_spearman, a401_spearman_err = [], [], [], [], [], []
bridge_results_rudnick, bridge_rudnick_err, bridge_pearson, bridge_pearson_err, bridge_spearman, bridge_spearman_err = [], [], [], [], [], []
a399_results_cb, a399_cb_err, a399_pearson_cb, a399_pearson_cb_err, a399_spearman_cb, a399_spearman_cb_err = [], [], [], [], [], []
a399trail_results_cb, a399trail_cb_err, a399trail_pearson_cb, a399trail_pearson_cb_err, a399trail_spearman_cb, a399trail_spearman_cb_err = [], [], [], [], [], []
a401_results_cb, a401_cb_err, a401_pearson_cb, a401_pearson_cb_err, a401_spearman_cb, a401_spearman_cb_err = [], [], [], [], [], []
bridge_results_cb, bridge_cb_err, bridge_pearson_cb, bridge_pearson_cb_err, bridge_spearman_cb, bridge_spearman_cb_err = [], [], [], [], [], []

for a339 in a399_cb:
    with open(a339) as f:
        try:
            lines = f.readlines()
            a399_pearson_cb.append(float(lines[-4].split()[0].replace('(','').replace(',','')))
            a399_pearson_cb_err.append(float(lines[-3].split()[-1]))
            a399_spearman_cb.append(float(lines[-2].split()[0].replace('(','').replace(',','')))
            a399_spearman_cb_err.append(float(lines[-1].split()[-1]))
            a399_results_cb.append(float(lines[3].replace('Linear regression slope is ', '').split()[0]))
            a399_cb_err.append(float(lines[3].replace('Linear regression slope is ', '').split()[-1]))
        except:
            pass

for a339 in a399trail_cb:
    with open(a339) as f:
        try:
            lines = f.readlines()
            a399trail_pearson_cb.append(float(lines[-4].split()[0].replace('(','').replace(',','')))
            a399trail_pearson_cb_err.append(float(lines[-3].split()[-1]))
            a399trail_spearman_cb.append(float(lines[-2].split()[0].replace('(','').replace(',','')))
            a399trail_spearman_cb_err.append(float(lines[-1].split()[-1]))
            a399trail_results_cb.append(float(lines[3].replace('Linear regression slope is ', '').split()[0]))
            a399trail_cb_err.append(float(lines[3].replace('Linear regression slope is ', '').split()[-1]))
        except:
            pass

for a401 in a401_cb:
    with open(a401) as f:
        try:
            lines = f.readlines()
            a401_pearson_cb.append(float(lines[-4].split()[0].replace('(','').replace(',','')))
            a401_pearson_cb_err.append(float(lines[-3].split()[-1]))
            a401_spearman_cb.append(float(lines[-2].split()[0].replace('(','').replace(',','')))
            a401_spearman_cb_err.append(float(lines[-1].split()[-1]))
            a401_results_cb.append(float(lines[3].replace('Linear regression slope is ', '').split()[0]))
            a401_cb_err.append(float(lines[3].replace('Linear regression slope is ', '').split()[-1]))
        except:
            pass


for bridge in bridge_cb:
    with open(bridge) as f:
        try:
            lines = f.readlines()
            bridge_pearson_cb.append(float(lines[-4].split()[0].replace('(','').replace(',','')))
            bridge_pearson_cb_err.append(float(lines[-3].split()[-1]))
            bridge_spearman_cb.append(float(lines[-2].split()[0].replace('(','').replace(',','')))
            bridge_spearman_cb_err.append(float(lines[-1].split()[-1]))
            bridge_results_cb.append(float(lines[3].replace('Linear regression slope is ', '').split()[0]))
            bridge_cb_err.append(float(lines[3].replace('Linear regression slope is ', '').split()[-1]))
        except:
            pass


for a339rudnick in a399_rudnick:
    with open(a339rudnick) as f:
        try:
            lines = f.readlines()
            a399_pearson.append(float(lines[-4].split()[0].replace('(','').replace(',','')))
            a399_pearson_err.append(float(lines[-3].split()[-1]))
            a399_spearman.append(float(lines[-2].split()[0].replace('(','').replace(',','')))
            a399_spearman_err.append(float(lines[-1].split()[-1]))
            a399_results_rudnick.append(float(lines[3].replace('Linear regression slope is ', '').split()[0]))
            a399_rudnick_err.append(float(lines[3].replace('Linear regression slope is ', '').split()[-1]))
        except:
            pass

for a339rudnick in a399trail_rudnick:
    with open(a339rudnick) as f:
        try:
            lines = f.readlines()
            a399trail_pearson.append(float(lines[-4].split()[0].replace('(','').replace(',','')))
            a399trail_pearson_err.append(float(lines[-3].split()[-1]))
            a399trail_spearman.append(float(lines[-2].split()[0].replace('(','').replace(',','')))
            a399trail_spearman_err.append(float(lines[-1].split()[-1]))
            a399trail_results_rudnick.append(float(lines[3].replace('Linear regression slope is ', '').split()[0]))
            a399trail_rudnick_err.append(float(lines[3].replace('Linear regression slope is ', '').split()[-1]))
        except:
            pass

for a401rudnick in a401_rudnick:
    with open(a401rudnick) as f:
        try:
            lines = f.readlines()
            a401_pearson.append(float(lines[-4].split()[0].replace('(','').replace(',','')))
            a401_pearson_err.append(float(lines[-1].split()[-1]))
            a401_spearman.append(float(lines[-2].split()[0].replace('(','').replace(',','')))
            a401_spearman_err.append(float(lines[-1].split()[-1]))
            a401_results_rudnick.append(float(lines[3].replace('Linear regression slope is ', '').split()[0]))
            a401_rudnick_err.append(float(lines[3].replace('Linear regression slope is ', '').split()[-1]))
        except:
            pass

for bridgerudnick in bridge_rudnick:
    with open(bridgerudnick) as f:
        try:
            lines = f.readlines()
            bridge_pearson.append(float(lines[-4].split()[0].replace('(','').replace(',','')))
            bridge_pearson_err.append(float(lines[-1].split()[-1]))
            bridge_spearman.append(float(lines[-2].split()[0].replace('(','').replace(',','')))
            bridge_spearman_err.append(float(lines[-1].split()[-1]))
            bridge_results_rudnick.append(float(lines[3].replace('Linear regression slope is ', '').split()[0]))
            bridge_rudnick_err.append(float(lines[3].replace('Linear regression slope is ', '').split()[-1]))
        except:
            pass


print('\nTotal:\n-----------------------------------------')
a399_tot = a399_results_rudnick+a399_results_cb
a399_err = a399_rudnick_err+a399_cb_err
a399_pearson_tot = a399_pearson+a399_pearson_cb
a399_pearson_err = a399_pearson_err + a399_pearson_cb_err
a399_spearman_tot = a399_spearman+a399_spearman_cb
a399_spearman_err = a399_spearman_err + a399_spearman_cb_err
print(f'A399\nSlope: {bootstrap(a399_tot, a399_err)}')
print(f'Pearson R: {bootstrap(a399_pearson_tot, a399_pearson_err)}')
print(f'Spearman R: {bootstrap(a399_spearman_tot, a399_spearman_err)}\n')

a399trail_tot = a399trail_results_rudnick+a399trail_results_cb
a399trail_err = a399trail_rudnick_err+a399trail_cb_err
a399trail_pearson_tot = a399trail_pearson+a399trail_pearson_cb
a399trail_pearson_err = a399trail_pearson_err + a399trail_pearson_cb_err
a399trail_spearman_tot = a399trail_spearman+a399trail_spearman_cb
a399trail_spearman_err = a399trail_spearman_err + a399trail_spearman_cb_err
print(f'A399 trail\nSlope: {bootstrap(a399trail_tot, a399trail_err)}')
print(f'Pearson R: {bootstrap(a399trail_pearson_tot, a399trail_pearson_err)}')
print(f'Spearman R: {bootstrap(a399trail_spearman_tot, a399trail_spearman_err)}\n')

a399_tot += a399trail_tot
a399_err += a399trail_err
a399_pearson_tot += a399trail_pearson_tot
a399_pearson_err += a399trail_pearson_err
a399_spearman_tot += a399trail_spearman_tot
a399_spearman_err += a399trail_spearman_err
print(f'A399 extended\nSlope: {bootstrap(a399_tot, a399_err)}')
print(f'Pearson R: {bootstrap(a399_pearson_tot, a399_pearson_err)}')
print(f'Spearman R: {bootstrap(a399_spearman_tot, a399_spearman_err)}\n')

a401_tot = a401_results_rudnick+a401_results_cb
a401_err = a401_rudnick_err+a401_cb_err
a401_pearson_tot = a401_pearson+a401_pearson_cb
a401_pearson_err = a401_pearson_err + a401_pearson_cb_err
a401_spearman_tot = a401_spearman+a401_spearman_cb
a401_spearman_err = a401_spearman_err + a401_spearman_cb_err
print(f'A401\nSlope: {bootstrap(a401_tot, a401_err)}')
print(f'Pearson R: {bootstrap(a401_pearson_tot, a401_pearson_err)}')
print(f'Spearman R: {bootstrap(a401_spearman_tot, a401_spearman_err)}\n')

bridge_tot = bridge_results_rudnick+bridge_results_cb
bridge_err = bridge_rudnick_err+bridge_cb_err
bridge_pearson_tot = bridge_pearson+bridge_pearson_cb
bridge_pearson_err = bridge_pearson_err + bridge_pearson_cb_err
bridge_spearman_tot = bridge_spearman+bridge_spearman_cb
bridge_spearman_err = bridge_spearman_err + bridge_spearman_cb_err
print(f'Bridge\nSlope: {bootstrap(bridge_tot, bridge_err)}')
print(f'Pearson R: {bootstrap(bridge_pearson_tot, bridge_pearson_err)}')
print(f'Spearman R: {bootstrap(bridge_spearman_tot, bridge_spearman_err)}')


# plt.hist(bridge_results_rudnick+bridge_results_cb)
# plt.title("Bridge")
# plt.show()
#
# plt.hist(a399_results_rudnick+a399_results_cb)
# plt.title("A399")
# plt.show()
#
# plt.hist(a401_results_rudnick+a401_results_cb)
# plt.title("A401")
# plt.show()