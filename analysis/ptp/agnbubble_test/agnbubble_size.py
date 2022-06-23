import matplotlib.pyplot as plt
from glob import glob
import numpy as np

def pixel_to_arcsec(v):
    return v*80/27

bridge_rudnick = sorted(glob(f'../ptp_results_bridgebubble/bridgeresults_rudnick_*.txt'))
bridge_cb = sorted(glob(f'../ptp_results_bridgebubble/bridgeresults_cb_*.txt'))

print(bridge_cb)

bridge_results_cb, bridge_cb_err, bridge_pearson_cb, bridge_pearson_cb_err, bridge_spearman_cb, bridge_spearman_cb_err = [], [], [], [], [], []
bridge_results_rudnick, bridge_rudnick_err, bridge_pearson, bridge_pearson_err, bridge_spearman, bridge_spearman_err = [], [], [], [], [], []

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

print(np.array(bridge_pearson_cb)+np.array(bridge_pearson))

pearson = np.array(bridge_pearson_cb)+np.array(bridge_pearson)
pearson_err = np.array(bridge_pearson_cb_err)+np.array(bridge_pearson_err)
spearman = np.array(bridge_spearman_cb)+np.array(bridge_spearman)
spearman_err = np.array(bridge_spearman_cb_err)+np.array(bridge_spearman_err)
slope = np.array(bridge_results_cb)+np.array(bridge_results_rudnick)
slope_err = np.array(bridge_cb_err)+np.array(bridge_rudnick_err)

pearson/=2
pearson_err/=2
spearman/=2
spearman_err/=2
slope/=2
slope_err/=2

scales = [int(pixel_to_arcsec(v)) for v in list(range(33, 33+len(bridge_results_rudnick)*4, 4))]

# plt.errorbar(scales, bridge_results_rudnick, yerr=bridge_rudnick_err, fmt='.', ecolor='red', elinewidth=0.4,
#             color='darkred', capsize=2, markersize=4)
# plt.xlabel('Cell resolution [arcsec]', size=14)
# plt.ylabel('Slope', size=14)
# plt.show()
#
# plt.errorbar(scales, bridge_spearman, yerr=bridge_spearman_err, fmt='.', ecolor='green', elinewidth=0.4,
#             color='darkgreen', capsize=2, markersize=4)
# plt.xlabel('Cell resolution [arcsec]', size=14)
# plt.ylabel('Spearman R', size=14)
# plt.show()
#
# plt.errorbar(scales, bridge_pearson, yerr=bridge_pearson_err, fmt='.', ecolor='blue', elinewidth=0.4,
#             color='darkblue', capsize=2, markersize=4)
# plt.xlabel('Cell resolution [arcsec]', size=14)
# plt.ylabel('Pearson R', size=14)
# plt.show()

plt.errorbar(scales, slope, yerr=slope_err, fmt='.', ecolor='red', elinewidth=0.4,
            color='darkred', capsize=2, markersize=4)
plt.xlabel('Cell resolution [arcsec]', size=14)
plt.ylabel('Slope', size=14)
plt.show()

plt.errorbar(scales, spearman, yerr=spearman_err, fmt='.', ecolor='green', elinewidth=0.4,
            color='darkgreen', capsize=2, markersize=4)
plt.xlabel('Cell resolution [arcsec]', size=14)
plt.ylabel('Spearman R', size=14)
plt.show()

plt.errorbar(scales, pearson, yerr=pearson_err, fmt='.', ecolor='blue', elinewidth=0.4,
            color='darkblue', capsize=2, markersize=4)
plt.xlabel('Cell resolution [arcsec]', size=14)
plt.ylabel('Pearson R', size=14)
plt.show()