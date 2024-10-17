"""
This script has been written to plot h5 solutions similar to was done in de Jong et al. 2024
The matrix plotting style is similar to LoSoTo (https://github.com/revoltek/losoto/tree/master/losoto)
"""


import matplotlib.pyplot as plt
import tables
import numpy as np

# desired font
font_name = "Serif"
# update Matplotlib font configuration
plt.rcParams['font.family'] = font_name

def make_utf8(s):
    """
    Convert input to utf8 instead of bytes

    :param inp: string input
    """
    return s.decode('utf-8') if isinstance(s, bytes) else s

def wrap_phase(phase):
    """
    Wrap phases (stolen from https://github.com/tikk3r/lofar-h5plot/blob/master/h5plot#L1528)
    """
    wphase = (phase + np.pi) % (2 * np.pi) - np.pi
    return wphase

def make_plot(h5s, stations, soltab, names=None, outputname=None):
    """
    Make plot from multiple h5s
    """

    # Set the dimensions of the grid
    rows = len(h5s)
    cols = len(stations)

    # Create a figure and a set of subplots
    fig, axs = plt.subplots(rows, cols, figsize=(int(2*cols*2), int(2*rows)))  # figsize is adjustable to your needs

    # Iterate over each subplot to customize
    for i, h5 in enumerate(h5s):
        print(h5, i)
        t = tables.open_file(h5)
        freqs = t.root.sol000._f_get_child(soltab).freq[:]
        time = t.root.sol000._f_get_child(soltab).time[:]

        timespan = abs(time.min() - time.max())/3600

        vals = t.root.sol000._f_get_child(soltab).val
        ants = [make_utf8(s) for s in t.root.sol000.phase000.ant[:]]
        axes = make_utf8(vals.attrs["AXES"]).split(',')

        try:
            vals = np.take(vals[:], indices=[0], axis=axes.index('dir'))
        except ValueError:
            pass

        try:
            idx = 0
            if names is not None:
                if 'RL' in names[i]:
                    idx = 1
            vals = np.take(vals, indices=[idx], axis=axes.index('pol'))
        except ValueError:
            pass

        for j, station in enumerate(stations):
            print(station, j)
            try:
                ref = np.take(vals, indices=[list(ants).index('CS001HBA0')], axis=axes.index('ant')).reshape(len(time), len(freqs))
            except:
                ref = np.take(vals, indices=[list(ants).index('ST001')], axis=axes.index('ant')).reshape(len(time), len(freqs))

            vals_im = np.take(vals, indices=[list(ants).index(make_utf8(station))], axis=axes.index('ant')).reshape(len(time), len(freqs))

            if 'amplitude' in soltab:
                vals_im = np.clip(vals_im, 0, 2)
                cmap = 'RdBu_r'
            if 'phase' in soltab:
                vals_im = wrap_phase(vals_im-ref)
                # vals_im = wrap_phase(vals_im)
                vmin, vmax = -np.pi, np.pi
                cmap = 'twilight_shifted'
            else:
                vmin, vmax = 0, 2
                vals_im = vals_im

            im = axs[i, j].imshow(vals_im.T, aspect='auto', origin='lower', vmin=vmin, vmax=vmax, cmap=cmap)

            if i == 0:
                axs[i, j].set_title(make_utf8(station), size=23)
            # if j == 2:
            #     if names is not None:
            #         bbox_props = dict(boxstyle='round', facecolor='white', edgecolor='white', alpha=0.75)
            #         axs[i, j].text(vals_im.shape[0] // 21, vals_im.shape[1]//2, names[i], color='black',
            #                       fontsize=17, ha='left', va='bottom', bbox=bbox_props)
            if j == 0:
                if names is not None:
                    bbox_props = dict(boxstyle='round', facecolor='white', edgecolor='white', alpha=0.75)
                    axs[i, j].text(vals_im.shape[0] // 22, vals_im.shape[1] //2, names[i], color='black',
                                  fontsize=15, ha='left', va='bottom', bbox=bbox_props)
                # if i%2==0:
                #     axs[i, j].set_ylabel('Freq. [MHz]', size=16)
                if i%2!=0:
                    y_ticks = np.divide(np.linspace(freqs.min(), freqs.max(), num=3), 1000000).astype(int)
                    axs[i, j].set_yticks(ticks=np.linspace(axs[i, j].get_ylim()[0], axs[i, j].get_ylim()[1], num=3).astype(int), labels=y_ticks, fontsize=18)
                else:
                    axs[i, j].set_yticks([])
            else:
                axs[i, j].set_yticks([])
            if i == rows-1:
                # axs[i, j].set_xlabel('Time [hrs]', size=16)
                if j%2==0:
                    x_ticks = np.linspace(0, 8, num=3).round(0).astype(int)
                    axs[i, j].set_xticks(ticks=np.linspace(axs[i, j].get_xlim()[0], axs[i, j].get_xlim()[1], num=3).astype(int), labels=x_ticks, fontsize=18)
                else:
                    axs[i, j].set_xticks([])

            else:
                axs[i, j].set_xticks([])


    # Adjust layout to make room for the colorbar
    fig.subplots_adjust(wspace=0, hspace=0)
    # fig.tight_layout(rect=[0.05, 0.05, 0.95, 0.95])

    # Create colorbar
    cbar_ax = fig.add_axes([0.91, 0.14, 0.02, 0.7])  # x, y, width, height
    cbar = fig.colorbar(im, ax=axs.ravel().tolist(), cax=cbar_ax)
    if 'phase' in soltab:
        cbar.set_ticks(ticks=[-3.1415, -1.57075, 0, 1.57075, 3.1415], labels=['$-\pi$', '$-\pi$/2', '0', '$\pi$/2', '$\pi$'],
                       size=22)
    else:
        cbar.set_ticks([0, 0.5, 1, 1.5, 2], labels=['0', '0.5', '1', '1.5', '$\geq 2$'], size=22)


    if 'phase' in soltab:
        cbar.set_label('Phase correction', fontsize=22)
    if 'amplitude' in soltab:
        cbar.set_label('Amplitude correction', fontsize=22)

    fig.text(0.05, 0.5, 'Frequency [MHz]', va='center', rotation='vertical', fontsize=23)  # Adjust position (0.04, 0.5) and fontsize as needed
    fig.text(0.5, 0.022, 'Time [hrs]', va='center', rotation='horizontal', fontsize=23)  # Adjust position (0.04, 0.5) and fontsize as needed

    # fig.tight_layout()  # Adjust the layout to not overlap
    plt.savefig(outputname, dpi=200)

def main():

    # UPDATE FOR OWN NEED

    stations = [b'CS002HBA1', b'RS208HBA', b'DE604HBA', b'PL611HBA', b'IE613HBA']
    h5s = ['/project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/delayselfcal_new/scalarcomplexgain4_skyselfcalcyle000_L686962_120_168MHz_averaged.ms.avg.h5',
           '/project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/delayselfcal_new/fulljones5_skyselfcalcyle000_L686962_120_168MHz_averaged.ms.avg.h5',
           '/project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/delayselfcal_new/fulljones5_skyselfcalcyle000_L686962_120_168MHz_averaged.ms.avg.h5',
           '/project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/delayselfcal_new/scalarcomplexgain6_skyselfcalcyle000_L686962_120_168MHz_averaged.ms.avg.h5',
           '/project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/delayselfcal_new/merged_skyselfcalcyle000_L686962_120_168MHz_averaged.ms.avg.h5',
           '/project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/delayselfcal_new/merged_skyselfcalcyle000_L686962_120_168MHz_averaged.ms.avg.h5']
    names = ['scalarcomplexgain 1',
             'fulljones (RR)', 'fulljones (RL)', 'scalarcomplexgain 2', 'merged solutions (RR)',
             'merged solutions (RL)']

    make_plot(h5s, stations, 'amplitude000', names, 'delay_amplitude_solutions.png')

    stations = [b'CS002HBA1', b'RS208HBA', b'DE604HBA', b'PL611HBA', b'IE613HBA']
    h5s = ['/project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/delayselfcal_new/scalarphasediff0_skyselfcalcyle000_L686962_120_168MHz_averaged.ms.avg.h5',
           '/project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/delayselfcal_new/scalarphase1_skyselfcalcyle000_L686962_120_168MHz_averaged.ms.avg.h5',
           '/project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/delayselfcal_new/scalarphase2_skyselfcalcyle000_L686962_120_168MHz_averaged.ms.avg.h5',
           '/project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/delayselfcal_new/scalarphase3_skyselfcalcyle000_L686962_120_168MHz_averaged.ms.avg.h5',
           '/project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/delayselfcal_new/scalarcomplexgain4_skyselfcalcyle000_L686962_120_168MHz_averaged.ms.avg.h5',
           '/project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/delayselfcal_new/fulljones5_skyselfcalcyle000_L686962_120_168MHz_averaged.ms.avg.h5',
           '/project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/delayselfcal_new/fulljones5_skyselfcalcyle000_L686962_120_168MHz_averaged.ms.avg.h5',
           '/project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/delayselfcal_new/scalarcomplexgain6_skyselfcalcyle000_L686962_120_168MHz_averaged.ms.avg.h5',
           '/project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/delayselfcal_new/merged_skyselfcalcyle000_L686962_120_168MHz_averaged.ms.avg.h5',
           '/project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/delayselfcal_new/merged_skyselfcalcyle000_L686962_120_168MHz_averaged.ms.avg.h5']
    names = ['scalarphasediff (RR-LL)', 'scalarphase 1', 'scalarphase 2', 'scalarphase 3',
             'scalarcomplexgain 1', 'fulljones (RR)', 'fulljones (RL)', 'scalarcomplexgain 2', 'merged solutions (RR)',
             'merged solutions (RL)']

    make_plot(h5s, stations, 'phase000', names, 'delay_phase_solutions.png')

    stations = [b'RS208HBA', b'RS503HBA', b'DE604HBA', b'PL611HBA', b'IE613HBA']
    h5s = [ '/project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/new_ddcal/allselfcals/P22459/merged_selfcalcyle011_flagged_L686962_P22459.ms.copy.phaseup.h5',
            '/project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/new_ddcal/allselfcals/P35307/merged_selfcalcyle011_flagged_L686962_P35307.ms.copy.phaseup.h5',
           '/project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/new_ddcal/allselfcals/P19951/merged_selfcalcyle011_flagged_L686962_P19951.ms.copy.phaseup.h5']

    names = ['Facet 10', 'Facet 11', 'Facet 20']

    make_plot(h5s, stations, 'phase000', names, 'dd1_phase_solutions.png')
    make_plot(h5s, stations, 'amplitude000', names, 'dd1_amplitude_solutions.png')

    stations = [b'CS032HBA0', b'CS103HBA0', b'RS208HBA', b'RS307HBA', b'RS509HBA']
    h5s = ['/project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/imaging/split_facets2/facet_0/1.2imaging/selfcaloutput/merged_selfcalcyle009_concat_L68.ms.avg.h5',
           '/project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/imaging/split_facets2/facet_26/1.2imaging/selfcaloutput/merged_selfcalcyle009_concat_L68.ms.avg.h5',
           '/project/lofarvwf/Share/jdejong/output/ELAIS/ALL_L/imaging/split_facets2/facet_29/1.2imaging/selfcaloutput/merged_selfcalcyle009_concat_L68.ms.avg.h5']
    # names = ['merged solutions 1',
    #          'merged solutions 4',
    #          'merged solutions 7',
    #          'merged solutions 10']

    names = ['Facet 13', 'Facet 21', 'Facet 23']

    make_plot(h5s, stations, 'phase000', names, 'dd_dutch1_phase_solutions.png')

    names = ['Facet 13', 'Facet 23', 'Facet 21']

    make_plot(h5s, stations, 'amplitude000', names, 'dd_dutch1_amplitude_solutions.png')



if __name__ == '__main__':
    main()