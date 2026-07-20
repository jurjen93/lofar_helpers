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
    wphase = (phase + np.pi) % (2 * np.pi) - np.pi
    return wphase

def make_plot(h5s, stations, soltab, names=None, outputname=None):
    """
    Make plot from multiple h5s
    """

    # Set the dimensions of the grid
    rows = len(names)
    cols = len(stations)

    # Create a figure and a set of subplots
    fig, axs = plt.subplots(rows, cols, figsize=(int(2*cols*2), int(2*rows)))  # figsize is adjustable to your needs

    # Iterate over each subplot to customize
    t = tables.open_file(h5s[0])
    for i, name in enumerate(names):
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
            vals = np.take(vals, indices=[i], axis=axes.index('pol'))
        except ValueError:
            pass

        for j, station in enumerate(stations):
            print(station, j)
            try:
                ref = np.take(vals, indices=[list(ants).index('CS001HBA0')], axis=axes.index('ant')).reshape(len(time), len(freqs))
            except:
                ref = np.take(vals, indices=[list(ants).index('ST001')], axis=axes.index('ant')).reshape(len(time), len(freqs))

            vals_im = np.take(vals, indices=[list(ants).index(make_utf8(station).split('\n')[0])], axis=axes.index('ant')).reshape(len(time), len(freqs))

            if 'amplitude' in soltab:
                vals_im = np.clip(vals_im, 0, 2)

            if 'phase' in soltab:
                vals_im = wrap_phase(vals_im-ref)
                vmin, vmax = -np.pi, np.pi
            else:
                vmin, vmax = 0, 1.5
                vals_im = vals_im

            if 'phase' in soltab:
                cmap = 'RdBu_r'
            else:
                cmap = 'Blues'
            im = axs[i, j].imshow(vals_im.T, aspect='auto', origin='lower', vmin=vmin, vmax=vmax, cmap=cmap)

            if i == 0:
                axs[i, j].set_title(make_utf8(station), size=26)
            # if j == 2:
            #     if names is not None:
            #         bbox_props = dict(boxstyle='round', facecolor='white', edgecolor='white', alpha=0.75)
            #         axs[i, j].text(vals_im.shape[0] // 21, vals_im.shape[1]//2, names[i], color='black',
            #                       fontsize=17, ha='left', va='bottom', bbox=bbox_props)
            if j == 0:
                if names is not None:
                    bbox_props = dict(boxstyle='round', facecolor='white', edgecolor='white', alpha=0.75)
                    axs[i, j].text(vals_im.shape[0] // 22, vals_im.shape[1] //2, names[i], color='black',
                                  fontsize=16, ha='left', va='bottom', bbox=bbox_props)
                # if i%2==0:
                #     axs[i, j].set_ylabel('Freq. [MHz]', size=16)
                if i%2!=0:
                    y_ticks = np.divide(np.linspace(freqs.min(), freqs.max(), num=3), 1000000).astype(int)
                    axs[i, j].set_yticks(ticks=np.linspace(axs[i, j].get_ylim()[0], axs[i, j].get_ylim()[1], num=3).astype(int), labels=y_ticks, fontsize=22)
                else:
                    axs[i, j].set_yticks([])
            else:
                axs[i, j].set_yticks([])
            if i == rows-1:
                # axs[i, j].set_xlabel('Time [hrs]', size=16)
                if j%2==0:
                    x_ticks = np.linspace(0, 8, num=3).round(0).astype(int)
                    axs[i, j].set_xticks(ticks=np.linspace(axs[i, j].get_xlim()[0], axs[i, j].get_xlim()[1], num=3).astype(int), labels=x_ticks, fontsize=22)
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
                       size=25)
    else:
        cbar.set_ticks([0, 0.5, 1, 1.5], labels=['0', '0.5', '1', '$\geq 1.5$'], size=22)


    if 'phase' in soltab:
        cbar.set_label('Phase correction', fontsize=24)
    if 'amplitude' in soltab:
        cbar.set_label('Amplitude correction', fontsize=24)

    fig.text(0.05, 0.5, 'Frequency [MHz]', va='center', rotation='vertical', fontsize=24)  # Adjust position (0.04, 0.5) and fontsize as needed
    fig.text(0.5, 0.03, 'Time [hrs]', va='center', rotation='horizontal', fontsize=24)  # Adjust position (0.04, 0.5) and fontsize as needed

    # fig.tight_layout()  # Adjust the layout to not overlap
    plt.savefig(outputname, dpi=200)

def main():
    stations = [b'CS002HBA0\n(Dutch core)', b'RS208HBA\n(Dutch remote)', b'DE604HBA\n(Germany)', b'SE607HBA\n(Sweden)', b'PL612HBA\n(Poland)']
    h5s = ['../merged_selfcalcycle011_linearfulljones_ILTJ174713.89+653235.8_delaycal.ms.copy.avg.h5']
    names = ['XX', 'XY', 'YX', 'YY']

    make_plot(h5s, stations, 'amplitude000', names, 'delay_amplitude_solutions.png')
    make_plot(h5s, stations, 'phase000', names, 'delay_phase_solutions.png')




if __name__ == '__main__':
    main()