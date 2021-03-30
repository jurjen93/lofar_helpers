"""
After one has created solution files from self calling on the extracted boxes,
one can use this script to merge solution files
Do this by importing the function merge_h5 from this script:
-------------------------
EXAMPLE:

from h5_merger import merge_h5
merge_h5(h5_out='test.h5',
        h5_files='*.h5',
        ms_files='*.ms',
        convert_tec=True)
"""

import os
import casacore.tables as ct
from glob import glob
from losoto.h5parm import h5parm
from losoto.lib_operations import reorderAxes
import numpy as np
from scipy.interpolate import interp1d

__author__ = "Jurjen de Jong"


# TODO: Weights keep on 1 -> future investigation
# TODO: 8) test with DPPP command to process h5 file to see if same image comes out in [applycal from runwscleanLBautoR.py]
# TODO: 9) Boxes -> Euclidean distance van de fluxes
# TODO: rotation000 rotates polarization -> check leakage

class MergeH5:
    """Merge multiple h5 files"""

    def __init__(self, h5_out, h5_files=None, ms_files=None, convert_tec=True):
        """
        :param h5_out: name of merged output h5 file
        :param files: h5 files to merge, can be both list and string
        :param ms_files: ms files to use, can be both list and string
        :param convert_tec: convert TEC to phase or not
        """

        self.file = h5_out

        if type(ms_files) == list:
            ms = ms_files
        elif type(ms_files) == str:
            ms = glob(ms_files)
        else:
            ms = glob('*dysco.sub.shift.avg.weights.ms.archive0.goodtimes')

        if type(h5_files) == list:
            self.h5_files = h5_files
        elif type(h5_files) == str:
            self.h5_files = glob(h5_files)
        else:
            max_num_h5 = sorted(glob('*ms.archive0*.avgsolsgrid*.h5'))[-1].split('.')[-2].split('_')[-1]
            self.h5_files = glob('*ms.archive0*.avgsolsgrid_%s.h5' % max_num_h5)

        t = ct.taql('SELECT CHAN_FREQ, CHAN_WIDTH FROM ' + ms[0] + '::SPECTRAL_WINDOW')
        self.ax_freq = t.getcol('CHAN_FREQ')[0]
        t.close()

        t = ct.table(ms[0])
        self.ax_time = sorted(np.unique(t.getcol('TIME')))
        t.close()

        self.from1to2pol = False  # Add polarization axis for tec or phase solutions
        self.convert_tec = convert_tec  # do not convert tec

        self.solaxnames = ['pol', 'dir', 'ant', 'freq', 'time']

    def get_values(self, ss, st, solset, soltab):

        if 'pol' in st.getAxesNames():
            print("polarization is in {solset}/{soltab}".format(solset=solset, soltab=soltab))
        else:
            print("polarization is not in {solset}/{soltab}".format(solset=solset, soltab=soltab))

        time_axes = st.getAxisValues('time')
        freq_axes = st.getAxisValues('freq')

        if self.ax_time[0] > time_axes[0] or self.ax_time[-1] < time_axes[-1]:
            print("Time axes between h5 and ms are not overlapping, so we reduce time_axis")
        if self.ax_freq[0] > freq_axes[0] or self.ax_freq[-1] < freq_axes[-1]:
            print("Frequency axes between h5 and ms are not overlapping, so we reduce frequency axis")
        if float(soltab[-3:]) > 0:
            print("Error: {soltab} does not end on 000".format(soltab=soltab))
        for av in self.axes_new:
            if av in st.getAxesNames() and st.getAxisLen(av) == 0:
                print("No {av} in {solset}/{soltab}".format(av=av, solset=solset, soltab=soltab))
        time_axes_new = [list(time_axes).index(t) for t in time_axes if t < self.ax_time[-1] and t > self.ax_time[0]]
        freq_axes_new = [list(freq_axes).index(f) for f in freq_axes if f < self.ax_freq[-1] and f > self.ax_freq[0]]
        print('Value shape before reduction {values}'.format(values=st.getValues()[0].shape))
        values = reorderAxes(st.getValues()[0], st.getAxesNames(), self.axes_current)[...,
                 freq_axes_new[0]:freq_axes_new[-1] + 1, time_axes_new[0]:time_axes_new[-1] + 1]
        if len(values.shape) == 5:
            time_axes = values[0, 0, 0, 0, :]
            freq_axes = values[0, 0, 0, :, 0]
        elif len(values.shape) == 4:
            time_axes = values[0, 0, 0, :]
            freq_axes = values[0, 0, :, 0]
        elif len(values.shape) == 3:
            time_axes = values[0, 0, :]
            freq_axes = values[0, :, 0]
        else:
            print('ERROR: missing axes?')
        return values, time_axes, freq_axes

    def sort_soltabs(self, soltabs):
        """
        Sort solution tables
        :param soltabs: solutions tables
        """
        soltabs = set(soltabs)
        if self.convert_tec:
            tp_phasetec = [li for li in soltabs if 'tec' in li or 'phase' in li]
            tp_amplitude = [li for li in soltabs if 'amplitude' in li]
            tp_rotation = [li for li in soltabs if 'rotation' in li]
            return [sorted(tp_amplitude, key=lambda x: float(x[-3:])),
                    sorted(tp_rotation, key=lambda x: float(x[-3:])),
                    sorted(sorted(tp_phasetec), key=lambda x: float(x[-3:]))]
        else:
            tp_phase = [li for li in soltabs if 'phase' in li]
            tp_tec = [li for li in soltabs if 'tec' in li]
            tp_amplitude = [li for li in soltabs if 'amplitude' in li]
            tp_rotation = [li for li in soltabs if 'rotation' in li]
            return [sorted(tp_phase, key=lambda x: float(x[-3:])),
                    sorted(tp_tec, key=lambda x: float(x[-3:])),
                    sorted(tp_amplitude, key=lambda x: float(x[-3:])),
                    sorted(tp_rotation, key=lambda x: float(x[-3:]))]

    @property
    def get_allkeys(self):
        """
        Get all solution sets, solutions tables, and ax names
        """
        self.all_soltabs, self.all_solsets, self.all_axes, self.antennas = [], [], [], []

        for h5_name in self.h5_files:
            h5 = h5parm(h5_name)
            for solset in h5.getSolsetNames():
                self.all_solsets += [solset]
                ss = h5.getSolset(solset)
                for soltab in ss.getSoltabNames():
                    self.all_soltabs += [soltab]
                    st = ss.getSoltab(soltab)
                    self.all_axes += ['/'.join([solset, soltab, an]) for an in st.getAxesNames()]
                    self.antennas += list(st.getAxisValues('ant'))
            h5.close()
        self.all_soltabs = self.sort_soltabs(self.all_soltabs)
        self.all_solsets = set(self.all_solsets)
        self.all_axes = set(self.all_axes)
        self.antennas = list(set(self.antennas))
        return self

    def get_clean_values(self, soltab, st):
        """
        Get default values, based on model h5 file
        :param soltab: solution table name
        :param st: solution table itself
        """
        if 'pol' in st.getAxesNames():
            self.polarizations = st.getAxisValues('pol')
        else:
            self.polarizations = []  # always five dimensions

        if 'amplitude' in soltab and 'pol' in st.getAxesNames():
            self.gains = np.ones(
                (len(self.polarizations), 1, len(self.antennas), len(self.ax_freq), len(self.ax_time)))
        else:
            self.gains = np.ones((1, len(self.antennas), len(self.ax_freq), len(self.ax_time)))

        if 'phase' in soltab and 'pol' in st.getAxesNames():
            self.phases = np.zeros(
                (max(len(self.polarizations), 1), 1, len(self.antennas), len(self.ax_freq), len(self.ax_time)))

        elif 'rotation' in soltab:
            self.phases = np.zeros((1, len(self.antennas), len(self.ax_freq), len(self.ax_time)))

        elif 'tec' in soltab:
            self.phases = np.zeros((1, 1, len(self.antennas), len(self.ax_freq), len(self.ax_time)))
        else:
            self.phases = np.zeros((1, len(self.antennas), len(self.ax_freq), len(self.ax_time)))

        self.directions = {}  # directions in a dictionary
        self.n = 0  # direction number back to zero

        return self

    @staticmethod
    def tecphase_conver(tec, freqs):
        """
        convert tec to phase
        :param tec: TEC
        :param freqs: frequencies
        """
        return -8.4479745e9 * tec / freqs

    @staticmethod
    def interp_along_axis(x, interp_from, interp_to, axis):
        """
        Interpolate along axis
        """
        interp_vals = interp1d(interp_from, x, axis=axis, kind='nearest', fill_value='extrapolate')
        new_vals = interp_vals(interp_to)
        return new_vals

    def get_model_h5(self, solset, soltab):
        """
        Get model (clean) h5 file
        :param solset: solution set name (sol000, sol001,..)
        :param soltab: solution table name
        """
        h5_N = 0

        if '000' in soltab:

            for h5_name_to_merge in self.h5_files:  # make template

                h5_to_merge = h5parm(h5_name_to_merge)

                if solset not in h5_to_merge.getSolsetNames():
                    h5_to_merge.close()
                    continue  # use other h5 file with solset
                else:
                    ss = h5_to_merge.getSolset(solset)
                if soltab not in ss.getSoltabNames():
                    h5_to_merge.close()
                    continue  # use other h5 file with soltab
                else:
                    st = ss.getSoltab(soltab)

                # if '/'.join([solset, soltab, 'pol']) in self.all_axes and 'pol' not in st.getAxesNames():
                # continue  # use h5 which has the polarizations as model file

                # make new clean values for new soltabs (ending on 000)
                if not self.convert_tec or (self.convert_tec and 'tec' not in soltab):
                    self.get_clean_values(soltab, st)
                    self.axes_new = [an for an in self.solaxnames if an in st.getAxesNames()]

                h5_to_merge.close()
                break

    def get_sol(self, solset, soltab):
        """
        Get solutions merged
        :param solset: solution set name
        :param soltab: solution table name
        """
        for h5_name in self.h5_files:

            h5 = h5parm(h5_name)
            if solset not in h5.getSolsetNames():
                h5.close()
                continue

            ss = h5.getSolset(solset)

            if soltab not in ss.getSoltabNames():
                h5.close()
                continue

            st = ss.getSoltab(soltab)

            self.axes_current = [an for an in self.solaxnames if an in st.getAxesNames()]

            values, time_axes, freq_axes = self.get_values(ss, st, solset, soltab)

            d = ss.getSou()
            source_coords = d[d.keys()[0]]
            d = 'Dir{:02d}'.format(self.n)

            if any([np.array_equal(source_coords, list(sv)) for sv in self.directions.values()]):
                # Direction already exists, add to the existing solutions.
                idx = list([list(l) for l in self.directions.values()]).index(list(source_coords))
            else:
                print('Adding new direction {:f},{:f}'.format(*source_coords))
                idx = self.n
                self.directions.update({d: source_coords})
                self.n += 1

            if st.getType() == 'tec':
                if self.convert_tec:
                    # Convert tec to phase.
                    if len(self.polarizations) > 0 and len(self.phases.shape) == 5:
                        valtmp = np.ones((len(self.polarizations),) + values.shape)
                        valtmp[0, ...] = values
                        valtmp[-1, ...] = values
                        values = valtmp
                        # -1 assumes the expected shape along the frequency axis.
                        if self.axes_new[-2] != 'freq':
                            print('Frequency axis is not on right position')
                        freqs = self.ax_freq.reshape(1, 1, 1, -1, 1)
                        tecphase = self.tecphase_conver(values, freqs)
                        tp = self.interp_along_axis(tecphase, time_axes, self.ax_time,
                                                    self.axes_new.index('time'))
                    elif len(self.phases.shape) == 4:
                        freqs = self.ax_freq.reshape(1, 1, -1, 1)
                        tecphase = self.tecphase_conver(values, freqs)
                        tp = self.interp_along_axis(tecphase, time_axes, self.ax_time,
                                                    self.axes_current.index('time'))
                    elif len(self.phases.shape) == 5:
                        freqs = self.ax_freq.reshape(1, 1, 1, -1, 1)
                        tecphase = self.tecphase_conver(values, freqs)
                        tp = self.interp_along_axis(tecphase, time_axes, self.ax_time,
                                                    self.axes_current.index('time'))
                    else:
                        print('ERROR: something went wrong with reshaping')

                    # Now add the phases to the total phase correction for this direction.
                    if len(self.phases.shape) == 4:
                        self.phases[idx, :, :, :] += tp[0, :, :, :]
                        phasetmp = np.zeros((2,) + self.phases.shape[:])
                        phasetmp[0, ...] = self.phases
                        phasetmp[-1, ...] = self.phases
                        self.phases = phasetmp
                        self.axes_new = ['pol'] + self.axes_new

                    elif len(self.phases.shape) == 5:
                        self.phases[:, idx, :, :, :] += tp[:, 0, :, :, :]

                else:
                    if 'dir' in self.axes_current:
                        values = values[0, :, 0, :]

                    tp = self.interp_along_axis(values, time_axes, self.ax_time, -1)
                    tp = tp.reshape((1, tp.shape[0], 1, tp.shape[1]))
                    # Now add the tecs to the total phase correction for this direction.
                    if 'dir' in self.axes_current:
                        self.phases[idx, :, :, :] += tp[0, ...]
                    else:
                        self.phases[idx, :, :] += tp


            elif st.getType() == 'phase' or st.getType() == 'rotation':
                if 'pol' in self.axes_current:
                    if st.getAxisLen('pol') == 4:
                        print("Add fulljones type with 4 polarizations")
                        print("Please check if this is correct")
                        if self.phases.shape[0] == 2:
                            phasetmp = np.zeros((4,) + self.phases.shape[1:])
                            phasetmp[0, ...] = self.phases[0, ...]
                            phasetmp[-1, ...] = self.phases[1, ...]
                            self.phases = phasetmp
                        elif len(self.phases.shape) < 5:
                            phasetmp = np.zeros((4,) + self.phases.shape)
                            phasetmp[0, ...] = self.phases
                            phasetmp[-1, ...] = self.phases
                            self.phases = phasetmp
                    elif st.getAxisLen('pol') == 2 and self.phases.shape[0] == 4:
                        print("Add to fulljones type with 4 polarizations")
                        phasetmp = np.zeros((4,) + values.shape[1:])
                        phasetmp[0, ...] = values[0, ...]
                        phasetmp[-1, ...] = values[1, ...]
                        values = phasetmp
                elif self.phases.shape[0] == 4 and len(self.phases.shape) == 5:
                    print("Add to fulljones type with 4 polarizations")
                    phasetmp = np.zeros((4,) + values.shape)
                    phasetmp[0, ...] = values
                    phasetmp[-1, ...] = values
                    values = phasetmp

                idxnan = np.where((~np.isfinite(values)))
                values[idxnan] = 0.0
                tp = self.interp_along_axis(values, time_axes, self.ax_time,
                                            self.axes_current.index('time'))
                if tp.shape[-2] == 1:
                    tptmp = tp
                    for ff in self.ax_freq[:-1]:
                        tp = np.append(tp, tptmp, axis=-2)
                else:
                    tp = self.interp_along_axis(tp, freq_axes, self.ax_freq,
                                                self.axes_current.index('freq'))

                if 'dir' in self.axes_current:
                    if 'pol' in self.axes_current and len(tp.shape) == len(self.phases.shape):
                        self.phases[:, idx, :, :, :] += tp[:, 0, ...]
                    elif 'pol' in self.axes_current and len(tp.shape) - len(self.phases.shape) == 1:
                        self.phases[idx, :, :, :] += tp[0, 0, ...]
                    else:
                        self.phases[idx, :, :, :] += tp[0, ...]
                else:
                    self.phases += tp

            elif st.getType() == 'amplitude':

                if 'pol' in self.axes_current:
                    if st.getAxisLen('pol') == 4:
                        print("Add fulljones type with 4 polarizations")
                        print("Please check if this is correct")
                        if self.gains.shape[0] == 2:
                            gaintmp = np.zeros((4,) + self.gains.shape[1:])
                            gaintmp[0, ...] = self.gains[0, ...]
                            gaintmp[-1, ...] = self.gains[1, ...]
                            self.gains = gaintmp
                        elif self.gains.shape < 5:
                            gaintmp = np.zeros((4,) + self.gains.shape)
                            gaintmp[0, ...] = self.gains
                            gaintmp[-1, ...] = self.gains
                            self.gains = gaintmp
                    elif st.getAxisLen('pol') == 2 and self.gains.shape[0] == 4:
                        print("Add to fulljones type with 4 polarizations")
                        gaintmp = np.zeros((4,) + values.shape[1:])
                        gaintmp[0, ...] = values[0, ...]
                        gaintmp[-1, ...] = values[1, ...]
                        values = gaintmp
                elif self.gains.shape[0] == 4 and len(self.gains.shape) == 5:
                    print("Add to fulljones type with 4 polarizations")
                    gaintmp = np.zeros((4,) + values.shape)
                    gaintmp[0, ...] = values
                    gaintmp[-1, ...] = values
                    gaintmp = values

                idxnan = np.where((~np.isfinite(values)))
                values[idxnan] = 1.0
                tp = self.interp_along_axis(values, time_axes, self.ax_time,
                                            self.axes_current.index('time'))

                if tp.shape[-2] == 1:
                    tptmp = tp
                    for ff in self.ax_freq[:-1]:
                        tp = np.append(tp, tptmp, axis=-2)
                else:
                    tp = self.interp_along_axis(tp, freq_axes, self.ax_freq,
                                                self.axes_current.index('freq'))

                if 'dir' in self.axes_current:
                    if 'pol' in self.axes_current:
                        self.gains[:, idx, :, :, :] *= tp[:, 0, ...]
                    else:
                        self.gains[idx, :, :, :] *= tp[0, ...]
                else:
                    self.gains *= tp

            h5.close()

        return self

    def DPPP_style(self, soltab):
        if 'pol' in self.axes_new and len(self.axes_new) == 5:
            DPPP_axes = ['time', 'freq', 'ant', 'dir', 'pol']
        elif 'pol' not in self.axes_new and len(self.axes_new) == 4:
            DPPP_axes = ['time', 'ant', 'dir', 'freq']
            if len(self.phases.shape) == 5:
                self.phases = self.phases[0]
        else:
            DPPP_axes = []

        if 'phase' in soltab or 'tec' in soltab or 'rotation' in soltab:
            self.phases = reorderAxes(self.phases, self.axes_new, DPPP_axes)
        elif 'amplitude' in soltab:
            self.gains = reorderAxes(self.gains, self.axes_new, DPPP_axes)

        return DPPP_axes

    def create_new_dataset(self, solset, soltab):
        """
        Create a new dataset in the h5 file
        :param solset: solution set name
        :param soltab: solution table name
        """
        if len(self.directions.keys()) == 0:
            return self

        self.h5_out = h5parm(self.file, readonly=False)
        if solset in self.h5_out.getSolsetNames():
            solsetout = self.h5_out.getSolset(solset)
        else:
            solsetout = self.h5_out.makeSolset(solset)
        sources = np.array(self.directions.items(), dtype=[('name', 'S128'), ('dir', '<f4', (2,))])
        solsetout.obj.source.append(sources)

        axes_vals = {'dir': self.directions.keys(), 'ant': self.antennas, 'freq': self.ax_freq, 'time': self.ax_time}
        DPPP_axes = self.DPPP_style(soltab)

        if 'pol' in self.axes_new:
            if len(self.polarizations) > 0:
                axes_vals.update({'pol': self.polarizations})
            else:
                axes_vals.update({'pol': ['XX', 'YY']})
            self.axes_new = DPPP_axes
        elif len(self.axes_new) == 4 and len(DPPP_axes) > 0:
            self.axes_new = DPPP_axes

        # right order vals
        axes_vals = [v[1] for v in sorted(axes_vals.items(), key=lambda pair: self.axes_new.index(pair[0]))]

        if 'phase' in soltab:
            weights = np.ones(self.phases.shape)
            print('Value shape after processing {values}'.format(values=weights.shape))
            solsetout.makeSoltab('phase', axesNames=self.axes_new, axesVals=axes_vals, vals=self.phases,
                                 weights=weights)
        if 'amplitude' in soltab:
            weights = np.ones(self.gains.shape)
            print('Value shape after processing {values}'.format(values=weights.shape))
            solsetout.makeSoltab('amplitude', axesNames=self.axes_new, axesVals=axes_vals, vals=self.gains,
                                 weights=weights)
        if 'tec' in soltab:
            if self.axes_new.index('freq') == 1:
                self.phases = self.phases[:, 0, :, :]
            elif self.axes_new.index('freq') == 3:
                self.phases = self.phases[:, :, :, 0]
            else:
                self.phases = self.phases[:, :, 0, :]
            weights = np.ones(self.phases.shape)
            print('Value shape after processing {values}'.format(values=weights.shape))
            solsetout.makeSoltab('tec', axesNames=['dir', 'ant', 'time'],
                                 axesVals=[self.ax_time, self.antennas, self.directions.keys()],
                                 vals=self.phases, weights=weights)

        print('{solset}/{soltab} DONE'.format(solset=solset, soltab=soltab))
        self.h5_out.close()
        return self


def merge_h5(h5_out=None, h5_files=None, ms_files=None, convert_tec=True):
    """
    Main function that uses the class MergeH5 to merge h5 files.
    :param h5_out (string): h5 file name out
    :param h5_files (string or list): h5 files to merge
    :param ms_files: ms files to use, can be both list and string
    :param convert_tec (boolean): convert TEC to phase or not
    """
    if h5_out in glob(h5_out):
        os.system('rm {}'.format(h5_out))
    merge = MergeH5(h5_out=h5_out, h5_files=h5_files, ms_files=ms_files, convert_tec=convert_tec)
    merge.get_allkeys
    for ss in merge.all_solsets:
        for st_group in merge.all_soltabs:
            if len(st_group) > 0:
                for st in st_group:
                    merge.get_model_h5(ss, st)
                    merge.get_sol(ss, st)
                if merge.convert_tec and (('phase' in st_group[0]) or ('tec' in st_group[0])):
                    merge.create_new_dataset(ss, 'phase')
                else:
                    merge.create_new_dataset(ss, st)


if __name__ == '__main__':
    # merge_h5(h5_out='test.h5',
    # h5_files='/net/ouderijn/data2/rvweeren/test3C320_7/tecandphase_phmin1_skyselfcalcyle0_L807090_concat.4s16ch.P1.calcorrected.3C320phaseshift.ms.h5',
    # ms_files='/net/ouderijn/data2/rvweeren/test3C320_7/L807090_concat.4s16ch.P1.calcorrected.3C320phaseshift.ms',
    # convert_tec=True)

    merge_h5('merged_solutions.h5', h5_files='*ms.archive0*.avgsolsgrid*_9.h5', convert_tec=True,
             ms_files='*dysco.sub.shift.avg.weights.ms.archive0.goodtimes')

    # test=True
    # if test:
    # from runwscleanLBauto import applycal, makeimage
    # applycal(glob('box_3.dysco.sub.shift.avg.weights.ms.archive0.goodtimes.avg')[0], ['merged_solutions.h5'])