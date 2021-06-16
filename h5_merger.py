"""
LAST UPDATE: 14-6-2021

IMPORTANT: PYTHON 3 VERSION WORKS BEST, PYTHON 2 RAISES ISSUES IN DIRECTION ORDERS

After one has created solution files from self calling on the extracted boxes,
one can use this script to merge solution files
Do this by importing the function merge_h5 from this script:
-------------------------
EXAMPLE:

from h5_merger import merge_h5
merge_h5(h5_out='test.h5',
        h5_tables='*.h5',
        ms_files='*.ms',
        convert_tec=True)

TODO:
- bug fix python 2 with opening tables
"""

import os
from casacore import tables as ct
from glob import glob
from losoto.h5parm import h5parm
from losoto.lib_operations import reorderAxes
import numpy as np
from scipy.interpolate import interp1d
import sys

__author__ = "Jurjen de Jong (jurjendejong@strw.leidenuniv.nl)"
__all__ = ['merge_h5', 'str2bool']

# TODO: Weights keep on 1 -> future investigation
# TODO: Boxes -> Euclidean distance van de fluxes
# TODO: test rotation (fulljones)
# TODO: test convert_tec==False


class MergeH5:
    """Merge multiple h5 tables"""

    def __init__(self, h5_out, h5_tables=None, ms_files=None, convert_tec=True, make_new_direction=True):
        """
        :param h5_out: name of merged output h5 table
        :param files: h5 tables to merge, can be both list and string
        :param ms_files: ms files to use, can be both list and string
        :param convert_tec: convert TEC to phase or not
        """

        self.file = h5_out

        if type(ms_files) == list:
            ms = ms_files
            print("WARNING: MS list given, while only one MS used.")
        elif type(ms_files) == str:
            ms = glob(ms_files)
        else:
            ms = []

        if type(h5_tables) == list:
            self.h5_tables = h5_tables
        elif type(h5_tables) == str:
            self.h5_tables = glob(h5_tables)
        else:
            max_num_h5 = sorted(glob('*ms.archive0*.avgsolsgrid*.h5'))[-1].split('.')[-2].split('_')[-1]
            self.h5_tables = glob('*ms.archive0*.avgsolsgrid_%s.h5' % max_num_h5)
        if len(ms)>0:# check if there is a valid ms file
            t = ct.taql('SELECT CHAN_FREQ, CHAN_WIDTH FROM ' + ms[0] + '::SPECTRAL_WINDOW')
            self.ax_freq = t.getcol('CHAN_FREQ')[0]
            t.close()

            t = ct.table(ms[0])
            self.ax_time = sorted(np.unique(t.getcol('TIME')))
            t.close()
        else:# if we dont have ms files, we use the time and frequency axis of the longest h5 table
            print('No MS file given, will use h5 table for frequency and time axis')
            self.ax_time = []
            self.ax_freq = []
            for h5_name in self.h5_tables:
                h5 = h5parm(h5_name)
                for solset in h5.getSolsetNames():
                    ss = h5.getSolset(solset)
                    for soltab in ss.getSoltabNames():
                        st = ss.getSoltab(soltab)
                        try:
                            if len(st.getAxisValues('time'))>len(self.ax_time):
                                self.ax_time = st.getAxisValues('time')
                            if len(st.getAxisValues('freq'))>len(self.ax_freq):
                                self.ax_freq = st.getAxisValues('freq')
                        except:
                            pass
                h5.close()

        self.convert_tec = convert_tec  # convert tec or not
        self.make_new_direction = make_new_direction

        self.solaxnames = ['pol', 'dir', 'ant', 'freq', 'time'] # standard solax order to do our manipulations

    def get_values(self, ss, st, solset, soltab):
        """
        Get the values from the h5 table to work with.
        Also do some checks on the time and frequency axis.
        :param ss: solution set
        :param st: solution table
        :param solset: solset name
        :param soltab: soltab name
        """
        if 'pol' in st.getAxesNames():
            print("polarization is in {solset}/{soltab}".format(solset=solset, soltab=soltab))
        else:
            print("polarization is not in {solset}/{soltab}".format(solset=solset, soltab=soltab))

        time_axes = st.getAxisValues('time')
        freq_axes = st.getAxisValues('freq')

        print('Value shape before --> {values}'.format(values=st.getValues()[0].shape))

        if self.ax_time[0] > time_axes[-1] or time_axes[0] > self.ax_time[-1] :
            print("WARNING: Time axes of h5 and MS are not overlapping.")
        if self.ax_freq[0] > freq_axes[-1] or freq_axes[0] > self.ax_freq[-1]:
            print("WARNING: Frequency axes of h5 and MS are not overlapping.")
        if float(soltab[-3:]) > 0:
            print("WARNING: {soltab} does not end on 000".format(soltab=soltab))
        for av in self.axes_new:
            if av in st.getAxesNames() and st.getAxisLen(av) == 0:
                print("No {av} in {solset}/{soltab}".format(av=av, solset=solset, soltab=soltab))
        values = reorderAxes(st.getValues()[0], st.getAxesNames(), self.axes_current)

        return values, time_axes, freq_axes

    def sort_soltabs(self, soltabs):
        """
        Sort solution tables.
        This is import to run the steps and add directions according to our algorithm.
        Dont touch if you dont have to.
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
        Get all solution sets, solutions tables, and ax names in lists.
        """
        self.all_soltabs, self.all_solsets, self.all_axes, self.antennas = [], [], [], []

        for h5_name in self.h5_tables:
            h5 = h5parm(h5_name)
            for solset in h5.getSolsetNames():
                self.all_solsets += [solset]
                ss = h5.getSolset(solset)
                for n, soltab in enumerate(ss.getSoltabNames()):
                    self.all_soltabs += [soltab]
                    st = ss.getSoltab(soltab)
                    self.all_axes += ['/'.join([solset, soltab, an]) for an in st.getAxesNames()]
                    if n==0:
                        self.antennas = st.getAxisValues('ant') #check if same for all h5
                    elif list(self.antennas)!=list(st.getAxisValues('ant')):
                        print('ERROR: antennas not the same')
                        sys.exit()
            h5.close()
        self.all_soltabs = self.sort_soltabs(self.all_soltabs)
        self.all_solsets = set(self.all_solsets)
        self.all_axes = set(self.all_axes)
        return self

    def get_clean_values(self, soltab, st):
        """
        Get default values, based on model h5 table
        :param soltab: solution table name
        :param st: solution table itself
        """

        if 'pol' in st.getAxesNames():
            self.polarizations = st.getAxisValues('pol')
        else:
            self.polarizations = []  # most cases you want to have 5 dimensions but no polarization is still an option

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
            self.phases = np.zeros((2, 1, len(self.antennas), len(self.ax_freq), len(self.ax_time)))
        else:
            self.phases = np.zeros((1, len(self.antennas), len(self.ax_freq), len(self.ax_time)))

        self.directions = {}  # directions in a dictionary
        self.n = 0  # direction number reset

        return self

    @staticmethod
    def tecphase_conver(tec, freqs):
        """
        convert tec to phase
        :param tec: TEC
        :param freqs: frequencies
        :return tec phase converted values
        """
        return -8.4479745e9 * tec / freqs

    @staticmethod
    def interp_along_axis(x, interp_from, interp_to, axis):
        """
        Interpolate along axis
        :param x: frequency or time axis. Must be equal to the length of interp_from.
        :param interp_from: interpolate from this axis.
        :param interp_to: interpolate to this axis
        :param axis: interpolation axis
        :return return the interpolated result
        """
        interp_vals = interp1d(interp_from, x, axis=axis, kind='nearest', fill_value='extrapolate')
        new_vals = interp_vals(interp_to)
        return new_vals

    def get_model_h5(self, solset, soltab):
        """
        Get model (clean) h5 table
        :param solset: solution set name (sol000, sol001,..)
        :param soltab: solution table name
        """

        if '000' in soltab:

            for h5_name_to_merge in self.h5_tables:  # make template

                h5_to_merge = h5parm(h5_name_to_merge)

                if solset not in h5_to_merge.getSolsetNames():
                    h5_to_merge.close()
                    continue  # use other h5 table with solset
                else:
                    ss = h5_to_merge.getSolset(solset)
                if soltab not in ss.getSoltabNames():
                    h5_to_merge.close()
                    continue  # use other h5 table with soltab
                else:
                    st = ss.getSoltab(soltab)

                # if '/'.join([solset, soltab, 'pol']) in self.all_axes and 'pol' not in st.getAxesNames():
                # continue  # use h5 which has the polarizations as model file

                # make new clean values for new soltabs (ending on 000)
                if not self.convert_tec or (self.convert_tec and 'tec' not in soltab):
                    self.get_clean_values(soltab, st)
                    self.axes_new = [an for an in self.solaxnames if an in st.getAxesNames()]
                elif 'tec' in soltab and self.convert_tec:
                    for st_group in self.all_soltabs:
                        if soltab in st_group and ('phase000' not in st_group and 'phase{n}'.format(n=soltab[-3:])):
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
        for h5_name in self.h5_tables:

            h5 = h5parm(h5_name)
            if solset not in h5.getSolsetNames():
                h5.close()
                continue

            ss = h5.getSolset(solset)

            if soltab not in ss.getSoltabNames():
                h5.close()
                continue

            st = ss.getSoltab(soltab)

            #current axes for reordering of axes
            self.axes_current = [an for an in self.solaxnames if an in st.getAxesNames()]

            #get values, time, and freq axis
            values, time_axes, freq_axes = self.get_values(ss, st, solset, soltab)
            
            #update current and new axes if missing pol axes
            if len(self.axes_current)==4 and ((len(self.phases.shape)==5 
               and (st.getType() in ['phase', 'rotation'] or (st.getType()=='tec' and self.convert_tec)))
               or st.getType()=='amplitude' and len(self.gains.shape)==5):
               self.axes_current = ['pol'] + self.axes_current
               if len(self.axes_new)==4:
                  self.axes_new = ['pol'] + self.axes_new
                  
            #get source coordinates
            d = ss.getSou()
            source_coords = d[list(d.keys())[0]]
            d = 'Dir{:02d}'.format(self.n)

            print('Merge new h5 table in {direction}'.format(direction=d))

            if not self.make_new_direction and self.n==1:
                idx = 0
            elif any([np.array_equal(source_coords, list(sv)) for sv in self.directions.values()]):
                # Direction already exists, add to the existing solutions.
                idx = list([list(l) for l in self.directions.values()]).index(list(source_coords))
            else:# new direction
                print('Adding new direction {:f},{:f}'.format(*source_coords))
                idx = self.n
                self.directions.update({d: source_coords})
                if self.make_new_direction:
                    self.n += 1
                if self.n>1:# for self.n==1 we dont have to do anything
                    if st.getType() in ['tec','phase','rotation']:
                        shape = list(self.phases.shape)
                        dir_index = len(self.phases.shape)-4
                        if dir_index<0:
                            print('ERROR: Missing axes')
                            sys.exit()
                        if self.n>shape[dir_index]:
                            shape[dir_index]=1
                            self.phases = np.append(self.phases, np.zeros(shape),
                                                     axis=dir_index)#add clean phase to merge with
                    elif st.getType()=='amplitude':
                        shape = list(self.gains.shape)
                        dir_index = len(self.gains.shape)-4
                        if dir_index<0:
                            print('ERROR: Missing axes')
                            sys.exit()
                        if self.n>shape[dir_index]:
                            shape[dir_index]=1
                            self.gains = np.append(self.gains, np.ones(shape),
                                                    axis=dir_index)#add clean gain to merge with
            if st.getType() == 'tec':
                if self.convert_tec:# Convert tec to phase.
                    if len(self.polarizations) > 0 and len(self.phases.shape) == 5:
                        valtmp = np.ones((len(self.polarizations),) + values.shape)
                        valtmp[0, ...] = values
                        valtmp[-1, ...] = values
                        values = valtmp
                        # -1 assumes the expected shape along the frequency axis.
                        if self.axes_new[-2] != 'freq':
                            print('WARNING: Frequency axis is not on right position')
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
                        print('ERROR: Something went wrong with reshaping. Shouldnt end up here..')
                        sys.exit()

                    # Make tp shape same as phases
                    if len(self.phases.shape)==5 and tp.shape[0]==1:
                        phasetmp = np.zeros(self.phases.shape)
                        phasetmp[0,...] = tp[0, ...]
                        phasetmp[1,...] = tp[0, ...]
                        tp = phasetmp

                    # Add phases together
                    if len(tp.shape) - len(self.phases.shape) == 1:
                        self.phases[idx, ...] += tp[0,0, ...]
                        phasetmp = np.zeros((2,) + self.phases.shape[:])
                        phasetmp[0, ...] = self.phases
                        phasetmp[-1, ...] = self.phases
                        self.phases = phasetmp
                        if 'pol' not in self.axes_new:
                            self.axes_new = ['pol'] + self.axes_new

                    elif len(self.phases.shape) - len(tp.shape) == 1:#probably never reaches here
                        self.phases[0, idx, ...] += tp[0, ...]
                        self.phases[1, idx, ...] += tp[0, ...]

                    elif len(self.phases.shape) == len(tp.shape):
                        if len(self.phases.shape)==5:
                           self.phases[:, idx, ...] += tp[:, 0, ...]
                        elif len(self.phases.shape)==4:
                           self.phases[idx, ...] += tp[0, ...]

                    elif len(self.phases.shape)==5 and len(tp.shape)==5:
                        if self.phases.shape[0]==2 and tp.shape[0]==1:
                           self.phases[0, idx, ...] += tp[0, 0, ...]
                           self.phases[1, idx, ...] += tp[1, 0, ...]

                else:
                    if 'dir' in self.axes_current:#this line is trivial and could be removed
                        values = values[0, :, 0, :]

                    tp = self.interp_along_axis(values, time_axes, self.ax_time, -1)
                    tp = tp.reshape((1, tp.shape[0], 1, tp.shape[1]))
                    # Now add the tecs to the total phase correction for this direction.
                    if 'dir' in self.axes_current:#this line is trivial and could be removed
                        self.phases[idx, ...] += tp[0, ...]
                    else:
                        self.phases[idx, :, :] += tp

            elif st.getType() == 'phase' or st.getType() == 'rotation':
                if 'pol' in self.axes_current and 'pol' in st.getAxesNames():
                    if st.getAxisLen('pol') == 4:
                        print("Add fulljones type with 4 polarizations")
                        print("WARNING: this part hasn't been tested yet. Please check if this is correct")
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
                        print("WARNING: this part hasn't been tested yet. Please check if this is correct")
                        phasetmp = np.zeros((4,) + values.shape[1:])
                        phasetmp[0, ...] = values[0, ...]
                        phasetmp[-1, ...] = values[1, ...]
                        values = phasetmp
                elif 'pol' in self.axes_current and 'pol' not in st.getAxesNames() and len(self.phases.shape)==5:
                    phasetmp = np.zeros((self.phases.shape[0],) + values.shape)
                    phasetmp[0, ...] = values
                    phasetmp[-1, ...] = values
                    values = phasetmp

                idxnan = np.where((~np.isfinite(values)))
                values[idxnan] = 0.0

                tp = self.interp_along_axis(values, time_axes, self.ax_time,
                                            self.axes_current.index('time'))

                if tp.shape[-2] == 1:
                    tptmp = tp
                    for _ in self.ax_freq[:-1]:
                        tp = np.append(tp, tptmp, axis=-2)
                else:
                    tp = self.interp_along_axis(tp, freq_axes, self.ax_freq,
                                                self.axes_current.index('freq'))

                if len(self.phases.shape)==5 and self.phases.shape[0]==1:
                    phasetmp = np.zeros((2,) +  self.phases.shape[1:])
                    phasetmp[0,...] = self.phases[0, ...]
                    phasetmp[-1,...] = self.phases[0, ...]
                    self.phases = phasetmp

                if len(tp.shape) == len(self.phases.shape):
                    if len(self.phases.shape)==5:
                       self.phases[:, idx, ...] += tp[:, 0, ...]
                    elif len(self.phases.shape)==4:
                       self.phases[idx,...] += tp[0,...]
                       phasetmp = np.zeros((2,) + self.phases.shape[:])
                       phasetmp[0, ...] = self.phases
                       phasetmp[-1, ...] = self.phases
                       self.phases = phasetmp
                       if 'pol' not in self.axes_new:
                           self.axes_new = ['pol'] + self.axes_new
                elif len(tp.shape) - len(self.phases.shape) == 1:
                    self.phases[idx, ...] += tp[0,0, ...]
                    phasetmp = np.zeros((2,) + self.phases.shape[:])
                    phasetmp[0, ...] = self.phases
                    phasetmp[-1, ...] = self.phases
                    self.phases = phasetmp
                    if 'pol' not in self.axes_new:
                        self.axes_new = ['pol'] + self.axes_new
                elif len(self.phases.shape) - len(tp.shape) == 1:
                    self.phases[0, idx, ...] += tp
                    self.phases[-1, idx,...] += tp

            elif st.getType() == 'amplitude':
                if 'pol' in self.axes_current and 'pol' in st.getAxesNames():
                    if st.getAxisLen('pol') == 4:
                        print("Add fulljones type with 4 polarizations")
                        print("WARNING: this part hasn't been tested yet. Please check if this is correct")
                        if self.gains.shape[0] == 2:
                            gaintmp = np.zeros((4,) + self.gains.shape[1:])
                            gaintmp[0, ...] = self.gains[0, ...]
                            gaintmp[-1, ...] = self.gains[1, ...]
                            self.gains = gaintmp
                        elif len(self.gains.shape) < 5:
                            gaintmp = np.zeros((4,) + self.gains.shape)
                            gaintmp[0, ...] = self.gains
                            gaintmp[-1, ...] = self.gains
                            self.gains = gaintmp
                    elif st.getAxisLen('pol') == 2 and self.gains.shape[0] == 4:
                        print("Add to fulljones type with 4 polarizations")
                        print("WARNING: this part hasn't been tested yet. Please check if this is correct")
                        gaintmp = np.zeros((4,) + values.shape[1:])
                        gaintmp[0, ...] = values[0, ...]
                        gaintmp[-1, ...] = values[1, ...]
                        values = gaintmp
                elif 'pol' in self.axes_current and 'pol' not in st.getAxesNames() and len(self.gains.shape)==5:
                    phasetmp = np.zeros((self.gains.shape[0],) + values.shape)
                    phasetmp[0, ...] = values
                    phasetmp[-1, ...] = values
                    values = phasetmp

                idxnan = np.where((~np.isfinite(values)))
                values[idxnan] = 1.0
                tp = self.interp_along_axis(values, time_axes, self.ax_time,
                                            self.axes_current.index('time'))

                if tp.shape[-2] == 1:
                    tptmp = tp
                    for _ in self.ax_freq[:-1]:
                        tp = np.append(tp, tptmp, axis=-2)
                else:
                    tp = self.interp_along_axis(tp, freq_axes, self.ax_freq,
                                                self.axes_current.index('freq'))

                if len(self.gains.shape)==5 and self.gains.shape[0]==1:
                    gaintmp = np.zeros((2,) + self.gains.shape[1:])
                    gaintmp[0,...] = self.gains[0, ...]
                    gaintmp[-1,...] = self.gains[0, ...]
                    self.gains = gaintmp

                if len(self.gains.shape)==5 and len(tp.shape)==5:
                    self.gains[:, idx, ...] *= tp[:, 0, ...]
                elif len(self.gains.shape)==4 and len(tp.shape)==4:
                    self.gains[idx, ...] *= tp[0, ...]
                    gaintmp = np.zeros((2,) + self.gains.shape)
                    gaintmp[0,...]=self.gains
                    gaintmp[-1,...]=self.gains
                    self.gains = gaintmp
                    if 'pol' not in self.axes_new:
                        self.axes_new = ['pol'] + self.axes_new
                elif len(self.gains.shape)==5 and len(tp.shape)==4:
                    self.gains[0, idx, ...] *= tp[0,...]
                    self.gains[-1, idx, ...] *= tp[0,...]
                elif len(self.gains.shape)==4 and len(tp.shape)==5:
                    gaintmp = np.zeros((2,) + self.gains.shape)
                    gaintmp[0,...] = self.gains
                    gaintmp[-1,...] = self.gains
                    self.gains = gaintmp
                    self.gains[:, idx, ...] *= tp[:, 0, ...]
                    if 'pol' not in self.axes_new:
                        self.axes_new = ['pol'] + self.axes_new

            h5.close()

        return self

    def DPPP_style(self, soltab):
        """
        Reorder the axes in DPPP style because that is needed in other LOFAR pipeline (parts)
        :param soltab: solution table
        """
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
        Create a new dataset in the h5 table
        :param solset: solution set name
        :param soltab: solution table name
        """
        if len(self.directions.keys()) == 0:# return if no directions
            return self

        self.h5_out = h5parm(self.file, readonly=False)
        if solset in self.h5_out.getSolsetNames():
            solsetout = self.h5_out.getSolset(solset)
        else:
            solsetout = self.h5_out.makeSolset(solset)

        sources = list({i: (np.round(j[0],4), np.round(j[1],4)) for i,j in self.directions.items()}.items())
        #validate if new source directions are not already existing
        current_sources = [source[0].decode('UTF-8') for source in solsetout.obj.source[:]]
        new_sources = [source for source in sources if source[0] not in current_sources]
        if len(new_sources)>0:
            solsetout.obj.source.append(new_sources)

        axes_vals = {'dir': list(self.directions.keys()),
                     'ant': self.antennas,
                     'freq': self.ax_freq,
                     'time': self.ax_time}

        DPPP_axes = self.DPPP_style(soltab)#reorder the axis to DPPP style

        if 'pol' in self.axes_new:
            if len(self.polarizations) > 0:
                axes_vals.update({'pol': self.polarizations})
            else:
                axes_vals.update({'pol': ['XX', 'YY']})#need to be updated for rotation where len(pol)==4
            self.axes_new = DPPP_axes
        elif len(self.axes_new) == 4 and len(DPPP_axes) > 0:
            self.axes_new = DPPP_axes

        # right order vals
        axes_vals = [v[1] for v in sorted(axes_vals.items(), key=lambda pair: self.axes_new.index(pair[0]))]

        #make new solution table
        if 'phase' in soltab:
            weights = np.ones(self.phases.shape)
            print('Value shape after --> {values}'.format(values=weights.shape))
            solsetout.makeSoltab('phase', axesNames=self.axes_new, axesVals=axes_vals, vals=self.phases,
                                 weights=weights)
        if 'amplitude' in soltab:
            weights = np.ones(self.gains.shape)
            print('Value shape after --> {values}'.format(values=weights.shape))
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
            print('Value shape after --> {values}'.format(values=weights.shape))
            solsetout.makeSoltab('tec', axesNames=['dir', 'ant', 'time'],
                                 axesVals=[self.ax_time, self.antennas, list(self.directions.keys())],
                                 vals=self.phases, weights=weights)

        print('DONE: {solset}/{soltab}'.format(solset=solset, soltab=soltab))
        self.h5_out.close()
        return self


def merge_h5(h5_out=None, h5_tables=None, ms_files=None, convert_tec=True, make_new_direction=True):
    """
    Main function that uses the class MergeH5 to merge h5 tables.
    :param h5_out (string): h5 table name out
    :param h5_tables (string or list): h5 tables to merge
    :param ms_files (string or list): ms files to use, can be both list and string
    :param convert_tec (boolean): convert TEC to phase or not
    """
    if h5_out in glob(h5_out):
        os.system('rm {}'.format(h5_out))
    merge = MergeH5(h5_out=h5_out, h5_tables=h5_tables, ms_files=ms_files, convert_tec=convert_tec, make_new_direction=make_new_direction)
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
        #try:#add amplitude and phase if not available in h5 table
        if 'amplitude000' not in [item for sublist in merge.all_soltabs for item in sublist]:
            merge.gains = np.ones((2, len(merge.directions.keys()), len(merge.antennas), len(merge.ax_freq), len(merge.ax_time)))
            merge.axes_new = ['time', 'freq', 'ant', 'dir', 'pol']
            merge.polarizations = ['XX', 'YY']
            merge.gains = reorderAxes(merge.gains, merge.solaxnames, merge.axes_new)
            merge.create_new_dataset(ss, 'amplitude')
            #if 'phase000' not in [item for sublist in merge.all_soltabs for item in sublist] and \
               #'tec000' not in [item for sublist in merge.all_soltabs for item in sublist]:
                #merge.phases = np.zeros((2, len(merge.directions.keys()), len(merge.antennas), len(merge.ax_freq), len(merge.ax_time)))
                #merge.axes_new = ['time', 'freq', 'ant', 'dir', 'pol']
                #merge.polarizations = ['XX', 'YY']
                #merge.phases = reorderAxes(merge.phases, merge.solaxnames, merge.axes_new)
                #merge.create_new_dataset(ss, 'phase')
        #except:#add try to except to be sure that adding extra phase and amplitude is not going to break the code
            #pass
    print('END: h5 solution file(s) merged')

if __name__ == '__main__':
    from argparse import ArgumentParser, ArgumentTypeError

    def str2bool(v):
        v = str(v)
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise ArgumentTypeError('Boolean value expected.')

    parser = ArgumentParser()
    parser.add_argument('-out', '--h5_out', type=str, help='h5 table name for output')
    parser.add_argument('-in', '--h5_tables', action='append', type=str, help='h5 tables to merge')
    parser.add_argument('-ms', '--ms_files', type=str, help='ms files')
    parser.add_argument('-ct', '--convert_tec', type=bool, default=True, help='convert tec to phase')
    parser.add_argument('-nd', '--make_new_direction', type=str2bool, nargs='?', const=True, default=True, help='make new directions')

    args = parser.parse_args()

    # make sure h5 tables in right format
    if '[' in args.h5_tables:
        h5tables = args.h5_tables.replace('[','').replace(']','').replace(' ','').split(',')
    elif ' ' in args.h5_tables:
        h5tables = args.h5_tables.split()
    else:
        h5tables = args.h5_tables

    merge_h5(h5_out=args.h5_out, h5_tables=h5tables, ms_files=args.ms_files, convert_tec=args.convert_tec, make_new_direction=args.make_new_direction)