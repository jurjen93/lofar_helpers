"""
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

Following input parameters are possible:
h5_out ---> the output name of the h5 table
h5_tables ---> h5 tables that have to be merged
ms_files ---> read time and frequency from measurement set
h5_time_freq ---> read time and frequency from h5
convert_tec ---> convert tec to phase
merge_all_in_one ---> merge all in one direction (default is False), if True it adds everything in one direction
lin2circ ---> convert linear to circular polarization (default is False)
circ2lin ---> convert circular to linear polarization (default is False)
"""

# TODO: test rotation (fulljones)
# TODO: test convert_tec==False  ---> now they are deleted if convert_tec==false
# TODO: test circ2lin and vice versa

__author__ = "Jurjen de Jong (jurjendejong@strw.leidenuniv.nl)"

import os
from casacore import tables as ct
from glob import glob
from losoto.h5parm import h5parm
from losoto.lib_operations import reorderAxes
from scipy.interpolate import interp1d
import sys
import re
import tables
from numpy import zeros, ones, round, unique, array_equal, append, where, isfinite, expand_dims, pi, array, all, complex128, exp, angle

__all__ = ['merge_h5', 'str2bool']

def remove_numbers(inp):
    return "".join(re.findall("[a-zA-z]+", inp))


class MergeH5:
    """Merge multiple h5 tables"""

    def __init__(self, h5_out, h5_tables=None, ms_files=None, h5_time_freq=None, convert_tec=True, merge_all_in_one=False, solset='sol000'):
        """
        :param h5_out: name of merged output h5 table
        :param files: h5 tables to merge, can be both list and string
        :param ms_files: read time and frequency from measurement set
        :param h5_time_freq: read time and frequency from h5
        :param convert_tec: convert TEC to phase or not
        :param merge_all_in_one: merge all in one direction
        :param solset: solset number
        """

        self.file = h5_out
        self.solset = solset # for now this is standard sol0000

        if type(ms_files) == list:
            ms = ms_files
            print("WARNING: MS list given. Only first MS used.")
        elif type(ms_files) == str:
            ms = glob(ms_files)
        else:
            ms = []

        if type(h5_tables) == list:
            self.h5_tables = h5_tables
        elif type(h5_tables) == str:
            self.h5_tables = glob(h5_tables)
        else:
            print('No h5 table given. We will use all h5 tables in current folder.')
            self.h5_tables = glob('*.h5')
        if h5_time_freq:
            if len(ms)>0:
                print('--h5_time and/or --h5_freq are given, so measurement sets will not be used.')
            T = tables.open_file(h5_time_freq)
            self.ax_time = T.root.sol000.phase000.time[:]
            self.ax_freq = T.root.sol000.phase000.freq[:]
            T.close()

        elif len(ms)>0: # if there are multiple ms files
            print('Will take the time and freq from the following measurement sets:\n'+'\n'.join(ms))
            self.ax_time = array([])
            self.ax_freq = array([])
            for m in ms:
                t = ct.taql('SELECT CHAN_FREQ, CHAN_WIDTH FROM ' + m + '::SPECTRAL_WINDOW')
                self.ax_freq = append(self.ax_freq, t.getcol('CHAN_FREQ')[0])
                t.close()

                t = ct.table(m)
                self.ax_time = append(self.ax_time, t.getcol('TIME'))
                t.close()
            self.ax_time = array(sorted(unique(self.ax_time)))
            self.ax_freq = array(sorted(unique(self.ax_freq)))

        else:  # if we dont have ms files, we use the time and frequency axis of the longest h5 table
            print('No MS file given, will use h5 table for frequency and time axis')
            self.ax_time = array([])
            self.ax_freq = array([])
            for h5_name in self.h5_tables:
                h5 = h5parm(h5_name)
                ss = h5.getSolset(self.solset)
                for soltab in ss.getSoltabNames():
                    st = ss.getSoltab(soltab)
                    try:
                        if len(st.getAxisValues('time')) > len(self.ax_time):
                            self.ax_time = st.getAxisValues('time')
                    except:
                        print('No time axis in {solset}/{soltab}'.format(solset=solset, soltab=soltab))
                    try:
                        if len(st.getAxisValues('freq')) > len(self.ax_freq):
                            self.ax_freq = st.getAxisValues('freq')
                    except:
                        print('No freq axis in {solset}/{soltab}'.format(solset=solset, soltab=soltab))
                h5.close()

        if len(self.ax_freq) == 0:
            sys.exit('Cannot read frequency axis from input MS set or input H5.')
        if len(self.ax_time) == 0:
            sys.exit('Cannot read time axis from input MS or input H5.')
        if not self.same_antennas:
            sys.exit('Antenna tables are not the same')

        self.convert_tec = convert_tec  # convert tec or not
        self.merge_all_in_one = merge_all_in_one

        self.solaxnames = ['pol', 'dir', 'ant', 'freq', 'time']  # standard solax order to do our manipulations

        self.directions = {}  # directions in a dictionary

    @property
    def same_antennas(self):
        H_ref = tables.open_file(self.h5_tables[0])
        antennas_ref = H_ref.root.sol000.antenna[:]
        H_ref.close()
        for h5_name in self.h5_tables:
            H = tables.open_file(h5_name)
            antennas = H.root.sol000.antenna[:]
            H.close()
            if not all(antennas_ref == antennas):
                print('Antennas from '+self.h5_tables[0]+':')
                print(antennas_ref)
                print('Antennas from '+h5_name+':')
                print(antennas)
                return False
        return True

    def get_and_check_values(self, st, solset, soltab):
        """
        Get the values from the h5 table to merge.
        Also do some checks on the time and frequency axis.
        :param st: solution table
        :param solset: solset name
        :param soltab: soltab name
        """
        if 'pol' in st.getAxesNames():
            print("polarization is in {solset}/{soltab}".format(solset=solset, soltab=soltab))
        else:
            print("polarization is not in {solset}/{soltab}".format(solset=solset, soltab=soltab))

        time_axes = st.getAxisValues('time')

        if 'freq' in st.getAxesNames():
            freq_axes = st.getAxisValues('freq')
        else:
            freq_axes = self.ax_freq

        print('Value shape before --> {values}'.format(values=st.getValues()[0].shape))

        if self.ax_time[0] > time_axes[-1] or time_axes[0] > self.ax_time[-1]:
            print("WARNING: Time axes of h5 and MS are not overlapping.")

        if self.ax_freq[0] > freq_axes[-1] or freq_axes[0] > self.ax_freq[-1]:
            print("WARNING: Frequency axes of h5 and MS are not overlapping.")
        if float(soltab[-3:]) > 0:
            print("WARNING: {soltab} does not end on 000".format(soltab=soltab))
        for av in self.axes_new:
            if av in st.getAxesNames() and st.getAxisLen(av) == 0:
                print("No {av} in {solset}/{soltab}".format(av=av, solset=solset, soltab=soltab))

        if len(st.getAxesNames())!=len(st.getValues()[0].shape):
            sys.exit('Axes ({axlen}) and Value dimensions ({vallen}) are not equal'.format(axlen=len(st.getAxesNames()), vallen=len(st.getValues()[0].shape)))

        if 'dir' in st.getAxesNames():
            values = reorderAxes(st.getValues()[0], st.getAxesNames(), self.axes_current)
        else:
            print('No direction axis, we will add this.')
            origin_values = st.getValues()[0]
            values = reorderAxes(origin_values.reshape(origin_values.shape+(1,)), st.getAxesNames()+['dir'], self.axes_current)

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
                    # sorted(tp_tec, key=lambda x: float(x[-3:])),
                    sorted(tp_amplitude, key=lambda x: float(x[-3:])),
                    sorted(tp_rotation, key=lambda x: float(x[-3:]))]

    @staticmethod
    def has_integer(input):
        try:
            for s in str(input):
                if s.isdigit():
                    return True
            return False
        except:
            return False

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
                    if n == 0:
                        self.antennas = st.getAxisValues('ant')  # check if same for all h5
                    elif list(self.antennas) != list(st.getAxisValues('ant')):
                        sys.exit('ERROR: antennas not the same')
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

        num_dir = max(len(self.directions), 1)  # number of directions already existing

        if 'pol' in st.getAxesNames():
            self.polarizations = st.getAxisValues('pol')
        else:
            self.polarizations = []  # most cases you want to have 5 dimensions but no polarization is still an option

        if 'amplitude' in soltab and 'pol' in st.getAxesNames():
            self.gains = ones(
                (len(self.polarizations), num_dir, len(self.antennas), len(self.ax_freq), len(self.ax_time)))
        else:
            self.gains = ones((1, len(self.antennas), len(self.ax_freq), len(self.ax_time)))

        if 'phase' in soltab and 'pol' in st.getAxesNames():
            self.phases = zeros(
                (max(len(self.polarizations), 1), num_dir, len(self.antennas), len(self.ax_freq), len(self.ax_time)))

        elif 'rotation' in soltab:
            self.phases = zeros((num_dir, len(self.antennas), len(self.ax_freq), len(self.ax_time)))

        elif 'tec' in soltab:
            self.phases = zeros((2, num_dir, len(self.antennas), len(self.ax_freq), len(self.ax_time)))
        else:
            self.phases = zeros((num_dir, len(self.antennas), len(self.ax_freq), len(self.ax_time)))

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
        :param solset: solution set name (sol000)
        :param soltab: solution table name
        """

        if '000' in soltab:

            for h5_name_to_merge in self.h5_tables:  # make template

                h5_to_merge = h5parm(h5_name_to_merge)

                if solset not in h5_to_merge.getSolsetNames():
                    h5_to_merge.close()
                    sys.exit(solset +' does not exist in '+h5_name_to_merge)
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

    @staticmethod
    def get_number_of_directions(st):
        """
        Get number of directions in solution table
        :param st: solution table
        """
        if 'dir' in st.getAxesNames():
            dir_index = st.getAxesNames().index('dir')
            return st.getValues()[0].shape[dir_index]
        else:
            return 1

    def merge_files(self, solset, soltab):
        """
        Merge the h5 files
        :param solset: solution set name
        :param soltab: solution table name
        """
        for h5_name in self.h5_tables:

            print(h5_name, len(self.h5_tables))

            h5 = h5parm(h5_name)
            if solset not in h5.getSolsetNames():
                h5.close()
                continue

            ss = h5.getSolset(solset)

            if soltab not in ss.getSoltabNames():
                h5.close()
                continue

            st = ss.getSoltab(soltab)

            # current axes for reordering of axes
            self.axes_current = [an for an in self.solaxnames if an in st.getAxesNames()]

            # add dir if missing
            if 'dir' not in self.axes_new:
                self.axes_new.insert(1, 'dir')
            if 'dir' not in self.axes_current:
                self.axes_current.insert(1, 'dir')

            init_dir_index = self.axes_current.index('dir')  # index of direction

            print('Solution table from {table}'.format(table=h5_name.split('/')[-1]))
            num_dirs = self.get_number_of_directions(st)  # number of directions
            print('This table has {numdirection} direction(s)'.format(numdirection=num_dirs))

            # get values, time, and freq axis
            table_values, time_axes, freq_axes = self.get_and_check_values(st, solset, soltab)

            for dir_idx in range(num_dirs):#loop over all directions

                shape = list(table_values.shape)
                shape[init_dir_index] = 1
                values = zeros(shape)

                if init_dir_index == 0:
                    values[0, ...] += table_values[dir_idx, ...]
                elif init_dir_index == 1:
                    values[:, 0, ...] += table_values[:, dir_idx, ...]
                elif init_dir_index == 2:
                    values[:, :, 0, ...] += table_values[:, :, dir_idx, ...]
                elif init_dir_index == 3:
                    values[:, :, :, 0, ...] += table_values[:, :, :, dir_idx, ...]
                elif init_dir_index == 4:
                    values[:, :, :, :, 0, ...] += table_values[:, :, :, :, dir_idx, ...]

                # update current and new axes if missing pol axes
                if len(self.axes_current) == 4 and ((len(self.phases.shape) == 5
                                                     and (st.getType() in ['phase', 'rotation'] or (
                                st.getType() == 'tec' and self.convert_tec)))
                                                    or st.getType() == 'amplitude' and len(self.gains.shape) == 5):
                    self.axes_current = ['pol'] + self.axes_current
                    if len(self.axes_new) == 4:
                        self.axes_new = ['pol'] + self.axes_new

                # get source coordinates
                d = ss.getSou()
                source_coords = d[list(d.keys())[dir_idx]]

                d = 'Dir{:02d}'.format(self.n)

                if self.merge_all_in_one and self.n == 1:
                    idx = 0
                    print('Merging direction {:f},{:f} with previous direction'.format(*source_coords))
                    if self.directions['Dir00'][0]>0 and self.directions['Dir00'][1]>0:
                        self.directions.update({'Dir00': source_coords}) # 0.0 coordinate bug
                        print('Adding new direction {:f},{:f}'.format(*source_coords))
                elif any([array_equal(source_coords, list(sv)) for sv in self.directions.values()]):
                    # Direction already exists, add to the existing solutions.
                    print('Direction {:f},{:f} already exists. Adding to this direction.'.format(*source_coords))
                    idx = list([list(l) for l in self.directions.values()]).index(list(source_coords))
                else:  # new direction
                    if source_coords[0] > 0 and source_coords[1] > 0:
                        print('Adding new direction {:f},{:f}'.format(*source_coords))
                    idx = self.n
                    self.directions.update({d: source_coords})
                    if not self.merge_all_in_one:
                        self.n += 1
                    if self.n > 1:  # for self.n==1 we dont have to do anything
                        if st.getType() in ['tec', 'phase', 'rotation']:
                            shape = list(self.phases.shape)
                            dir_index = len(self.phases.shape) - 4
                            if dir_index < 0:
                                sys.exit('ERROR: Missing axes')
                            if self.n > shape[dir_index]:
                                shape[dir_index] = 1
                                self.phases = append(self.phases, zeros(shape),
                                                        axis=dir_index)  # add clean phase to merge with
                        elif st.getType() == 'amplitude':
                            shape = list(self.gains.shape)
                            dir_index = len(self.gains.shape) - 4
                            if dir_index < 0:
                                sys.exit('ERROR: Missing axes')
                            if self.n > shape[dir_index]:
                                shape[dir_index] = 1
                                self.gains = append(self.gains, ones(shape),
                                                       axis=dir_index)  # add clean gain to merge with
                if st.getType() == 'tec':

                    # add frequencies
                    if 'freq' not in st.getAxesNames() and len(st.getAxesNames())==3:
                        ax = self.axes_new.index('freq') - len(self.axes_new)
                        values = expand_dims(values, axis=ax)
                        if 'pol' not in st.getAxesNames():
                            self.axes_current = ['dir', 'ant', 'freq', 'time']
                        else:
                            self.axes_current = ['pol', 'dir', 'ant', 'freq', 'time']
                        valuestmp = values
                        for _ in range(len(self.ax_freq)-1):
                            values = append(values, valuestmp, axis=-2)

                    if self.convert_tec:  # Convert tec to phase.
                        if len(self.polarizations) > 0 and len(self.phases.shape) == 5 and 'pol' not in self.axes_current:
                            valtmp = ones((len(self.polarizations),) + values.shape)
                            valtmp[0, ...] = values
                            valtmp[-1, ...] = values
                            values = valtmp

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
                            sys.exit('ERROR: Something went wrong with reshaping. Shouldnt end up here..')
                        # Make tp shape same as phases

                        if len(self.phases.shape) == 5 and tp.shape[0] == 1:
                            phasetmp = zeros(self.phases.shape)
                            phasetmp[0, ...] = tp[0, ...]
                            phasetmp[1, ...] = tp[0, ...]
                            tp = phasetmp

                        # Add phases together
                        if len(tp.shape) - len(self.phases.shape) == 1:
                            self.phases[idx, ...] += tp[0, 0, ...]
                            phasetmp = zeros((2,) + self.phases.shape[:])
                            phasetmp[0, ...] = self.phases
                            phasetmp[-1, ...] = self.phases
                            self.phases = phasetmp
                            if 'pol' not in self.axes_new:
                                self.axes_new = ['pol'] + self.axes_new

                        elif len(self.phases.shape) - len(tp.shape) == 1:  # probably never reaches here
                            self.phases[0, idx, ...] += tp[0, ...]
                            self.phases[1, idx, ...] += tp[0, ...]

                        elif len(self.phases.shape) == len(tp.shape):
                            if len(self.phases.shape) == 5:
                                self.phases[:, idx, ...] += tp[:, 0, ...]
                            elif len(self.phases.shape) == 4:
                                self.phases[idx, ...] += tp[0, ...]

                        elif len(self.phases.shape) == 5 and len(tp.shape) == 5:
                            if self.phases.shape[0] == 2 and tp.shape[0] == 1:
                                self.phases[0, idx, ...] += tp[0, 0, ...]
                                self.phases[1, idx, ...] += tp[1, 0, ...]

                    else:
                        if 'dir' in self.axes_current:  # this line is trivial and could maybe be removed
                            values = values[0, :, 0, :]

                        tp = self.interp_along_axis(values, time_axes, self.ax_time, -1)
                        tp = tp.reshape((1, tp.shape[0], 1, tp.shape[1]))
                        # Now add the tecs to the total phase correction for this direction.
                        if 'dir' in self.axes_current:  # this line is trivial and could be removed
                            self.phases[idx, ...] += tp[0, ...]
                        else:
                            self.phases[idx, :, :] += tp

                elif st.getType() == 'phase' or st.getType() == 'rotation':
                    if 'pol' in self.axes_current and 'pol' in st.getAxesNames():
                        if st.getAxisLen('pol') == 4:
                            print("Add fulljones type with 4 polarizations")
                            print("WARNING: this part hasn't been properly tested yet. Please check if output is correct.")
                            if self.phases.shape[0] == 2:
                                phasetmp = zeros((4,) + self.phases.shape[1:])
                                phasetmp[0, ...] = self.phases[0, ...]
                                phasetmp[-1, ...] = self.phases[1, ...]
                                self.phases = phasetmp
                            elif len(self.phases.shape) < 5:
                                phasetmp = zeros((4,) + self.phases.shape)
                                phasetmp[0, ...] = self.phases
                                phasetmp[-1, ...] = self.phases
                                self.phases = phasetmp
                        elif st.getAxisLen('pol') == 2 and self.phases.shape[0] == 4:
                            print("Add to fulljones type with 4 polarizations")
                            print("WARNING: this part hasn't been properly tested yet. Please check if output is correct.")
                            phasetmp = zeros((4,) + values.shape[1:])
                            phasetmp[0, ...] = values[0, ...]
                            phasetmp[-1, ...] = values[1, ...]
                            values = phasetmp
                    elif 'pol' in self.axes_current and 'pol' not in st.getAxesNames() and len(self.phases.shape) == 5:
                        phasetmp = zeros((self.phases.shape[0],) + values.shape)
                        phasetmp[0, ...] = values
                        phasetmp[-1, ...] = values
                        values = phasetmp

                    idxnan = where((~isfinite(values)))
                    values[idxnan] = 0.0

                    tp = self.interp_along_axis(values, time_axes, self.ax_time,
                                                self.axes_current.index('time'))

                    if tp.shape[-2] == 1:
                        tptmp = tp
                        for _ in self.ax_freq[:-1]:
                            tp = append(tp, tptmp, axis=-2)
                    else:
                        tp = self.interp_along_axis(tp, freq_axes, self.ax_freq,
                                                    self.axes_current.index('freq'))

                    if len(self.phases.shape) == 5 and self.phases.shape[0] == 1:
                        phasetmp = zeros((2,) + self.phases.shape[1:])
                        phasetmp[0, ...] = self.phases[0, ...]
                        phasetmp[-1, ...] = self.phases[0, ...]
                        self.phases = phasetmp

                    if len(tp.shape) == len(self.phases.shape):
                        if len(self.phases.shape) == 5:
                            self.phases[:, idx, ...] += tp[:, 0, ...]
                        elif len(self.phases.shape) == 4:
                            self.phases[idx, ...] += tp[0, ...]
                            phasetmp = zeros((2,) + self.phases.shape[:])
                            phasetmp[0, ...] = self.phases
                            phasetmp[-1, ...] = self.phases
                            self.phases = phasetmp
                            if 'pol' not in self.axes_new:
                                self.axes_new = ['pol'] + self.axes_new
                    elif len(tp.shape) - len(self.phases.shape) == 1:
                        self.phases[idx, ...] += tp[0, 0, ...]
                        phasetmp = zeros((2,) + self.phases.shape[:])
                        phasetmp[0, ...] = self.phases
                        phasetmp[-1, ...] = self.phases
                        self.phases = phasetmp
                        if 'pol' not in self.axes_new:
                            self.axes_new = ['pol'] + self.axes_new
                    elif len(self.phases.shape) - len(tp.shape) == 1:
                        self.phases[0, idx, ...] += tp
                        self.phases[-1, idx, ...] += tp

                elif st.getType() == 'amplitude':
                    if 'pol' in self.axes_current and 'pol' in st.getAxesNames():
                        if st.getAxisLen('pol') == 4:
                            print("Add fulljones type with 4 polarizations")
                            print("WARNING: this part hasn't been properly tested yet. Please check if output is correct.")
                            if self.gains.shape[0] == 2:
                                gaintmp = zeros((4,) + self.gains.shape[1:])
                                gaintmp[0, ...] = self.gains[0, ...]
                                gaintmp[-1, ...] = self.gains[1, ...]
                                self.gains = gaintmp
                            elif len(self.gains.shape) < 5:
                                gaintmp = zeros((4,) + self.gains.shape)
                                gaintmp[0, ...] = self.gains
                                gaintmp[-1, ...] = self.gains
                                self.gains = gaintmp
                        elif st.getAxisLen('pol') == 2 and self.gains.shape[0] == 4:
                            print("Add to fulljones type with 4 polarizations")
                            print("WARNING: this part hasn't been properly tested yet. Please check if output is correct.")
                            gaintmp = zeros((4,) + values.shape[1:])
                            gaintmp[0, ...] = values[0, ...]
                            gaintmp[-1, ...] = values[1, ...]
                            values = gaintmp
                    elif 'pol' in self.axes_current and 'pol' not in st.getAxesNames() and len(self.gains.shape) == 5:
                        phasetmp = zeros((self.gains.shape[0],) + values.shape)
                        phasetmp[0, ...] = values
                        phasetmp[-1, ...] = values
                        values = phasetmp

                    idxnan = where((~isfinite(values)))
                    values[idxnan] = 1.0
                    tp = self.interp_along_axis(values, time_axes, self.ax_time,
                                                self.axes_current.index('time'))

                    if tp.shape[-2] == 1:
                        tptmp = tp
                        for _ in self.ax_freq[:-1]:
                            tp = append(tp, tptmp, axis=-2)
                    else:
                        tp = self.interp_along_axis(tp, freq_axes, self.ax_freq,
                                                    self.axes_current.index('freq'))

                    if len(self.gains.shape) == 5 and self.gains.shape[0] == 1:
                        gaintmp = zeros((2,) + self.gains.shape[1:])
                        gaintmp[0, ...] = self.gains[0, ...]
                        gaintmp[-1, ...] = self.gains[0, ...]
                        self.gains = gaintmp

                    if len(self.gains.shape) == 5 and len(tp.shape) == 5:
                        self.gains[:, idx, ...] *= tp[:, 0, ...]
                    elif len(self.gains.shape) == 4 and len(tp.shape) == 4:
                        self.gains[idx, ...] *= tp[0, ...]
                        gaintmp = zeros((2,) + self.gains.shape)
                        gaintmp[0, ...] = self.gains
                        gaintmp[-1, ...] = self.gains
                        self.gains = gaintmp
                        if 'pol' not in self.axes_new:
                            self.axes_new = ['pol'] + self.axes_new
                    elif len(self.gains.shape) == 5 and len(tp.shape) == 4:
                        self.gains[0, idx, ...] *= tp[0, ...]
                        self.gains[-1, idx, ...] *= tp[0, ...]
                    elif len(self.gains.shape) == 4 and len(tp.shape) == 5:
                        gaintmp = zeros((2,) + self.gains.shape)
                        gaintmp[0, ...] = self.gains
                        gaintmp[-1, ...] = self.gains
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

    def order_directions(self):
        """
        This method will be called when the user is using python 2, as there is a bug in the direction that we can resolve
        with this extra step.
        """
        import h5py
        from numpy import sort

        T = h5py.File(self.file, 'r+')
        for ss in T.keys():
            T[ss+'/source'][:] = sort(T[ss+'/source'][:])
            if '00' in ss: # for some reason it sorts randomly ascending or descending, this extra step is a fix for that
                T[ss + '/source'][:] = sort(T[ss + '/source'][:])
            for st in T[ss].keys():
                if self.has_integer(st):
                    for ax in T['/'.join([ss, st])]:
                        if 'dir' == ax:
                            H = tables.open_file(self.file, 'r+')
                            if ss=='sol000':
                                if st=='phase000': # the following if-elif-else statements are ugly but there is a bug in h5py which isn't resolved yet.
                                    H.root.sol000.phase000.dir[:] = [c[0] for c in H.root.sol000.source[:]]
                                elif st=='tec000':
                                    H.root.sol000.tec000.dir[:] = [c[0] for c in H.root.sol000.source[:]]
                                elif st=='amplitude000':
                                    H.root.sol000.amplitude000.dir[:] = [c[0] for c in H.root.sol000.source[:]]
                            else:
                                print('ERROR: {ss}!=sol000 and you are using Python 2.\n'
                                      'This means your direction table could be wrongly ordered.\n'
                                      'Best solution is to use Python 3.'.format(ss=ss))
                            H.close()
        T.close()
        return self

    def reduce_memory_source(self):
        """
        We need to store the data in 136 bytes per directions.
        Python 3 saves it automatically in more than that number.
        """
        # tables.file._open_files.close_all()
        T = tables.open_file(self.file, 'r+')
        new_source = array(T.root.sol000.source[:], dtype=[('name', 'S128'), ('dir', '<f4', (2,))])
        T.root.sol000.source._f_remove()
        T.create_table(T.root.sol000, 'source', new_source, "Source names and directions")
        T.close()
        return self

    @staticmethod
    def keep_new_sources(current_sources, new_sources):
        """
        Remove sources from new_sources that are already in current_sources
        :param current_sources: current sources that we need to compare with new_sources
        :param new_sources: new sources to be add
        :return ---> New unique sources
        """
        current_sources_dir = [source[0].decode('UTF-8') for source in current_sources]
        current_sources_coor = [source[1] for source in current_sources]
        new_sources = [source for source in new_sources if source[0] not in current_sources_dir]
        del_index = []
        for coor in current_sources_coor:
            for n, source in enumerate(new_sources):
                if round(coor[0], 4) == round(new_sources[1][0], 4) and round(coor[1], 4) == round(new_sources[1][1], 4):
                    del_index.append(n)
        return [source for i, source in enumerate(new_sources) if i not in del_index]

    def create_new_dataset(self, solset, soltab):
        """
        Create a new dataset in the h5 table
        :param solset: solution set name
        :param soltab: solution table name
        """
        if len(self.directions.keys()) == 0:  # return if no directions
            return self

        self.h5_out = h5parm(self.file, readonly=False)
        if solset in self.h5_out.getSolsetNames():
            solsetout = self.h5_out.getSolset(solset)
        else:
            solsetout = self.h5_out.makeSolset(solset)

        # validate if new source directions are not already existing
        # current_sources_dir = [source[0].decode('UTF-8') for source in solsetout.obj.source[:]]
        # current_sources_coor = [source[1] for source in solsetout.obj.source[:]]
        # new_sources = [source for source in sources if
        #                (source[0] not in current_sources) and (float(source[1][0])>0 and float(source[1][1])>0)]
        # new_sources = [source for source in sources if (source[0] not in current_sources_dir) and (source[1] not in current_sources_coor)]

        new_sources = self.keep_new_sources(solsetout.obj.source[:], list(self.directions.items()))

        # antennas = st.getAxisvalues('ant')
        # ss_antennas = ss.obj.antennas.read()
        # antennasout = solsetout.getAnt()
        # antennatable=solsetout.obj._f_get_child('antenna')
        # antennatable.append(ss_antennas)

        if len(new_sources) > 0:
            solsetout.obj.source.append(new_sources)

        axes_vals = {'dir': list(self.directions.keys()),
                     'ant': self.antennas,
                     'freq': self.ax_freq,
                     'time': self.ax_time}

        DPPP_axes = self.DPPP_style(soltab)  # reorder the axis to DPPP style

        if 'pol' in self.axes_new:
            if len(self.polarizations) > 1:
                axes_vals.update({'pol': self.polarizations})
            else:
                axes_vals.update({'pol': ['XX', 'YY']})  # need to be updated for rotation where len(pol)==4
            self.axes_new = DPPP_axes
        elif len(self.axes_new) == 4 and len(DPPP_axes) > 0:
            self.axes_new = DPPP_axes

        # right order vals
        axes_vals = [v[1] for v in sorted(axes_vals.items(), key=lambda pair: self.axes_new.index(pair[0]))]

        # make new solution table
        if 'phase' in soltab:
            weights = ones(self.phases.shape)
            print('Value shape after --> {values}'.format(values=weights.shape))
            solsetout.makeSoltab('phase', axesNames=self.axes_new, axesVals=axes_vals, vals=self.phases,
                                 weights=weights)
        if 'amplitude' in soltab:
            weights = ones(self.gains.shape)
            print('Value shape after --> {values}'.format(values=weights.shape))
            solsetout.makeSoltab('amplitude', axesNames=self.axes_new, axesVals=axes_vals, vals=self.gains,
                                 weights=weights)
        if 'tec' in soltab:
            print('ADD TEC')
            if self.axes_new.index('freq') == 1:
                self.phases = self.phases[:, 0, :, :]
            elif self.axes_new.index('freq') == 3:
                self.phases = self.phases[:, :, :, 0]
            else:
                self.phases = self.phases[:, :, 0, :]
            weights = ones(self.phases.shape)
            print('Value shape after --> {values}'.format(values=weights.shape))
            solsetout.makeSoltab('tec', axesNames=['dir', 'ant', 'time'],
                                 axesVals=[self.ax_time, self.antennas, list(self.directions.keys())],
                                 vals=self.phases, weights=weights)

        print('DONE: {solset}/{soltab}'.format(solset=solset, soltab=soltab))
        self.h5_out.close()
        return self

    def add_directions(self, add_directions=None):
        """
        Add default directions (phase all zeros, amplitude all ones)
        :param add_directions: list with directions
        """
        if not add_directions:
            return self

        h5 = h5parm(self.file, readonly=True)
        filetemp = self.file.replace('.h5','_temp.h5')
        h5_temp = h5parm(filetemp, readonly=False)
        solset = h5.getSolset(self.solset)
        solsettemp = h5_temp.makeSolset(self.solset)
        if type(add_directions[0])==list:
            sources = list([source[1] for source in solset.obj.source[:]]) + add_directions
        else:
            sources = list([source[1] for source in solset.obj.source[:]]) + [add_directions]
        sources = [(bytes('Dir' + str(n).zfill(2), 'utf-8'), list(ns)) for n, ns in enumerate(sources)]
        if len(sources) > 0:
            solsettemp.obj.source.append(sources)

        for st in h5.getSolset(self.solset).getSoltabNames():
            solutiontable = h5.getSolset(self.solset).getSoltab(st)
            axes = solutiontable.getValues()[1]
            values = solutiontable.getValues()[0]
            axes['dir'] = [ns[0] for ns in sources]
            dir_index = solutiontable.getAxesNames().index('dir')
            new_shape = list(values.shape)
            last_idx = new_shape[dir_index]
            new_idx = last_idx+len(add_directions)-1
            new_shape[dir_index] = new_idx

            if 'phase' in st:
                values_new = zeros(tuple(new_shape))
            elif 'amplitude' in st:
                values_new = ones(tuple(new_shape))
            else:
                values_new = zeros(tuple(new_shape))

            if dir_index == 0:
                values_new[0:last_idx, ...] = values
            elif dir_index == 1:
                values_new[:, 0:last_idx, ...] = values
            elif dir_index == 2:
                values_new[:, :, 0:last_idx, ...] = values
            elif dir_index == 3:
                values_new[:, :, :, 0:last_idx, ...] = values
            elif dir_index == 4:
                values_new[:, :, :, :, 0:last_idx, ...] = values

            weights = ones(values_new.shape)
            solsettemp.makeSoltab(remove_numbers(st), axesNames=list(axes.keys()), axesVals=list(axes.values()), vals=values_new,
                             weights=weights)

            print('Default directions added for '+self.solset+'/'+st)
            print('Shape change: '+str(values.shape)+' ---> '+str(values_new.shape))

        h5.close()
        h5_temp.close()

        os.system('rm '+self.file +' && mv '+filetemp+' '+self.file)

        return self

    def remove_pol(self, single=False):
        """
        Reduce table to one single polarization
        :param single: if single==True we leave a single pole such that values have shape=(..., 1), if False we remove pol-axis entirely
        """
        T = tables.open_file(self.file, 'r+')
        if single:
            newpol = array([b'I'], dtype='|S2')
        try:
            if T.root.sol000.phase000.val[0,0,0,0,0] == T.root.sol000.phase000.val[0,0,0,0,-1] and \
                T.root.sol000.phase000.val[-1, 0, 0, 0, 0] == T.root.sol000.phase000.val[-1, 0, 0, 0, -1]:
                if single:
                    print('Phase has same values for XX and YY polarization.\nReducing into one Polarization I.')
                else:
                    print('Phase has same values for XX and YY polarization.\nRemoving Polarization.')
                if single:
                    newphase = T.root.sol000.phase000.val[:, :, :, :, 0:1]
                else:
                    newphase = T.root.sol000.phase000.val[:, :, :, :, 0]
                T.root.sol000.phase000.val._f_remove()
                T.root.sol000.phase000.pol._f_remove()
                T.create_array(T.root.sol000.phase000, 'val', newphase)
                if single:
                    T.root.sol000.phase000.val.attrs['AXES'] = b'time,freq,ant,dir,pol'
                    T.create_array(T.root.sol000.phase000, 'pol', newpol)
                else:
                    T.root.sol000.phase000.val.attrs['AXES'] = b'time,freq,ant,dir'
            else:
                sys.exit('ERROR: Phase has not the same values for XX and YY polarization.\nERROR: No polarization reduction will be done.')
        except:
            pass

        try:
            if T.root.sol000.amplitude000.val[0, 0, 0, 0, 0] == T.root.sol000.amplitude000.val[0, 0, 0, 0, -1] and \
                T.root.sol000.amplitude000.val[-1, 0, 0, 0, 0] == T.root.sol000.amplitude000.val[-1, 0, 0, 0, -1]:
                if single:
                    print('Amplitude has same values for XX and YY polarization.\nReducing into one Polarization I.')
                else:
                    print('Amplitude has same values for XX and YY polarization.\nRemoving Polarization.')
                if single:
                    newampl = T.root.sol000.amplitude000.val[:, :, :, :, 0:1]
                else:
                    newampl = T.root.sol000.amplitude000.val[:, :, :, :, 0]
                T.root.sol000.amplitude000.val._f_remove()
                T.root.sol000.amplitude000.pol._f_remove()
                T.create_array(T.root.sol000.amplitude000, 'val', newampl)
                if single:
                    T.root.sol000.amplitude000.val.attrs['AXES'] = b'time,freq,ant,dir,pol'
                    T.create_array(T.root.sol000.amplitude000, 'pol', newpol)
                else:
                    T.root.sol000.amplitude000.val.attrs['AXES'] = b'time,freq,ant,dir'
            else:
                sys.exit('ERROR: Amplitude has not the same values for XX and YY polarization.\nERROR: No polarization reduction will be done.')
        except:
            pass

        T.close()
        return self

    def add_antennas(self):
        "Add antennas to output table"
        print('Add antenna table.')
        T = tables.open_file(self.h5_tables[0])
        antennas = T.root.sol000.antenna[:]
        T.close()
        H = tables.open_file(self.file, 'r+')
        H.root.sol000.antenna._f_remove()
        H.create_table(H.root.sol000, 'antenna', antennas, "Antenna names and positions")
        H.close()
        return self


def make_h5_name(h5_name):
    if '.h5' != h5_name[-3:]:
        h5_name += '.h5'
    return h5_name

def change_solset(h5, solset_in, solset_out, delete=True, overwrite=True):
    """
    This function is to change the solset numbers.

    1) Copy solset_in to solset_out (overwriting if overwrite==True)
    2) Delete solset_in if delete==True
    """
    H = tables.open_file(h5, 'r+')
    if solset_in=='sol002':
        H.root.sol002._f_copy(H.root, newname=solset_out, overwrite=overwrite, recursive=True)
        print('Succesfully copied ' + solset_in + ' to ' + solset_out)
        if delete:
            H.root.sol002._f_remove(recursive=True)
            print('Removed ' + solset_in)
        H.close()
    elif solset_in=='sol001':
        H.root.sol001._f_copy(H.root, newname=solset_out, overwrite=overwrite, recursive=True)
        print('Succesfully copied ' + solset_in + ' to ' + solset_out)
        if delete:
            H.root.sol001._f_remove(recursive=True)
            print('Removed ' + solset_in)
        H.close()
    elif solset_in=='sol000':
        H.root.sol000._f_copy(H.root, newname=solset_out, overwrite=overwrite, recursive=True)
        print('Succesfully copied ' + solset_in + ' to ' + solset_out)
        if delete:
            H.root.sol000._f_remove(recursive=True)
            print('Removed ' + solset_in)
        H.close()

    else:
        H.close()
        sys.exit('Code is not advanced enough to change the solset name.\n '
                 'Please contact jurjendejong@strw.leidenuniv.nl.')

class PolChange:
    """
    This Python class helps to convert polarization from linear to circular or vice versa.
    """
    def __init__(self, h5_in, h5_out):
        """
        :param h5_in: h5 input name
        :param h5_out: h5 output name
        """
        self.h5_in = h5parm(h5_in, readonly=True)
        self.h5_out = h5parm(h5_out, readonly=False)
        self.axes_names = ['time', 'freq', 'ant', 'dir', 'pol']

    @staticmethod
    def lin2circ(G):
        """Convert linear polarization to circular polarization"""
        RR = (G[..., 0] + G[..., -1]).astype(complex128)
        LL = (G[..., 0] + G[..., -1]).astype(complex128)
        RL = (G[..., 0] - G[..., -1]).astype(complex128)
        LR = (G[..., 0] - G[..., -1]).astype(complex128)

        if G.shape[-1] == 4:
            RR += 1j * (G[..., 2] - G[..., 1])
            LL += 1j * (G[..., 1] - G[..., 2])
            RL += 1j * (G[..., 2] + G[..., 1])
            LR -= 1j * (G[..., 2] + G[..., 1])

        RR /= 2
        LL /= 2
        RL /= 2
        LR /= 2

        G_new = zeros(G.shape[0:-1] + (4,)).astype(complex128)
        G_new[..., 0] += RR
        G_new[..., 1] += RL
        G_new[..., 2] += LR
        G_new[..., 3] += LL
        return G_new

    @staticmethod
    def circ2lin(G):
        """Convert circular polarization to linear polarization"""
        XX = (G[..., 0] + G[..., -1]).astype(complex128)
        YY = (G[..., 0] + G[..., -1]).astype(complex128)
        XY = 1j * (G[..., 0] - G[..., -1]).astype(complex128)
        YX = 1j * (G[..., -1] - G[..., 0]).astype(complex128)

        if G.shape[-1] == 4:
            XX += (G[..., 2] + G[..., 1]).astype(complex128)
            YY -= (G[..., 1] + G[..., 2]).astype(complex128)
            XY += 1j * (G[..., 2] - G[..., 1])
            YX += 1j * (G[..., 2] - G[..., 1])

        XX /= 2
        YY /= 2
        XY /= 2
        YX /= 2

        G_new = zeros(G.shape[0:-1] + (4,)).astype(complex128)
        G_new[..., 0] += XX
        G_new[..., 1] += XY
        G_new[..., 2] += YX
        G_new[..., 3] += YY
        return G_new

    @staticmethod
    def add_polarization(values, dim_pol):
        """
        Add extra polarization if there is no polarization
        :param values: values which need to get a polarization
        :param dim_pol: number of dimensions
        """
        values_temp = ones(values.shape+(dim_pol,))
        for i in range(dim_pol):
            values_temp[..., i] = values

        return values_temp


    def make_template(self, soltab):
        """
        Make template of the Gain matrix with only ones
        :param soltab: solution table (phase, amplitude)
        """
        self.G, self.axes_vals = array([]), {}
        for ss in self.h5_in.getSolsetNames():
            for st in self.h5_in.getSolset(ss).getSoltabNames():
                solutiontable = self.h5_in.getSolset(ss).getSoltab(st)
                if soltab in st:
                    try:
                        if 'pol' in solutiontable.getAxesNames():
                            values = reorderAxes(solutiontable.getValues()[0], solutiontable.getAxesNames(), self.axes_names)
                            self.G = ones(values.shape).astype(complex128)
                        else:
                            values = reorderAxes(solutiontable.getValues()[0], solutiontable.getAxesNames(), self.axes_names[0:-1])
                            self.G = ones(values.shape+(2,)).astype(complex128)
                    except:
                        sys.exit('ERROR:\nReceived '+str(solutiontable.getAxesNames())+', but expect at least [time, freq, ant, dir] or [time, freq, ant, dir, pol]')

                    self.axes_vals = {'time': solutiontable.getAxisValues('time'),
                                 'freq': solutiontable.getAxisValues('freq'),
                                 'ant': solutiontable.getAxisValues('ant'),
                                 'dir': solutiontable.getAxisValues('dir'),
                                 'pol': ['XX', 'XY', 'YX', 'YY']}
                    break

        print('Shape of input {shape}'.format(shape=self.G.shape))
        return self

    def add_tec(self, solutiontable):
        """
        :param solutiontable: the solution table for the TEC
        """
        tec_axes_names = [ax for ax in self.axes_names if solutiontable.getAxesNames()]
        tec = reorderAxes(solutiontable.getValues()[0], solutiontable.getAxesNames(), tec_axes_names)
        if 'freq' in solutiontable.getAxesNames():
            axes_vals_tec = {'time': solutiontable.getAxisValues('time'),
                             'freq': solutiontable.getAxisValues('freq'),
                             'ant': solutiontable.getAxisValues('ant'),
                             'dir': solutiontable.getAxisValues('dir')}
        else:
            axes_vals_tec = {'dir': solutiontable.getAxisValues('dir'),
                             'ant': solutiontable.getAxisValues('ant'),
                             'time': solutiontable.getAxisValues('time')}
        if 'pol' in solutiontable.getAxesNames():
            if tec.shape[-1] == 2:
                axes_vals_tec.update({'pol': ['XX', 'YY']})
            elif tec.shape[-1] == 4:
                axes_vals_tec.update({'pol': ['XX', 'XY', 'YX', 'YY']})
        axes_vals_tec = [v[1] for v in
                         sorted(axes_vals_tec.items(), key=lambda pair: self.axes_names.index(pair[0]))]
        self.solsetout.makeSoltab('tec', axesNames=tec_axes_names, axesVals=axes_vals_tec, vals=tec, weights=ones(tec.shape))

    def make_new_gains(self, lin2circ, circ2lin):
        """
        :param lin2circ: boolean for linear to circular conversion
        :param circ2lin: boolean for circular to linear conversion
        """
        for ss in self.h5_in.getSolsetNames():

            self.solsetout = self.h5_out.makeSolset(ss)
            solsetin = self.h5_in.getSolset(ss)

            self.solsetout.obj.source.append(solsetin.obj.source[:])

            for st in self.h5_in.getSolset(ss).getSoltabNames():
                solutiontable = self.h5_in.getSolset(ss).getSoltab(st)
                print('Reading {st} from {ss}'.format(ss=ss, st=st))
                if 'phase' in st:
                    if 'pol' in solutiontable.getAxesNames():
                        values = reorderAxes(solutiontable.getValues()[0], solutiontable.getAxesNames(), self.axes_names)
                        self.G *= exp(values * 1j)
                    else:
                        values = reorderAxes(solutiontable.getValues()[0], solutiontable.getAxesNames(), self.axes_names[0:-1])
                        self.G *= exp(self.add_polarization(values, 2) * 1j)

                elif 'amplitude' in st:
                    if 'pol' in solutiontable.getAxesNames():
                        values = reorderAxes(solutiontable.getValues()[0], solutiontable.getAxesNames(), self.axes_names)
                        self.G *= values * 1j
                    else:
                        values = reorderAxes(solutiontable.getValues()[0], solutiontable.getAxesNames(), self.axes_names[0:-1])
                        self.G *= self.add_polarization(values, 2)

                elif 'tec' in st:
                    self.add_tec(solutiontable)
                else:
                    print("Didn't include {st} in this version yet".format(st=st))
                    print("Let me (Jurjen) know if you need to include this.")

            if lin2circ:
                print('Convert linear polarization to circular polarization')
                G_new = self.lin2circ(self.G)
            elif circ2lin:
                print('Convert circular polarization to linear polarization')
                G_new = self.circ2lin(self.G)
            else:
                sys.exit('ERROR: No conversion given')
            print('Shape of output for amplitude and phase: {shape}'.format(shape=G_new.shape))

            phase = angle(G_new)
            amplitude = abs(G_new)

            self.axes_vals = [v[1] for v in sorted(self.axes_vals.items(), key=lambda pair: self.axes_names.index(pair[0]))]

            self.solsetout.makeSoltab('phase', axesNames=self.axes_names, axesVals=self.axes_vals, vals=phase, weights=ones(phase.shape))
            print('Created new phase solutions')

            self.solsetout.makeSoltab('amplitude', axesNames=self.axes_names, axesVals=self.axes_vals, vals=amplitude, weights=ones(amplitude.shape))
            print('Created new amplitude solutions')

        self.h5_in.close()
        self.h5_out.close()

        return self


def merge_h5(h5_out=None, h5_tables=None, ms_files=None, h5_time_freq=None, convert_tec=True, merge_all_in_one=False,
             lin2circ=False, circ2lin=False, add_directions=None, single_pol=None, no_pol=None, use_solset='sol000'):
    """
    Main function that uses the class MergeH5 to merge h5 tables.
    :param h5_out (string): h5 table name out
    :param h5_tables (string or list): h5 tables to merge
    :param ms_files (string or list): ms files to take freq and time axis from
    :param h5_time_freq (str or list): h5 file to take freq and time axis from
    :param convert_tec (boolean): convert TEC to phase or not
    :param merge_all_in_one: merge all in one direction
    :param lin2circ: boolean for linear to circular conversion
    :param circ2lin: boolean for circular to linear conversion
    :param add_directions: add default directions by giving a list of directions (coordinates)
    :param single_pol: only one polarization
    :param use_solset: use specific solset number
    """

    h5_out = make_h5_name(h5_out)

    #if alternative solset number is given, we will make a temp h5 file that has the alternative solset number because the code runs on sol000
    if use_solset!='sol000':
        for h5_ind, h5 in enumerate(h5_tables):
            temph5 = h5.replace('.h5', 'temp.h5')
            os.system('cp '+h5+' '+temph5)
            change_solset(temph5, use_solset, 'sol000')
            h5_tables[h5_ind] = temph5

    if h5_out.split('/')[-1] in [f.split('/')[-1] for f in glob(h5_out)]:
        os.system('rm {}'.format(h5_out))
    merge = MergeH5(h5_out=h5_out, h5_tables=h5_tables, ms_files=ms_files, convert_tec=convert_tec,
                    merge_all_in_one=merge_all_in_one, h5_time_freq=h5_time_freq)

    merge.get_allkeys()

    for st_group in merge.all_soltabs:
        if len(st_group) > 0:
            for st in st_group:
                merge.get_model_h5('sol000', st)
                merge.merge_files('sol000', st)
            if merge.convert_tec and (('phase' in st_group[0]) or ('tec' in st_group[0])):
                merge.create_new_dataset('sol000', 'phase')
            else:
                merge.create_new_dataset('sol000', st)
    # try:#add amplitude and phase if not available in h5 table
    if 'amplitude000' not in [item for sublist in merge.all_soltabs for item in sublist]:
        merge.gains = ones(
            (2, len(merge.directions.keys()), len(merge.antennas), len(merge.ax_freq), len(merge.ax_time)))
        merge.axes_new = ['time', 'freq', 'ant', 'dir', 'pol']
        merge.polarizations = ['XX', 'YY']
        merge.gains = reorderAxes(merge.gains, merge.solaxnames, merge.axes_new)
        merge.create_new_dataset('sol000', 'amplitude')
        # if 'phase000' not in [item for sublist in merge.all_soltabs for item in sublist] and \
        # 'tec000' not in [item for sublist in merge.all_soltabs for item in sublist]:
        # merge.phases = zeros((2, len(merge.directions.keys()), len(merge.antennas), len(merge.ax_freq), len(merge.ax_time)))
        # merge.axes_new = ['time', 'freq', 'ant', 'dir', 'pol']
        # merge.polarizations = ['XX', 'YY']
        # merge.phases = reorderAxes(merge.phases, merge.solaxnames, merge.axes_new)
        # merge.create_new_dataset(ss, 'phase')
    # except:#add try to except to be sure that adding extra phase and amplitude is not going to break the code
    # pass

    #Add antennas
    merge.add_antennas()
    print('END: h5 solution file(s) merged')

    if add_directions:
        merge.add_directions(add_directions)

    if lin2circ and circ2lin:
        sys.exit('Both polarization conversions are given, please choose 1.')
    elif lin2circ or circ2lin:

        print('THIS FEATURE IS NOT PROPERLY TESTED YET, PLEASE GIVE FEEDBACK IF THE RESULT IS NOT AS EXPECTED (jurjendejong@strw.leidenuniv.nl)')

        if lin2circ:
            h5_polchange = h5_out[0:-3]+'_circ.h5'
            print('Polarization will be converted from linear to circular')
        else:
            h5_polchange = h5_out[0:-3]+'_lin.h5'
            print('Polarization will be converted from circular to linear')

        Pol = PolChange(h5_in=h5_out, h5_out=h5_polchange)

        Pol.make_template('phase')
        if len(Pol.G.shape) > 1:
            Pol.make_template('amplitude')

        Pol.make_new_gains(lin2circ, circ2lin)
        print('{file} has been created'.format(file=h5_polchange))

    if sys.version_info.major == 2:
        print('You are using python 2. For this version we need to do an extra reordering step.')
        merge.order_directions()

    #Check table source size
    T = tables.open_file(merge.file)
    first_source = T.root.sol000.source[:][0]
    T.close()
    if first_source.nbytes>140:
        print('The source table memory size is too big. We will change the dtype to reduce size (probably a Python 3 issue).')
        merge.reduce_memory_source()

    #remove polarization axis if double
    if single_pol:
        print('Make a single polarization')
        merge.remove_pol(single=True)
    elif no_pol:
        print('Remove polarization')
        merge.remove_pol()

    if use_solset!='sol000':
        for h5 in h5_tables:
            if 'temp.h5' in h5:
                os.system('rm '+h5)



if __name__ == '__main__':
    from argparse import ArgumentParser

    def str2bool(v):
        v = str(v)
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            return True

    parser = ArgumentParser()
    parser.add_argument('-out', '--h5_out', type=str, help='h5 table name for output', required=True)
    parser.add_argument('-in', '--h5_tables', type=str, nargs='+', help='h5 tables to merge', required=True)
    parser.add_argument('-ms', '--ms_time_freq', type=str, help='ms files to use time and frequency arrays from')
    parser.add_argument('--h5_time_freq', type=str, help='h5 file to use time and frequency arrays from')
    parser.add_argument('-ct', '--convert_tec', type=str2bool, nargs='?', const=True, default=True, help='convert tec to phase')
    parser.add_argument('--merge_all_in_one', action='store_true', help='merge all solutions in one direction')
    parser.add_argument('--lin2circ', action='store_true', help='transform linear polarization to circular')
    parser.add_argument('--circ2lin', action='store_true', help='transform circular polarization to linear')
    parser.add_argument('--add_direction', default=None, help='add direction with amplitude 1 and phase 0 [ex: --add_direction [0.73,0.12]')
    parser.add_argument('--single_pol', action='store_true', default=None, help='Return only a single polarization axis if both polarizations are the same.')
    parser.add_argument('--no_pol', action='store_true', default=None, help='Remove polarization axis if both polarizations are the same.')
    parser.add_argument('--combine_h5', action='store_true', default=None, help='Combine h5 with different time axis into 1.')
    parser.add_argument('--usesolset', type=str, default='sol000', help='Choose a solset to merge from your input h5 files.')
    args = parser.parse_args()

    # check if solset name is accepted
    if 'sol' not in args.usesolset or sum(c.isdigit() for c in args.usesolset)!=3:
        sys.exit(args.usesolse+' not an accepted name. Only sol000, sol001, sol002, ... are accepted names for solsets.')

    # make sure h5 tables in right format
    if '[' in args.h5_tables:
        h5tables = args.h5_tables.replace('[', '').replace(']', '').replace(' ', '').split(',')
    elif ' ' in args.h5_tables:
        h5tables = args.h5_tables.split()
    else:
        h5tables = args.h5_tables

    if type(h5tables)==str:
        h5tables = glob(h5tables)
    elif type(h5tables)==list and len(h5tables)==1:
        h5tables = glob(h5tables[0])
    elif type(h5tables)==list:
        h5tablestemp=[]
        for h5 in h5tables:
            h5tablestemp+=glob(h5)

    if args.add_direction:
        add_directions = args.add_direction.replace('[','').replace(']','').split(',')
        add_directions = [float(add_directions[0]), float(add_directions[1])]
        if add_directions[0]>pi*6 or add_directions[1]>pi*6:
            sys.exit('ERROR: Please give values in radian')
    else:
        add_directions = None


    merge_h5(h5_out=args.h5_out,
             h5_tables=h5tables,
             ms_files=args.ms_time_freq,
             h5_time_freq=args.h5_time_freq,
             convert_tec=args.convert_tec,
             merge_all_in_one=args.merge_all_in_one,
             lin2circ=args.lin2circ,
             circ2lin=args.circ2lin,
             add_directions=add_directions,
             single_pol=args.single_pol,
             no_pol=args.no_pol,
             use_solset=args.usesolset)