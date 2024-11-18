"""Subtraction script, developed for de Jong et al. (2024) --> https://arxiv.org/pdf/2407.13247"""

__author__ = "Jurjen de Jong"

import os
import random
import re
import shutil
import sys
from argparse import ArgumentParser
from glob import glob
from itertools import repeat

import numpy as np
import pandas as pd
import pyregion
import tables
from astropy.io import fits
from astropy.wcs import WCS
from casacore.tables import table, taql


def fast_copy(filein, dest):
    """
    Fast copy of folders

    Args:
        filein: file in
        dest: out destination
    """
    available_cores = os.cpu_count()

    # Prepare the destination path
    if os.path.isdir(filein):
        os.system(f'cp -r {filein} {dest}')
        return

        #TODO: Figure out how to do below for MS
        if os.path.exists(dest):
            dest = os.path.join(dest, os.path.basename(filein))

        cmd = f"mkdir -p {dest} && "

        # Copy files in parallel using the available cores
        cmd += f"find {filein} -type f | parallel -j {available_cores} cp -rp {{}} {dest}"

        # Run the final command
        os.system(cmd)

    else:
        # regular copy
        os.system(f'cp {filein} {dest}')


def add_trailing_zeros(s, digitsize=4):
    """
     Repeat the zero character and add it to front

     :param s: string
     :param digitsize: number of digits (default 4)
     :return: trailing zeros + number --> example: 0001, 0021, ...
     """
    padded_string = "".join(repeat("0", digitsize)) + s
    return padded_string[-digitsize:]


def unlink(symlink_path):
    """
    Replaces a symbolic link with the actual data it points to.

    :param symlink_path: Path to the symbolic link to be replaced
    """
    try:
        # Check if the provided path is a symbolic link
        if os.path.islink(symlink_path):
            # Get the actual path the symlink points to
            target_path = os.readlink(symlink_path)

            # Remove the symlink
            os.unlink(symlink_path)
            print(f"Symlink '{symlink_path}' removed.")

            # Copy the data from the target path to the symlink location
            if os.path.isdir(target_path):
                shutil.copytree(target_path, symlink_path)
                print(f"Directory '{target_path}' copied to '{symlink_path}'.")
            else:
                shutil.copy2(target_path, symlink_path)
                print(f"File '{target_path}' copied to '{symlink_path}'.")
        else:
            print(f"'{symlink_path}' is not a symbolic link.")
    except FileNotFoundError:
        print(f"The symlink '{symlink_path}' does not exist.")
    except OSError as e:
        print(f"Error: {e} - Could not replace the symlink with data.")


def get_largest_divider(inp, max=1000):
    """
    Get largest divider

    :param inp: input number
    :param max: max divider

    :return: largest divider from inp bound by max
    """
    for r in range(max)[::-1]:
        if inp % r == 0:
            return r
    sys.exit("ERROR: code should not arrive here.")


def parse_history(ms, hist_item):
    """
    Grep specific history item from MS

    :param ms: measurement set
    :param hist_item: history item

    :return: parsed string
    """
    hist = os.popen('taql "SELECT * FROM ' + ms + '::HISTORY" | grep ' + hist_item).read().split(' ')
    for item in hist:
        if hist_item in item and len(hist_item) <= len(item):
            return item
    print('WARNING:' + hist_item + ' not found')
    return None


def make_utf8(inp):
    """
    Convert input to utf8 instead of bytes

    :param inp: string input
    :return: input in utf-8 format
    """

    try:
        inp = inp.decode('utf8')
        return inp
    except (UnicodeDecodeError, AttributeError):
        return inp


def get_time_preavg_factor(ms: str = None):
    """
    Get time pre-averaging factor (given by demixer.timestep)

    :param ms: measurement set

    :return: averaging integer
    """
    parse_str = "demixer.timestep="
    parsed_history = parse_history(ms, parse_str)
    avg_num = re.findall(r'\d+', parsed_history.replace(parse_str, ''))[0]
    if avg_num.isdigit():
        factor = int(float(avg_num))
        if factor != 1:
            print("WARNING: " + ms + " time has been pre-averaged with factor " + str(
                factor) + ". This might cause stronger time smearing effects in your final image.")
        return factor
    elif type(avg_num)==float:
        factor = float(avg_num)
        print("WARNING: parsed factor in " + ms + " is not a digit but a float")
        return factor
    else:
        print("WARNING: parsed factor in " + ms + " is not a float or digit")
        return None


def split_facet_h5(h5parm: str = None, dirname: str = None):
    """
    Split multi-facet h5parm

    :param h5parm: multi-facet h5parm
    :param dirname: direction name
    """

    outputh5 = f'{h5parm}.{dirname}.h5'
    os.system(f'cp {h5parm} {outputh5}')

    with tables.open_file(outputh5, 'r+') as outh5:

        axes = make_utf8(outh5.root.sol000.phase000.val.attrs["AXES"])

        dir_axis = axes.split(',').index('dir')

        sources = outh5.root.sol000.source[:]
        dirs = [make_utf8(dir) for dir in sources['name']]
        dir_idx = dirs.index(dirname)

        def get_data(soltab, axis):
            return np.take(outh5.root.sol000._f_get_child(soltab)._f_get_child(axis)[:], indices=[dir_idx], axis=dir_axis)

        phase_w = get_data('phase000', 'weight')
        amplitude_w = get_data('amplitude000', 'weight')
        phase_v = get_data('phase000', 'val')
        amplitude_v = get_data('amplitude000', 'val')
        new_dirs = np.array([outh5.root.sol000.source[:][dir_idx]])
        # new_dirs['name'][0] = bytes('Dir' + str(0).zfill(2), 'utf-8')
        dirs = np.array([outh5.root.sol000.phase000.dir[:][dir_idx]])

        outh5.remove_node("/sol000/phase000", "val", recursive=True)
        outh5.remove_node("/sol000/phase000", "weight", recursive=True)
        outh5.remove_node("/sol000/phase000", "dir", recursive=True)
        outh5.remove_node("/sol000/amplitude000", "val", recursive=True)
        outh5.remove_node("/sol000/amplitude000", "weight", recursive=True)
        outh5.remove_node("/sol000/amplitude000", "dir", recursive=True)
        outh5.remove_node("/sol000", "source", recursive=True)

        outh5.create_array('/sol000/phase000', 'val', phase_v)
        outh5.create_array('/sol000/phase000', 'weight', phase_w)
        outh5.create_array('/sol000/phase000', 'dir', dirs)
        outh5.create_array('/sol000/amplitude000', 'val', amplitude_v)
        outh5.create_array('/sol000/amplitude000', 'weight', amplitude_w)
        outh5.create_array('/sol000/amplitude000', 'dir', dirs)
        outh5.create_table('/sol000', 'source', new_dirs, title='Source names and directions')

        outh5.root.sol000.phase000.val.attrs['AXES'] = bytes(axes, 'utf-8')
        outh5.root.sol000.phase000.weight.attrs['AXES'] = bytes(axes, 'utf-8')
        outh5.root.sol000.amplitude000.val.attrs['AXES'] = bytes(axes, 'utf-8')
        outh5.root.sol000.amplitude000.weight.attrs['AXES'] = bytes(axes, 'utf-8')

    return outputh5


def isfulljones(h5: str = None):
    """
    Verify if file is fulljones

    :param h5: h5 file
    """
    T = tables.open_file(h5)
    soltab = list(T.root.sol000._v_groups.keys())[0]
    if 'pol' in make_utf8(T.root.sol000._f_get_child(soltab).val.attrs["AXES"]):
        if T.root.sol000._f_get_child(soltab).pol[:].shape[0] == 4:
            T.close()
            return True
    T.close()
    return False


class SubtractWSClean:
    def __init__(self, mslist: list = None, region: str = None, localnorth: bool = True, inverse: bool = False):
        """
        Subtract image with WSClean

        :param mslist: measurement set list
        :param region: region file to mask
        :param model_image: model image
        """

        # list with MS
        self.mslist = mslist

        # wsclean model image
        model_images = glob('*-model*.fits')
        hdu = fits.open(model_images[0])
        self.imshape = (hdu[0].header['NAXIS1'], hdu[0].header['NAXIS2'])

        if len(glob('*-????-model-pb.fits')) >= 1:
            self.model_images = glob('*-????-model-pb.fits')
        elif len(glob('*-????-model.fits')) >= 1:
            self.model_images = glob('*-????-model.fits')
        elif len(glob('*-model*.fits')) >= 1:
            self.model_images = glob('*-model*.fits')
        f = fits.open(self.model_images[0])
        history = str(f[0].header['HISTORY']).replace('\n', '').split()
        if '-taper-gaussian' in history:
            self.scale = history[history.index('-taper-gaussian') + 1]
        elif '-scale' in history:
            self.scale = history[history.index('-scale') + 1]
        else:
            self.scale = ''

        # region file to mask
        if region is None:
            self.region = None
        elif localnorth:
            self.region = self.box_to_localnorth(region)
        else:
            self.region = pyregion.open(region)

        self.inverse = inverse

    def clean_model_images(self):
        """
        Delete model images that do not match with MS
        """

        freqs = []
        for ms in self.mslist:
            with table(ms + "::SPECTRAL_WINDOW", ack=False) as t:
                freqs += list(t.getcol("CHAN_FREQ")[0])
        self.fmax_ms = max(freqs)
        self.fmin_ms = min(freqs)
        model_images = glob('*-model*.fits')
        for modim in model_images:
            fts = fits.open(modim)[0]
            fdelt, fcent = fts.header['CDELT3'] / 2, fts.header['CRVAL3']
            fmin, fmax = fcent - fdelt, fcent + fdelt
            if fmin > self.fmax_ms or fmax < self.fmin_ms:
                print(modim + ' does not overlap with MS bandwidth --> DELETE')
                os.system('rm ' + modim)
            else:
                print(modim + ' overlaps with MS bandwidth --> KEEP')

        # rename and resort model images --> remove trailing zeros when only 1 model image, otherwise renumber model images
        if len(glob('*-????-model.fits')) > 1:
            for n, modim in enumerate(sorted(glob('*-????-model.fits'))):
                os.system('mv ' + modim + ' ' + re.sub(r'\d{4}', add_trailing_zeros(str(n), 4), modim))
        elif len(glob('*-????-model.fits')) == 1:
            for n, modim in enumerate(sorted(glob('*-????-model.fits'))):
                os.system('mv ' + modim + ' ' + re.sub(r'\-\d{4}', '', glob('*-????-model.fits')[0]))
        if len(glob('*-????-model-pb.fits')) > 1:
            for n, modim in enumerate(sorted(glob('*-????-model-pb.fits'))):
                os.system('mv ' + modim + ' ' + re.sub(r'\d{4}', add_trailing_zeros(str(n), 4), modim))
        elif len(glob('*-????-model-pb.fits')) == 1:
            for n, modim in enumerate(sorted(glob('*-????-model-pb.fits'))):
                os.system('mv ' + modim + ' ' + re.sub(r'\-\d{4}', '', glob('*-????-model-pb.fits')[0]))

        # select correct model images
        if len(glob('*-????-model-pb.fits')) >= 1:
            self.model_images = glob('*-????-model-pb.fits')
        elif len(glob('*-????-model.fits')) >= 1:
            self.model_images = glob('*-????-model.fits')
        elif len(glob('*-model-pb.fits')) >= 1:
            self.model_images = glob('*-model-pb.fits')
        elif len(glob('*-model.fits')) >= 1:
            self.model_images = glob('*-model.fits')
        else:
            sys.exit("ERROR: No model images found.")

        return self

    def box_to_localnorth(self, region: str = None):
        """
        Adjust box for local north
        :param regionfile: region file
        :param fitsfile: fits file
        """
        r = pyregion.open(region)

        hduflat = self.flat_model_image(self.model_images[0])
        CRVAL1 = hduflat.header['CRVAL1']  # RA

        if len(r[:]) > 1:
            print('Only one region can be specified, your file contains', len(r[:]))
            sys.exit()

        if r[0].name != 'box':
            print('Only box region supported')
            return pyregion.open(region)

        ra = r[0].coord_list[0]

        r[0].coord_list[4] = CRVAL1 - ra  # rotate box
        print('Angle adjusted box', CRVAL1 - ra)

        if os.path.isfile('adjustedbox.reg'):
            os.system('rm -rf adjustedbox.reg')

        r.write("adjustedbox.reg")
        return pyregion.open("adjustedbox.reg")

    @staticmethod
    def flat_model_image(fitsfile: str = None):
        """
        Flatten a fits file so that it becomes a 2D image. Return new header and data
        (taken from sub_sources_outside_region.py)
        """
        hdu = fits.open(fitsfile)
        naxis = hdu[0].header['NAXIS']
        if naxis < 2:
            raise sys.exit('Can\'t make map from this')
        if naxis == 2:
            return fits.PrimaryHDU(header=hdu[0].header, data=hdu[0].data)

        w = WCS(hdu[0].header)
        wn = WCS(naxis=2)

        wn.wcs.crpix[0] = w.wcs.crpix[0]
        wn.wcs.crpix[1] = w.wcs.crpix[1]
        wn.wcs.cdelt = w.wcs.cdelt[0:2]
        wn.wcs.crval = w.wcs.crval[0:2]
        wn.wcs.ctype[0] = w.wcs.ctype[0]
        wn.wcs.ctype[1] = w.wcs.ctype[1]

        header = wn.to_header()
        header["NAXIS"] = 2
        copy = ('EQUINOX', 'EPOCH', 'BMAJ', 'BMIN', 'BPA', 'RESTFRQ', 'TELESCOP', 'OBSERVER')
        for k in copy:
            r = hdu[0].header.get(k)
            if r is not None:
                header[k] = r

        slice = []
        for i in range(naxis, 0, -1):
            if i <= 2:
                slice.append(np.s_[:], )
            else:
                slice.append(0)

        hdu = fits.PrimaryHDU(header=header, data=hdu[0].data[tuple(slice)])
        return hdu

    @staticmethod
    def invert_mask(mask):
        """
        invert mask (True=False and vice versa)
        :param mask: mask
        :return: inverted mask
        """
        return np.invert(mask)

    def mask_region(self, region_cube: bool = False):
        """
        :param region_cube: if region_cube make cube, otherwise 2D (flatten)
        """

        for fits_model in self.model_images:

            print('Mask ' + fits_model)


            with fits.open(fits_model) as hdu:

                b = not self.inverse

                if region_cube:
                    manualmask = self.region.get_mask(hdu=hdu[0], shape=self.imshape)
                    for i in range(hdu[0].header['NAXIS4']):
                        hdu[0].data[i][0][np.where(manualmask == b)] = 0.0
                    hdu.writeto(fits_model, overwrite=True)

                else:
                    hduflat = self.flat_model_image(fits_model)
                    manualmask = self.region.get_mask(hdu=hduflat)
                    hdu[0].data[0][0][np.where(manualmask == b)] = 0.0
                    hdu.writeto(fits_model, overwrite=True)

        return self

    def subtract_col(self, out_column: str = None):

        """
        Subtract column in MeasurementSet

        :param out_column: out column name
        """

        for ms in self.mslist:
            with table(ms, readonly=False, ack=False) as ts:
                colnames = ts.colnames()

                if "MODEL_DATA" not in colnames:
                    sys.exit(
                        f"ERROR: MODEL_DATA does not exist in {ms}.\nThis is most likely due to a failed predict step.")

            if out_column not in colnames:
                # Creating the column with DP3 ensures we can directly compress the data
                os.system(f"DP3 msin={ms} msout=. msout.datacolumn={out_column} steps=[] msout.storagemanager=dysco"
                          f"msout.storagemanager.databitrate=6 msout.storagemanager.weightbitrate=6")
            else:
                print(out_column, ' already exists')

            if 'SUBTRACT_DATA' in colnames:
                colmn = 'SUBTRACT_DATA'
            else:
                colmn = 'DATA'

            if self.inverse:
                sign = '+'
            else:
                sign = '-'

            print(f'Subtracting --> {colmn} {sign} MODEL_DATA for {ms}')

            # subtraction or addition
            taql(f"UPDATE {ms} SET {out_column}={colmn}{sign}MODEL_DATA")

            # remove MODEL_DATA to save memory
            taql(f"ALTER TABLE {ms} DROP COLUMN MODEL_DATA")

        return self

    def predict(self, h5parm: str = None, facet_regions: str = None):
        """
        Predict image by reading history from model images

        :param h5parm: h5 solutions (optional)
        :param facet_regions: facet regions (if h5 solutions given)
        """

        f = fits.open(self.model_images[0])
        comparse = str(f[0].header['HISTORY']).replace('\n', '').split()
        command = ['wsclean', '-predict', f'-name {self.model_images[0].split("-")[0]}', '-mem 50',
                   '-parallel-gridding 2']

        for n, argument in enumerate(comparse):
            if argument in ['-gridder', '-padding',
                            '-idg-mode', '-beam-aterm-update', '-pol', '-scale']:
                if ' '.join(comparse[n:n + 2]) == '-gridder wgridder-apply-primary-beam':
                    command.append('-gridder wgridder')
                    command.append('-apply-primary-beam')
                else:
                    command.append(' '.join(comparse[n:n + 2]))
            elif argument in ['-size']:
                command.append(' '.join(comparse[n:n + 3]))
            elif argument in ['-use-differential-lofar-beam', '-grid-with-beam',
                              '-use-idg', '-log-time', '-gap-channel-division',
                              '-apply-primary-beam']:
                if argument not in command:
                    command.append(argument)
            if argument == '-taper-gaussian':
                self.scale = comparse[n + 1]
            elif argument == '-scale' and '-taper-gaussian' not in comparse:
                self.scale = comparse[n + 1]

        if len(self.model_images) > 1:
            command += ['-channels-out ' + str(len(self.model_images))]

        freqboundary = []
        for modim in sorted(self.model_images)[:-1]:
            fts = fits.open(modim)[0]
            fdelt, fcent = fts.header['CDELT3'] / 2, fts.header['CRVAL3']

            freqboundary.append(str(int(fcent + fdelt)))
            # fts.close()

        if len(freqboundary) > 0:
            command += ['-channel-division-frequencies ' + ','.join(freqboundary)]

        if h5parm is not None:
            command += [f'-apply-facet-solutions {h5parm} amplitude000,phase000',
                        f'-facet-regions {facet_regions}', '-apply-facet-beam',
                        f'-facet-beam-update {comparse[comparse.index("-facet-beam-update") + 1]}']

        command += [' '.join(self.mslist)]

        # run
        print('\n'.join(command))
        predict_cmd = open("predict.cmd", "w")
        predict_cmd.write('\n'.join(command))
        predict_cmd.close()

        os.system(' '.join(command) + ' > predict.log')

        return self


    def run_DP3(self, phaseshift: str = None, freqavg: str = None,
                timeres: str = None, concat: bool = None, speedup_facet_subtract: bool = None,
                applybeam: bool = None, applycal_h5: str = None, dirname: str = None):
        """
        Run DP3 command

        :param phaseshift: do phase shift to specific center
        :param freqavg: frequency averaging
        :param timeres: time resolution in seconds
        :param concat: concat the measurement sets
        :param applybeam: apply beam in phaseshifted phase center (or otherwise center of field)
        :param applycal_h5: applycal solution file
        :param dirname: direction name from h5parm
        """

        steps = []

        command = ['DP3',
                   'msin.missingdata=True',
                   'msin.datacolumn=SUBTRACT_DATA' if not self.inverse else 'msin.datacolumn=DATA',
                   'msin.orderms=False',
                   'msout.storagemanager=dysco',
                   'msout.storagemanager.databitrate=6',
                   'msout.storagemanager.weightbitrate=6']

        # 1) PHASESHIFT
        if phaseshift is not None:
            phasecenter = phaseshift.replace('[', '').replace(']', '').split(',')
            phasecenter = f'[{phasecenter[0]},{phasecenter[1]}]'
            steps.append('ps')
            command += ['ps.type=phaseshifter',
                        'ps.phasecenter=' + phasecenter]

        # For facet subtraction we can do first averaging before applying beam and solutions (saves compute time)
        if speedup_facet_subtract:
            # 2) AVERAGING
            if freqavg is not None or timeres is not None:
                steps.append('avg')
                command += ['avg.type=averager']
                if freqavg is not None:
                    if str(freqavg).isdigit() or not str(freqavg)[-1].isalpha():
                        command += [f'avg.freqstep={int(freqavg)}']
                    else:
                        command += [f'avg.freqresolution={freqavg}']
                if timeres is not None:
                    command += [f'avg.timeresolution={timeres}']

            # 3) APPLY BEAM
            if applybeam:
                steps.append('beam1')
                command += ['beam1.type=applybeam',
                            'beam1.updateweights=True']
                if applycal_h5 is not None and phaseshift is not None and dirname is not None:
                    H = tables.open_file(applycal_h5)
                    sources = H.root.sol000.source[:]
                    H.close()
                    dirs = [make_utf8(dir) for dir in sources['name']]
                    dir_idx = dirs.index(dirname)
                    ra, dec = sources['dir'][dir_idx]
                    dir = str(f"[{round(ra,5)},{round(dec,5)}]")
                    command += ['beam1.direction='+dir]
                else:
                    command += ['beam1.direction=[]']

            # 4) APPLYCAL
            if applycal_h5 is not None:
                # add fulljones solutions apply
                if isfulljones(applycal_h5):
                    steps.append('ac')
                    command += ['ac.type=applycal',
                                'ac.parmdb=' + applycal_h5,
                                'ac.correction=fulljones',
                                'ac.soltab=[amplitude000,phase000]']
                    if phaseshift is not None and dirname is not None:
                        command += ['ac.direction=' + dirname]
                # add non-fulljones solutions apply
                else:
                    ac_count = 0
                    T = tables.open_file(applycal_h5)
                    for corr in T.root.sol000._v_groups.keys():
                        if 'phase' in corr or 'amp' in corr:
                            command += [f'ac{ac_count}.type=applycal',
                                        f'ac{ac_count}.parmdb={applycal_h5}',
                                        f'ac{ac_count}.correction={corr}']
                            if phaseshift is not None and dirname is not None:
                                command += [f'ac{ac_count}.direction=' + dirname]
                            steps.append(f'ac{ac_count}')
                            ac_count += 1
                    T.close()

            # 5) APPLY BEAM (OPTIONAL IF APPLY BEAM USED FOR APPLYCAL)
            if applybeam and applycal_h5 is not None and phaseshift is not None and dirname is not None:
                steps.append('beam2')
                command += ['beam2.type=applybeam',
                            'beam2.updateweights=True',
                            'beam2.direction=[]']

        # Most accurate version, but more expensive than the --speedup_facet_subtract approach
        else:
            # 2) APPLY BEAM
            if applybeam:
                steps.append('beam1')
                command += ['beam1.type=applybeam',
                            'beam1.updateweights=True']
                if applycal_h5 is not None and phaseshift is not None and dirname is not None:
                    H = tables.open_file(applycal_h5)
                    sources = H.root.sol000.source[:]
                    H.close()
                    dirs = [make_utf8(dir) for dir in sources['name']]
                    dir_idx = dirs.index(dirname)
                    ra, dec = sources['dir'][dir_idx]
                    dir = str(f"[{round(ra,5)},{round(dec,5)}]")
                    command += ['beam1.direction='+dir]
                else:
                    command += ['beam1.direction=[]']

            # 3) APPLYCAL
            if applycal_h5 is not None:
                # add fulljones solutions apply
                if isfulljones(applycal_h5):
                    steps.append('ac')
                    command += ['ac.type=applycal',
                                'ac.parmdb=' + applycal_h5,
                                'ac.correction=fulljones',
                                'ac.soltab=[amplitude000,phase000]']
                    if phaseshift is not None and dirname is not None:
                        command += ['ac.direction=' + dirname]
                # add non-fulljones solutions apply
                else:
                    ac_count = 0
                    T = tables.open_file(applycal_h5)
                    for corr in T.root.sol000._v_groups.keys():
                        if 'phase' in corr or 'amp' in corr:
                            command += [f'ac{ac_count}.type=applycal',
                                        f'ac{ac_count}.parmdb={applycal_h5}',
                                        f'ac{ac_count}.correction={corr}']
                            if phaseshift is not None and dirname is not None:
                                command += [f'ac{ac_count}.direction=' + dirname]
                            steps.append(f'ac{ac_count}')
                            ac_count += 1
                    T.close()

            # 4) APPLY BEAM (OPTIONAL IF APPLY BEAM USED FOR APPLYCAL)
            if applybeam and applycal_h5 is not None and phaseshift is not None and dirname is not None:
                steps.append('beam2')
                command += ['beam2.type=applybeam',
                            'beam2.updateweights=True',
                            'beam2.direction=[]']

            # 5) AVERAGING
            if freqavg is not None or timeres is not None:
                steps.append('avg')
                command += ['avg.type=averager']
                if freqavg is not None:
                    if str(freqavg).isdigit() or not str(freqavg)[-1].isalpha():
                        command += [f'avg.freqstep={int(freqavg)}']
                    else:
                        command += [f'avg.freqresolution={freqavg}']
                if timeres is not None:
                    command += [f'avg.timeresolution={timeres}']

        command += ['steps=' + str(steps).replace(" ", "").replace("\'", "")]

        if concat:
            msout = ['subtract_concat.ms']
            command += [f'msin={",".join(self.mslist)}',
                        'msout=subtract_concat.ms']

            print('\n'.join(command))
            dp3_cmd = open("dp3.cmd", "w")
            dp3_cmd.write('\n'.join(command))
            dp3_cmd.close()

            print(f"Make subtract_concat.ms")

            os.system(' '.join(command) + " > dp3.subtract.log")
        else:
            msout = []
            for n, ms in enumerate(self.mslist):
                command += [f'msin={ms}', f'msout=sub{self.scale}_{ms}']

                print('\n'.join(command))
                dp3_cmd = open("dp3.cmd", "w")
                dp3_cmd.write('\n'.join(command))
                dp3_cmd.close()

                print(f"Make sub{self.scale}_{ms}")

                os.system(' '.join(command + [f'msin={ms}', f'msout=sub{self.scale}_{ms}']) + f" > dp3.sub{n}.log")
            msout.append(f'sub{self.scale}_{ms}')

        return msout


def parse_args():
    """
    Command line argument parser
    """
    parser = ArgumentParser(description='Subtract or predict region with WSClean')
    parser.add_argument('--mslist', nargs='+', help='MeasurementSet(s)', required=True)
    parser.add_argument('--facets_predict', type=str, help='Multi-facet region file for prediction')
    parser.add_argument('--h5parm_predict', type=str, help='Multi-dir h5 solution file corresponding to --facets_predict')
    parser.add_argument('--region', type=str, help='Region file to mask for subtraction or predict back to data when --inverse')
    parser.add_argument('--model_image_folder', type=str, help='Directory with model images (if not given model images from run are selected)')
    parser.add_argument('--model_images', nargs='+', help='Instead of --model_image_folder, you can also specify the model images to use')
    parser.add_argument('--no_local_north', action='store_true', help='Do not move box to local north')
    parser.add_argument('--use_region_cube', action='store_true', help='Use region cube')
    parser.add_argument('--phasecenter', type=str, help='Phaseshift to given point (example: --phaseshift 16h06m07.61855,55d21m35.4166)')
    parser.add_argument('--freqavg', type=str, help='Frequency averaging')
    parser.add_argument('--timeres', type=str, help='Time resolution averaging in seconds')
    parser.add_argument('--concat', action='store_true', help='Concat MeasurementSets')
    parser.add_argument('--applybeam', action='store_true', help='Apply beam in phaseshift center or center of field')
    parser.add_argument('--applycal', action='store_true', help='Applycal solutions from facet after subtraction and phaseshifting')
    parser.add_argument('--applycal_h5', type=str, help='Applycal solution file')
    parser.add_argument('--forwidefield', action='store_true', help='Will search for the polygon_info.csv file to extract information from')
    parser.add_argument('--skip_predict', action='store_true', help='Skip predict and do only subtract')
    parser.add_argument('--even_time_avg', action='store_true', help='Only allow even time averaging (in case of combining observations with different averaging factors) and --forwidefield')
    parser.add_argument('--inverse', action='store_true', help='Instead of subtracting, you predict and add model data from a single facet')
    parser.add_argument('--scratch_toil', action='store_true', help='Experts only: Run on scratch when using toil')
    parser.add_argument('--speedup_facet_subtract', action='store_true', help='DP3 speedup for facet subtraction by performing averaging earlier (may introduce accuracy issues)')

    return parser.parse_args()


def main():
    """
    Main function
    """

    args = parse_args()

    if not args.skip_predict:

        # copy model images from model_image_folder
        if args.model_image_folder is not None:
            if len(glob(args.model_image_folder + '/*-????-model-pb.fits')) > 1:
                os.system('cp ' + args.model_image_folder + '/*-????-model-pb.fits .')
                if len(glob(args.model_image_folder + '/*-????-model.fits')) > 1:
                    os.system('cp ' + args.model_image_folder + '/*-????-model.fits .')
            elif len(glob(args.model_image_folder + '/*-????-model.fits')) > 1:
                os.system('cp ' + args.model_image_folder + '/*-????-model.fits .')
            elif len(glob(args.model_image_folder + '/*-model-pb.fits')) > 1:
                os.system('cp ' + args.model_image_folder + '/*-model-pb.fits .')
                if len(glob(args.model_image_folder + '/*-model.fits')) > 1:
                    os.system('cp ' + args.model_image_folder + '/*-model.fits .')
            elif len(glob(args.model_image_folder + '/*-model.fits')) > 1:
                os.system('cp ' + args.model_image_folder + '/*-model.fits .')
            else:
                sys.exit("ERROR: missing model images in folder " + args.model_image_folder)
        elif args.model_images is not None:
            for model in args.model_images:
                os.system('cp ' + model + ' .')

    # --forwidefield --> will read averaging and phasecenter from polygon_info.csv
    if args.forwidefield:
        if os.path.isfile('polygon_info.csv'):
            polygon_info = pd.read_csv('polygon_info.csv')
        elif os.path.isfile('../polygon_info.csv'):
            polygon_info = pd.read_csv('../polygon_info.csv')
        elif os.path.isfile('../../polygon_info.csv'):
            polygon_info = pd.read_csv('../../polygon_info.csv')
        elif os.path.isfile('../../../polygon_info.csv'):
            polygon_info = pd.read_csv('../../../polygon_info.csv')
        else:
            sys.exit('ERROR: using --forwidefield option needs polygon_info.csv file to read polygon information from')

        with table(args.mslist[0] + "::SPECTRAL_WINDOW", ack=False) as t:
            channum = len(t.getcol("CHAN_FREQ")[0])

        with table(args.mslist[0], ack=False) as t:
            time = np.unique(t.getcol("TIME"))
            dtime = abs(time[1] - time[0])

        polygon = polygon_info.loc[polygon_info.polygon_file == args.region.split('/')[-1]]
        try:
            phasecenter = polygon['poly_center'].values[0]
        except AttributeError:
            print('WARNING: no poly center in polygon_info.csv, use dir instead.')
            phasecenter = polygon['dir'].values[0]
        except KeyError:
            print('WARNING: no poly center in polygon_info.csv, use dir instead.')
            phasecenter = polygon['dir'].values[0]

        avg = int(polygon['avg'].values[0])

        # take only averaging factors that are channum%avg==0
        freqavg = get_largest_divider(channum, avg + 1)

        try:
            # if there is pre averaging done on the ms, we need to take this into account
            timeres = int(avg // get_time_preavg_factor(args.mslist[0]) * dtime)

        except:
            timeres = int(avg * dtime)

        # in case of widefield imaging and stacking multiple nights, you might want to have only even time resolutions
        if args.even_time_avg and timeres % 2 != 0:
            timeres -= 1
            print(f"Time resolution from {timeres + 1} to {timeres}")

        dirname = polygon['dir_name'].values[0]

    else:
        phasecenter = args.phasecenter
        freqavg = args.freqavg
        timeres = args.timeres
        dirname = None

    if args.scratch_toil:
        hashvalue = random.getrandbits(128)
        hasfolder = "%032x" % hashvalue
        hasfolder = hasfolder[0:10]
        absolute_path = os.path.abspath('/tmp')
        runpath = absolute_path+'/'+hasfolder

        # mkdir and copy files
        os.system(f'mkdir -p {runpath}')
        for dataset in args.mslist: fast_copy(dataset, runpath)
        command = [f'cp *-model*.fits {runpath}']
        if args.region is not None:
            command += [f'cp {args.region} {runpath}']
        # when running with scratch + toil, the next commands are to clean up the tmp* files
        command += ['rm *-model*.fits', f'rm -rf {args.model_image_folder}']
        if not args.applybeam and not args.applycal:
            command += [f'rm -rf {dataset}' for dataset in args.mslist]
        os.system('&&'.join(command))
        outpath = os.getcwd()
        os.chdir(runpath)

        # replace symlinks with data to correct
        for ms in args.mslist:
            unlink(ms.split('/')[-1])


    # set subtract object
    subpred = SubtractWSClean(mslist=args.mslist if not args.scratch_toil else [ms.split('/')[-1] for ms in args.mslist],
                             region=args.region if not args.scratch_toil or args.region is None else args.region.split('/')[-1],
                             localnorth=not args.no_local_north,
                             inverse=args.inverse)

    if not args.skip_predict:

        # clean model images
        subpred.clean_model_images()

        # mask
        print('############## MASK REGION ##############')
        if args.region is not None:
            subpred.mask_region(region_cube=args.use_region_cube)

        if args.scratch_toil:
            os.system(f'cp {outpath}/{args.h5parm_predict.split("/")[-1]} {runpath}')

        if args.inverse:
            faceth5 = split_facet_h5(h5parm=args.h5parm_predict if not args.scratch_toil else args.h5parm_predict.split('/')[-1],
                                            dirname=dirname)
            # predict with 1 facet
            print('############## PREDICT ##############')
            subpred.predict(h5parm=faceth5,
                           facet_regions=args.region if not args.scratch_toil else args.region.split('/')[-1])

        else:
            # predict with multiple facets
            print('############## PREDICT ##############')
            if args.scratch_toil:
                os.system(f'cp {outpath}/{args.facets_predict.split("/")[-1]} {runpath}')
            subpred.predict(h5parm=args.h5parm_predict if not args.scratch_toil else args.h5parm_predict.split('/')[-1],
                           facet_regions=args.facets_predict if not args.scratch_toil else args.facets_predict.split('/')[-1])

    # subtract or add (if inverse=True)
    print('############## SUBTRACT ##############')
    subpred.subtract_col(out_column='SUBTRACT_DATA' if not args.inverse else "DATA")

    # DP3 steps
    if args.phasecenter is not None or \
            args.freqavg is not None or \
            args.timeres is not None or \
            args.concat or \
            args.applybeam or \
            args.applycal:
        print('############## RUN DP3 ##############')
        if args.applycal_h5 is not None:
            applycalh5 = args.applycal_h5
        elif args.applycal and args.applycal_h5 is None and args.h5parm_predict is not None:
            applycalh5 = args.h5parm_predict
        elif args.applycal and not args.applycal_h5:
            sys.exit("ERROR: need a solution file for applycal (give with --applycal_h5)")
        else:
            applycalh5 = None

        if args.scratch_toil:
            os.system(f'cp {outpath}/{applycalh5.split("/")[-1]} {runpath}')

        # run DP3
        msout = subpred.run_DP3(phaseshift=phasecenter, freqavg=freqavg, timeres=timeres,
                       concat=args.concat, applybeam=args.applybeam, speedup_facet_subtract=args.speedup_facet_subtract,
                       applycal_h5=applycalh5 if not args.scratch_toil else applycalh5.split('/')[-1], dirname=dirname)

        if args.scratch_toil and args.forwidefield:
            # copy averaged MS back to output folder
            for ms in msout: fast_copy(ms, f'{outpath}/{dirname.replace("Dir","facet_")}-{ms.split("/")[-1]}')
            # clean up scratch directory (for big MS)
            os.system(f'cp *.log {outpath} && rm -rf *.ms')
            os.chdir(outpath)

        elif args.forwidefield:
            for ms in msout: os.system(f"mv {ms} {dirname.replace('Dir','facet_')}-{ms.split('/')[-1]}")

    elif args.scratch_toil:
        # copy back the subtracted MS to the output path
        for ms in subpred.mslist: fast_copy(ms, f'{outpath}/subfov_{ms.split("/")[-1]}')
        os.system(f'cp *.log {outpath}')
        os.chdir(outpath)
    else:
        for ms in subpred.mslist: os.system(f'mv {ms} subfov_{ms.split("/")[-1]}')

if __name__ == "__main__":
    main()
