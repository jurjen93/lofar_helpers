import numpy as np
import sys
import os
import casacore.tables as ct
import pyregion
from astropy.io import fits
from astropy.wcs import WCS
from glob import glob
import tables
from itertools import repeat
import re
import pandas as pd


def add_trailing_zeros(s, digitsize=4):
    """
     Repeat the zero character and add it to front

     :param s: string
     :param digitsize: number of digits (default 4)
     :return: trailing zeros + number --> example: 0001, 0021, ...
     """
    padded_string = "".join(repeat("0", digitsize)) + s
    return padded_string[-digitsize:]

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

def isfloat(num):
    """
    Check if value is a float
    """
    try:
        float(num)
        return True
    except ValueError:
        return False

def parse_history(ms, hist_item):
    """
    Grep specific history item from MS

    :param ms: measurement set
    :param hist_item: history item

    :return: parsed string
    """
    hist = os.popen('taql "SELECT * FROM '+ms+'::HISTORY" | grep '+hist_item).read().split(' ')
    for item in hist:
        if hist_item in item and len(hist_item)<=len(item):
            return item
    print('WARNING:' + hist_item + ' not found')
    return None


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
        if factor!=1:
            print("WARNING: " + ms + " time has been pre-averaged with factor "+str(factor)+". This might cause time smearing effects.")
        return factor
    elif isfloat(avg_num):
        factor = float(avg_num)
        print("WARNING: parsed factor in " + ms + " is not a digit but a float")
        return factor
    else:
        print("WARNING: parsed factor in " + ms + " is not a float or digit")
        return None

class SubtractWSClean:
    def __init__(self, mslist: list = None, region: str = None, localnorth: bool = True, onlyprint: bool = False):
        """
        Subtract image with WSClean

        :param mslist: measurement set list
        :param region: region file
        :param model_image: model image
        :param onlyprint: print only the commands (for testing purposes)
        """

        # list with MS
        self.mslist = mslist

        # wsclean model image
        model_images = glob('*-model*.fits')
        hdu = fits.open(model_images[0])
        self.imshape = (hdu[0].header['NAXIS1'], hdu[0].header['NAXIS2'])

        self.onlyprint = onlyprint

        self.model_images = glob("*-model*.fits")
        f = fits.open(self.model_images[0])
        history = str(f[0].header['HISTORY']).replace('\n', '').split()
        if '-taper-gaussian' in history:
            self.scale = history[history.index('-taper-gaussian') + 1]
        elif '-scale' in history:
            self.scale = history[history.index('-scale') + 1]
        else:
            self.scale = ''

        # region file to mask
        if localnorth:
            self.region = self.box_to_localnorth(region)
        else:
            self.region = pyregion.open(region)

    def clean_model_images(self):
        """
        Delete model images that do not match with MS
        """

        freqs = []
        for ms in self.mslist:
            t = ct.table(ms + "::SPECTRAL_WINDOW")
            freqs += list(t.getcol("CHAN_FREQ")[0])
            t.close()
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
        if len(glob('*-????-model-pb.fits')) >=1:
            self.model_images = glob('*-????-model-pb.fits')
        elif len(glob('*-????-model.fits')) >=1:
            self.model_images = glob('*-????-model.fits')
        else:
            self.model_images = glob('*-model-*.fits')

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
        (taken from sub-sources-outside-region.py)
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

            if not self.onlyprint:

                hdu = fits.open(fits_model)

                if region_cube:
                    manualmask = self.region.get_mask(hdu=hdu[0], shape=self.imshape)
                    for i in range(hdu[0].header['NAXIS4']):
                        hdu[0].data[i][0][np.where(manualmask == True)] = 0.0
                    hdu.writeto(fits_model, overwrite=True)

                else:
                    hduflat = self.flat_model_image(fits_model)
                    manualmask = self.region.get_mask(hdu=hduflat)
                    hdu[0].data[0][0][np.where(manualmask == True)] = 0.0
                    hdu.writeto(fits_model, overwrite=True)

                # hdu.close()

        return self


    def subtract_col(self, out_column: str = None):

        """
        Subtract column in Measurement Set
        :param out_column: out column name
        """

        for ms in self.mslist:
            print('Subtract ' + ms)
            ts = ct.table(ms, readonly=False)
            colnames = ts.colnames()

            if "MODEL_DATA" not in colnames:
                sys.exit(f"ERROR: MODEL_DATA does not exist in {ms}.\nThis is most likely due to a failed predict step.")

            if not self.onlyprint:
                if out_column not in colnames:
                    # get column description from DATA
                    desc = ts.getcoldesc('DATA')
                    # create output column
                    desc['name'] = out_column
                    # create template for output column
                    ts.addcols(desc)

                else:
                    print(out_column, ' already exists')

            # get number of rows
            nrows = ts.nrows()
            # make sure every slice has the same size
            best_slice = get_largest_divider(nrows, 1000)
            for c in range(0, nrows, best_slice):
                if 'CORRECTED_DATA' in colnames:
                    if c==0:
                        print('SUBTRACT --> CORRECTED_DATA - MODEL_DATA')
                    if not self.onlyprint:
                        data = ts.getcol('CORRECTED_DATA', startrow=c, nrow=best_slice)
                else:
                    if c==0:
                        print('SUBTRACT --> DATA - MODEL_DATA')
                    if not self.onlyprint:
                        data = ts.getcol('DATA', startrow=c, nrow=best_slice)
                if not self.onlyprint:
                    model = ts.getcol('MODEL_DATA', startrow=c, nrow=best_slice)
                    ts.putcol(out_column, data - model, startrow=c, nrow=best_slice)
            ts.close()

        return self

    def predict(self, h5parm: str = None, facet_regions: str = None):
        """
        Predict image

        :param h5parm: h5 solutions (optional)
        :param facet_regions: facet regions (if h5 solutions given)
        """

        f = fits.open(self.model_images[0])
        comparse = str(f[0].header['HISTORY']).replace('\n', '').split()
        command = ['wsclean', '-predict', f'-name {self.model_images[0].split("-")[0]}']

        for n, argument in enumerate(comparse):
            if argument in ['-gridder', '-padding', '-parallel-gridding',
                            '-idg-mode', '-beam-aterm-update', '-pol', '-scale']:
                command.append(' '.join(comparse[n:n + 2]))
            elif argument in ['-size']:
                command.append(' '.join(comparse[n:n + 3]))
            elif argument in ['-use-differential-lofar-beam', '-grid-with-beam',
                              '-use-idg', '-log-time', '-gap-channel-division',
                              '-apply-primary-beam']:
                command.append(argument)
            if argument == '-taper-gaussian':
                self.scale = comparse[n + 1]
            elif argument == '-scale' and '-taper-gaussian' not in comparse:
                self.scale = comparse[n + 1]

        if len(self.model_images)>1:
            command += ['-channels-out ' + str(len(self.model_images))]

        freqboundary = []
        for modim in sorted(self.model_images)[:-1]:
            fts = fits.open(modim)[0]
            fdelt, fcent = fts.header['CDELT3'] / 2, fts.header['CRVAL3']

            freqboundary.append(str(int(fcent + fdelt)))
            # fts.close()

        if len(freqboundary)>0:
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

        if not self.onlyprint:
            os.system(' '.join(command) + ' > log_predict.log')

        return self

    @staticmethod
    def isfulljones(h5: str = None):
        """
        Verify if file is fulljones

        :param h5: h5 file
        """
        T = tables.open_file(h5)
        soltab = list(T.root.sol000._v_groups.keys())[0]
        if 'pol' in T.root.sol000._f_get_child(soltab).val.attrs["AXES"].decode('utf8'):
            if T.root.sol000._f_get_child(soltab).pol[:].shape[0] == 4:
                T.close()
                return True
        T.close()
        return False

    def run_DP3(self, phaseshift: str = None, freqavg: str = None,
                timeavg: str = None, concat: bool = None,
                applybeam: bool = None, applycal_h5: str = None, dirname: str = None):
        """
        Run DP3 command

        :param phaseshift: do phase shift to specific center
        :param freqavg: frequency averaging
        :param timeavg: time averaging
        :param concat: concat the measurement sets
        :param applybeam: apply beam in phaseshifted phase center (or otherwise center of field)
        :param applycal_h5: applycal solution file
        """

        steps = []

        command = ['DP3',
                   'msin.missingdata=True',
                   'msin.datacolumn=SUBTRACT_DATA',
                   'msin.orderms=False',
                   'msout.storagemanager=dysco']

        # 1) PHASESHIFT
        if phaseshift is not None:
            phasecenter = phaseshift.replace('[', '').replace(']', '').split(',')
            phasecenter = f'[{phasecenter[0]},{phasecenter[1]}]'
            steps.append('ps')
            command += ['ps.type=phaseshifter',
                        'ps.phasecenter=' + phasecenter]

        # 2) APPLY BEAM
        if applybeam:
            steps.append('beam')
            command += ['beam.type=applybeam',
                        'beam.direction=[]',
                        'beam.updateweights=True']

        # 3) APPLYCAL
        if applycal_h5 is not None:
            # add fulljones solutions apply
            if self.isfulljones(applycal_h5):
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
                    command += [f'ac{ac_count}.type=applycal',
                                f'ac{ac_count}.parmdb={applycal_h5}',
                                f'ac{ac_count}.correction={corr}']
                    if phaseshift is not None and dirname is not None:
                        command += [f'ac{ac_count}.direction=' + dirname]
                    steps.append(f'ac{ac_count}')
                    ac_count += 1

        # 4) AVERAGING
        if freqavg is not None or timeavg is not None:
            steps.append('avg')
            command += ['avg.type=averager']
            if freqavg is not None:
                if str(freqavg).isdigit() or not str(freqavg)[-1].isalpha():
                    command += [f'avg.freqstep={int(freqavg)}']
                else:
                    command += [f'avg.freqresolution={freqavg}']
            if timeavg is not None:
                if str(timeavg).isdigit():
                    command += [f'avg.timestep={int(timeavg)}']
                else:
                    command += [f'avg.timeresolution={timeavg}']

        command += ['steps=' + str(steps).replace(" ", "").replace("\'", "")]

        if concat:
            command += [f'msin={",".join(self.mslist)}',
                        'msout=subtract_concat.ms']

            print('\n'.join(command))
            dp3_cmd = open("dp3.cmd", "w")
            dp3_cmd.write('\n'.join(command))
            dp3_cmd.close()

            if not self.onlyprint:
                os.system(' '.join(command) + " > dp3.subtract.log")
        else:
            for n, ms in enumerate(self.mslist):
                command+=[f'msin={ms}', f'msout=sub{self.scale}_{ms}']

                print('\n'.join(command))
                dp3_cmd = open("dp3.cmd", "w")
                dp3_cmd.write('\n'.join(command))
                dp3_cmd.close()

                if not self.onlyprint:
                    os.system(' '.join(command + [f'msin={ms}', f'msout=sub{self.scale}_{ms}']) + f" > dp3.sub{n}.log")

        return self


if __name__ == "__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser(description='Subtract region with WSClean')
    parser.add_argument('--mslist', nargs='+', help='measurement sets', required=True)
    parser.add_argument('--region', type=str, help='region file', required=True)
    parser.add_argument('--output_name', type=str, help='name of output files (default is model image name)')
    parser.add_argument('--model_image_folder', type=str,
                        help='folder where model images are stored (if not given script takes model images from run folder)')
    parser.add_argument('--no_local_north', action='store_true', help='do not move box to local north')
    parser.add_argument('--use_region_cube', action='store_true', help='use region cube')
    parser.add_argument('--h5parm_predict', type=str, help='h5 solution file')
    parser.add_argument('--facets_predict', type=str, help='facet region file with all facets to apply solutions')
    parser.add_argument('--phasecenter', type=str,
                        help='phaseshift to given point (example: --phaseshift 16h06m07.61855,55d21m35.4166)')
    parser.add_argument('--freqavg', type=str, help='frequency averaging')
    parser.add_argument('--timeavg', type=str, help='time averaging')
    parser.add_argument('--concat', action='store_true', help='concat MS')
    parser.add_argument('--applybeam', action='store_true', help='apply beam in phaseshift center or center of field')
    parser.add_argument('--applycal', action='store_true', help='applycal after subtraction and phaseshifting')
    parser.add_argument('--applycal_h5', type=str, help='applycal solution file')
    parser.add_argument('--print_only_commands', action='store_true', help='only print commands for testing purposes')
    parser.add_argument('--forwidefield', action='store_true',
                        help='will search for the polygon_info.csv file to extract information from')
    parser.add_argument('--skip_predict', action='store_true',
                        help='skip predict and do only subtract')
    args = parser.parse_args()

    if not args.skip_predict:

        # copy model images from model_image_folder
        if args.model_image_folder is not None:
            if len(glob(args.model_image_folder + '/*-????-model-pb.fits'))>1:
                os.system('cp ' + args.model_image_folder + '/*-????-model-pb.fits .')
            elif len(glob(args.model_image_folder + '/*-????-model.fits'))>1:
                os.system('cp ' + args.model_image_folder + '/*-????-model.fits .')
            elif len(glob(args.model_image_folder + '/*-model-pb.fits'))>1:
                os.system('cp ' + args.model_image_folder + '/*-model-pb.fits .')
            elif len(glob(args.model_image_folder + '/*-model.fits'))>1:
                os.system('cp ' + args.model_image_folder + '/*-model.fits .')
            else:
                sys.exit("ERROR: missing model images in folder "+args.model_image_folder)

        # remove MFS images if in folder
        if len(glob("*-????-model*.fits"))>1 and len(glob("*MFS-model*.fits"))>1:
            for mfs in glob("*MFS-model*.fits"):
                os.system('rm '+mfs)

        # rename model images
        if args.output_name is not None:
            model_images = glob('*-model*.fits')
            oldname = model_images[0].split("-")[0]
            for model in model_images:
                os.system('mv ' + model + ' ' + model.replace(oldname, args.output_name))

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

        t = ct.table(args.mslist[0] + "::SPECTRAL_WINDOW")
        channum = len(t.getcol("CHAN_FREQ")[0])
        t.close()

        polygon = polygon_info.loc[polygon_info.polygon_file == args.region.split('/')[-1]]
        try:
            phasecenter = polygon['poly_center'].values[0]
        except AttributeError:
            print('WARNING: no poly center in polygon_info.csv, use dir instead.')
            phasecenter = polygon['dir'].values[0]
        except KeyError:
            print('WARNING: no poly center in polygon_info.csv, use dir instead.')
            phasecenter = polygon['dir'].values[0]

        # take only averaging factors that are channum%avg==0
        avg = get_largest_divider(channum, int(polygon['avg'].values[0]))

        freqavg = int(avg)
        try:
            # if there is pre averaging done on the ms, we need to take this into account
            timeavg = int(freqavg/get_time_preavg_factor(args.mslist[0]))
        except:
            timeavg = int(freqavg)
        dirname = polygon['dir_name'].values[0]

    else:
        phasecenter = args.phasecenter
        freqavg = args.freqavg
        timeavg = args.timeavg
        dirname = None

    object = SubtractWSClean(mslist=args.mslist, region=args.region, localnorth=not args.no_local_north,
                             onlyprint=args.print_only_commands)

    if not args.skip_predict:

        # clean model images
        object.clean_model_images()

        # mask
        print('############## MASK REGION ##############')
        object.mask_region(region_cube=args.use_region_cube)

        # predict
        print('############## PREDICT ##############')
        object.predict(h5parm=args.h5parm_predict, facet_regions=args.facets_predict)

    # subtract
    print('############## SUBTRACT ##############')
    object.subtract_col(out_column='SUBTRACT_DATA')

    # extra DP3 step
    if args.phasecenter is not None or \
            args.freqavg is not None or \
            args.timeavg is not None or \
            args.concat is not None or \
            args.applybeam is not None or \
            args.applycal is not None:
        print('############## RUN DP3 ##############')
        if args.applycal_h5 is not None:
            applycalh5 = args.applycal_h5
        elif args.applycal and args.applycal_h5 is None and args.h5parm_predict is not None:
            applycalh5 = args.h5parm_predict
        elif args.applycal and not args.applycal_h5:
            sys.exit("ERROR: need a solution file for applycal (give with --applycal_h5)")
        else:
            applycalh5 = None

        object.run_DP3(phaseshift=phasecenter, freqavg=freqavg, timeavg=timeavg,
                       concat=args.concat, applybeam=args.applybeam, applycal_h5=applycalh5, dirname=dirname)
