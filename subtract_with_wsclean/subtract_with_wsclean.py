import numpy as np
import sys
import pyrap.tables as pt
import os
import casacore.tables as pt
import pyregion
from astropy.io import fits
from astropy.wcs import WCS
from glob import glob

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
        self.model_images = glob('*-model*.fits')
        hdu = fits.open(self.model_images[0])
        self.imshape = (hdu[0].header['NAXIS1'], hdu[0].header['NAXIS2'])

        # region file to mask
        if localnorth:
            self.region = self.box_to_localnorth(region)
        else:
            self.region = pyregion.open(region)

        self.onlyprint = onlyprint

        self.scale=''

    def box_to_localnorth(self, region):
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
    def flat_model_image(fitsfile):
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

    def mask_region(self, region_cube:bool = False):
        """
        :param region_cube: if region_cube make cube, otherwise 2D (flatten)
        """

        for fits_model in self.model_images:

            print('Mask '+fits_model)

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

                hdu.close()

        return self


    def subtract_col(self, out_column):

        """
        Subtract column in Measurement Set
        :param out_column: out column name
        """

        for ms in self.mslist:
            print('Subtract '+ms)
            if not self.onlyprint:
                ts = pt.table(ms, readonly=False)
                colnames = ts.colnames()
                if out_column not in colnames:
                    desc = ts.getcoldesc('DATA')
                    desc['name'] = out_column
                    ts.addcols(desc)
                    ts.close()  # to write results

                else:
                    print(out_column, ' already exists')
                    ts.close()

        for ms in self.mslist:
            ts = pt.table(ms, readonly=False)
            colnames = ts.colnames()
            if 'CORRECTED_DATA' in colnames:
                print('SUBTRACT --> CORRECTED_DATA - MODEL_DATA')
                if not self.onlyprint:
                    data = ts.getcol('CORRECTED_DATA')
            else:
                print('SUBTRACT --> DATA - MODEL_DATA')
                if not self.onlyprint:
                    data = ts.getcol('DATA')
            if not self.onlyprint:
                model = ts.getcol('MODEL_DATA')
                ts.putcol(out_column, data - model)
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
            if argument in ['-channels-out', '-gridder',
                            '-padding', '-parallel-gridding',
                            '-idg-mode', '-beam-aterm-update', '-pol', '-minuv-l']:
                command.append(' '.join(comparse[n:n+2]))
            elif argument in ['-size']:
                command.append(' '.join(comparse[n:n+3]))
            elif argument in ['-use-differential-lofar-beam', '-grid-with-beam',
                              '-use-idg', '-log-time', '-gap-channel-division',
                              '-apply-primary-beam']:
                command.append(argument)
            if argument=='-taper-gaussian':
                self.scale=comparse[n+1]
            elif argument=='-scale' and '-taper-gaussian' not in comparse:
                self.scale=comparse[n+1]

        if h5parm is not None:
            command+=[f'-apply-facet-solutions {h5parm} amplitude000,phase000',
                      f' -facet-regions {facet_regions}', '-apply-facet-beam',
                      f'-facet-beam-update {comparse[comparse.index("-facet-beam-update")+1]}']

        command += [' '.join(self.mslist)]

        #run
        print('\n'.join(command))
        if not self.onlyprint:
            os.system(' '.join(command)+' > log_predict.txt')

        return self


    def run_DP3(self, phaseshift: str = None, freqavg: str = None,
                timeavg: str = None, concat: bool = None,
                applybeam: bool = None, applycal_h5: str = None):
        """
        Run DP3 command

        :param phaseshift: do phase shift to specific center
        :param freqavg: frequency averaging
        :param timeavg: time averaging
        :param concat: concat the measurement sets
        :param applybeam: apply beam in phaseshifted phase center (or otherwise center of field)
        :param applycal_h5: applycal solution file
        """

        steps=[]

        command = ['DP3',
                   'msin.missingdata=True',
                   'msin.datacolumn=SUBTRACT_DATA',
                   'msin.orderms=False',
                   'msout.storagemanager=dysco',
                   'msout.writefullreslag=False']

        #1) PHASESHIFT
        if phaseshift is not None:
            phasecenter = phaseshift.replace('[','').replace(']','').split(',')
            phasecenter = f'[{phasecenter[0]},{phasecenter[1]}]'
            steps.append('ps')
            command += ['ps.type=phaseshifter',
                        'ps.phasecenter='+phasecenter]

        #2) APPLY BEAM
        if applybeam:
            steps.append('beam')
            command += ['beam.type=applybeam',
                        'beam.direction=[]',
                        'beam.updateweights=True']

        #3) APPLYCAL
        if applycal_h5 is not None:
            steps.append('ac')
            command += ['ac.type=applycal',
                        'ac.parmdb='+applycal_h5,
                        'ac.correction=fulljones',
                        'ac.soltab=[amplitude000,phase000]']
            if phaseshift is not None:
                command += ['ac.direction='+phasecenter]

        #4) AVERAGING
        if freqavg is not None or timeavg is not None:
            steps.append('avg')
            command += ['avg.type=averager']
            if freqavg is not None:
                if freqavg.isdigit():
                    command+=[f'avg.freqstep={freqavg}']
                else:
                    command += [f'avg.freqresolution={freqavg}']
            if timeavg is not None:
                if timeavg.isdigit():
                    command+=[f'avg.timestep={timeavg}']
                else:
                    command+=[f'avg.timeresolution={timeavg}']

        command+=['steps='+str(steps).replace(" ","").replace("\'","")]

        if concat:
            command+=[f'msin={",".join(self.mslist)}',
                      'msout=subtract_concat.ms']
            print('\n'.join(command))
            if not self.onlyprint:
                os.system(' '.join(command))
        else:
            for ms in self.mslist:
                print('\n'.join(command+[f'msin={ms}', f'msout=sub{self.scale}_{ms}']))
                if not self.onlyprint:
                    os.system(' '.join(command+[f'msin={ms}', f'msout=sub{self.scale}_{ms}']))

        return self


if __name__ == "__main__":

    from argparse import ArgumentParser
    parser = ArgumentParser(description='Subtract region with WSClean')
    parser.add_argument('--mslist', nargs='+', help='measurement sets', required=True)
    parser.add_argument('--region', type=str, help='region file', required=True)
    parser.add_argument('--output_name', type=str, help='name of output files (default is model image name)')
    parser.add_argument('--model_image_folder', type=str, help='folder where model images are stored (if not given script takes model images from run folder)')
    parser.add_argument('--no_local_north', action='store_true', help='do not move box to local north')
    parser.add_argument('--use_region_cube', action='store_true', help='use region cube')
    parser.add_argument('--h5parm_predict', type=str, help='h5 solution file')
    parser.add_argument('--facets_predict', type=str, help='facet region file with all facets to apply solutions')
    parser.add_argument('--phasecenter', type=str, help='phaseshift to given point (example: --phaseshift 16h06m07.61855,55d21m35.4166)')
    parser.add_argument('--freqavg', type=str, help='frequency averaging')
    parser.add_argument('--timeavg', type=str, help='time averaging')
    parser.add_argument('--concat', action='store_true', help='concat MS')
    parser.add_argument('--applybeam', action='store_true', help='apply beam in phaseshift center or center of field')
    parser.add_argument('--applycal', action='store_true', help='applycal after subtraction and phaseshifting')
    parser.add_argument('--applycal_h5', type=str, help='applycal solution file')
    parser.add_argument('--print_only_commands', action='store_true', help='only print commands for testing purposes')
    parser.add_argument('--forwidefield', action='store_true', help='will search for the polygon_info.csv file')
    args = parser.parse_args()

    if args.model_image_folder is not None:
        os.system('cp '+args.model_image_folder+'/*-model*.fits .')

    if len(glob('*-model*.fits'))==0:
        sys.exit("ERROR: missing model images in folder.\nPlease copy model images to run folder or give --model_image_folder.")

    if args.output_name is not None:
        model_images = glob('*-model*.fits')
        oldname = model_images[0].split("-")[0]
        for model in model_images:
            os.system('mv '+model+' '+model.replace(oldname, args.output_name))

    #--forwidefield --> will read averaging and phasecenter from polygon_info.csv
    if args.forwidefield:
        import pandas as pd
        if os.path.isfile('polygon_info.csv'):
            polygon_info = pd.read_csv('polygon_info.csv')
        elif os.path.isfile('../polygon_info.csv'):
            polygon_info = pd.read_csv('../polygon_info.csv')
        else:
            sys.exit('ERROR: using --forwidefield option needs polygon_info.csv file to read polygon information from')

        polygon = polygon_info[polygon_info.polygon_file==args.region].reset_index().to_dict()['dir']
        phasecenter = polygon['dir'][0]
        freqavg = polygon['avg'][0]
        timeavg = polygon['avg'][0]
    else:
        phasecenter = args.phasecenter
        freqavg = args.freqavg
        timeavg = args.timeavg


    object = SubtractWSClean(mslist=args.mslist, region=args.region, localnorth=not args.no_local_north,
                             onlyprint=args.print_only_commands)

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
            applycalh5=None

        object.run_DP3(phaseshift=phasecenter, freqavg=freqavg, timeavg=timeavg,
                       concat=args.concat, applybeam=args.applybeam, applycal_h5=applycalh5)
        