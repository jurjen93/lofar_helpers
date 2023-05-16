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
    def __init__(self, mslist: list = None, region: str = None):
        """
        Subtract image with WSClean

        :param mslist: measurement set list
        :param region: region file
        :param model_image: model image
        """

        # list with MS
        self.mslist = mslist

        # region file to mask
        self.region = pyregion.open(region)

        # wsclean model image
        self.model_images = glob('*-model.fits')
        hdu = fits.open(self.model_images[0])
        self.imshape = (hdu[0].header['NAXIS1'], hdu[0].header['NAXIS2'])

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

    def mask_region(self, region_cube:bool = False, mask_fits:str = 'mask.fits'):
        """
        :param region_cube: if region_cube make cube, otherwise 2D (flatten)
        :param mask_fits: name of fits file
        """

        for fits_model in self.model_images:

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

        return self

    def subtract_col(self, out_column):

        """
        Subtract column in Measurement Set
        :param out_column: out column name
        """

        for ms in self.mslist:
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
                data = ts.getcol('CORRECTED_DATA')
            else:
                data = ts.getcol('DATA')
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
        command = ['wsclean',
                   '-gridder wgridder',
                   f'-size {self.imshape[0]} {self.imshape[1]}',
                   '-channels-out 6',
                   '-padding 1.2',
                   '-predict',
                   '-parallel-gridding 5',
                    f'-name {self.model_images[0].split("-")[0]}']

        if h5parm is not None:
            command+=[f'-apply-facet-solutions {h5parm} amplitude000,phase000',
                      f' -facet-regions {facet_regions}']

        command += [' '.join(self.mslist)]

        #run
        print('\n'.join(command))
        os.system(' '.join(command)+' > log_predict.txt')

        return self

    def run_DP3(self, phaseshift: bool = False, freqavg: str = None, timeavg: str = None, concat: bool = False):
        """
        Run DP3 command

        :param phaseshift: do phase shift to specific center
        :param freqavg: frequency averaging
        :param timeavg: time averaging
        :param concat: concat the measurement sets
        """

        steps=[]

        command = ['DP3',
                   'msin.missingdata=True',
                   'msin.datacolumn=SUBTRACT_DATA',
                   'msin.orderms=False',
                   'msout.storagemanager=dysco',
                   'msout.writefullreslag=False']

        if phaseshift is not None:
            steps.append('ps')
            command += ['ps.type=phaseshifter', f'ps.phasecenter={phaseshift}']

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

        command+=[f'steps={steps}']

        if concat:
            command+=[f'msin={",".join(self.mslist)}',
                      'msout=subtract_concat.ms']
            print('\n'.join(command))
            os.system(' '.join(command))
        else:
            for ms in self.mslist:
                command+=[f'msin={ms}', f'msout=sub_{ms}']
                print('\n'.join(command))
                os.system(' '.join(command))

        return self

if __name__ == "__main__":

    from argparse import ArgumentParser
    parser = ArgumentParser(description='Subtract region with WSClean')
    parser.add_argument('--mslist', nargs='+', help='measurement sets', required=True)
    parser.add_argument('--region', type=str, help='region file', required=True)
    parser.add_argument('--use_region_cube', action='store_true', help='use region cube')
    parser.add_argument('--mask_fits', type=str, help='name for mask_fits', default='mask.fits')
    parser.add_argument('--h5parm', type=str, help='h5 solution file', default=None)
    parser.add_argument('--facet_region', type=str, help='facet region file', default=None)
    parser.add_argument('--phaseshift', type=str, help='phaseshift to given point (example: --phaseshift 16h06m07.61855,55d21m35.4166)', default=None)
    parser.add_argument('--freqavg', type=str, help='frequency averaging', default=None)
    parser.add_argument('--timeavg', type=str, help='time averaging', default=None)
    parser.add_argument('--concat', action='store_true', help='concat MS', default=None)
    args = parser.parse_args()

    if len(glob('*-model.fits'))==0:
        sys.exit('ERROR: missing *-model.fits images from WSClean.\nCopy these images to run folder.')

    object = SubtractWSClean(mslist=args.mslist, region=args.region)

    # mask
    print('############## MASK REGION ##############')
    object.mask_region(region_cube=args.use_region_cube, mask_fits=args.mask_fits)

    # predict
    print('############## PREDICT ##############')
    object.predict(h5parm=args.h5parm, facet_regions=args.facet_region)

    # subtract
    print('############## SUBTRACT ##############')
    object.subtract_col(out_column='SUBTRACT_DATA')

    # extra DP3 step
    if args.phaseshift is not None or \
        args.freqavg is not None or \
        args.timeavg is not None or \
        args.concat is not None:
        print('############## RUN DP3 ##############')
        object.run_DP3(phaseshift = args.phaseshift, freqavg = args.freqavg, timeavg = args.timeavg, concat = args.concat)