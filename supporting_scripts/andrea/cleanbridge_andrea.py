#!/usr/bin/env python

import matplotlib

matplotlib.use('Agg')
import os, sys
import numpy as np
import pyrap.tables as pt
import os.path
import pyregion
import argparse
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)


def cleanup(keepunmasked=False):

    os.system('rm -rf *-00*-*.fits')  # remove channel maps
    os.system('rm -rf *-dirty.fits')  # remove all dirty images

    if keepunmasked == False:
        os.system('rm -rf *_compact-*.fits')  # remove non-masked images
        os.system('rm -rf *_subROBUST*.fits')  # remove non-masked images
        os.system('rm -rf *_ROBUST*fits')  # remove non-masked images

    return


def getimsize(boxfile, cellsize=1.5):
    """
    find imsize need to image a DS9 boxfile region
    """
    r = pyregion.open(boxfile)

    xs = np.ceil((r[0].coord_list[2]) * 1.6 * 3600. / cellsize)
    ys = np.ceil((r[0].coord_list[3]) * 1.6 * 3600. / cellsize)

    imsize = np.ceil(xs)  # // Round up decimals to an integer
    if (imsize % 2 == 1):
        imsize = imsize + 1
    return np.int(imsize)


def compute_uvmin(redshift, sourceLLS=1.0):
    '''
    sourceLLS in units of Mpc#
    taper output for WSClean in arcsec
    '''
    oneradinmpc = cosmo.angular_diameter_distance(redshift) / (360. / (2. * np.pi))
    scalebarlengthdeg = sourceLLS / oneradinmpc.value

    return 1. / (scalebarlengthdeg * np.pi / 180.)


def compute_taper(redshift, taperscale):
    '''
    taperscale in units of kpc#
    '''
    oneradinmpc = cosmo.angular_diameter_distance(redshift) / (360. / (2. * np.pi))
    taper = 1e-3 * taperscale / (oneradinmpc.value)

    return taper * 3600


def adjustniter_for_taper(taper, niter):
    if taper < 5:
        return np.int(niter)
    if taper >= 5 and taper < 15:
        return np.int(niter / 2)
    if taper >= 15:
        return np.int(niter / 4)


def makeimage(mslist, imageout, pixsize, imsize, channelsout=6, niter=15000, robust=-0.5, minuv=80, mgain=0.8,
              uvtaper=None, predict=True, deepmultiscale=False, column=None, instrument='lofarHBA', h5=None, facet=None):

    # some setup
    automaskval = 2.5

    if uvtaper != None:
        if uvtaper < 0:
            print('Not supported uvtaper', uvtaper)
        else:
            print(imsize, pixsize, uvtaper)
            imsizein = np.int(np.float(imsize) * (pixsize / (uvtaper / 5.)))
            pixsizein = np.int(uvtaper / 5.)
            if np.float(pixsizein) < pixsize:  # to deal with rounding issues which cause a 0arcsec pixelsize
                pixsizein = pixsize
                imsizein = imsize
    else:
        imsizein = imsize
        pixsizein = pixsize

    if imsizein < 511:  # otherwise images get too small for multiscales
        imsizein = 512

    # few simple checks to make sure we have useful data
    msliststring = ' '.join(map(str, mslist))
    os.system('rm -f ' + imageout + '-*.fits')
    imcol = 'CORRECTED_DATA'
    t = pt.table(mslist[0], readonly=True)  # just test for first ms in mslist
    colnames = t.colnames()
    if 'CORRECTED_DATA' not in colnames:  # check which column to image
        imcol = 'DATA'
    t.close()

    if column != None:
        imcol = column

    # build wsclean command
    cmd = 'wsclean -apply-facet-solutions '+h5+' amplitude000,phase000 -facet-regions ' + facet + ' '
    cmd += '-no-update-model-required -minuv-l ' + str(minuv) + ' '
    cmd += '-size ' + str(imsizein) + ' ' + str(imsizein) + ' -reorder '
    cmd += '-weight briggs ' + str(robust) + ' -weighting-rank-filter 3 -clean-border 1 '
    cmd += '-mgain ' + str(mgain) + ' -fit-beam -data-column ' + imcol + ' -padding 1.4 '
    cmd += '-multiscale -multiscale-max-scales 10 '
    cmd += '-join-channels '

    if channelsout != 0:
        cmd += '-channels-out ' + str(channelsout) + ' '

    # cmd += '-parallel-deconvolution ' + str(np.int(imsizein)/2) + ' '


    cmd += '-auto-mask ' + str(automaskval) + ' -auto-threshold 0.5 '

    if uvtaper != None:
        cmd += '-taper-gaussian ' + str(uvtaper) + 'arcsec '

    if instrument == 'lofarHBA':

        baselineav = 2.5e3 * 60000. * 2. * np.pi * np.float(pixsizein) / (24. * 60. * 60 * np.float(imsizein))

        # limit baseline averaging to 10, fixes prob
        if baselineav > 10.0:
            baselineav = 10.

        baselineav = str(baselineav)
        cmd += '-pol i '
        cmd += '-fit-spectral-pol 3 '  # -beam-shape 6arcsec 6arcsec 0deg '
        cmd += '-baseline-averaging ' + baselineav + ' '
    elif (instrument == 'ugmrt3') or (instrument == 'ugmrt4') or (instrument == 'gmrt325') or (instrument == 'gmrt610'):
        cmd += '-pol RR '
    elif instrument == 'jvlaC':
        cmd += '-pol i '

    cmd += '-name ' + imageout + ' -scale ' + str(pixsizein) + 'arcsec '

    print('WSCLEAN: ', cmd + '-niter ' + str(niter) + ' ' + msliststring)
    # logging.info(cmd + '-niter ' + str(niter) + ' ' + msliststring)
    os.system(cmd + '-niter ' + str(niter) + ' ' + msliststring)

    if deepmultiscale:
        # predict first to fill MODEL_DATA so we can continue with clean
        cmdp = 'wsclean -apply-facet-solutions '+h5+' amplitude000,phase000 -facet-regions ' + facet + ' -size '
        cmdp += str(imsizein) + ' ' + str(imsizein) + ' -channels-out ' + str(channelsout) + ' -padding 1.4 -predict '

        cmdp += '-name ' + imageout + ' -scale ' + str(pixsizein) + 'arcsec ' + msliststring
        print('PREDICT STEP for continue: ', cmdp)
        os.system(cmdp)

        # NOW continue cleaning
        cmd += '-niter ' + str(niter / 15) + ' -multiscale -continue ' + msliststring
        print('WSCLEAN continue: ', cmd)
        os.system(cmd)

    if predict:
        cmd = 'wsclean -size '
        cmd += str(imsizein) + ' ' + str(imsizein) + ' '

        if channelsout != 0:
            cmd += '-channels-out ' + str(channelsout) + ' '

        cmd += '-padding 1.4 -predict '

        if (instrument == 'ugmrt3') or (instrument == 'ugmrt4'):
            cmd += '-pol RR '

        cmd += '-name ' + imageout + ' -scale ' + str(pixsizein) + 'arcsec ' + msliststring
        print('PREDICT STEP: ', cmd)
        os.system(cmd)


def subtractcompact(mslist, imageout, pixsize, imsize, minuv, channelsout=6, niter=15000, robust=-0.5,
                    outcolumn='DIFFUSE_SUB', instrument='lofarHBA', predictcompact=None, h5=None, facet=None):
    if predictcompact == 'skip':  # it assumes that the DIFFUSE_SUB column is already present
        return

    if predictcompact != 'start':  # if predictcompact == 'start' skips this block

        # some setup

        makeimage(mslist, imageout + '_compact', pixsize, imsize, channelsout=channelsout, niter=niter, robust=robust,
                  minuv=minuv, predict=False, instrument=instrument, h5=h5, facet=facet)

        # re-image with mask
        if predictcompact == 'stop':
            makeimage(mslist, imageout + '_compactmask', pixsize, imsize, channelsout=channelsout, niter=niter,
                      robust=robust, minuv=minuv,
                      predict=False, deepmultiscale=False,
                      instrument=instrument, h5=h5, facet=facet)  # note that here the predict is False, it will be done later in the if predictcompact == 'start' block
            print('Run stopped before the predict of the compact image. Now edit the model and then re-run with --predictcompact=start')
            return
        else:
            makeimage(mslist, imageout + '_compactmask', pixsize, imsize, channelsout=channelsout, niter=niter,
                      robust=robust, minuv=minuv,
                      predict=True, deepmultiscale=False,
                      instrument=instrument, h5=h5, facet=facet)  # note that here the predict is True

    elif predictcompact == 'start':  # now do the predict

        cmd = 'wsclean -size '
        cmd += str(imsize) + ' ' + str(imsize) + ' '

        if channelsout != 0:
            cmd += '-channels-out ' + str(channelsout) + ' '

        cmd += '-padding 1.4 -predict '

        if (instrument == 'ugmrt3') or (instrument == 'ugmrt4'):
            cmd += '-pol RR '

        msliststring = ' '.join(map(str, mslist))

        cmd += '-name ' + imageout + '_compactmask' + ' -scale ' + str(pixsize) + 'arcsec ' + msliststring
        print('PREDICT STEP: ', cmd)
        os.system(cmd)

    # now subtract the columns

    for ms in mslist:
        ts = pt.table(ms, readonly=False)
        colnames = ts.colnames()
        if outcolumn not in colnames:
            desc = ts.getcoldesc('DATA')
            desc['name'] = outcolumn
            ts.addcols(desc)
            ts.close()  # to write results

        else:
            print(outcolumn, ' already exists')
            ts.close()

    for ms in mslist:
        ts = pt.table(ms, readonly=False)
        colnames = ts.colnames()
        if 'CORRECTED_DATA' in colnames:
            data = ts.getcol('CORRECTED_DATA')
        else:
            data = ts.getcol('DATA')
        model = ts.getcol('MODEL_DATA')
        ts.putcol(outcolumn, data - model)
        ts.close()
    return


parser = argparse.ArgumentParser(
    description='Make images from extraction run. Requires working version of the DR2-pipeline software and WSClean (Oct 2018 or newer')
parser.add_argument('--boxfile', help='optional boxfile to set imsize automatically', type=str)
parser.add_argument('--imsize', help='image size, you can take it from selfcal.log', type=int, default=1500)
parser.add_argument('--niter', help='niter, default=25000', default=25000, type=int)
parser.add_argument('--channelsout', help='channelsout (set it to 0 to switch off joint deconvolution), default=6',
                    default=6, type=int)
parser.add_argument('--minuv', help='inner uv-cut for image in lambda, default=80', default=80., type=float)
parser.add_argument('--mgain', help='major iteration gain, default=0.8', default=0.8, type=float)
parser.add_argument('--tapers_arcsec', help='space-separated values of tapers in arcsec units (set it to -1 to skip arcsec-tapered images), default=10 15 30 60', default=[10,15,30,60], nargs='+', type=int)
parser.add_argument('--skipstandardimg', help='skip the standard robust -0.5 and the robust -1.25 images',
                    default=False, action='store_true')
parser.add_argument('--pixelscale', help='pixels size in arcsec, default=1.5', default=1.5, type=float)
parser.add_argument('--sourceLLS', help='size in Mpc of diffuse emission for uvcut, default=0.4', default=0.4,
                    type=float)
parser.add_argument('--minuvforsub',
                    help='if provided (in lambda) it overwrites the uvcut corresponding to sourceLLS, default=None',
                    default=None, type=float)
parser.add_argument('--predictcompact',
                    help='use "stop" to stop the run before the predict of compact image to allow manual editing of the model, then re-run the same command using "start". Use "skip" if you already have done the prediction',
                    choices=['stop', 'start', 'skip'])
parser.add_argument('--instrument',
                    help='adjust some parameters for [lofarHBA|ugmrt3|ugmrt4|gmrt325|gmrt610|jvlaC], default=lofarHBA',
                    default='lofarHBA', type=str)
parser.add_argument('--z', help='redshift of cluster, not required if --nodosub is used', default=-1., type=float)
parser.add_argument('-i', '--imagename', help='imagename, default=image', required=True, type=str)
parser.add_argument('--maskthreshold', help='threshold for MakeMask.py, default=3.0', default=3., type=float)
parser.add_argument('--keepunmasked', help='keep unmasked images and masks (for debugging purposes)', default=False,
                    action='store_true')
parser.add_argument('--h5', help='h5 image', required=True, type=str)
parser.add_argument('--facetregions', help='facet regions', required=True, type=str)
parser.add_argument('--test', action='store_true', help='use this for tests to cut time')
parser.add_argument('--ms', nargs='*', help='msfile(s)')

args = vars(parser.parse_args())

makemask = '/net/para10/data1/shimwell/software/killmsddf/new-install/DDFacet/SkyModel/MakeMask.py'

if args['test']:
    import pyrap.tables as pt
    for n, ms in enumerate(args['ms']):
        # pt.taql(
        #     'SELECT FROM {MSIN} WHERE TIME IN (SELECT DISTINCT TIME FROM {MSIN} OFFSET {time[0]} LIMIT {time[1]}) GIVING {MSOUT} AS PLAIN'.format(
        #         MSIN=ms, MSOUT=ms+'.test', time=[0, 300]))
        args['ms'][n] += '.test'

minuv = args['minuv']
pixsize = args['pixelscale']
niter = args['niter']
mslist = sorted(args['ms'])
imageout = args['imagename']

if args['instrument'] == 'lofarHBA':
    pass
elif args['instrument'] == 'ugmrt3':
    args['imsize'] = 6000
    pixsize = 1.25
elif args['instrument'] == 'ugmrt4':
    args['imsize'] = 3000
    pixsize = 1.0
elif args['instrument'] == 'gmrt325':
    args['channelsout'] = 0
    args['imsize'] = 6000
    pixsize = 1.25
elif args['instrument'] == 'gmrt610':
    args['channelsout'] = 0
    args['imsize'] = 3000
    pixsize = 1.0
elif args['instrument'] == 'jvlaC':
    print('Instrument not implemented yet')  # for the vla robust 0.0 is better, set it as default?
    sys.exit()
else:
    print('Instrument not available.')
    sys.exit()

if args['boxfile'] == None and args['imsize'] == None:
    print('Incomplete input detected, either boxfile or imsize is required')
    sys.exit()
if args['boxfile'] != None and args['imsize'] != None:
    print('Wrong input detected, both boxfile and imsize are set')
    sys.exit()

if args['boxfile'] != None:
    imsize = str(getimsize(args['boxfile'], args['pixelscale']))
if args['imsize'] != None:
    imsize = str(args['imsize'])


if args['minuvforsub'] == None:
    if args['z'] < 0:
        print('You need provide a redshift, none was given')
        sys.exit()
    else:
        minuv_forsub = compute_uvmin(args['z'], sourceLLS=args['sourceLLS'])
else:
    minuv_forsub = args['minuvforsub']


subtractcompact(mslist, imageout, pixsize, imsize, minuv_forsub, channelsout=args['channelsout'],
                niter=np.int(niter / 1.25), robust=-0.5, outcolumn='DIFFUSE_SUB', instrument=args['instrument'],
                predictcompact='stop', h5=args['h5'], facet=args['facetregions'])

subtractcompact(mslist, imageout, pixsize, imsize, minuv_forsub, channelsout=args['channelsout'],
                niter=np.int(niter / 1.25), robust=-0.5, outcolumn='DIFFUSE_SUB', instrument=args['instrument'],
                predictcompact='start', h5=args['h5'], facet=args['facetregions'])


if args['skipstandardimg'] == False:

    #  -----------------------------------------------------------------------------
    #  --- make the standard image robust -0.5 image, compact source subtracted ----
    #  -----------------------------------------------------------------------------
    makeimage(mslist, imageout + '_subROBUST-0.5', pixsize, imsize, channelsout=args['channelsout'],
              niter=np.int(niter / 1.5), robust=-0.5, minuv=minuv, mgain=0.8, column='DIFFUSE_SUB', predict=False,
              instrument=args['instrument'], h5=args['h5'], facet=args['facetregions'])


    # re-image with mask
    makeimage(mslist, imageout + '_masksubROBUST-0.5', pixsize, imsize, channelsout=args['channelsout'],
              niter=np.int(niter / 1.5), robust=-0.5, minuv=minuv, mgain=0.8, predict=False, column='DIFFUSE_SUB',
              deepmultiscale=False, instrument=args['instrument'], h5=args['h5'], facet=args['facetregions'])

for t in args['tapers_arcsec']:

    if t == -1:  # to skip the imaging of arcsec-tapered images
        break

    newniter = adjustniter_for_taper(float(t), niter)

    # make the tapered image without mask, compact source subtracted
    makeimage(mslist, imageout + '_subROBUST-0.5TAPER' + str(t), pixsize, imsize, channelsout=args['channelsout'],
              niter=newniter, robust=-0.5, minuv=minuv, mgain=args['mgain'], predict=False, column='DIFFUSE_SUB',
              uvtaper=float(t), instrument=args['instrument'], h5=args['h5'], facet=args['facetregions'])


    # re-image with mask
    makeimage(mslist, imageout + '_masksubROBUST-0.5TAPER' + str(t), pixsize, imsize, channelsout=args['channelsout'],
              niter=newniter, robust=-0.5, minuv=minuv, mgain=args['mgain'], predict=False, deepmultiscale=False,
              column='DIFFUSE_SUB', uvtaper=float(t), instrument=args['instrument'], h5=args['h5'], facet=args['facetregions'])

if args['skipstandardimg'] == False:
    #  --------------------------------------------------
    #  --- make the standard image robust -0.5 image ----
    #  --------------------------------------------------
    makeimage(mslist, imageout + '_ROBUST-0.5', pixsize, imsize, channelsout=args['channelsout'], niter=niter,
              robust=-0.5, minuv=minuv, mgain=0.8, predict=False, instrument=args['instrument'], h5=args['h5'], facet=args['facetregions'])


    # re-image with mask
    makeimage(mslist, imageout + '_maskROBUST-0.5', pixsize, imsize, channelsout=args['channelsout'], niter=niter,
              robust=-0.5, minuv=minuv, mgain=0.8,
              predict=False, deepmultiscale=False, instrument=args['instrument'], h5=args['h5'], facet=args['facetregions'])

    #  --------------------------------------------------
    #  --- make the high-res image robust -1.25 image ----
    #  --------------------------------------------------
    makeimage(mslist, imageout + '_maskROBUST-1.25', pixsize, imsize, channelsout=args['channelsout'], niter=niter,
              robust=-1.25, minuv=minuv, mgain=0.8,
              predict=False, deepmultiscale=False, instrument=args['instrument'], h5=args['h5'], facet=args['facetregions'])

for t in args['tapers_arcsec']:

    if t == -1:  # to skip the imaging of arcsec-tapered images
        break

    if t >= 60:  # to skip images with taper > 60arcsec and with point sources
        continue

    newniter = adjustniter_for_taper(float(t), niter)

    # make the tapered image without mask
    makeimage(mslist, imageout + '_ROBUST-0.5TAPER' + str(t), pixsize, imsize, channelsout=args['channelsout'],
              niter=newniter, robust=-0.5, minuv=minuv, mgain=args['mgain'], predict=False, uvtaper=float(t),
              instrument=args['instrument'], h5=args['h5'], facet=args['facetregions'])

    # re-image with mask
    makeimage(mslist, imageout + '_maskROBUST-0.5TAPER' + str(t), pixsize, imsize, channelsout=args['channelsout'],
              niter=newniter, robust=-0.5, minuv=minuv, mgain=args['mgain'], predict=False, deepmultiscale=False,
              uvtaper=float(t), instrument=args['instrument'], h5=args['h5'], facet=args['facetregions'])

cleanup(keepunmasked=args['keepunmasked'])