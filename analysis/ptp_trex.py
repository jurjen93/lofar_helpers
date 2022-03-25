#!/usr/bin/python
# Filename:PT-REX.py

###########################
"""
Author: Alessandro Ignesti
Ver 2.0
Contacts: https://aignesti.github.io/

Point-to-point TRend EXtractor [PT-REX]

##########################
REQUIRMENTS:
-Python 3.*
-CASA 6.0 Python modular version (https://casa.nrao.edu/casadocs/casa-6.1.0/usingcasa/obtaining-and-installing)
[ CASA 6.0 or lower is strongly advised to use the imviewer task in parallel with PT-REX]
-matplotlib.pyplot as plt
-astropy
-numpy
-scipy
-bces (https://github.com/rsnemmen/BCES)
##########################

DATA PREPARATION:
The input maps have to re-named as:
-radio map: [name]-radio.fits
-X-rays instensity (counts) map: [name]-img.fits
-X-rays background map: [name]-bmap.fits
-X-rays exposure map: [name]-emap.fits

The region and mask box files must be in CASA region format (.crtf).
##########################

USAGE:

Run PT-REX.py in the same folder of the images and the region files
$ python3 PT-REX.py

      __
      /oo\
     |    |
 ^^  (vvvv)   ^^
 \\  /\__/\  //
  \\/      \//
   /        \
  |          |    ^
  /          \___/ |
 (            )     |
  \----------/     /
    //    \\_____/
   W       W


"""
###########################
import sys

sys.path.append('PTREX/linmix/linmix')

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from astropy.io import fits
import numpy as np
import os
from scipy.stats import norm
import scipy
import scipy.stats
from scipy.optimize import curve_fit
import bces.bces
import linmix
import gc
from casatasks import *
from casatools import *
import readline

readline.parse_and_bind("tab: complete")

casalog.filter('SEVERE')
os.system('mv *.log PTREX/logs')


def fit_func(x, a, b):
    return b * x ** a


def fit_alg(key, x, ex, y, ey, k_gaus, n_chain_lin):  # k,kerr,A,Aerr,
    gc.collect()
    k, kerr, A, Aerr = 0., 0., 0., 0.
    if key == 'LS':
        i_g = [1.0, np.min(y)]

        z, pcov = curve_fit(fit_func, x, y, sigma=ey, p0=i_g, maxfev=1000)
        perr = np.sqrt(np.diag(pcov))
        k, kerr = z[0], perr[0]
        A, Aerr = z[1], perr[1]

    if key == 'BCES_bisec':
        a, b, aerr, berr, covab = bces.bces.bces(np.log10(x), np.divide(ex, x) * 0.434, np.log10(y),
                                                 np.divide(ey, y) * 0.434, np.zeros_like(
                f_x))  # np.cov(np.divide(ex,x)*0.434,np.divide(ey,y)*0.434)
        k, kerr, A, Aerr = a[2], aerr[2], 10. ** b[2], np.log(10.) * (10. ** b[2]) * berr[2]

    if key == 'BCES_ort':
        a, b, aerr, berr, covab = bces.bces.bces(np.log10(x), np.divide(ex, x) * 0.434, np.log10(y),
                                                 np.divide(ey, y) * 0.434, np.zeros_like(f_x))

        k, kerr, A, Aerr = a[3], aerr[3], 10. ** b[3], np.log(10.) * (10. ** b[3]) * berr[3]

    if key == 'LinMix':
        lm = linmix.LinMix(np.log10(x), np.log10(y), xsig=np.array(ex) / np.array(x) * 0.434,
                           ysig=np.array(ey) / np.array(y) * 0.434, K=k_gaus,
                           nchains=n_chain_lin)  # , delta=delta, K=K, nchains=nchains)xycov=xycov
        lm.run_mcmc(silent=True)  # , silent=silent)miniter=10, maxiter=100
        k, kerr = lm.chain['beta'].mean(), lm.chain['beta'].std()
        Am, Aerrm = lm.chain['alpha'].mean(), lm.chain['alpha'].std()
        A = 10. ** Am
        Aerr = np.log(10.) * Aerrm * (10. ** Am)
    return k, kerr, A, Aerr


def start():
    img = []
    bmap = []
    emap = []
    lame = False
    for file in os.listdir("./"):
        if file.endswith("-radio.fits"):
            map_radio = file
        if file.endswith("-img.fits"):
            img.append(file)
        if file.endswith("-bmap.fits"):
            bmap.append(file)
        if file.endswith("-emap.fits"):
            emap.append(file)

    if len(bmap) == 0:
        radio = fits.open(map_radio)
        data = radio[0].data
        header = radio[0].header
        bmap = data * 0
        fits.writeto('dummy_b.fits', bmap, header, overwrite=True)
        bmap = ['dummy_b.fits']
        print('No background map was found. X-ray background is set to 0.')
    if len(emap) == 0:
        lame = True
        radio = fits.open(map_radio)
        header = radio[0].header
        data = radio[0].data
        emap = data / data
        fits.writeto('dummy_e.fits', emap, header, overwrite=True)
        emap = ['dummy_e.fits']
        print('No exposure map was found. X-ray exposure is set to 1.')
    img.sort()
    bmap.sort()
    emap.sort()

    head = imhead(imagename=map_radio)

    if head['restoringbeam']['major']['unit'] == 'arcsec':
        bmaj, bmin, scala = head['restoringbeam']['major']['value'], head['restoringbeam']['minor']['value'], \
                            head['incr'][1] * 206264.8
        beam_area = 3.14 * bmaj * bmin / 4.0 / 0.693147181 / scala ** 2
        beam_area_arcsec = 3.14 * bmaj * bmin / 4.0 / 0.693147181
    if head['restoringbeam']['major']['unit'] == 'deg':
        bmaj, bmin, scala = head['restoringbeam']['major']['value'] * 3600, head['restoringbeam']['minor'][
            'value'] * 3600, head['incr'][1] * 206264.8
        beam_area = 3.14 * bmaj * bmin / 4.0 / 0.693147181 / scala ** 2
        beam_area_arcsec = 3.14 * bmaj * bmin / 4.0 / 0.693147181
    camp = np.sqrt(bmaj * bmin) / scala
    print('-----------------------------------------')
    print('Radio map: ', map_radio)
    print('X-rays map(s):', img)
    print('X-rays background map(s):', bmap)
    print('X-rays exposure map(s):', emap)
    print('Beam radio: ', bmaj, 'X', bmin, ' arcsec')
    print('Beam area [px]: ', beam_area)
    print('Beam area [arcsec]: ', beam_area_arcsec)
    print('1 px= ', scala, ' arcsec')
    print('Beam sampling factor: ', camp)
    print('------------------------------------------')
    cal_e = np.float32(input('Calibration error of the radio image [es. 0.05]: '))
    stat_fit = input('Fitting method [LS | BCES_ort | BCES_bisec | LinMix]: ')
    if stat_fit == 'LinMix':
        k_gaus = int(input('Number of Gaussian for posterior [es. 2]: '))
        n_chain_lin = int(input('Number of MCMC chains [es. 100]: '))
    else:
        k_gaus, n_chain_lin = 0, 0
    return (img, bmap, emap, map_radio, bmaj, bmin, scala, beam_area, beam_area_arcsec, camp, cal_e, stat_fit, k_gaus,
            n_chain_lin, lame)


logo = mpimg.imread('PTREX/logo2.png')
imgplot = plt.imshow(logo)
plt.axis('off')
plt.show(block=False)

print('\n---Point-to-point TRend EXtractor---')
img = []
bmap = []
emap = []
img, bmap, emap, map_radio, bmaj, bmin, scala, beam_area, beam_area_arcsec, camp, cal_e, stat_fit, k_gaus, n_chain_lin, lame = start()

extra = ', linewidth=1, linestyle=-, color=magenta\n'

while True:
    print('-------------------------')
    print('TASK LIST')
    print('-Create a new mask.image from a mask.reg [1]')
    print('-Create a new J2000 mesh.crtf for the radio map [2]')
    print('-Single mesh analysis [3]')
    print('-Monte-Carlo analysis [4]')
    print('UTILITIES:')
    print('--Reload images and configuration [r]')
    print('--quit [q]')
    print('INPUTS:')
    print('Calibration error: ', cal_e)
    print('Fitting method: ', stat_fit)
    print('-------------------------')
    print('   \n')
    task = input('Task to execute: ')
    fit = []

    if task == '1':
        prjct = input('New mask name: ')
        msk_r = input('Mask regions (.crtf): ')
        os.system('rm -r *_temp.image')
        importfits(fitsimage=map_radio, imagename='radio_temp.image', overwrite=True)
        immath(imagename='radio_temp.image', region=msk_r, outfile='mask_temp.image', expr='IM0/IM0')
        imregrid(imagename='mask_temp.image', template='radio_temp.image', decimate=10, interpolation='nearest',
                 output=prjct + '.image', overwrite=True)

        exportfits(imagename=prjct + '.image', fitsimage=prjct + '.fits', overwrite=True)
        os.system('rm -r *_temp.image')
        print('New mask created: ', prjct, '.fits')

    if task == '2':
        flusso = 0.0
        prjct = input('New mesh name: ')
        region = input('Working region [.crtf]: ')
        reg = imstat(imagename=map_radio, region=region)
        xtr, ytr, xbl, ybl = reg['trc'][0], reg['trc'][1], reg['blc'][0], reg['blc'][1]
        msk = input('Mask: ')
        deltax = np.float32(input('Box x-size [pixel]: '))
        deltay = np.float32(input('Box y-size [pixel]: '))
        thresh = np.float32(input('Radio flux density threshold [Jy/beam]: '))
        if beam_area > deltax * deltay:
            choice = (input('WARNING: the box si smaller than the beam. Do you want to choose a new box size? (y/n) '))
            if choice == 'y':
                print('Smallest allowed size: ', np.sqrt(beam_area))
                deltax = np.float32(input('Box x-size [pixel]: '))
                deltay = np.float32(input('box y-size [pixel]: '))
            if choice == 'n':
                print('Proceeding anyway')
        area_arcsec = deltax * deltay * scala ** 2

        g = open(prjct + '.crtf', "w")

        x = np.arange(xbl, xtr, deltax)
        y = np.arange(ybl, ytr, deltay)
        g.write('#CRTFv0 CASA Region Text Format version 0\n')
        print('Raw mesh: ', len(x), 'x', len(y))
        print('Total box: ', len(x) * len(y))
        print('Box size: ', deltax, 'x', deltay, 'pix= ', deltax * deltay, ' pix^2= ', area_arcsec, 'arcsec^2')
        print('-Adapting mesh-')
        i = 0

        for k in range(0, len(x) - 1):

            for j in range(0, len(y) - 1):

                box = 'box [[' + str(x[k + 1]) + 'pix,' + str(y[j + 1]) + 'pix], [' + str(x[k]) + 'pix, ' + str(
                    y[j]) + 'pix]]'
                try:
                    mystat = imstat(imagename=map_radio, region=box, listit=False, verbose=False)
                    # flusso=np.float(mystat['flux'])

                    if np.float32(mystat['flux']) / area_arcsec > thresh / beam_area_arcsec:
                        mystat_m = imstat(imagename=msk, region=box, listit=False, verbose=False)['sum']
                        if not mystat_m > 0.0:
                            g.write('box [[' + mystat['trcf'][0:12] + ', ' + mystat['trcf'][14:27] + '], [' + mystat[
                                                                                                                  'blcf'][
                                                                                                              0:12] + ', ' +
                                    mystat['blcf'][14:27] + ']] coord=J2000' + extra)
                except:
                    continue

                s = ((np.int(i * 20 / ((len(x) - 1) * (len(y) - 1))) * '#') + (
                            np.int(20 - i * 20 / ((len(x) - 1) * (len(y) - 1))) * ' ') + (' ') + str(
                    np.int(i * 100 / ((len(x) - 1) * (len(y) - 1) - 1))) + (' %'))
                sys.stdout.write("\rAdapting - Status: [" + s + "]")
                sys.stdout.flush()
                i = i + 1

        g.close()
        gc.collect()
        print('\nNew mesh created: ', prjct + '.crtf')

    if task == '3':
        prjct = input('Project name: ')
        grid = (input('Mesh name: '))
        griglia = [line.rstrip('\n') for line in open(grid)]

        rms = np.float32(input('RMS of the radio map [Jy/beam]: '))
        p1 = open(prjct + '.dat', "w")
        f_x = []
        e_x = []
        f_r = []
        e_r = []

        top = len(griglia)
        for k in range(1, top):
            try:
                mystat_R = imstat(imagename=map_radio, region=str(griglia[k]), listit=False, verbose=False)
                area_arcsec = (mystat_R['trc'][0] - mystat_R['blc'][0]) * (
                            mystat_R['trc'][1] - mystat_R['blc'][1]) * scala ** 2
                flusso_R = np.float32(mystat_R['flux']) / area_arcsec

                error_r = np.sqrt((cal_e * flusso_R) ** 2 + (
                            rms * np.sqrt(np.float32(mystat_R['npts']) / beam_area) / beam_area_arcsec) ** 2)

                ##################################
                l_img = []
                l_bmap = []
                l_emap = []
                for t in range(0, len(img)):
                    mystat_x = imstat(imagename=img[t], region=str(griglia[k]), listit=False, verbose=False)['sum']
                    l_img.append(np.float32(mystat_x))
                    mystat_x_bkg = imstat(imagename=bmap[t], region=str(griglia[k]), listit=False, verbose=False)['sum']
                    l_bmap.append(np.float32(mystat_x_bkg))
                    mystat_x_exp = imstat(imagename=emap[t], region=str(griglia[k]), listit=False, verbose=False)[
                        'mean']
                    l_emap.append(np.float32(mystat_x_exp))

                flusso_x = (sum(l_img) - sum(l_bmap)) / sum(l_emap) / area_arcsec
                error_x = (np.sqrt(sum(l_img) + sum(l_bmap))) / sum(l_emap) / area_arcsec
                if flusso_x[0] > 0 and flusso_R[0] > 0 and error_x[0] < flusso_x[0] and error_r[0] < flusso_R[0]:
                    f_x.append(flusso_x[0])
                    f_r.append(flusso_R[0])
                    e_x.append(error_x[0])
                    e_r.append(error_r[0])
                if flusso_x < 0 or flusso_R < 0 or error_x > flusso_x:
                    continue
            except:
                continue

            ################
            p1.write(('{} {} {} {}\n'.format(float(flusso_x), float(error_x), float(flusso_R), float(error_r))))

            s = ((np.int(k * 20 / top) * '#') + ((np.int(20 - k * 20 / top)) * ' ') + (' ') + str(
                np.int(k * 100 / (top - 1))) + (' %'))
            sys.stdout.write("\rAnalysys - Status: [" + s + "]")
            sys.stdout.flush()
        # print '\nVIB trovati: ',VIB
        print('\nData file created: ', prjct, '.dat')
        p1.close()
        pers = scipy.stats.pearsonr(np.log10(f_x), np.log10(f_r))
        sper = scipy.stats.spearmanr(np.log10(f_x), np.log10(f_r))
        fit = fit_alg(stat_fit, f_x, e_x, f_r, e_r, k_gaus, n_chain_lin)
        print('Best-fit slope k: ', round(fit[0], 2), '+-', round(fit[1], 2))
        print('Best-fit normalization A: ', round(fit[2], 2), '+-', round(fit[3], 2))
        print('Person coeff: ', str(round(pers[0], 2)), ' P-val: ', str('{:0.3e}'.format(pers[1])))
        print('Spearman coeff: ', str(round(sper[0], 2)), ' P-val: ', str('{:0.3e}'.format(sper[1])))
        x2 = np.linspace(0.3 * np.min(f_x), 1.5 * np.max(f_x), 100)
        fit3 = fit_func(x2, fit[0], fit[2])

        ########PLOT###################
        plt.clf()
        plt.rc('font', size=12)
        plt.rc('axes', labelsize=15)
        plt.rc('legend', fontsize=15)
        # k_mock=[fit[0]-fit[1],fit[0]+fit[1]]

        # A_mock=[np.mean(f_r)/(np.mean(f_x)**k_mock[0]),np.mean(f_r)/(np.mean(f_x)**k_mock[1])]
        # fit_low=fit_func(x2,k_mock[0],A_mock[0])
        # fit_up=fit_func(x2,k_mock[1],A_mock[1])
        # plt.fill_between(x2,fit_low,fit_up,facecolor='blue',alpha=0.3)

        plt.errorbar(f_x, f_r, xerr=e_x, yerr=e_r, linewidth=0, elinewidth=1, color='black', capsize=2)
        plt.plot(x2, fit3, color='blue', linewidth=1.0,
                 label='k$_{SM}$=' + str(round(fit[0], 2)) + '$\pm$' + str(round(fit[1], 2)))
        ###rand k

        term1, term2 = 0.0, 0.0
        conf_Y = []
        for i in range(0, len(f_x)):
            term1 = term1 + ((f_r[i] - fit_func(f_x[i], fit[0], fit[2])) ** 2) / (len(f_x) - 2.)
            term2 = term2 + (f_x[i] - np.mean(f_x)) ** 2
        for i in x2:
            conf_Y.append(1.96 * np.sqrt(term1) * np.sqrt((1. / len(f_x)) + ((i - np.mean(f_x)) ** 2) / term2))

        plt.fill_between(x2, fit3 - conf_Y, fit3 + conf_Y, color='blue', alpha=0.4, linewidth=0)
        # plt.plot(x2,fit_func(x2,fit[0]+fit[1],fit[2]+fit[3]),linewidth=0.5,color='blue',alpha=0.6)
        # plt.plot(x2,fit_func(x2,fit[0]-fit[1],fit[2]-fit[3]),linewidth=0.5,color='blue',alpha=0.6)
        plt.xscale("log")
        plt.yscale("log")
        plt.xlim(0.5 * np.min(f_x), 1.5 * np.max(f_x))
        plt.ylim(0.5 * np.min(f_r), 1.5 * np.max(f_r))
        plt.ylabel('$I_{R}$ [Jy arcsec$^{-2}$]')
        if lame == True:

            plt.xlabel('$I_{X}$ [Arbitrary units arcsec$^{-2}$]')
        else:
            plt.xlabel('$I_{X}$ [ph cm$^{-2}$ s$^{-1}$ arcsec$^{-2}$]')
        # plt.title('k= '+str(a[2])+' st.dev= '+str(aerr[2])+' corr_p= '+str(pers[0]))
        plt.legend()
        plt.show(block=False)
        plt.savefig('plot_' + prjct + '.pdf')
        gc.collect()
    # gc.collect()
    if task == '4':
        prjct = input('Project name: ')
        niter = input('Number of MC iterations: ')
        niter = int(niter)
        region = input('Working region [.crtf]: ')
        reg = imstat(imagename=map_radio, region=region)
        xtr, ytr, xbl, ybl = reg['trc'][0], reg['trc'][1], reg['blc'][0], reg['blc'][1]
        xc, yc, width, height = (xtr + xbl) / 2., (ytr + ybl) / 2., xtr - xbl, ytr - ybl
        msk = input('Mask: ')
        deltax = np.float32(input('Box x-size [pixel]: '))
        deltay = np.float32(input('Box y-size [pixel]: '))
        area_arcsec = deltax * deltay * scala ** 2
        rms = np.float32(input('RMS of the radio map [Jy/beam]: '))
        thresh = np.float32(input('Radio flux density threshold [Jy/beam]: '))
        if beam_area > deltax * deltay:
            choice = (
                input('WARNING: the box size is smaller than the beam. Do you want to choose a new box size? (y/n) '))
            if choice == 'y':
                print('Smallest allowed size: ', np.sqrt(beam_area))
                deltax = np.float32(input('Box x-size [pixel]: '))
                deltay = np.float32(input('box y-size [pixel]: '))
            if choice == 'n':
                print('Proceeding anyway')
        area_arcsec = deltax * deltay * scala ** 2
        rand_x = np.random.randint(xc - width / 4.0, xc + width / 4.0, niter)
        rand_y = np.random.randint(yc - height / 4.0, yc + height / 4.0, niter)
        # slope=[]

        i = 0
        nbox = []

        for i in range(0, niter):

            f_x = []
            e_x = []
            f_r = []
            e_r = []
            fit = []

            cont = 0.0
            x_grid = np.arange(int(rand_x[i] - 1.25 * width), int(rand_x[i] + 1.25 * width), deltax)
            y_grid = np.arange(int(rand_y[i] - 1.25 * height), int(rand_y[i] + 1.25 * height), deltay)
            # g = open('grdi_'+str(i)+'.crtf',"w")
            # g.write('#CRTFv0 CASA Region Text Format version 0\n')
            for k in range(0, len(x_grid) - 1):
                for j in range(0, len(y_grid) - 1):
                    box = 'box [[' + str(x_grid[k + 1]) + 'pix,' + str(y_grid[j + 1]) + 'pix], [' + str(
                        x_grid[k]) + 'pix, ' + str(y_grid[j]) + 'pix]]'

                    try:
                        mystat_R = imstat(imagename=map_radio, region=box, listit=False, verbose=False)
                        # flusso=np.np.float32(mystat_R['flux'])
                        if mystat_R['flux'] / area_arcsec > thresh / beam_area_arcsec:
                            mystat_m = np.float16(imstat(imagename=msk, region=box, listit=False, verbose=False)['sum'])
                            if not mystat_m > 0:
                                cont = cont + 1
                                box_j2000 = 'box [[' + mystat_R['trcf'][0:12] + ', ' + mystat_R['trcf'][
                                                                                       14:27] + '], [' + mystat_R[
                                                                                                             'blcf'][
                                                                                                         0:12] + ', ' + \
                                            mystat_R['blcf'][14:27] + ']] coord=J2000'
                                # g.write(box_j2000+extra)
                                flusso_R = np.float32(mystat_R['flux']) / area_arcsec

                                error_r = np.sqrt((cal_e * flusso_R) ** 2 + (rms * np.sqrt(
                                    np.float32(mystat_R['npts']) / beam_area) / beam_area_arcsec) ** 2)

                                l_img = []
                                l_bmap = []
                                l_emap = []
                                for t in range(0, len(img)):
                                    mystat_x = np.float32(
                                        imstat(imagename=img[t], region=box_j2000, listit=False, verbose=False)['sum'])
                                    l_img.append(mystat_x)
                                    mystat_x_bkg = np.float32(
                                        imstat(imagename=bmap[t], region=box_j2000, listit=False, verbose=False)['sum'])
                                    l_bmap.append(mystat_x_bkg)
                                    mystat_x_exp = np.float32(
                                        imstat(imagename=emap[t], region=box_j2000, listit=False, verbose=False)[
                                            'mean'])
                                    l_emap.append(mystat_x_exp)

                                flusso_x = (sum(l_img) - sum(l_bmap)) / sum(l_emap) / area_arcsec
                                error_x = (np.sqrt(sum(l_img) + sum(l_bmap))) / sum(l_emap) / area_arcsec
                                if flusso_x > 0 and flusso_R > 0 and error_x < flusso_x:
                                    f_x.append(flusso_x[0])
                                    f_r.append(flusso_R[0])
                                    e_x.append(error_x[0])
                                    e_r.append(error_r[0])
                                if flusso_x < 0 or flusso_R < 0 or error_x > flusso_x:
                                    continue
                        # del mystat_R,mystat_m,mystat_x,mystat_x_bkg,mystat_x_exp
                    except:
                        continue
            # g.close()
            nbox.append(cont)
            fit_mc = fit_alg(stat_fit, f_x, e_x, f_r, e_r, k_gaus, n_chain_lin)
            k_boot = np.random.normal(fit_mc[0], fit_mc[1],
                                      1)  # BOOTSTRAPING del valore da usare nella distribuzione!!!!
            # g.close()
            gc.collect()
            with open(prjct + '.dat', "a") as myfile:
                myfile.write(str(float(k_boot)) + '\n')
            # del f_x,f_r,e_x,e_r,l_img,l_emap,l_bmap,x,y,box,box_j2000

            s = ((np.int(i * 20 / (niter)) * '#') + (np.int(20 - i * 20 / (niter)) * ' ') + (' ') + str(
                np.int(i * 100 / (niter - 1))) + (' %'))
            sys.stdout.write("\rMC analysis - Status: [" + s + "]")
            sys.stdout.flush()
        # print i,z[1]
        # gc.collect()

        # g.close()
        slope = np.genfromtxt(prjct + '.dat')
        print('\nMean slope: ', np.mean(slope))
        print('st. dev: ', np.std(slope))
        print('Boxes: ', np.mean(nbox), ' (', np.std(nbox), ')')

        ########PLOT###################
        # sc3=(input('Pulire plot? (s/n) '))
        # if sc3=='s':
        #	plt.cla()

        plt.clf()
        plt.rc('font', size=12)
        plt.rc('axes', labelsize=15)
        plt.rc('legend', fontsize=15)
        plt.hist(slope, color='blue', alpha=0.6)

        plt.title('$k_{MC}=$' + str(round(np.mean(slope), 2)) + '$\pm$' + str(round(np.std(slope), 2)))
        # plt.text(0.1,0.9,'media: '+str(np.mean(slope))+' st.dev: '+str(np.std(slope)))
        plt.xlabel('$k_{SM}')
        plt.ylabel('Frequency')
        plt.show(block=False)
        plt.savefig('plot_' + prjct + '.pdf')
    # gc.collect()

    if task == 'r':
        img = []
        bmap = []
        emap = []
        img, bmap, emap, map_radio, bmaj, bmin, scala, beam_area, beam_area_arcsec, camp, cal_e, stat_fit, k_gaus, n_chain_lin, lame = start()
    # plt.close()
    # logo=mpimg.imread('casa6/logo2.png')
    # imgplot = plt.imshow(logo)
    # plt.axis('off')
    # plt.show(block=False)
    if task == 'q':
        print('Good bye!')
        break