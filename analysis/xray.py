from astropy.io import fits

hdu = fits.open('../fits/mosaic_a399_a401.fits')
hdu_bkg = fits.open('../fits/mosaic_a399_a401_bkg.fits')
hdu_exp = fits.open('../fits/mosaic_a399_a401_exp.fits')

data = hdu[0].data
data_bkg = hdu_bkg[0].data
data_exp = hdu_exp[0].data

print(data.shape)
print(data_bkg.shape)
print(data_exp.shape)

findat=(data - data_bkg)/data_exp