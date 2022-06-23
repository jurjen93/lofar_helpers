import matplotlib.pyplot as plt
from astropy.visualization import make_lupton_rgb
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename

# Read in the three images downloaded from here:
g_name = get_pkg_data_filename('/home/jurjen/Documents/Python/lofar_helpers/analysis/t044.5800+13.3000.g.fits')
r_name = get_pkg_data_filename('/home/jurjen/Documents/Python/lofar_helpers/analysis/t044.5800+13.3000.r.fits')
i_name = get_pkg_data_filename('/home/jurjen/Documents/Python/lofar_helpers/analysis/t044.5800+13.3000.i.fits')

g = fits.open(g_name)[0].data * 1.5
r = fits.open(r_name)[0].data * 2
i = fits.open(i_name)[0].data * 1


rgb_default = make_lupton_rgb(i, r, g, minimum=300, stretch=400, Q=30, filename="test.jpeg")
plt.imshow(rgb_default, origin='lower')
plt.show()

# table = getimages([crd[0]], [crd[1]], filters="rgi")
#
# fnames = []
#
# for row in table:
#     if k==1:
#         ra = row['ra']+0.1
#     else:
#         ra = row['ra']
#     dec = row['dec']
#     projcell = row['projcell']
#     subcell = row['subcell']
#     filter = row['filter']
#
#     # create a name for the image -- could also include the projection cell or other info
#     fname = "t{:08.4f}{:+07.4f}.{}.fits".format(ra, dec, filter)
#     fnames.append(fname)
#
#     # if '/home/jurjen/Documents/Python/lofar_helpers/analysis/'+fname not in glob('/home/jurjen/Documents/Python/lofar_helpers/analysis/*'):
#     url = row["url"]
#     print("%11.6f %10.6f skycell.%4.4d.%3.3d %s" % (ra, dec, projcell, subcell, fname))
#     r = requests.get(url)
#     open(fname, "wb").write(r.content)
#     # else:
#     #     print(fname+' already exists')


# for fname in fnames:
#     if '.g.' in fname:
#         self.reproject_map('/home/jurjen/Documents/Python/lofar_helpers/analysis/'+fname, '/home/jurjen/Documents/Python/lofar_helpers/analysis/'+fname.replace('.fits', '_reproject.fits'))
#     elif '.r.' in fname:
#         self.reproject_map('/home/jurjen/Documents/Python/lofar_helpers/analysis/'+fname, '/home/jurjen/Documents/Python/lofar_helpers/analysis/'+fname.replace('.fits', '_reproject.fits'))
#     elif '.i.' in fname:
#         self.reproject_map('/home/jurjen/Documents/Python/lofar_helpers/analysis/'+fname, '/home/jurjen/Documents/Python/lofar_helpers/analysis/'+fname.replace('.fits', '_reproject.fits'))
#
# for fname in fnames:
#     if '.g.' in fname:
#         g = fits.open('/home/jurjen/Documents/Python/lofar_helpers/analysis/'+fname.replace('.fits', '_reproject.fits'))[0].data * 1.5
#         gcut = Cutout2D(
#             data=g,
#             position=(s[0], s[1]),
#             size=(s[3], s[2]),
#             wcs=self.wcs,
#             mode='partial'
#         )
#         g = gcut.data
#         del gcut
#     elif '.r.' in fname:
#         r = fits.open('/home/jurjen/Documents/Python/lofar_helpers/analysis/'+fname.replace('.fits', '_reproject.fits'))[0].data * 2
#         rcut = Cutout2D(
#             data=r,
#             position=(s[0], s[1]),
#             size=(s[3], s[2]),
#             wcs=self.wcs,
#             mode='partial'
#         )
#         r = rcut.data
#         del rcut
#     elif '.i.' in fname:
#         i = fits.open('/home/jurjen/Documents/Python/lofar_helpers/analysis/'+fname.replace('.fits', '_reproject.fits'))[0].data * 1
#         icut = Cutout2D(
#             data=i,
#             position=(s[0], s[1]),
#             size=(s[3], s[2]),
#             wcs=self.wcs,
#             mode='partial'
#         )
#         i = icut.data
#         del icut
#
# rgb_default = make_lupton_rgb(i, r, g, filename="test.jpeg")
# plt.imshow(rgb_default, origin='lower')
# plt.show()
