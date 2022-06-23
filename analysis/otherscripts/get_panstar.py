"""ps1bulk.py: Get PS1 stack FITS cutout images at a list of positions"""

import numpy as np
from astropy.table import Table
import requests
import time
from io import StringIO

ps1filename = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
fitscut = "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi"


def getimages(tra, tdec, size=5000, filters="grizy", format="fits", imagetypes="stack"):
    """Query ps1filenames.py service for multiple positions to get a list of images
    This adds a url column to the table to retrieve the cutout.

    tra, tdec = list of positions in degrees
    size = image size in pixels (0.25 arcsec/pixel)
    filters = string with filters to include
    format = data format (options are "fits", "jpg", or "png")
    imagetypes = list of any of the acceptable image types.  Default is stack;
        other common choices include warp (single-epoch images), stack.wt (weight image),
        stack.mask, stack.exp (exposure time), stack.num (number of exposures),
        warp.wt, and warp.mask.  This parameter can be a list of strings or a
        comma-separated string.

    Returns an astropy table with the results
    """

    if format not in ("jpg", "png", "fits"):
        raise ValueError("format must be one of jpg, png, fits")
    # if imagetypes is a list, convert to a comma-separated string
    if not isinstance(imagetypes, str):
        imagetypes = ",".join(imagetypes)
    # put the positions in an in-memory file object
    cbuf = StringIO()
    cbuf.write('\n'.join(["{} {}".format(ra, dec) for (ra, dec) in zip(tra, tdec)]))
    cbuf.seek(0)
    # use requests.post to pass in positions as a file
    r = requests.post(ps1filename, data=dict(filters=filters, type=imagetypes),
                      files=dict(file=cbuf))
    r.raise_for_status()
    tab = Table.read(r.text, format="ascii")

    urlbase = "{}?size={}&format={}".format(fitscut, size, format)
    tab["url"] = ["{}&ra={}&dec={}&red={}".format(urlbase, ra, dec, filename)
                  for (filename, ra, dec) in zip(tab["filename"], tab["ra"], tab["dec"])]
    return tab


if __name__ == "__main__":
    t0 = time.time()

    # create a test set of image positions
    tdec = [13.3]
    tra = [44.6]

    # get the PS1 info for those positions
    table = getimages(tra, tdec, filters="rgbi")
    print("{:.1f} s: got list of {} images for {} positions".format(time.time() - t0, len(table), len(tra)))

    # extract cutout for each position/filter combination
    for row in table:
        ra = row['ra']
        dec = row['dec']
        projcell = row['projcell']
        subcell = row['subcell']
        filter = row['filter']

        # create a name for the image -- could also include the projection cell or other info
        fname = "t{:08.4f}{:+07.4f}.{}.fits".format(ra, dec, filter)

        url = row["url"]
        print("%11.6f %10.6f skycell.%4.4d.%3.3d %s" % (ra, dec, projcell, subcell, fname))
        r = requests.get(url)
        open(fname, "wb").write(r.content)
    print("{:.1f} s: retrieved {} FITS files for {} positions".format(time.time() - t0, len(table), len(tra)))