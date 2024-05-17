'''
This is a simple script to get a small cutout from a big image
'''

from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS

position = (548, 2070)   ## pixel position
size = (120,120)         ## pixel size

hdu = fits.open("G0.0+000-VE-avg-I.fits")[0]
wcs = WCS(hdu.header)
cutout = Cutout2D(hdu.data, position=position, size=size, wcs=wcs)
hdu.data = cutout.data
hdu.header.update(cutout.wcs.to_header())
hdu.writeto("G1.6.fits", overwrite=True)
