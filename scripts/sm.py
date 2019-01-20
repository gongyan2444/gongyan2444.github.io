#####################################
## This is designed to convolve   ###
## single-dish data to coarser    ###
## angular resolutions and match  ###
## the grid as the other file.    ###
## contact: gongyan2444@gmail.com ###
## python sm.py -h ! help info    ###
#####################################
from astropy.io import fits 
import scipy.ndimage as ndimage
from pylab import *
import sys
import numpy as np
from FITS_tools.hcongrid import hcongrid  
import argparse

parser=argparse.ArgumentParser(
    description='''Usage: python sm.py infile tfile outname reso''',
    epilog="""End.""")
parser.add_argument('infile', type=str, help='the name of input file')
parser.add_argument('tfile', type=str, help='fits with target header')
parser.add_argument('outname', type=str, help='the name of output file')
parser.add_argument('reso', type=float, help='angular resolution (arcsec)')
args=parser.parse_args()




####input parameters####
infile  = sys.argv[1]       # input
tfile   = sys.argv[2]       # target header
outname = sys.argv[3]       # output
reso    = float(sys.argv[4])# in arcsec

#############################
#### Main code ##############
#############################
hd      = fits.getheader(infile) 
dat     = fits.getdata(infile)
projhd  = fits.getheader(tfile)


cornel = np.sqrt(reso**2.-(hd['BMAJ']*3600.)**2.)
sigma = cornel/(hd['CDELT2']*3600.)/(2.*np.sqrt(2.*np.log(2.)))

### Gaussian smoothing the image ######
gimage=ndimage.filters.gaussian_filter(dat, sigma, order=0, output=None, mode='reflect', cval=0.0)

### project image  ####
new_image = hcongrid(gimage, hd, projhd)

fits.writeto(outname, new_image, header=projhd,overwrite=True)




