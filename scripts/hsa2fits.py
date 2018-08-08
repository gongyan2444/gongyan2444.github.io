from herschel.hifi.dp.tools import polarPair
#
# load observation into session
obsid = 1342268465
obs = getObservation(obsid, useHsa=True)
#
# extract wideband spectra (WBS) in H and V polarizations
# box 001 product 0001
wbs_h = obs.refs["level2"].product.refs["WBS-H-LSB"].product.refs["box_001"].product["0001"]
wbs_v = obs.refs["level2"].product.refs["WBS-V-LSB"].product.refs["box_001"].product["0001"]
#
# create polar pairs
pp_wbs = polarPair(ds1=wbs_h, ds2=wbs_v, tolerance=100)
#
# stitch spectra
st_wbs = stitch(ds=pp_wbs, variant = 'midPoints', unit = 'kHz')
#
# export to class file in FITS format
hiClass(data = st_wbs, fileName = '/Users/ygong/465.fits')

# ------------------------------------------------------------------------------ #

# extract wideband spectra (WBS) in H and V polarizations
# box 001 product 0002
wbs_h = obs.refs["level2"].product.refs["WBS-H-LSB"].product.refs["box_001"].product["0001"]
wbs_v = obs.refs["level2"].product.refs["WBS-V-LSB"].product.refs["box_001"].product["0001"]
#
# create polar pairs
pp_wbs = polarPair(ds1=wbs_h, ds2=wbs_v, tolerance=100)
#
# average H and V polarizations
#wbs_av = pp_wbs.avg()
#
# stitch spectra
st_wbs = stitch(ds=pp_wbs, variant = 'midPoints', unit = 'kHz')
#
# export to class file in FITS format
hiClass(data = st_wbs, fileName = '/Users/ygong/try2.fits')
#HiClassTask()(dataset = st_wbs, fileName = '/aux/pc20001a/somewhere/1342190101-wbs-0002.fits')

# ------------------------------------------------------------------------------ #

# extract wideband spectra (WBS) in H and V polarizations
# box 001 product 0003
wbs_h = obs.refs["level2"].product.refs["WBS-H-USB"].product.refs["box_001"].product["0003"]
wbs_v = obs.refs["level2"].product.refs["WBS-V-USB"].product.refs["box_001"].product["0003"]
#
# create polar pairs
pp_wbs = polarPair(ds1=wbs_h, ds2=wbs_v, tolerance=100)
#
# average H and V polarizations
#wbs_av = pp_wbs.avg()
#
# stitch spectra
st_wbs = stitch(ds=pp_wbs, variant = 'midPoints', unit = 'kHz')
#
# export to class file in FITS format
#HiClassTask()(dataset = st_wbs, fileName = '/aux/pc20001a/somewhere/1342190101-wbs-0003.fits')

# ------------------------------------------------------------------------------ #
