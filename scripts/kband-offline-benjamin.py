# now with two-phases

# OfflinePipeline2

import os
import glob
import itertools
from collections import OrderedDict
from Entities import EntitiesError, ObsEntity


project = '08-21'

FITS_DIR = '/daten/mbfits/mbfits-2021-12'
CLASS_DIR = '/homes/ygong/'
FEBES = ['S14mm-WFFTS', 'S14mm-OPTOCBE']
GAIN = (0.95400, 3.1900e-3, -5.4200e-5)
TAPERS = {
    # BB: edge_db, edge_order
    1: (-15.7, 3),  # 22250 MHz
    2: (-18.2, 3),  # 23770 MHz
    }
# TCAL_VALUES = {
#     # BB: (tcal1, .., tcaln) with n: number of pixels (incl. polarizations)
#     1: (11.61, 11.89, 1., 1.),  # 1L, 1R, 2L, 2R (latter two are dummies)
#     2: (9.8, 10.15, 1., 1.),
#    }
TCAL_VALUES = {  # April 2019
    # BB: (tcal1, .., tcaln) with n: number of pixels (incl. polarizations)
    1: (12.87, 13.09, 1., 1.),  # 1L, 1R, 2L, 2R (latter two are dummies)
    2: (10.66, 10.81, 1., 1.),
    }


proj_files = []
for febe in FEBES:
    proj_files.extend(glob.glob(os.path.join(
        FITS_DIR, 'EFFBG_2021-??-??_????_{:s}_{:s}.fits'.format(project, febe)
        )))

scan_dict = OrderedDict()
scan_dict['.'] = [
    int(g.split('_')[2])
    for g in proj_files
    ]


error_list = []

for sub_rawdir, scan_list in scan_dict.items():

    setFitsDir(os.path.join(FITS_DIR, sub_rawdir))

    for scan, febe in itertools.product(scan_list, FEBES):

        scan_info = getScanInfo(scan)
        basebands = sorted(scan_info.getBasebandList(febe))
        subscans = scan_info.getSubscanNumList()

        for bband in basebands:

            reverseArrayOrientation = False

            # read bandwidth
            obs = ObsEntity()
            mbf = internal.getDatasetForScan(scan)
            obs.fill(mbf, 1, febe, bband)
            bw = int(obs.BANDWID / 1.e6 + 0.5)
            if bw != 300:
                continue

            setClassName(os.path.join(
                CLASS_DIR,
                "<project>_bb{}.100m".format(bband)
                ))

            for subscan in subscans:
                try:
                    reduceSubscan(
                        scan, subscan, baseband=bband,
                        regridMethod='gaussian', edgePercent=5,
                        applyAatm=True,
                        reverseArrayOrientation=reverseArrayOrientation,
                        gainCurve=GAIN, taperParameters=TAPERS[bband],
                        tcal=TCAL_VALUES[bband],
                        )
                except (IndexError, EntitiesError, TypeError):
                    error_list.append(
                        (sub_rawdir, scan, bband, subscan)
                        )
                    continue

print('error_list', error_list)
