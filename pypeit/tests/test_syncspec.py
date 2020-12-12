import os

import pytest

from pypeit.scripts import coadd_1dspec
from pypeit.tests.tstutils import dev_suite_required, cooked_required

#def test_lrisr():
#    dpath = '/home/xavier/scratch/FRB190714/keck_lris_red_C/Science'
#    file1 = 'spec1d_LR.20200128.55304-frb1907B_LRISr_2020Jan28T152144.726.fits'
#    file2 = 'spec1d_LR.20200128.56019-frb1907B_LRISr_2020Jan28T153339.168.fits'
#    file3 = 'spec1d_LR.20200128.56734-frb1907B_LRISr_2020Jan28T154534.733.fits'
#    file4 = 'spec1d_LR.20200128.57449-frb1907B_LRISr_2020Jan28T155729.779.fits'
#
#    coadd_1dspec.coadd1d_filelist([os.path.join(dpath, file1),
#                                   os.path.join(dpath, file2),
#                                   os.path.join(dpath, file3),
#                                   os.path.join(dpath, file4)],
#                                  'tst_lris_files',
#                                  1)

@cooked_required
def test_kastb():
    dpath = os.path.join(os.environ['PYPEIT_DEV'], 'Cooked', 'Science')
    file1 = 'spec1d_b27-J1217p3905_KASTb_2015May20T045733.560.fits' #SPAT0176-SLIT0000-DET01
    file2 = 'spec1d_b28-J1217p3905_KASTb_2015May20T051801.470.fits' #SPAT0175-SLIT0000-DET01

    outfiles = coadd_1dspec.coadd1d_filelist([os.path.join(dpath, file1),
                                             os.path.join(dpath, file2)],
                                             'tst_coadd_files', 1)
    
    for ofile in outfiles:
        os.remove(ofile)

