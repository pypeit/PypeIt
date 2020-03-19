import os

from pypeit.scripts import coadd_1dspec

dpath = '/home/xavier/local/Python/PypeIt-development-suite/REDUX_OUT/Shane_Kast_blue/600_4310_d55/shane_kast_blue_A'
file1 = 'Science/spec1d_b27-J1217p3905_KASTb_2015May20T045733.560.fits' #SPAT0176-SLIT0000-DET01
file2 = 'Science/spec1d_b28-J1217p3905_KASTb_2015May20T051801.470.fits' #SPAT0175-SLIT0000-DET01

coadd_1dspec.coadd1d_filelist([os.path.join(dpath, file1),
                            os.path.join(dpath, file2)],
                           'tst_coadd_files',
                           1)
