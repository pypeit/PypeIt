""" Simple script to copy over key spec1d files from the DevSuite"""

import os
import shutil

# shane_kast_blue
sfile = os.path.join(os.getenv('PYPEIT_DEV'), 'REDUX_OUT', 'shane_kast_blue',
                     '600_4310_d55', 'shane_kast_blue_A', 'Science',
                     'spec1d_b27-J1217p3905_KASTb_2015May20T045733.560.fits')
shutil.copyfile(sfile, 'spec1d_b27.fits')
sfile = os.path.join(os.getenv('PYPEIT_DEV'), 'REDUX_OUT', 'shane_kast_blue',
                     '600_4310_d55', 'shane_kast_blue_A', 'Science',
                     'spec1d_b28-J1217p3905_KASTb_2015May20T051801.470.fits')
shutil.copyfile(sfile, 'spec1d_b28.fits')

# shane kast red
sfile = os.path.join(os.getenv('PYPEIT_DEV'), 'REDUX_OUT', 'shane_kast_red',
                     '600_7500_d55_ret', 'Science',
                     'spec1d_r153-J0025-0312_KASTr_2015Jan23T025323.850.fits')
shutil.copyfile(sfile, 'spec1d_r153-J0025-0312_KASTr_2015Jan23T025323.850.fits')

# Gemini GNIRS
basefiles = ['spec1d_cN20170331S0216-pisco_GNIRS_2017Mar31T085412.181.fits',
             'spec1d_cN20170331S0217-pisco_GNIRS_2017Mar31T085933.097.fits',
             'spec1d_cN20170331S0220-pisco_GNIRS_2017Mar31T091538.731.fits',
             'spec1d_cN20170331S0221-pisco_GNIRS_2017Mar31T092059.597.fits']
gpath = os.path.join(os.getenv('PYPEIT_DEV'), 'REDUX_OUT', 'gemini_gnirs',
                     '32_SB_SXD', 'Science')

for basefile in basefiles:
    shutil.copyfile(os.path.join(gpath, basefile), basefile)

