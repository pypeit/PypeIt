""" Generate the wavelength templates for MDM/Modspec """

import os

from pypeit.core.wavecal import templates
import numpy as np

def mdm_modspec_echelle(overwrite=True):
    """mdm_modspec_echelle Template from the 1200/5100 grism [INSERT UNITS HERE]

    Args:
        overwrite: bool, optional
          Overwrite the existing file? [Default: False]
    """

    binspec = 1
    outroot = 'mdm_modspec_echelle_Ar.fits'

    # PypeIt fits
    ##IF could find dev files/folder system, would set up properly and somewhat like this:
    ## wpath = os.path.join(templates.template_path, 'mdm_modspec_echelle')
    ##basefiles = ['MasterWaveCalib_ArI_5100.fits','MasterWaveCalib_NeI_5100.fits','MasterWaveCalib_XeI_5100.fits'] 
    ##wfiles = [os.path.join(wpath, basefile) for basefile in basefiles]
    
    ##but since I cannot find them, for now:
    #filePath = r'C:\Users\thela\Desktop\coolege\Fall2022\upTo\Summer\Space_ze_Final_Frontier\python_scripts\testingSite\folderForReductionTest\testHolyGrail5\__template'
    #basefiles = ['MasterWaveCalib_ArI_5100.fits','MasterWaveCalib_NeI_5100.fits','MasterWaveCalib_XeI_5100.fits']
    #wfiles = [os.path.join(filePath, basefile) for basefile in basefiles] ## these are the MasterWaveCalib files

    filePath = r'C:\Users\thela\Desktop\coolege\Fall2022\upTo\Summer\Space_ze_Final_Frontier\python_scripts\testingSite\folderForReductionTest\testHolyGrail6\ArI\mdm_modspec_echelle_A'
    wfiles = os.path.join(filePath, 'wvcalib.fits')

    # Snippets
    slits = [150] ## spatial id
    wv_cuts = []
    assert len(wv_cuts) == len(slits)-1
    
    
    templates.build_template(wfiles, slits, wv_cuts, binspec, outroot,
                             det_cut=None, chk=True,
                             normalize=True, lowredux=False,
                             subtract_conti=True, overwrite=overwrite,
                             shift_wave=True)
    

    
if __name__ == '__main__':
    mdm_modspec_echelle(overwrite=True)

