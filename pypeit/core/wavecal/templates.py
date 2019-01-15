""" Module to generate templates for the PypeIt full_template wavelength calibration routine"""

import os
import numpy as np
import pdb

from pkg_resources import resource_filename
from scipy.io import readsav

from astropy.table import Table
from astropy import units

from linetools import utils as ltu

from pypeit import utils
from pypeit.core.wave import airtovac


# Data Model
# FITS table
#  wave -- Wavelength values
#  flux -- Arc spectrum flux values
#
# Meta must include BINNING of the template with 1=native

# LRIS
def build_LRISb_300(outroot='keck_lris_blue_300_d680.fits',
                    outpath=resource_filename('pypeit', 'data/arc_lines/reid_arxiv')):
    # Load xidl file
    xidl_file = os.path.join(os.getenv('XIDL_DIR'),
                             'Spec/Longslit/calib/linelists/lris_blue_300.sav')
    xidl_dict = readsav(xidl_file)
    nspec = xidl_dict['archive_arc'].shape[0]
    npix = xidl_dict['archive_arc'].shape[1]
    # This is the best one (well-centered)
    slit = 15
    calib = xidl_dict['calib'][slit]
    # Generate the wavelengths
    wv_air = cheby_val(calib['FFIT'], np.arange(npix),
                       calib['NRM'], calib['NORD'])
    wv_vac = airtovac(wv_air * units.AA)
    # Flip to blue to red
    wv_vac = wv_vac[::-1]
    spec = xidl_dict['archive_arc'][slit][::-1]
    # Table
    tbl = Table()
    tbl['wave'] = wv_vac.value
    tbl['flux'] = spec
    tbl.meta['INSTR'] = 'keck_lris_blue_300_d680'
    tbl.meta['BINSPEC'] = 1
    # Write
    outfile = os.path.join(outpath, outroot)
    tbl.write(outfile, overwrite=True)
    print("Wrote: {}".format(outfile))


def build_LRISb_600(outroot='keck_lris_blue_600_d560.fits',
              outpath=resource_filename('pypeit', 'data/arc_lines/reid_arxiv')):
    # Load from the Dev suite
    wfile = os.path.join(os.getenv('PYPEIT_DEV'), 'REDUX_OUT/Keck_LRIS_blue/multi_600_4000_d560',
                         'MF_keck_lris_blue', 'MasterWaveCalib_A_1_01.json')
    wv_dict = ltu.loadjson(wfile)
    # Hard-coded parameters (could be wrong in a future version..
    good = [0,7]
    lcut = 4500.
    # Load and splice
    yvals = []
    lvals = []
    for kk, slit in enumerate(good):
        iwv_calib = wv_dict[str(slit)]
        x = np.arange(iwv_calib['nspec'])
        tmpwv = utils.func_val(iwv_calib['fitc'], x/iwv_calib['xnorm'], iwv_calib['function'],
                                           minx=iwv_calib['fmin'], maxx=iwv_calib['fmax'])
        #
        if kk == 0:
            gdi = tmpwv < lcut
        else:
            gdi = tmpwv > lcut
        # Save
        yvals.append(np.array(iwv_calib['spec'])[gdi])
        lvals.append(tmpwv[gdi])
    nwspec = np.concatenate(yvals)
    nwwv = np.concatenate(lvals)
    # Generate the table
    tbl = Table()
    tbl['wave'] = nwwv
    tbl['flux'] = nwspec
    tbl.meta['INSTR'] = 'keck_lris_blue_600_d560'
    tbl.meta['BINNING'] = 2
    # Write
    outfile = os.path.join(outpath, outroot)
    tbl.write(outfile, overwrite=True)
    print("Wrote: {}".format(outfile))

## Low-Redux routines

def fcheby(xnrm,order):
    leg = np.zeros((len(xnrm),order))
    leg[:,0] = 1.
    if order >= 2:
        leg[:,1] = xnrm
    # For loop
    for j in range(2,order):
        leg[:,j] = 2.0 * xnrm * leg[:,j-1] - leg[:,j-2]
    # Return
    return leg

def cheby_val(coeff, x, nrm, order):
    #
    xnrm = 2. * (x - nrm[0])/nrm[1]
    # Matrix first
    leg = fcheby(xnrm, order)
    # Dot
    return np.dot(leg, coeff)

# Command line execution
if __name__ == '__main__':
    # LRISb
    build_LRISb_300()
    build_LRISb_600()

