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

outpath=resource_filename('pypeit', 'data/arc_lines/reid_arxiv')

def build_template(in_file, slits, wv_cuts, binspec, outroot, lowredux=True):
    # Load xidl file
    # Grab it
    # Load and splice
    yvals = []
    lvals = []
    for kk, slit in enumerate(slits):
        if lowredux:
            wv_vac, spec = xidl_arcspec(in_file, slit)
        else:
            wv_vac, spec = pypeit_arcspec(in_file, slit)
        # Cut
        if len(slits) > 1:
            if kk == 0:
                llow = 0.
                lhi = wv_cuts[0]
            elif kk == len(slits)-1:
                llow = wv_cuts[kk-1]
                lhi = 1e9
            else:
                llow = wv_cuts[kk-1]
                lhi = wv_cuts[kk]
            #
            gdi = (wv_vac > llow) & (wv_vac < lhi)
        else:
            gdi = np.arange(spec.size).astype(int)
        # Append
        yvals.append(spec[gdi])
        lvals.append(wv_vac[gdi])
    # Concatenate
    nwspec = np.concatenate(yvals)
    nwwv = np.concatenate(lvals)
    # Generate the table
    write_template(nwwv, nwspec, binspec, outpath, outroot)


def pypeit_arcspec(in_file, slit):
    wv_dict = ltu.loadjson(in_file)
    iwv_calib = wv_dict[str(slit)]
    x = np.arange(iwv_calib['nspec'])
    wv_vac = utils.func_val(iwv_calib['fitc'], x/iwv_calib['xnorm'], iwv_calib['function'],
                           minx=iwv_calib['fmin'], maxx=iwv_calib['fmax'])
    # Return
    return wv_vac, np.array(iwv_calib['spec'])

def write_template(nwwv, nwspec, binspec, outpath, outroot):
    tbl = Table()
    tbl['wave'] = nwwv
    tbl['flux'] = nwspec
    tbl.meta['INSTR'] = 'keck_lris_blue_600_d560'
    tbl.meta['BINSPEC'] = binspec
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


def xidl_arcspec(xidl_file, slit):
    xidl_dict = readsav(xidl_file)
    nspec = xidl_dict['archive_arc'].shape[0]
    npix = xidl_dict['archive_arc'].shape[1]
    # This is the best one (well-centered)
    calib = xidl_dict['calib'][slit]
    # Generate the wavelengths
    wv_air = cheby_val(calib['FFIT'], np.arange(npix),
                       calib['NRM'], calib['NORD'])
    wv_vac = airtovac(wv_air * units.AA)
    # Flip to blue to red
    wv_vac = wv_vac[::-1]
    spec = xidl_dict['archive_arc'][slit][::-1]

    # Return
    return wv_vac.value, spec

def main(flg):

    # Keck LRISb
    if flg & (2**0): # B300, all lamps
        binspec = 1
        slits = [15]
        xidl_file = os.path.join(os.getenv('XIDL_DIR'),
                             'Spec/Longslit/calib/linelists/lris_blue_300.sav')
        outroot = 'keck_lris_blue_300_d680.fits'
        build_template(xidl_file, slits, None, binspec, outroot, lowredux=True)

    if flg & (2**1): # B400, all lamps I think)
        binspec = 2
        outroot='keck_lris_blue_400_d560.fits'
        slits = [19,14]
        lcut = [5500.]
        xidl_file = os.path.join(os.getenv('XIDL_DIR'),
                                 'Spec/Longslit/calib/linelists/lris_blue_400_d560.sav')
        build_template(xidl_file, slits, lcut, binspec, outroot, lowredux=True)

    if flg & (2**2): # B600, all lamps
        binspec = 2
        outroot='keck_lris_blue_600_d560.fits'
        slits = [0,7]
        lcut = [4500.]
        wfile = os.path.join(os.getenv('PYPEIT_DEV'), 'REDUX_OUT/Keck_LRIS_blue/multi_600_4000_d560',
                         'MF_keck_lris_blue', 'MasterWaveCalib_A_1_01.json')
        build_template(wfile, slits, lcut, binspec, outroot, lowredux=False)

# Command line execution
if __name__ == '__main__':

    # LRISb
    flg = 0
    #flg += 2**0  # LRISb 300, all lamps
    flg += 2**1  # LRISb 400, all lamps
    #flg += 2**2  # LRISb 600, all lamps

    main(flg)

