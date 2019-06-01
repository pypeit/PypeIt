""" Module to generate templates for the PypeIt full_template wavelength calibration routine"""

import os
import numpy as np
from IPython import embed

from pkg_resources import resource_filename
from scipy.io import readsav

from astropy.table import Table
from astropy import units

from linetools import utils as ltu

from pypeit import utils
from pypeit.core.wave import airtovac
from pypeit.core.wavecal import waveio
from pypeit.core.wavecal import autoid
from pypeit.core.wavecal import fitting

from pypeit import debugger

# Data Model
# FITS table
#  wave -- Wavelength values
#  flux -- Arc spectrum flux values
#
# Meta must include BINNING of the template with 1=native
template_path = os.path.join(os.getenv('PYPEIT_DEV'), 'dev_algorithms/wavelengths/template_files/')

outpath=resource_filename('pypeit', 'data/arc_lines/reid_arxiv')

def build_template(in_files, slits, wv_cuts, binspec, outroot,
                   lowredux=True, ifiles=None, det_cut=None, chk=False):
    """
    Generate a full_template for a given instrument

    Args:
        in_files (list or str):
        slits:
        wv_cuts:
        binspec:
        outroot:
        lowredux:
        ifiles:
        det_cut:
        chk:

    Returns:

    """
    # Load xidl file
    # Grab it
    # Load and splice
    yvals = []
    lvals = []
    if not isinstance(in_files,list):
        in_files = [in_files]
        ifiles = [0]*len(slits)
    for kk, slit in enumerate(slits):
        in_file = in_files[ifiles[kk]]
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
    # Check
    if chk:
        debugger.plot1d(nwwv, nwspec)
    # Generate the table
    write_template(nwwv, nwspec, binspec, outpath, outroot, det_cut=det_cut)


def pypeit_arcspec(in_file, slit):
    wv_dict = ltu.loadjson(in_file)
    iwv_calib = wv_dict[str(slit)]
    x = np.arange(iwv_calib['nspec'])
    wv_vac = utils.func_val(iwv_calib['fitc'], x/iwv_calib['xnorm'], iwv_calib['function'],
                           minx=iwv_calib['fmin'], maxx=iwv_calib['fmax'])
    # Return
    return wv_vac, np.array(iwv_calib['spec'])


def write_template(nwwv, nwspec, binspec, outpath, outroot, det_cut=None):
    tbl = Table()
    tbl['wave'] = nwwv
    tbl['flux'] = nwspec
    tbl.meta['BINSPEC'] = binspec
    # Detector snippets??
    if det_cut is not None:
        tbl['det'] = 0
        for dets, wcuts in zip(det_cut['dets'], det_cut['wcuts']):
            gdwv = (tbl['wave'] > wcuts[0]) & (tbl['wave'] < wcuts[1])
            deti = np.sum([2**ii for ii in dets])
            tbl['det'][gdwv] += deti
    # Write
    outfile = os.path.join(outpath, outroot)
    tbl.write(outfile, overwrite=True)
    print("Wrote: {}".format(outfile))


#####################################################################################################
#####################################################################################################
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

def poly_val(coeff, x, nrm):
    """
    IDL style function for polynomial

    Args:
        coeff (np.ndarray):  Polynomial coefficients
        x (np.ndarray):  x array
        nrm (np.ndarray): Normalization terms

    Returns:
        np.ndarray:  Same shape as x

    """
    #
    xnrm = 2. * (x - nrm[0])/nrm[1]
    #
    n = len(coeff)-1
    y = coeff[n]
    #for i=n-1,0,-1 do y = TEMPORARY(y) * x + c[i]
    for ii in range(n-1,-1,-1):
        y = y*xnrm + coeff[ii]
    return y


def xidl_arcspec(xidl_file, slit):
    xidl_dict = readsav(xidl_file)
    if xidl_dict['archive_arc'].ndim == 2:
        nspec = xidl_dict['archive_arc'].shape[0]
        npix = xidl_dict['archive_arc'].shape[1]
    else:
        npix = xidl_dict['archive_arc'].shape[0]
    # This is the best one (well-centered)
    calib = xidl_dict['calib'][slit]
    # Generate the wavelengths
    if calib['FUNC'] == b'CHEBY':
        wv_air = cheby_val(calib['FFIT'], np.arange(npix),
                       calib['NRM'], calib['NORD'])
    elif calib['FUNC'] == b'POLY':
        wv_air = poly_val(calib['FFIT'], np.arange(npix), calib['NRM'])

    wv_vac = airtovac(wv_air * units.AA)
    if xidl_dict['archive_arc'].ndim == 2:
        spec = xidl_dict['archive_arc'][slit]
    else:
        spec = xidl_dict['archive_arc']
    # Flip to blue to red?
    if wv_vac[1] < wv_vac[0]:
        wv_vac = wv_vac[::-1]
        spec = spec[::-1]
    # Return
    return wv_vac.value, spec

def main(flg):

    # Keck LRISb
    if flg & (2**0): # B300, all lamps
        binspec = 1
        slits = [15]
        xidl_file = os.path.join(template_path, 'Keck_LRIS', 'B300', 'lris_blue_300.sav')
        outroot = 'keck_lris_blue_300_d680.fits'
        build_template(xidl_file, slits, None, binspec, outroot, lowredux=True)

    if flg & (2**1): # B400, all lamps I think)
        binspec = 2
        outroot='keck_lris_blue_400_d560.fits'
        slits = [19,14]
        lcut = [5500.]
        xidl_file = os.path.join(template_path, 'Keck_LRIS', 'B400', 'lris_blue_400_d560.sav')
        build_template(xidl_file, slits, lcut, binspec, outroot, lowredux=True)

    if flg & (2**2): # B600, all lamps
        binspec = 2
        outroot='keck_lris_blue_600_d560.fits'
        slits = [0,7]
        lcut = [4500.]
        wfile = os.path.join(template_path, 'Keck_LRIS', 'B600', 'MasterWaveCalib_A_1_01.json')
        build_template(wfile, slits, lcut, binspec, outroot, lowredux=False)

    if flg & (2**3): # B1200, all lamps?
        binspec = 2
        outroot='keck_lris_blue_1200_d460.fits'
        slits = [19,44]
        lcut = [3700.]
        xidl_file = os.path.join(template_path, 'Keck_LRIS', 'B1200', 'lris_blue_1200.sav')
        build_template(xidl_file, slits, lcut, binspec, outroot, lowredux=True)

    # Shane Kastb
    if flg & (2**4):  # 452/3306
        binspec = 1
        slits = [0]
        xidl_file = os.path.join(template_path, 'Shane_Kast', '452_3306', 'kast_452_3306.sav')
        outroot = 'shane_kast_blue_452.fits'
        build_template(xidl_file, slits, None, binspec, outroot, lowredux=True)

    if flg & (2**5):  # 600/4310
        binspec = 1
        slits = [0,3]
        lcut = [4550.]
        xidl_file = os.path.join(template_path, 'Shane_Kast', '600_4310', 'kast_600_4310.sav')
        outroot = 'shane_kast_blue_600.fits'
        build_template(xidl_file, slits, lcut, binspec, outroot, lowredux=True)

    if flg & (2**6):  # 830/3460
        binspec = 1
        slits = [0]
        xidl_file = os.path.join(template_path, 'Shane_Kast', '830_3460', 'kast_830_3460.sav')
        outroot = 'shane_kast_blue_830.fits'
        build_template(xidl_file, slits, None, binspec, outroot, lowredux=True)

    # Keck/DEIMOS
    if flg & (2**7):  # 600ZD :: Might not go red enough
        binspec = 1
        slits = [0,1]
        lcut = [7192.]
        xidl_file = os.path.join(template_path, 'Keck_DEIMOS', '600ZD', 'deimos_600.sav')
        outroot = 'keck_deimos_600.fits'
        build_template(xidl_file, slits, lcut, binspec, outroot, lowredux=True)

    if flg & (2**8):  # 830G
        binspec = 1
        outroot='keck_deimos_830G.fits'
        # 3-12 = blue  6508 -- 8410
        # 7-24 = blue  8497 -- 9925 (no lines after XeI)
        ifiles = [0, 0, 1]
        slits = [12, 14, 24]
        lcut = [8400., 8480]
        wfile1 = os.path.join(template_path, 'Keck_DEIMOS', '830G_M_8600', 'MasterWaveCalib_A_1_03.json')
        wfile2 = os.path.join(template_path, 'Keck_DEIMOS', '830G_M_8600', 'MasterWaveCalib_A_1_07.json')
        # det_dict
        det_cut = {}
        det_cut['dets'] = [[1,2,3,4], [5,6,7,8]]
        det_cut['wcuts'] = [[0,9000.], [8200,1e9]]  # Significant overlap is fine
        #
        build_template([wfile1,wfile2], slits, lcut, binspec, outroot, lowredux=False,
                       ifiles=ifiles, det_cut=det_cut)

    if flg & (2**9):  # 1200
        binspec = 1
        outroot='keck_deimos_1200G.fits'
        # 3-3 = blue  6268.23 -- 7540
        # 3-14 = red   6508 -- 7730
        # 7-3 = blue  7589 -- 8821
        # 7-17 = red  8000 - 9230
        # 7c-0 = red  9120 -- 9950
        ifiles = [0, 0, 1, 1, 2]
        slits = [3, 14, 3, 17, 0]
        lcut = [7450., 7730., 8170, 9120]
        wfile1 = os.path.join(template_path, 'Keck_DEIMOS', '1200G', 'MasterWaveCalib_A_1_03.json')
        wfile2 = os.path.join(template_path, 'Keck_DEIMOS', '1200G', 'MasterWaveCalib_A_1_07.json')
        wfile3 = os.path.join(template_path, 'Keck_DEIMOS', '1200G', 'MasterWaveCalib_A_1_07c.json')
        # det_dict
        det_cut = None
        #det_cut = {}
        #det_cut['dets'] = [[1,2,3,4], [5,6,7,8]]
        #det_cut['wcuts'] = [[0,9000.], [8200,1e9]]  # Significant overlap is fine
        #
        build_template([wfile1,wfile2,wfile3], slits, lcut, binspec, outroot, lowredux=False,
                       ifiles=ifiles, det_cut=det_cut, chk=True)

    # Keck/LRISr
    if flg & (2**10): # R400
        binspec = 2
        outroot='keck_lris_red_400.fits'
        slits = [7]  # Quite blue, but not the bluest
        lcut = []
        wfile = os.path.join(template_path, 'Keck_LRIS', 'R400', 'MasterWaveCalib_A_1_01.json')
        build_template(wfile, slits, lcut, binspec, outroot, lowredux=False)

    if flg & (2**11):  # R1200
        # slits = [2-3]  # 7726 -- 9250
        # slits = [1-4]  # 9250 -- 9925
        binspec = 1
        outroot='keck_lris_red_1200_9000.fits'
        ifiles = [0, 1]
        slits = [3, 7]
        lcut = [9250.]
        wfile1 = os.path.join(template_path, 'Keck_LRIS', 'R1200_9000', 'MasterWaveCalib_A_1_02.json')  # Original Dev
        wfile2 = os.path.join(template_path, 'Keck_LRIS', 'R1200_9000', 'MasterWaveCalib_A_1_01.json')  # Dev suite 2x1
        build_template([wfile1,wfile2], slits, lcut, binspec, outroot, lowredux=False,
                       ifiles=ifiles)

    if flg & (2**12):  # R600/5000
        # slits = [1-4]  # 5080 -- 7820
        # slits = [1-7]  # 7820 -- 9170
        binspec = 2
        outroot='keck_lris_red_600_5000.fits'
        slits = [4, 7]
        lcut = [7820.]
        wfile = os.path.join(template_path, 'Keck_LRIS', 'R600_5000', 'MasterWaveCalib_B_1_01.json')
        build_template(wfile, slits, lcut, binspec, outroot, lowredux=False)

    if flg & (2**13):  # Magellan/MagE
        # Load
        mase_path = os.path.join(os.getenv('XIDL_DIR'), 'Magellan', 'MAGE', 'mase', 'Calib')
        sav_file = os.path.join(mase_path, 'MagE_wvguess_jfh.idl')
        mase_dict = readsav(sav_file)
        mase_sol = Table(mase_dict['all_arcfit'])
        # Do it
        all_wave = np.zeros((2048, 15))
        all_flux = np.zeros_like(all_wave)
        for order in np.arange(15):
            all_flux[:,order] = mase_dict['sv_aspec'][order]
            # Build the wavelengths
            wv_air = cheby_val(mase_sol['FFIT'][order], np.arange(2048), mase_sol['NRM'][order],
                                         mase_sol['NORD'][order])
            all_wave[:,order] = airtovac(wv_air * units.AA).value
        # Write
        tbl = Table()
        tbl['wave'] = all_wave.T
        tbl['flux'] = all_flux.T
        tbl['order'] = np.arange(20, 5, -1, dtype=int)
        tbl.meta['BINSPEC'] = 1
        # Write
        outroot='magellan_mage.fits'
        outfile = os.path.join(template_path, outroot)
        tbl.write(outfile, overwrite=True)
        print("Wrote: {}".format(outfile))

    if flg & (2**14):  # Magellan/MagE Plots
        outpath = os.path.join(resource_filename('pypeit', 'data'), 'arc_lines', 'plots')
        new_mage_file = os.path.join(resource_filename('pypeit', 'data'), 'arc_lines', 'reid_arxiv',
                                     'magellan_mage.fits')
        # Load
        mage_wave = Table.read(new_mage_file)
        llist = waveio.load_line_lists(['ThAr_MagE'])
        #
        for kk in range(mage_wave['wave'].shape[1]):
            wv = mage_wave['wave'][:, kk]
            fx = mage_wave['flux'][:, kk]
            order = 20 - kk
            # Reidentify
            detections, spec_cont_sub, patt_dict = autoid.reidentify(fx, fx, wv, llist, 1)
            # Fit
            final_fit = fitting.fit_slit(fx, patt_dict, detections, llist)
            # Output
            outfile=os.path.join(outpath, 'MagE_order{:2d}_IDs.pdf'.format(order))
            autoid.arc_fit_qa(final_fit, outfile=outfile, ids_only=True)
            print("Wrote: {}".format(outfile))
            autoid.arc_fit_qa(final_fit, outfile=os.path.join(outpath, 'MagE_order{:2d}_full.pdf'.format(order)))

    if flg & (2**15):  # VLT/X-Shooter reid_arxiv
        # VIS
        reid_path = os.path.join(resource_filename('pypeit', 'data'), 'arc_lines', 'reid_arxiv')
        old_file = os.path.join(reid_path, 'vlt_xshooter_vis1x1.json')
        # Load
        odict, par = waveio.load_reid_arxiv(old_file)

        # Do it
        orders = odict['fit2d']['orders'][::-1].astype(int)  # Flipped
        all_wave = np.zeros((odict['0']['nspec'], orders.size))
        all_flux = np.zeros_like(all_wave)
        for kk,order in enumerate(orders):
            all_flux[:,kk] = odict[str(kk)]['spec']
            all_wave[:,kk] = airtovac(odict[str(kk)]['wave_soln'] * units.AA).value
        # Write
        tbl = Table()
        tbl['wave'] = all_wave.T
        tbl['flux'] = all_flux.T
        tbl['order'] = orders
        tbl.meta['BINSPEC'] = 1
        # Write
        outroot='vlt_xshooter_vis1x1.fits'
        outfile = os.path.join(reid_path, outroot)
        tbl.write(outfile, overwrite=True)
        print("Wrote: {}".format(outfile))

    if flg & (2**16):  # VLT/X-Shooter line list
        line_path = os.path.join(resource_filename('pypeit', 'data'), 'arc_lines', 'lists')
        old_file = os.path.join(line_path, 'ThAr_XSHOOTER_VIS_air_lines.dat')
        # Load
        air_list = waveio.load_line_list(old_file)
        # Vacuum
        vac_wv = airtovac(air_list['wave']*units.AA).value
        vac_list = air_list.copy()
        vac_list['wave'] = vac_wv
        # Write
        new_file = os.path.join(line_path, 'ThAr_XSHOOTER_VIS_lines.dat')
        vac_list.write(new_file, format='ascii.fixed_width', overwrite=True)
        print("Wrote: {}".format(new_file))


# Command line execution
if __name__ == '__main__':
    flg = 0

    # Keck/LRISb
    #flg += 2**0  # LRISb 300, all lamps
    #flg += 2**1  # LRISb 400, all lamps
    #flg += 2**2  # LRISb 600, all lamps
    #flg += 2**3  # LRISb 1200, all lamps?

    # Shane/Kastb
    #flg += 2**4  # Kastb 452/3306 -- Not yet tested
    #flg += 2**5  # Kastb 600/4310
    #flg += 2**6  # Kastb 830/3460 -- Not yet tested

    # Keck/DEIMOS
    #flg += 2**7  # 600
    #flg += 2**8  # 830G
    #flg += 2**9  # 1200

    # Keck/LRISr
    #flg += 2**10  # R400
    #flg += 2**11  # R1200
    #flg += 2**12  # R600/5000

    # Shane/Kastr
    #  Need several arcs to proceed this way

    # MagE
    #flg += 2**13
    #flg += 2**14  # Plots

    # VLT/X-Shooter
    flg += 2**15  # Convert JSON to FITS
    #flg += 2**16  # Convert JSON to FITS

    main(flg)

