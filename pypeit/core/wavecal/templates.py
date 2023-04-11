""" Module to generate templates for the PypeIt full_template wavelength calibration routine

.. include:: ../include/links.rst
"""

import os
import pathlib

import numpy as np
from IPython import embed

from matplotlib import pyplot as plt

from scipy.io import readsav
from scipy.interpolate import interp1d

from astropy.table import Table
from astropy import units

import linetools.utils

from pypeit import msgs
from pypeit import utils
from pypeit import io
from pypeit import wavecalib
from pypeit.core import arc
from pypeit.core.wave import airtovac
from pypeit.core.wavecal import waveio
from pypeit.core.wavecal import wvutils
from pypeit.core.wavecal import autoid
from pypeit.core.wavecal import wv_fitting
from pypeit.core import fitting
from pypeit import data

from pypeit.spectrographs.util import load_spectrograph

# Data Model
# FITS table
#  wave -- Wavelength values
#  flux -- Arc spectrum flux values
#
# Meta must include BINNING of the template with 1=native
if os.getenv('PYPEIT_DEV') is not None:
    template_path = (
        pathlib.Path(os.getenv('PYPEIT_DEV')) / 
        'dev_algorithms' / 'wavelengths' / 'template_files'
    )
else:
    # print("You may wish to set the PYPEIT_DEV environment variable")
    pass

# TODO: Rename the "Master" files here?


def build_template(in_files, slits, wv_cuts, binspec, outroot, outdir=None,
                   normalize=False, subtract_conti=False, wvspec=None,
                   lowredux=False, ifiles=None, det_cut=None, chk=False,
                   miny=None, overwrite=True, ascii_tbl=False, in_vac=True,
                   shift_wave=False, binning=None, micron=False):
    """
    Generate a full_template for a given instrument

    Args:
        in_files (list or str):
            Wavelength solution files, XIDL or PypeIt
            If PypeIt, they can be a mix of WaveCalib JSON and FITS files
        slits (list):
            Slits in the archive files to use
        wv_cuts (list):
            Wavelengths to cut each slit at. The elements of the list
            correspond to the wavelengths where two spectra are stitched together.
        binspec (int):
            Spectral binning of the archived spectrum
        outroot (str):
            Name of output archive
        outdir (str):
            Name of output directory
        lowredux (bool, optional):
            If true, in_files are from LowRedux
        wvspec (ndarray, optional):
            Manually input the wavelength values
        ifiles (list, optional):
            Ordering of the in_files.  Default is np.arange(len(in_files))
        det_cut (dict, optional):
            Cut the detector into pieces.  Important for long detectors with wavelengths on one side
        chk (bool, optional):
            Show a plot or two
        miny (float):
            Impose a minimum value
        normalize (bool, optional):
            If provided multiple in_files, normalize each
            snippet to have the same maximum amplitude.
        subtract_conti (bool, optional):
            Subtract the continuum for the final archive
        ascii_tbl (bool, optional):
            Table is a simple ASCII 2 column wave,flux table
        in_vac (bool, optional):
            True if input wavelengths are in vacuum
        shift_wave (bool, optional):
            Shift wavelengths when splicing to sync up precisely (Recommended)
            Requires PypeIt file (old JSON works for now)
        binning (list, optional):
            Allows for multiple binnings for input files
        micron (bool, optional):
            If True, assume the in_files archive has a solution in microns, so convert this to Angstroms
            which is the PypeIt convention. Default=False
    """
    if outdir is None:
        outdir = data.Paths.reid_arxiv
    if ifiles is None:
        ifiles = np.arange(len(in_files))
    if binning is None:
        binning = [None]*len(ifiles)
    # Load xidl file
    # Grab it
    # Load and splice
    yvals = []
    lvals = []
    if not isinstance(in_files, list):
        in_files = [in_files]
        ifiles = [0]*len(slits)
    # Loop on the files
    for kk, slit in enumerate(slits):
        # Load up
        if wvspec is None:
            in_file = in_files[ifiles[kk]]
            if lowredux:
                wv_vac, spec = xidl_arcspec(in_file, slit)
            elif ascii_tbl:
                wv_vac, spec = read_ascii(in_file, in_vac=in_vac)
            else:
                wv_vac, spec, pypeitFit = pypeit_arcspec(in_file, slit, binspec, binning[kk])
            if micron:
                wv_vac = 1e4*wv_vac
        else:
            wv_vac, spec = wvspec['wv_vac'], wvspec['spec']
        # Diagnostics
        msgs.info("wvmin, wvmax of {}: {}, {}".format(in_file, wv_vac.min(), wv_vac.max()))
        # Cut
        if len(slits) > 1:
            wvmin, wvmax = grab_wvlim(kk, wv_cuts, len(slits))
            # Default
            gdi = (wv_vac > wvmin) & (wv_vac < wvmax)
            if shift_wave:
                if len(lvals) > 0:
                    # Find pixel closet to end
                    ipix = np.argmin(np.abs(lvals[-1][-1] - wv_vac))
                    # Difference between the two spectra at this pixel
                    dwv_specs = wv_vac[ipix] - lvals[-1][-1]
                    # Delta wv per pix
                    dwv_snipp = wv_vac - np.roll(wv_vac, 1)
                    dwv_snipp[0] = dwv_snipp[1]
                    # Delta pix -- approximate but should be pretty good
                    dpix = dwv_specs / dwv_snipp[ipix]
                    # Calculate new wavelengths
                    npix = wv_vac.size
                    # Rebin spec?
                    if binning is not None and binning[kk] != binspec:
                        npix_orig = spec.size
                        x_orig = np.arange(npix_orig) / float(npix_orig - 1)
                        x = np.arange(npix) / float(npix - 1)
                        spec = (interp1d(x_orig, spec, axis=0,
                                         bounds_error=False, fill_value='extrapolate'))(x)
                    # Evaluate
                    new_wave = pypeitFit.eval((-dpix + np.arange(npix)) / (npix - 1))
                    # Range
                    iend = np.argmin(np.abs(new_wave - wvmax))
                    # Interpolate
                    f = interp1d(wv_vac, spec)
                    spec = f(new_wave[ipix + 1:iend])
                    wv_vac = new_wave[ipix+1:iend]
                    # Over-write gdi
                    gdi = np.ones_like(wv_vac, dtype=bool)
        else:
            gdi = np.arange(spec.size).astype(int)
        # Append
        yvals.append(spec[gdi])
        lvals.append(wv_vac[gdi])
    # Continuum
    if subtract_conti:
        for kk,spec in enumerate(yvals):
            _, _, _, _, spec_cont_sub = wvutils.arc_lines_from_spec(spec)
            yvals[kk] = spec_cont_sub
    # Normalize?
    if normalize:
        norm_val = 10000.
        # Max values
        maxs = []
        for kk,spec in enumerate(yvals):
            mx = np.max(spec)
            spec = spec * norm_val / mx
            yvals[kk] = spec
    # Concatenate
    nwspec = np.concatenate(yvals)
    nwwv = np.concatenate(lvals)
    # Min y?
    if miny is not None:
        nwspec = np.maximum(nwspec, miny)
    # Check
    if chk:
        plt.clf()
        ax = plt.gca()
        ax.plot(nwwv, nwspec)
        plt.show()
    # Generate the table
    wvutils.write_template(nwwv, nwspec, binspec, outdir, outroot, det_cut=det_cut, overwrite=overwrite)

def grab_wvlim(kk, wv_cuts, nslits):
    """
    Set the wavelength range to cut on

    Args:
        kk (int):
        wv_cuts (list):
        nslits (int):

    Returns:
        tuple: wv_min, wv_max (float, float)

    """
    if kk == 0:
        llow = 0.
        lhi = wv_cuts[0]
    elif kk == nslits - 1:
        llow = wv_cuts[kk - 1]
        lhi = 1e9
    else:
        llow = wv_cuts[kk - 1]
        lhi = wv_cuts[kk]
    #
    return llow, lhi


def pypeit_arcspec(in_file, slit, binspec, binning=None):
    """
    Load up the arc spectrum from an input JSON file

    Args:
        in_file (str):
            File containing the arc spectrum and or fit
        slit (int):
            slit index

    Returns:
        tuple: np.ndarray, np.ndarray, PypeItFit:  wave, flux, pypeitFitting

    """
    if '.json' in in_file:
        # Force any possible pathlib.Path object to string before `loadjson`
        wv_dict = linetools.utils.loadjson(str(in_file))
        iwv_calib = wv_dict[str(slit)]
        pypeitFitting = fitting.PypeItFit(fitc=np.array(iwv_calib['fitc']),
                                          func=iwv_calib['function'],
                                          minx=iwv_calib['fmin'], maxx=iwv_calib['fmax'])

        if binning is not None and binning != binspec:
            raise NotImplementedError('Cannot rebin yet')
        else:
            x = np.arange(len(iwv_calib['spec']))

        wv_vac = pypeitFitting.eval(x / iwv_calib['xnorm'])
        #wv_vac = utils.func_val(iwv_calib['fitc'], x/iwv_calib['xnorm'], iwv_calib['function'],
        #                   minx=iwv_calib['fmin'], maxx=iwv_calib['fmax'])
        flux = np.array(iwv_calib['spec']).flatten()
    elif '.fits' in in_file:
        wvcalib = wavecalib.WaveCalib.from_file(in_file, chk_version=False)
        idx = np.where(wvcalib.spat_ids == slit)[0][0]
        flux = wvcalib.arc_spectra[:,idx]
        #
        npix = flux.size
        if binning is not None and binning != binspec:
            npix = int(npix * binning / binspec)
            x = np.arange(npix) / (npix - 1)
        else:
            x = np.arange(npix) / (npix - 1)
        # Evaluate
        wv_vac = wvcalib.wv_fits[idx].pypeitfit.eval(x)
        pypeitFitting = wvcalib.wv_fits[idx].pypeitfit
    else:
        raise IOError("Bad in_file {}".format(in_file))

    # Return
    return wv_vac, flux, pypeitFitting


def pypeit_identify_record(iwv_calib, binspec, specname, gratname, dispangl, outdir=None):
    """From within PypeIt, generate a template file if the user manually identifies an arc spectrum

    Parameters
    ----------

    iwv_calib : dict
        Wavelength calibration returned by final_fit
    binspec : int
        Spectral binning
    specname : str
        Name of instrument
    gratname : str
        Name of grating
    dispangl : str
        Dispersion angle
    outdir : str, None
        Output directory
    """
    x = np.arange(len(iwv_calib['spec']))
    wv_vac = utils.func_val(iwv_calib['fitc'], x / iwv_calib['xnorm'], iwv_calib['function'],
                            minx=iwv_calib['fmin'], maxx=iwv_calib['fmax'])
    wvspec = dict(wv_vac=wv_vac, spec=np.array(iwv_calib['spec']))
    # Derive an output file
    cntr = 1
    extstr = ""
    while True:
        outroot = '{0:s}_{1:s}_{2:s}{3:s}.fits'.format(specname, gratname, dispangl, extstr)
        if (data.Paths.reid_arxiv / outroot).exists():
            extstr = "_{0:02d}".format(cntr)
        else:
            break
        cntr += 1
    slits = [0]
    lcut = [3200.]
    build_template("", slits, lcut, binspec, outroot, outdir=outdir, wvspec=wvspec, lowredux=False, overwrite=False)
    # Return
    return outroot

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


def read_ascii(tbl_file, in_vac=True):
    """

    The columns need to be wave, flux
    And the data should be monoonically increasing in wavelength

    Args:
        tbl_file (str):
            file of the table
        in_vac (bool, optional):
            If True, wavelenghts are already in vacuum

    Returns:
        tuple: np.ndarray, np.ndarray  of wavelength, flux

    """
    arc_spec = Table.read(tbl_file, format='ascii')

    # Wavelengths
    wv_vac = arc_spec['wave']
    if not in_vac:
        wv_vac = airtovac(wv_vac * units.AA)

    # Return
    return wv_vac.value, arc_spec['flux']


def xidl_arcspec(xidl_file, slit):
    """
    Read an XIDL format solution for a Longslit

    Note:  These are in air

    Parameters
    ----------
    xidl_file : str
       XIDL file
    slit : int
       The slit in question.

    Returns
    -------
    wave : np.ndarray
        Wavelengths in vacuum for that slit
    spec : np.ndarray
        Arc spectrum for that slit

    """
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


def xidl_hires(xidl_file, specbin=1):
    """
    Read an XIDL format solution for Keck/HIRES

    Note:  They used air

    Args:
        xidl_file (str):
            Keck/HIRES save file

    Returns:

    """
    xidl_dict = readsav(xidl_file)
    order_vec = xidl_dict['guess_ordr']
    norders = order_vec.size
    nspec = xidl_dict['sv_aspec'].shape[1]

    # Wavelengths
    wave = np.zeros((norders, specbin*nspec))
    spec = np.zeros((norders, specbin*nspec))

    calib = xidl_dict['all_arcfit']
    order_mask = np.ones(norders, dtype=bool)

    # Here we go on the fits
    for kk in range(norders):
        # Generate the wavelengths
        if calib['FUNC'][kk] == b'CHEBY':
            log10_wv_air = cheby_val(calib['FFIT'][kk], 
                               np.arange(nspec),
                        calib['NRM'][kk], calib['NORD'][kk])
        elif calib['FUNC'][kk] == b'POLY':
            log10_wv_air = poly_val(calib['FFIT'][kk], 
                              np.arange(nspec),
                              calib['NRM'][kk])
        else:
            order_mask[kk]=False
            continue

        wv_vac = airtovac(10**log10_wv_air * units.AA).value
        ispec = xidl_dict['sv_aspec'][kk,:]
        # Flip to blue to red?
        if wv_vac[1] < wv_vac[0]:
            wv_vac = wv_vac[::-1]
            ispec = ispec[::-1]
        # Fill
        if specbin != 1:
            wave[kk,:] = arc.resize_spec(wv_vac, nspec*specbin)
            spec[kk,:] = arc.resize_spec(ispec, nspec*specbin)
        else:
            wave[kk,:] = wv_vac
            spec[kk,:] = ispec
    # Return

    return order_vec[order_mask], wave[order_mask,:], spec[order_mask,:]


def main(flg):

    # Keck LRISb
    if flg & (2**0): # B300, all lamps
        binspec = 1
        slits = [15]
        xidl_file = template_path / 'Keck_LRIS' / 'B300' / 'lris_blue_300.sav'
        outroot = 'keck_lris_blue_300_d680.fits'
        build_template(xidl_file, slits, None, binspec, outroot, lowredux=True)

    if flg & (2**1): # B400, all lamps I think)
        binspec = 2
        outroot='keck_lris_blue_400_d560.fits'
        slits = [19,14]
        lcut = [5500.]
        xidl_file = template_path / 'Keck_LRIS' / 'B400' / 'lris_blue_400_d560.sav'
        build_template(xidl_file, slits, lcut, binspec, outroot, lowredux=True)

    if flg & (2**2): # B600, all lamps
        binspec = 2
        outroot='keck_lris_blue_600_d560.fits'
        slits = [0,7]
        lcut = [4500.]
        wfile = template_path / 'Keck_LRIS' / 'B600' / 'MasterWaveCalib_A_1_01.json'
        build_template(wfile, slits, lcut, binspec, outroot, lowredux=False)

    if flg & (2**3): # B1200, all lamps?
        binspec = 2
        outroot='keck_lris_blue_1200_d460.fits'
        slits = [19,44]
        lcut = [3700.]
        xidl_file = template_path / 'Keck_LRIS' / 'B1200' /'lris_blue_1200.sav'
        build_template(xidl_file, slits, lcut, binspec, outroot, lowredux=True)

    # ###############################################3
    # Keck/LRISr
    if flg & (2**10): # R400
        binspec = 2
        outroot='keck_lris_red_400.fits'
        slits = [7]  # Quite blue, but not the bluest
        lcut = []
        wfile = template_path / 'Keck_LRIS' / 'R400' / 'MasterWaveCalib_A_1_01.json'
        build_template(wfile, slits, lcut, binspec, outroot, lowredux=False)

    if flg & (2**11):  # R1200
        # slits = [2-3]  # 7726 -- 9250
        # slits = [1-4]  # 9250 -- 9925
        binspec = 1
        outroot='keck_lris_red_1200_9000.fits'
        ifiles = [0, 1]
        slits = [3, 7]
        lcut = [9250.]
        wfile1 = template_path / 'Keck_LRIS' / 'R1200_9000' / 'MasterWaveCalib_A_1_02.json'  # Original Dev
        wfile2 = template_path / 'Keck_LRIS' / 'R1200_9000' / 'MasterWaveCalib_A_1_01.json'  # Dev suite 2x1
        build_template([wfile1,wfile2], slits, lcut, binspec, outroot, lowredux=False,
                       ifiles=ifiles)

    if flg & (2**12):  # R600/5000
        # slits = [1-4]  # 5080 -- 7820
        # slits = [1-7]  # 7820 -- 9170
        binspec = 2
        outroot='keck_lris_red_600_5000.fits'
        slits = [4, 7]
        lcut = [7820.]
        wfile = template_path / 'Keck_LRIS' / 'R600_5000' /'MasterWaveCalib_B_1_01.json'
        build_template(wfile, slits, lcut, binspec, outroot, lowredux=False)

    if flg & (2**27):  # R600/7500
        # slits = [1-10]  # 5000 -- 7840
        # slits = [1-4]  # 7840 -- 9230
        binspec = 2
        outroot='keck_lris_red_600_7500.fits'
        slits = [10, 4]
        lcut = [7840.]
        wfile = template_path / 'Keck_LRIS' / 'R600_7500' / 'MasterWaveCalib_I_1_01.json'
        build_template(wfile, slits, lcut, binspec, outroot, lowredux=False,
                       chk=True, normalize=True, subtract_conti=True)

    # ##################################
    # Magellan/MagE
    if flg & (2**13):
        # Load
        mase_path = pathlib.Path(os.getenv('XIDL_DIR')) / 'Magellan' / 'MAGE' / 'mase' / 'Calib'
        sav_file = mase_path / 'MagE_wvguess_jfh.idl'
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
        outfile = template_path / outroot
        tbl.write(outfile, overwrite=True)
        print("Wrote: {}".format(outfile))

    if flg & (2**14):  # Magellan/MagE Plots
        new_mage_file = data.Paths.reid_arxiv / 'magellan_mage.fits'
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
            final_fit = wv_fitting.fit_slit(fx, patt_dict, detections, llist)
            # Output
            outfile=data.Paths.arc_plot / f'MagE_order{order:2d}_IDs.pdf'
            autoid.arc_fit_qa(final_fit, outfile=outfile, ids_only=True)
            print(f"Wrote: {outfile}")
            autoid.arc_fit_qa(
                final_fit,
                outfile=data.Paths.arc_plot / f'MagE_order{order:2d}_full.pdf'
            )

    if flg & (2**15):  # VLT/X-Shooter reid_arxiv
        # VIS
        for iroot, iout in zip(['vlt_xshooter_vis1x1.json', 'vlt_xshooter_nir.json'],
            ['vlt_xshooter_vis1x1.fits', 'vlt_xshooter_nir.fits']):
            # Load
            old_file = data.Paths.reid_arxiv / iroot
            odict, par = waveio.load_reid_arxiv(old_file)

            # Do it
            orders = odict['fit2d']['orders'][::-1].astype(int)  # Flipped
            all_wave = np.zeros((odict['0']['nspec'], orders.size))
            all_flux = np.zeros_like(all_wave)
            for kk,order in enumerate(orders):
                all_flux[:,kk] = odict[str(kk)]['spec']
                if 'nir' in iroot:
                    all_wave[:,kk] = odict[str(kk)]['wave_soln']
                else:
                    all_wave[:,kk] = airtovac(odict[str(kk)]['wave_soln'] * units.AA).value
            # Write
            tbl = Table()
            tbl['wave'] = all_wave.T
            tbl['flux'] = all_flux.T
            tbl['order'] = orders
            tbl.meta['BINSPEC'] = 1
            # Write
            outfile = data.Paths.reid_arxiv / iout
            tbl.write(outfile, overwrite=True)
            print("Wrote: {}".format(outfile))

    if flg & (2**16):  # VLT/X-Shooter line list
        old_file = data.get_linelist_filepath('ThAr_XSHOOTER_VIS_air_lines.dat')
        # Load
        air_list = waveio.load_line_list(old_file)
        # Vacuum
        vac_wv = airtovac(air_list['wave']*units.AA).value
        vac_list = air_list.copy()
        vac_list['wave'] = vac_wv
        # Write
        new_file = data.get_linelist_filepath('ThAr_XSHOOTER_VIS_lines.dat')
        vac_list.write(new_file, format='ascii.fixed_width', overwrite=True)
        print("Wrote: {}".format(new_file))

    if flg & (2**17):  # NIRES
        iroot = 'keck_nires.json'
        iout = 'keck_nires.fits'
        # Load
        old_file = data.Paths.reid_arxiv / iroot
        odict, par = waveio.load_reid_arxiv(old_file)

        # Do it
        orders = odict['fit2d']['orders'][::-1].astype(int)  # Flipped
        all_wave = np.zeros((odict['0']['nspec'], orders.size))
        all_flux = np.zeros_like(all_wave)
        for kk,order in enumerate(orders):
            all_flux[:,kk] = odict[str(kk)]['spec']
            if 'nir' in iroot:
                all_wave[:,kk] = odict[str(kk)]['wave_soln']
            else:
                all_wave[:,kk] = airtovac(odict[str(kk)]['wave_soln'] * units.AA).value
        # Write
        tbl = Table()
        tbl['wave'] = all_wave.T
        tbl['flux'] = all_flux.T
        tbl['order'] = orders
        tbl.meta['BINSPEC'] = 1
        # Write
        outfile = data.Paths.reid_arxiv / iout
        tbl.write(outfile, overwrite=True)
        print("Wrote: {}".format(outfile))


    if flg & (2**18):  # Gemini/GNIRS
        iroot = 'gemini_gnirs.json'
        iout = 'gemini_gnirs.fits'
        # Load
        old_file = data.Paths.reid_arxiv / iroot
        odict, par = waveio.load_reid_arxiv(old_file)

        # Do it
        orders = odict['fit2d']['orders'][::-1].astype(int)  # Flipped
        all_wave = np.zeros((odict['0']['nspec'], orders.size))
        all_flux = np.zeros_like(all_wave)
        for kk,order in enumerate(orders):
            all_flux[:,kk] = odict[str(kk)]['spec']
            if 'nir' in iroot:
                all_wave[:,kk] = odict[str(kk)]['wave_soln']
            else:
                all_wave[:,kk] = airtovac(odict[str(kk)]['wave_soln'] * units.AA).value
        # Write
        tbl = Table()
        tbl['wave'] = all_wave.T
        tbl['flux'] = all_flux.T
        tbl['order'] = orders
        tbl.meta['BINSPEC'] = 1
        # Write
        outfile = data.Paths.reid_arxiv / iout
        tbl.write(outfile, overwrite=True)
        print("Wrote: {}".format(outfile))


    if flg & (2**23):  # WHT/ISIS
        iroot = 'wht_isis_blue_1200_4800.json'
        outroot = 'wht_isis_blue_1200_4800.fits'
        wfile = template_path / 'WHT_ISIS' / '1200B' / iroot
        binspec = 2
        slits = [0]
        lcut = [3200.]
        build_template(wfile, slits, lcut, binspec, outroot, lowredux=False)

    if flg & (2**24):  # Magellan/FIRE
        iroot = 'magellan_fire_echelle.json'
        iout = 'magellan_fire_echelle.fits'
        # Load
        old_file = data.Paths.reid_arxiv / iroot
        odict, par = waveio.load_reid_arxiv(old_file)

        # Do it
        orders = odict['fit2d']['orders'][::-1].astype(int)  # Flipped
        all_wave = np.zeros((odict['0']['nspec'], orders.size))
        all_flux = np.zeros_like(all_wave)
        for kk,order in enumerate(orders):
            all_flux[:,kk] = odict[str(kk)]['spec']
            if 'nir' in iroot:
                all_wave[:,kk] = odict[str(kk)]['wave_soln']
            else:
                all_wave[:,kk] = airtovac(odict[str(kk)]['wave_soln'] * units.AA).value

        # Write
        tbl = Table()
        tbl['wave'] = all_wave.T
        tbl['flux'] = all_flux.T
        tbl['order'] = orders
        tbl.meta['BINSPEC'] = 1
        # Write
        outfile = data.Paths.reid_arxiv / iout
        tbl.write(outfile, overwrite=True)
        print("Wrote: {}".format(outfile))

    if flg & (2**25): # FIRE longslit
        binspec = 1
        outroot = 'magellan_fire_long.fits'
        xidl_file = pathlib.Path(os.getenv('FIRE_DIR')) / 'LowDispersion' / 'NeNeAr_archive_fit.fits'
        spec_file = pathlib.Path(os.getenv('FIRE_DIR')) / 'LowDispersion' / 'NeNeAr2.sav'
        fire_sol = Table.read(xidl_file)
        wave = cheby_val(fire_sol['FFIT'].data[0], np.arange(2048), fire_sol['NRM'].data[0], fire_sol['NORD'].data[0])
        wv_vac = airtovac(wave * units.AA)
        xidl_dict = readsav(spec_file)
        flux = xidl_dict['arc1d']
        wvutils.write_template(wv_vac.value, flux, binspec, data.Paths.reid_arxiv, outroot, det_cut=None)

    # Gemini/Flamingos2
    if flg & (2**26):
        iroot = ['Flamingos2_JH_JH.json','Flamingos2_HK_HK.json']
        outroot=['Flamingos2_JH_JH.fits','Flamingos2_HK_HK.fits']
        binspec = 1
        slits = [0]
        lcut = []
        for ii in range(len(iroot)):
            wfile = data.Paths.reid_arxiv / iroot[ii]
            build_template(wfile, slits, lcut, binspec, outroot[ii], lowredux=False)


    # MDM/OSMOS -- MDM4K
    if flg & (2 ** 28):
        # ArI 4159 -- 6800
        wfile = template_path / 'MDM_OSMOS' / 'MasterWaveCalib_MDM4K_01.json'
        outroot = 'mdm_osmos_mdm4k.fits'
        binspec = 1
        slits = [0]
        lcut = [3200.]
        build_template(wfile, slits, lcut, binspec, outroot, lowredux=False,
                       chk=True, subtract_conti=True)

    # MMT/MMIRS
    if flg & (2**32):
        iroot = ['mmt_mmirs_HK_zJ.json','mmt_mmirs_J_zJ.json','mmt_mmirs_K3000_Kspec.json']
        outroot=['mmt_mmirs_HK_zJ.fits','mmt_mmirs_J_zJ.fits','mmt_mmirs_K3000_Kspec.fits']
        binspec = 1
        slits = [1020,1020,1020]
        lcut = []
        for ii in range(len(iroot)):
            wfile = data.Paths.reid_arxiv / iroot[ii]
            build_template(wfile, slits, lcut, binspec, outroot[ii], lowredux=False)
    # LBT/MODS
    if flg & (2**33):
        iroot = ['lbt_mods1r_red.json','lbt_mods2r_red.json']
        outroot=['lbt_mods1r_red.fits','lbt_mods2r_red.fits']
        binspec = 1
        slits = [[1557],[1573]]
        lcut = []
        for ii in range(len(iroot)):
            wfile = data.Paths.reid_arxiv / iroot[ii]
            build_template(wfile, slits[ii], lcut, binspec, outroot[ii], lowredux=False)
    # P200 Triplespec
    if flg & (2**34):
        iroot = 'p200_triplespec_MasterWaveCalib.fits'
        iout = 'p200_triplespec.fits'
        # Load
        old_file = data.Paths.reid_arxiv / iroot
        par = io.fits_open(old_file)
        pyp_spec = par[0].header['PYP_SPEC']
        spectrograph  = load_spectrograph(pyp_spec)
        orders = spectrograph.orders

        # Do it
        all_wave = np.zeros((par[2].data['spec'].size, orders.size))
        all_flux = np.zeros_like(all_wave)
        for kk, order in enumerate(orders):
            all_flux[:, kk] = par[2*kk+2].data['spec']
            all_wave[:, kk] = par[2*kk+2].data['wave_soln']
        # Write
        tbl = Table()
        tbl['wave'] = all_wave.T
        tbl['flux'] = all_flux.T
        tbl['order'] = orders
        tbl.meta['BINSPEC'] = 1
        # Write
        outfile = data.Paths.reid_arxiv / iout
        tbl.write(outfile, overwrite=True)
        print("Wrote: {}".format(outfile))


# Command line execution
if __name__ == '__main__':
    flg = 0

    # TODO : There must be a better way to index these solutions...
    # Keck/LRISb
    #flg += 2**0  # LRISb 300, all lamps
    #flg += 2**1  # LRISb 400, all lamps
    #flg += 2**2  # LRISb 600, all lamps
    #flg += 2**3  # LRISb 1200, all lamps?

    # Keck/LRISr
    #flg += 2**10  # R400
    #flg += 2**11  # R1200

    # MagE
    #flg += 2**13
    #flg += 2**14  # Plots

    # VLT/X-Shooter
    #flg += 2**15  # Convert JSON to FITS
    #flg += 2**16  # Line list

    # Keck/NIRES
    #flg += 2**17  # Convert JSON to FITS

    # Gemini/GNIRS
    #flg += 2**18  # Convert JSON to FITS
    #flg += 2**19  # Convert JSON to FITS

    # Gemini/GMOS
    #flg += 2**20  # Hamamatsu R400 Convert JSON to FITS
    #flg += 2**21  # E2V Convert JSON to FITS
    #flg += 2**22  # Hamamatsu B600 XIDL

    # WHT/ISIS
    #flg += 2**23  # Convert JSON to FITS

    # Magellan/FIRE
    #flg += 2**24  # FIRE Echelle
    #flg += 2**25  # FIRElongslit

    # Gemini Flamingos2
    #flg += 2**26  # Longslit

    # Keck/LRIS r
    #flg += 2**27  # R600/7500

    # MDM/OSMMOS
    #flg += 2**28

    # MMT MMIRS
    #flg += 2**32

    # LBT MODS
    #flg += 2**33

    # P200 Triplespec
    #flg += 2**34

    main(flg)

