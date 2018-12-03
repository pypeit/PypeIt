""" Module for loading PypeIt files
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import os

import numpy as np

from astropy import units
from astropy.time import Time
from astropy.io import fits
from astropy.table import Table

from linetools.spectra.xspectrum1d import XSpectrum1D

from pypeit import msgs
from pypeit import specobjs
from pypeit import debugger
from pypeit.core import parse

def load_extraction(name, frametype='<None>', wave=True):
    msgs.info('Loading a pre-existing {0} extraction frame:'.format(frametype)
                + msgs.newline() + name)
    props_savas = dict({"ORDWN":"ordwnum"})
    props = dict({})
    props_allow = props_savas.keys()
    infile = fits.open(name)
    sciext = np.array(infile[0].data, dtype=np.float)
    sciwav = -999999.9*np.ones((sciext.shape[0],sciext.shape[1]))
    hdr = infile[0].header
    norders = hdr["NUMORDS"]
    pxsz    = hdr["PIXSIZE"]
    props = dict({})
    for o in range(norders):
        hdrname = "CDELT{0:03d}".format(o+1)
        cdelt = hdr[hdrname]
        hdrname = "CRVAL{0:03d}".format(o+1)
        crval = hdr[hdrname]
        hdrname = "CLINV{0:03d}".format(o+1)
        clinv = hdr[hdrname]
        hdrname = "CRPIX{0:03d}".format(o+1)
        crpix = hdr[hdrname]
        hdrname = "CNPIX{0:03d}".format(o+1)
        cnpix = hdr[hdrname]
        sciwav[:cnpix,o] = 10.0**(crval + cdelt*(np.arange(cnpix)-crpix))
        #sciwav[:cnpix,o] = clinv * 10.0**(cdelt*(np.arange(cnpix)-crpix))
        #sciwav[:cnpix,o] = clinv * (1.0 + pxsz/299792.458)**np.arange(cnpix)
    for k in props_allow:
        prsav = np.zeros(norders)
        try:
            for o in range(norders):
                hdrname = "{0:s}{1:03d}".format(k,o+1)
                prsav[o] = hdr[hdrname]
            props[props_savas[k]] = prsav.copy()
        except:
            pass
    del infile, hdr, prsav
    if wave is True:
        return sciext, sciwav, props
    else:
        return sciext, props


def load_ordloc(fname):
    # Load the files
    mstrace_bname, mstrace_bext = os.path.splitext(fname)
    lname = mstrace_bname+"_ltrace"+mstrace_bext
    rname = mstrace_bname+"_rtrace"+mstrace_bext
    # Load the order locations
    ltrace = np.array(fits.getdata(lname, 0),dtype=np.float)
    msgs.info("Loaded left order locations for frame:"+msgs.newline()+fname)
    rtrace = np.array(fits.getdata(rname, 0),dtype=np.float)
    msgs.info("Loaded right order locations for frame:"+msgs.newline()+fname)
    return ltrace, rtrace


def load_specobj(fname):
    """ Load a spec1d file into a list of SpecObjExp objects
    Parameters
    ----------
    fname : str

    Returns
    -------
    specObjs : list of SpecObjExp
    head0
    """
    speckeys = ['WAVE', 'SKY', 'MASK', 'FLAM', 'FLAM_IVAR', 'FLAM_SIG', 'COUNTS_IVAR', 'COUNTS']
    #
    specObjs = []
    hdulist = fits.open(fname)
    head0 = hdulist[0].header
    for hdu in hdulist:
        if hdu.name == 'PRIMARY':
            continue
        # Parse name
        idx = hdu.name
        objp = idx.split('-')
        if objp[-2][0:3] == 'DET':
            det = int(objp[-2][3:])
        else:
            det = int(objp[-2][1:])
        # Load data
        spec = Table(hdu.data)
        shape = (len(spec), 1024)  # 2nd number is dummy
        # Init
        #specobj = specobjs.SpecObj(shape, 'dum_config', int(objp[-1][1:]),
        #                           int(objp[-2][1:]), [float(objp[1][1:])/10000.]*2, 0.5,
        #                           float(objp[0][1:])/1000., 'unknown')
        # New and wrong
        try:
            specobj = specobjs.SpecObj(shape, None, None, idx = idx)
        except:
            debugger.set_trace()
            msgs.error("BUG ME")
        # TODO -- Figure out if this is a default
        # Add trace
        try:
            specobj.trace_spat = spec['TRACE']
        except:
            # KLUDGE!
            specobj.trace_spat = np.arange(len(spec['BOX_WAVE']))
        # Add spectrum
        if 'BOX_COUNTS' in spec.keys():
            for skey in speckeys:
                try:
                    specobj.boxcar[skey] = spec['BOX_{:s}'.format(skey)].data
                except KeyError:
                    pass
            # Add units on wave
            specobj.boxcar['WAVE'] = specobj.boxcar['WAVE'] * units.AA

        if 'OPT_COUNTS' in spec.keys():
            for skey in speckeys:
                try:
                    specobj.optimal[skey] = spec['OPT_{:s}'.format(skey)].data
                except KeyError:
                    pass
            # Add units on wave
            specobj.optimal['WAVE'] = specobj.optimal['WAVE'] * units.AA
        # Append
        specObjs.append(specobj)
    # Return
    return specObjs, head0


def load_tilts(fname):
    # Load the files
    msarc_bname, msarc_bext = os.path.splitext(fname)
    tname = msarc_bname+"_tilts"+msarc_bext
    sname = msarc_bname+"_satmask"+msarc_bext
    # Load the order locations
    tilts = np.array(fits.getdata(tname, 0),dtype=np.float)
    msgs.info("Loaded order tilts for frame:"+msgs.newline()+fname)
    satmask = np.array(fits.getdata(sname, 0),dtype=np.float)
    msgs.info("Loaded saturation mask for frame:"+msgs.newline()+fname)
    return tilts, satmask


def load_1dspec(fname, exten=None, extract='OPT', objname=None, flux=False):
    """
    Parameters
    ----------
    fname : str
      Name of the file
    exten : int, optional
      Extension of the spectrum
      If not given, all spectra in the file are loaded
    extract : str, optional
      Extraction type ('opt', 'box')
    objname : str, optional
      Identify extension based on input object name
    flux : bool, optional
      Return fluxed spectra?

    Returns
    -------
    spec : XSpectrum1D

    """
    # Keywords for Table
    rsp_kwargs = {}
    rsp_kwargs['wave_tag'] = '{:s}_WAVE'.format(extract)
    if flux:
        rsp_kwargs['flux_tag'] = '{:s}_FLAM'.format(extract)
        rsp_kwargs['sig_tag'] = '{:s}_FLAM_SIG'.format(extract)
    else:
        rsp_kwargs['flux_tag'] = '{:s}_COUNTS'.format(extract)
        rsp_kwargs['sig_tag'] = '{:s}_COUNTS_SIG'.format(extract)

    # Identify extension from objname?
    if objname is not None:
        hdulist = fits.open(fname)
        hdu_names = [hdu.name for hdu in hdulist]
        exten = hdu_names.index(objname)
        if exten < 0:
            msgs.error("Bad input object name: {:s}".format(objname))

    # Load
    spec = XSpectrum1D.from_file(fname, exten=exten, **rsp_kwargs)
    # Return
    return spec

def load_std_trace(spec1dfile, det):

    sdet = parse.get_dnum(det, prefix=False)
    hdulist_1d = fits.open(spec1dfile)
    det_nm = 'DET{:s}'.format(sdet)
    for hdu in hdulist_1d:
        if det_nm in hdu.name:
            tbl = Table(hdu.data)
            # TODO what is the data model for echelle standards? This routine needs to distinguish between echelle and longslit
            from IPython import embed
            embed()
            trace = tbl['TRACE']

    return trace


def waveids(fname):
    infile = fits.open(fname)
    pixels=[]
    msgs.info("Loading fitted arc lines")
    try:
        o = 1
        while True:
            pixels.append(infile[o].data.astype(np.float))
            o+=1
    except:
        pass
    return pixels

