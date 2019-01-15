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


def load_specobjs(fname):
    """ Load a spec1d file into a list of SpecObjExp objects
    Parameters
    ----------
    fname : str

    Returns
    -------
    specObjs : list of SpecObjExp
    head0
    """
    sobjs = specobjs.SpecObjs()
    speckeys = ['WAVE', 'SKY', 'MASK', 'FLAM', 'FLAM_IVAR', 'FLAM_SIG', 'COUNTS_IVAR', 'COUNTS']
    # sobjs_keys gives correspondence between header cards and sobjs attribute name
    sobjs_key = specobjs.SpecObj.sobjs_key()
    hdulist = fits.open(fname)
    head0 = hdulist[0].header
    #pypeline = head0['PYPELINE']
    # Is this an Echelle reduction?
    #if 'Echelle' in pypeline:
    #    echelle = True
    #else:
    #    echelle = False

    for hdu in hdulist:
        if hdu.name == 'PRIMARY':
            continue
        # Parse name
        idx = hdu.name
        specobj = specobjs.SpecObj(None, None, None, idx = idx)
        # Assign specobj attributes from header cards
        for attr, hdrcard in sobjs_key.items():
            try:
                value = hdu.header[hdrcard]
            except:
                continue
            setattr(specobj, attr, value)
        # Load data
        spec = Table(hdu.data)
        shape = (len(spec), 1024)  # 2nd number is dummy
        specobj.shape = shape
        specobj.trace_spat = spec['TRACE']
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
        sobjs.add_sobj(specobj)

    # Return
    return sobjs, head0

def load_spec_order(fname,objid=None,order=None,extract='OPT',flux=True):
    """Loading single order spectrum from a PypeIt 1D specctrum fits file.
        it will be called by ech_load_spec
    Parameters:
        fname (str) : The file name of your spec1d file
        objid (str) : The id of the object you want to load. (default is the first object)
        order (int) : which order you want to load (default is None, loading all orders)
        extract (str) : 'OPT' or 'BOX'
        flux (bool) : default is True, loading fluxed spectra
    returns:
        spectrum_out : XSpectrum1D
    """
    if objid is None:
        objid = 0
    if order is None:
        msgs.error('Please specify which order you want to load')

    # read extension name into a list
    primary_header = fits.getheader(fname, 0)
    nspec = primary_header['NSPEC']
    extnames = [primary_header['EXT0001']] * nspec
    for kk in range(nspec):
        extnames[kk] = primary_header['EXT' + '{0:04}'.format(kk + 1)]
    extnameroot = extnames[0]

    # Figure out which extension is the required data
    ordername = '{0:04}'.format(order)
    extname = extnameroot.replace('OBJ0000', objid)
    extname = extname.replace('ORDER0000', 'ORDER' + ordername)
    try:
        exten = extnames.index(extname) + 1
        msgs.info("Loading extension {:s} of spectrum {:s}".format(extname, fname))
    except:
        msgs.error("Spectrum {:s} does not contain {:s} extension".format(fname, extname))

    spectrum = load.load_1dspec(fname, exten=exten, extract=extract, flux=flux)
    # Polish a bit -- Deal with NAN, inf, and *very* large values that will exceed
    #   the floating point precision of float32 for var which is sig**2 (i.e. 1e38)
    bad_flux = np.any([np.isnan(spectrum.flux), np.isinf(spectrum.flux),
                       np.abs(spectrum.flux) > 1e30,
                       spectrum.sig ** 2 > 1e10,
                       ], axis=0)
    # Sometimes Echelle spectra have zero wavelength
    bad_wave = spectrum.wavelength < 1000.0*units.AA
    bad_all = bad_flux + bad_wave
    ## trim bad part
    wave_out,flux_out,sig_out = spectrum.wavelength[~bad_all],spectrum.flux[~bad_all],spectrum.sig[~bad_all]
    spectrum_out = XSpectrum1D.from_tuple((wave_out,flux_out,sig_out), verbose=False)
    #if np.sum(bad_flux):
    #    msgs.warn("There are some bad flux values in this spectrum.  Will zero them out and mask them (not ideal)")
    #    spectrum.data['flux'][spectrum.select][bad_flux] = 0.
    #    spectrum.data['sig'][spectrum.select][bad_flux] = 0.

    return spectrum_out

# JFH This routine is deprecated.
def ech_load_spec(files,objid=None,order=None,extract='OPT',flux=True):
    """Loading Echelle spectra from a list of PypeIt 1D spectrum fits files
    Parameters:
        files (str) : The list of file names of your spec1d file
        objid (str) : The id (one per fits file) of the object you want to load. (default is the first object)
        order (int) : which order you want to load (default is None, loading all orders)
        extract (str) : 'OPT' or 'BOX'
        flux (bool) : default is True, loading fluxed spectra
    returns:
        spectrum_out : XSpectrum1D
    """

    nfiles = len(files)
    if objid is None:
        objid = ['OBJ0000'] * nfiles
    elif len(objid) == 1:
        objid = objid * nfiles
    elif len(objid) != nfiles:
        msgs.error('The length of objid should be either 1 or equal to the number of spectra files.')

    fname = files[0]
    ext_final = fits.getheader(fname, -1)
    norder = ext_final['ECHORDER'] + 1
    msgs.info('spectrum {:s} has {:d} orders'.format(fname, norder))
    if norder <= 1:
        msgs.error('The number of orders have to be greater than one for echelle. Longslit data?')

    # Load spectra
    spectra_list = []
    for ii, fname in enumerate(files):

        if order is None:
            msgs.info('Loading all orders into a gaint spectra')
            for iord in range(norder):
                spectrum = load_spec_order(fname,objid=objid[ii],order=iord,extract=extract,flux=flux)
                # Append
                spectra_list.append(spectrum)
        elif order >= norder:
            msgs.error('order number cannot greater than the total number of orders')
        else:
            spectrum = load_spec_order(fname, objid=objid[ii], order=order, extract=extract, flux=flux)
            # Append
            spectra_list.append(spectrum)
    # Join into one XSpectrum1D object
    spectra = collate(spectra_list)
    # Return
    return spectra

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

