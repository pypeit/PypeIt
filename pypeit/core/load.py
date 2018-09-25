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


def get_utc(headarr):
    """
    Find and return the UTC for a file based on the headers read from
    all extensions.

    The value returned is the first UT or UTC keyword found any in any
    header object.

    Args:
        headarr (list):
            List of :obj:`astropy.io.fits.Header` objects to search
            through for a UTC of the observation.
    Returns:
        object: The value of the header keyword.
    """
    for h in headarr:
        if h == 'None':
            continue
        if 'UTC' in h.keys():
            return h['UTC']
        elif 'UT' in h.keys():
            return h['UT']
    return None


def convert_time(in_time, spectrograph):
    """
    Convert the time read from a file header to hours for all
    spectrographs.

    Args:
        in_time (str):
            The time read from the file header
        spectrograph
            (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            The spectrograph used to collect the data.  The class is
            used to set the time unit used in the file header.

    Returns:
        float: The time in hours.
    """
    # Convert seconds to hours
    if spectrograph.timeunit == 's':
        return float(in_time)/3600.0
    
    # Convert minutes to hours
    if spectrograph.timeunit == 'm':
        return float(in_time)/60.0

    # Convert from an astropy.Time format
    if spectrograph.timeunit in Time.FORMATS.keys():
        ival = float(in_time) if spectrograph.timeunit == 'mjd' else in_time
        tval = Time(ival, scale='tt', format=spectrograph.timeunit)
        # Put MJD in hours
        return tval.mjd * 24.0
        
    msgs.error('Bad time unit')


def create_fitstbl(file_list, spectrograph, strict=True):
    """
    Create a table of relevant fits file metadata used during the
    reduction.

    The content of the fits table is dictated by the header keywords
    specified for the provided spectrograph.  It is expected that this
    table can be used to set the frame type of each file.

    The metadata is validated using checks specified by the provided
    spectrograph class.

    .. todo::
        This should get abstracted to be independent of the
        spectrograph, with all the spectrograph dependent keywords be
        passed into the function.  The validation of the table would
        then happen in PypeItSetup.

    Args:
        file_list (list):
            The list of files to include in the table.
        spectrograph
            (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            The spectrograph used to collect the data save to each file.
            The class is used to provide the header keyword data to
            include in the table and specify any validation checks.
            
    Returns:
        :class:`astropy.table.Table`: A table with the relevant metadata
        for the provided fits files.
    """

    # Get the header keywords specific to the provided spectrograph.
    # head_keys is a nested dictionary listing the header keywords for
    # each extension in the file.  The top-level dictionary key is just
    # the 0-indexed number of the extension
    head_keys = spectrograph.header_keys()

    # The table is declared based on this input dictionary:
    # The directory, filename, instrument, utc are always included
    default_keys = [ 'directory', 'filename', 'instrume', 'utc' ]
    fitsdict = {k:[] for k in default_keys}

    # Add columns to the output table for each keyword.  The format is
    # extension.keyword.
    ext = {}
    for i in head_keys.keys():
        for k in head_keys[i].keys():
            if k in fitsdict.keys():
                raise ValueError('Keywords are not unique across all extensions!')
            ext[k] = i
            fitsdict[k] = []

    # TODO: Stopgap for required keys in fitstbl used by other parts of
    # the code.  Need to decide how to handle these.
    required_columns = ['time', 'date', 'target']
    required_for_ABBA = ['ra', 'dec'] # and for standards?
    added_by_fsort = ['frametype', 'framebit']
    if any([ c not in fitsdict.keys() for c in required_columns]):
        msgs.warn('Columns are missing.')

    # Number of files to read
    numfiles = len(file_list)

    # Loop on files
    for i in range(numfiles):

        # Read the fits headers
        headarr = spectrograph.get_headarr(file_list[i], strict=strict)

        # Check that the header is valid
        # TODO: The check_headers function currently always passes!
        # Needs to be implemented for each instrument.
        # spectrograph.check_headers() should raise an exception with an
        # appropriate message.
        try:
            # TODO: Move this into spectrograph.validate_fitstbl()
            spectrograph.check_headers(headarr)
        except Exception as e:
            msgs.warn('Reading of headers from file:' + msgs.newline() + file_list[i]
                      + msgs.newline() + 'failed with the following exception' + msgs.newline()
                      + e.__repr__() + msgs.newline() +
                      'Please check that the file was taken with the provided instrument:'
                      + msgs.newline() + '{0}'.format(spectrograph.spectrograph) + msgs.newline()
                      + 'Then either change the instrument or remove/skip the file.'
                      + msgs.newline()+ 'Continuing by ignoring this file...')
            numfiles -= 1
            continue

        # Add the directory, file name, and instrument to the table
        d,f = os.path.split(file_list[i])
        fitsdict['directory'].append(d)
        fitsdict['filename'].append(f)
        fitsdict['instrume'].append(spectrograph.spectrograph)

        # Add the time of the observation
        utc = get_utc(headarr)
        fitsdict['utc'].append('None' if utc is None else utc)
        if utc is None:
            msgs.warn("UTC is not listed as a header keyword in file:" + msgs.newline()
                      + file_list[i])

        # TODO: Read binning-dependent detector properties here? (maybe read speed too)

        # Now get the rest of the keywords
        for k in fitsdict.keys():
            if k in default_keys:
                continue

            # Try to read the header item
            try:
                value = headarr[ext[k]][head_keys[ext[k]][k]]
            except KeyError:
                # Keyword not found in header
                msgs.warn("{:s} keyword not in header. Setting to None".format(k))
                value = 'None'

            # Convert the time to hours
            # TODO: Done here or as a method in Spectrograph?
            if k == 'time' and value != 'None':
                value = convert_time(value, spectrograph)

            # Set the value
            vtype = type(value)
            if np.issubdtype(vtype, str):
                value = value.strip()
            if np.issubdtype(vtype, np.integer) or np.issubdtype(vtype, np.floating) \
                    or np.issubdtype(vtype, str) or np.issubdtype(vtype, np.bool_):
                fitsdict[k].append(value)
            else:
                msgs.bug("I didn't expect a useful header ({0:s}) to contain type {1:s}".format(k,
                             vtype).replace('<type ','').replace('>',''))

        msgs.info("Successfully loaded headers for file:" + msgs.newline() + file_list[i])

    # Report
    msgs.info("Headers loaded for {0:d} files successfully".format(numfiles))
    if numfiles != len(file_list):
        msgs.warn("Headers were not loaded for {0:d} files".format(len(file_list) - numfiles))
    if numfiles == 0:
        msgs.error("The headers could not be read from the input data files." + msgs.newline() +
                   "Please check that the settings file matches the data.")

    # Return an astropy.table.Table
    return Table(fitsdict)


def load_extraction(name, frametype='<None>', wave=True):
    msgs.info("Loading a pre-existing {0:s} extraction frame:".format(frametype)+msgs.newline()+name)
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
        objp = hdu.name.split('-')
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
            specobj = specobjs.SpecObj(shape, [float(objp[1][4:])/10000.]*2,
                                       int(objp[0][4:]),
                                       config='dummy_config',
                                       slitid=1, det=det,
                                       spat_pixpos=100)  # DUMMY
        except:
            msgs.error("BUG ME")
            debugger.set_trace()
        # Add trace
        specobj.trace = spec['TRACE']
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

