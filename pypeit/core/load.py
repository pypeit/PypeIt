""" Module for loading PypeIt files
"""
import os

import numpy as np

from astropy import units
from astropy.time import Time
from astropy.io import fits
from astropy.table import Table

from linetools.spectra.xspectrum1d import XSpectrum1D
from linetools.spectra.utils import collate
import linetools.utils


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


def load_specobjs(fname,order=None):
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
    speckeys = ['WAVE', 'WAVE_GRID_MASK', 'WAVE_GRID','WAVE_GRID_MIN','WAVE_GRID_MAX', 'SKY', 'MASK', 'FLAM', 'FLAM_IVAR', 'FLAM_SIG',
                'COUNTS_IVAR', 'COUNTS', 'COUNTS_SIG']
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
        objp = idx.split('-')
        if objp[-2][:5] == 'ORDER':
            iord = int(objp[-2][5:])
            msgs.info('Loading Echelle data.')
        else:
            msgs.info('Loading longslit data.')
            iord = int(-1)
        if (order is not None) and (iord !=order):
            continue
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


def load_ext_to_array(hdulist, ext_id, ex_value='OPT', flux_value=True):
    '''
    It will be called by load_1dspec_to_array.
    Load one-d spectra from ext_id in the hdulist
    Args:
        hdulist: FITS HDU list
        ext_id: extension name, i.e., 'SPAT1073-SLIT0001-DET03', 'OBJID0001-ORDER0003', 'OBJID0001-ORDER0002-DET01'
        ex_value: 'OPT' or 'BOX'
        flux_value: if True load fluxed data, else load unfluxed data
    Returns:
         wave, flux, ivar, mask
    '''

    if (ex_value != 'OPT') and (ex_value != 'BOX'):
        msgs.error('{:} is not recognized. Please change to either BOX or OPT.'.format(ex_value))

    # get the order/slit information
    ntrace0 = np.size(hdulist)-1
    idx_names = []
    for ii in range(ntrace0):
        idx_names.append(hdulist[ii+1].name) # idx name

    # Initialize ext
    ext = None
    for indx in (idx_names):
        if ext_id in indx:
            ext = indx
    if ext is None:
        msgs.error('Can not find extension {:}.'.format(ext_id))
    else:
        hdu_iexp = hdulist[ext]

    wave = hdu_iexp.data['{:}_WAVE'.format(ex_value)]
    mask = hdu_iexp.data['{:}_MASK'.format(ex_value)]
    if flux_value:
        flux = hdu_iexp.data['{:}_FLAM'.format(ex_value)]
        ivar = hdu_iexp.data['{:}_FLAM_IVAR'.format(ex_value)]
    else:
        msgs.warn('Loading unfluxed spectra')
        flux = hdu_iexp.data['{:}_COUNTS'.format(ex_value)]
        ivar = hdu_iexp.data['{:}_COUNTS_IVAR'.format(ex_value)]

    return wave, flux, ivar, mask

def load_1dspec_to_array(fnames, gdobj=None, order=None, ex_value='OPT', flux_value=True):
    '''
    Load the spectra from the 1d fits file into arrays.
    If Echelle, you need to specify which order you want to load.
    It can NOT load all orders for Echelle data.
    Args:
        fnames: 1D spectra fits file(s)
        gdobj: extension name (longslit/multislit) or objID (Echelle)
        order: None or an int number
        ex_value: 'OPT' or 'BOX'
        flux_value: True or False
    Returns:
        waves, fluxes, ivars, masks
        The shapes of all returns are exactly the same.
        Case 1: np.size(fnames)=np.size(gdobj)=1, order=None for Longslit or order=N (an int number) for Echelle
            Longslit/single order for a single fits file, they are 1D arrays with the size equal to Nspec
        Case 2: np.size(fnames)=np.size(gdobj)>1, order=None for Longslit or order=N (an int number) for Echelle
            Longslit/single order for a list of fits files, 2D array, the shapes are Nspec by Nexp
        Case 3: np.size(fnames)=np.size(gdobj)=1, order=None
            All Echelle orders for a single fits file, 2D array, the shapes are Nspec by Norders
        Case 4: np.size(fnames)=np.size(gdobj)>1, order=None
            All Echelle orders for a list of fits files, 3D array, the shapres are Nspec by Norders by Nexp
    '''

    # read in the first fits file
    nexp = np.size(fnames)
    if nexp == 1:
        fname0 = fnames
    else:
        fname0 = fnames[0]
    hdulist = fits.open(fname0)
    npix = hdulist[0].header['NPIX']
    pypeline = hdulist[0].header['PYPELINE']

    # get the order/slit information
    ntrace0 = np.size(hdulist)-1
    idx_names = []
    idx_objids = []
    idx_orders = []
    for ii in range(ntrace0):
        idx_names.append(hdulist[ii+1].name) # idx name
        idx_objids.append(hdulist[ii+1].name.split('-')[0])
        idx_orders.append(int(hdulist[ii+1].name.split('-')[1][5:])) # slit ID or order ID

    if pypeline == "Echelle":
        order_vec = np.unique(idx_orders)
        norder = np.size(order_vec)
    else:
        order_vec = None
        norder = 1

    ## Loading data from a single fits file
    if nexp == 1:
        # initialize arrays
        if (order is None) and (pypeline == "Echelle"):
            waves = np.zeros((npix, norder))
            fluxes = np.zeros_like(waves)
            ivars = np.zeros_like(waves)
            masks = np.zeros_like(waves, dtype=bool)

            for ii, iord in enumerate(order_vec):
                ext_id = gdobj+'-ORDER{:04d}'.format(iord)
                wave_iord, flux_iord, ivar_iord, mask_iord = load_ext_to_array(hdulist, ext_id, ex_value=ex_value,
                                                                               flux_value=flux_value)
                waves[:,ii] = wave_iord
                fluxes[:,ii] = flux_iord
                ivars[:,ii] = ivar_iord
                masks[:,ii] = mask_iord
        else:
            if pypeline == "Echelle":
                ext_id = gdobj+'-ORDER{:04d}'.format(order)
            else:
                ext_id = gdobj
            waves, fluxes, ivars, masks = load_ext_to_array(hdulist, ext_id, ex_value=ex_value, flux_value=flux_value)

    ## Loading data from a list of fits files
    else:
        # initialize arrays
        if (order is None) and (pypeline == "Echelle"):
            # store all orders into one single array
            waves = np.zeros((npix, norder, nexp))
        else:
            # store a specific order or longslit
            waves = np.zeros((npix, nexp))
        fluxes = np.zeros_like(waves)
        ivars = np.zeros_like(waves)
        masks = np.zeros_like(waves,dtype=bool)

        for iexp in range(nexp):
            hdulist_iexp = fits.open(fnames[iexp])
            if (order is None) and (pypeline == "Echelle"):
                for ii, iord in enumerate(order_vec):
                    ext_id = gdobj[iexp]+'-ORDER{:04d}'.format(iord)
                    wave_iord, flux_iord, ivar_iord, mask_iord = load_ext_to_array(hdulist_iexp, ext_id, ex_value=ex_value,
                                                                                   flux_value=flux_value)
                    waves[:,ii,iexp] = wave_iord
                    fluxes[:,ii,iexp] = flux_iord
                    ivars[:,ii,iexp] = ivar_iord
                    masks[:,ii,iexp] = mask_iord
            else:
                if pypeline == "Echelle":
                    ext_id = gdobj[iexp]+'-ORDER{:04d}'.format(order)
                else:
                    ext_id = gdobj[iexp]
                wave, flux, ivar, mask = load_ext_to_array(hdulist_iexp, ext_id, ex_value=ex_value, flux_value=flux_value)
                waves[:, iexp] = wave
                fluxes[:, iexp] = flux
                ivars[:, iexp] = ivar
                masks[:, iexp] = mask

    return waves,fluxes,ivars,masks

def load_spec_order(fname,norder, objid=None,order=None,extract='OPT',flux=True):
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

    # Figure out which extension is the required data
    #extnameroot = extnames[0]
    #ordername = '{0:04}'.format(order)
    #extname = extnameroot.replace('OBJ0001', objid)
    #extname = extname.replace('ORDER0000', 'ORDER' + ordername)
    extnames_array = np.reshape(np.array(extnames),(norder,int(nspec/norder)))
    extnames_good = extnames_array[:,int(objid[3:])-1]
    extname = extnames_good[order]

    try:
        exten = extnames.index(extname) + 1
        msgs.info("Loading extension {:s} of spectrum {:s}".format(extname, fname))
    except:
        msgs.error("Spectrum {:s} does not contain {:s} extension".format(fname, extname))

    spectrum = load_1dspec(fname, exten=exten, extract=extract, flux=flux)
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
    ext_first = fits.getheader(fname, 1)
    ext_final = fits.getheader(fname, -1)
    norder = abs(ext_final['ECHORDER'] - ext_first['ECHORDER']) + 1
    msgs.info('spectrum {:s} has {:d} orders'.format(fname, norder))
    if norder <= 1:
        msgs.error('The number of orders have to be greater than one for echelle. Longslit data?')

    # Load spectra
    spectra_list = []
    for ii, fname in enumerate(files):

        if order is None:
            msgs.info('Loading all orders into a gaint spectra')
            for iord in range(norder):
                spectrum = load_spec_order(fname, norder, objid=objid[ii],order=iord,extract=extract,flux=flux)
                # Append
                spectra_list.append(spectrum)
        elif order >= norder:
            msgs.error('order number cannot greater than the total number of orders')
        else:
            spectrum = load_spec_order(fname,norder, objid=objid[ii], order=order, extract=extract, flux=flux)
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

    # Identify extension from objname?
    if objname is not None:
        hdulist = fits.open(fname)
        hdu_names = [hdu.name for hdu in hdulist]
        exten = hdu_names.index(objname)
        if exten < 0:
            msgs.error("Bad input object name: {:s}".format(objname))

    # Keywords for Table
    rsp_kwargs = {}
    if flux:
        rsp_kwargs['flux_tag'] = '{:s}_FLAM'.format(extract)
        rsp_kwargs['sig_tag'] = '{:s}_FLAM_SIG'.format(extract)
    else:
        rsp_kwargs['flux_tag'] = '{:s}_COUNTS'.format(extract)
        rsp_kwargs['sig_tag'] = '{:s}_COUNTS_SIG'.format(extract)

    # Use the WAVE_GRID (for 2d coadds) if it exists, otherwise use WAVE
    rsp_kwargs['wave_tag'] = '{:s}_WAVE_GRID'.format(extract)
    # Load
    try:
        spec = XSpectrum1D.from_file(fname, exten=exten, **rsp_kwargs)
    except ValueError:
        rsp_kwargs['wave_tag'] = '{:s}_WAVE'.format(extract)
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
            trace = tbl['TRACE']

    return trace



def load_sens_dict(filename):
    """
    Load a full (all slit) wv_calib dict

    Includes converting the JSON lists of particular items into ndarray

    Fills self.wv_calib and self.par

    Args:
        filename (str): Master file

    Returns:
        dict or None: self.wv_calib

    """


    # Does the master file exist?
    if not os.path.isfile(filename):
        msgs.warn("No sensfunc file found with filename {:s}".format(filename))
        return None
    else:
        msgs.info("Loading sensfunc from file {:s}".format(filename))
        sens_dict = linetools.utils.loadjson(filename)
        # Recast a few items as arrays
        for key in sens_dict.keys():
            try:
                int(key)
            except ValueError:
                continue
            else:
                for tkey in sens_dict[key].keys():
                    if isinstance(sens_dict[key][tkey], list):
                        sens_dict[key][tkey] = np.array(sens_dict[key][tkey])

        return sens_dict



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

