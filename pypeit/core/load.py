""" Module for loading PypeIt files
"""
import os
import warnings

import numpy as np

from astropy import units
from astropy.time import Time
from astropy.io import fits
from astropy.table import Table

from linetools.spectra.xspectrum1d import XSpectrum1D
from linetools.spectra.utils import collate
import linetools.utils

from pypeit import msgs
from IPython import embed
from pypeit.core import parse



# TODO I don't think we need this routine
def load_ext_to_array(hdulist, ext_id, ex_value='OPT', flux_value=True, nmaskedge=None):
    '''
    It will be called by load_1dspec_to_array.
    Load one-d spectra from ext_id in the hdulist

    Args:
        hdulist: FITS HDU list
        ext_id: extension name, i.e., 'SPAT1073-SLIT0001-DET03', 'OBJID0001-ORDER0003', 'OBJID0001-ORDER0002-DET01'
        ex_value: 'OPT' or 'BOX'
        flux_value: if True load fluxed data, else load unfluxed data

    Returns:
        tuple: Returns wave, flux, ivar, mask
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

    # Mask Edges
    if nmaskedge is not None:
        mask[:int(nmaskedge)] = False
        mask[-int(nmaskedge):] = False

    if flux_value:
        flux = hdu_iexp.data['{:}_FLAM'.format(ex_value)]
        ivar = hdu_iexp.data['{:}_FLAM_IVAR'.format(ex_value)]
    else:
        msgs.warn('Loading unfluxed spectra')
        flux = hdu_iexp.data['{:}_COUNTS'.format(ex_value)]
        ivar = hdu_iexp.data['{:}_COUNTS_IVAR'.format(ex_value)]

    return wave, flux, ivar, mask

# TODO merge this with unpack orders
def load_1dspec_to_array(fnames, gdobj=None, order=None, ex_value='OPT', flux_value=True, nmaskedge=None):
    '''
    Load the spectra from the 1d fits file into arrays.
    If Echelle, you need to specify which order you want to load.
    It can NOT load all orders for Echelle data.

    Args:
        fnames (list): 1D spectra fits file(s)
        gdobj (list): extension name (longslit/multislit) or objID (Echelle)
        order (None or int): order number
        ex_value (str): 'OPT' or 'BOX'
        flux_value (bool): if True it will load fluxed spectra, otherwise load counts

    Returns:
        tuple: Returns the following:
            - waves (ndarray): wavelength array of your spectra, see
              below for the shape information of this array.
            - fluxes (ndarray): flux array of your spectra
            - ivars (ndarray): ivars of your spectra
            - masks (ndarray, bool): mask array of your spectra

        The shapes of all returns are exactly the same.
            - Case 1: np.size(fnames)=np.size(gdobj)=1, order=None for
              Longslit or order=N (an int number) for Echelle
              Longslit/single order for a single fits file, they are 1D
              arrays with the size equal to Nspec
            - Case 2: np.size(fnames)=np.size(gdobj)>1, order=None for
              Longslit or order=N (an int number) for Echelle
              Longslit/single order for a list of fits files, 2D array,
              the shapes are Nspec by Nexp
            - Case 3: np.size(fnames)=np.size(gdobj)=1, order=None All
              Echelle orders for a single fits file, 2D array, the
              shapes are Nspec by Norders
            - Case 4: np.size(fnames)=np.size(gdobj)>1, order=None All
              Echelle orders for a list of fits files, 3D array, the
              shapres are Nspec by Norders by Nexp
    '''

    # read in the first fits file
    if isinstance(fnames, (list, np.ndarray)):
        nexp = np.size(fnames)
        fname0 = fnames[0]
    elif isinstance(fnames, str):
        nexp = 1
        fname0 = fnames

    hdulist = fits.open(fname0)
    header = hdulist[0].header
    npix = header['NPIX']
    pypeline = header['PYPELINE']

    # get the order/slit information
    ntrace0 = np.size(hdulist)-1
    idx_orders = []
    for ii in range(ntrace0):
        idx_orders.append(int(hdulist[ii+1].name.split('-')[1][5:])) # slit ID or order ID


    if pypeline == "Echelle":
        ## np.unique automatically sort the returned array which is not what I want!!!
        ## order_vec = np.unique(idx_orders)
        dum, order_vec_idx = np.unique(idx_orders, return_index=True)
        order_vec = np.array(idx_orders)[np.sort(order_vec_idx)]
        norder = np.size(order_vec)
    else:
        norder = 1

    #TODO This is unneccessarily complicated. The nexp=1 case does the same operations as the nexp > 1 case. Refactor
    # this so that it just does the same set of operations once and then reshapes the array at the end to give you what
    # you want. Let's merge this with unpack orders

    ## Loading data from a single fits file
    if nexp == 1:
        # initialize arrays
        if (order is None) and (pypeline == "Echelle"):
            waves = np.zeros((npix, norder,nexp))
            fluxes = np.zeros_like(waves)
            ivars = np.zeros_like(waves)
            masks = np.zeros_like(waves, dtype=bool)

            for ii, iord in enumerate(order_vec):
                ext_id = gdobj[0]+'-ORDER{:04d}'.format(iord)
                wave_iord, flux_iord, ivar_iord, mask_iord = load_ext_to_array(hdulist, ext_id, ex_value=ex_value,
                                                                               flux_value=flux_value, nmaskedge=nmaskedge)
                waves[:,ii,0] = wave_iord
                fluxes[:,ii,0] = flux_iord
                ivars[:,ii,0] = ivar_iord
                masks[:,ii,0] = mask_iord
        else:
            if pypeline == "Echelle":
                ext_id = gdobj[0]+'-ORDER{:04d}'.format(order)
            else:
                ext_id = gdobj[0]
            waves, fluxes, ivars, masks = load_ext_to_array(hdulist, ext_id, ex_value=ex_value, flux_value=flux_value,
                                                            nmaskedge=nmaskedge)

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

            # ToDo: The following part can be removed if all data are reduced using the leatest pipeline
            if pypeline == "Echelle":
                ntrace = np.size(hdulist_iexp) - 1
                idx_orders = []
                for ii in range(ntrace):
                    idx_orders.append(int(hdulist_iexp[ii + 1].name.split('-')[1][5:]))  # slit ID or order ID
                    dum, order_vec_idx = np.unique(idx_orders, return_index=True)
                    order_vec = np.array(idx_orders)[np.sort(order_vec_idx)]
            # ToDo: The above part can be removed if all data are reduced using the leatest pipeline

            if (order is None) and (pypeline == "Echelle"):
                for ii, iord in enumerate(order_vec):
                    ext_id = gdobj[iexp]+'-ORDER{:04d}'.format(iord)
                    wave_iord, flux_iord, ivar_iord, mask_iord = load_ext_to_array(hdulist_iexp, ext_id, ex_value=ex_value,
                                                                                   nmaskedge = nmaskedge, flux_value=flux_value)
                    waves[:,ii,iexp] = wave_iord
                    fluxes[:,ii,iexp] = flux_iord
                    ivars[:,ii,iexp] = ivar_iord
                    masks[:,ii,iexp] = mask_iord
            else:
                if pypeline == "Echelle":
                    ext_id = gdobj[iexp]+'-ORDER{:04d}'.format(order)
                else:
                    ext_id = gdobj[iexp]
                wave, flux, ivar, mask = load_ext_to_array(hdulist_iexp, ext_id, ex_value=ex_value, flux_value=flux_value,
                                                           nmaskedge=nmaskedge)
                waves[:, iexp] = wave
                fluxes[:, iexp] = flux
                ivars[:, iexp] = ivar
                masks[:, iexp] = mask

    return waves, fluxes, ivars, masks, header

def load_spec_order(fname,norder, objid=None, order=None, extract='OPT', flux=True):
    """
    Loading single order spectrum from a PypeIt 1D specctrum fits file.
    it will be called by ech_load_spec

    Args:
        fname (str) : The file name of your spec1d file
        objid (str) : The id of the object you want to load. (default is the first object)
        order (int) : which order you want to load (default is None, loading all orders)
        extract (str) : 'OPT' or 'BOX'
        flux (bool) : default is True, loading fluxed spectra

    Returns:
        XSpectrum1D: spectrum_out
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
    """
    Loading Echelle spectra from a list of PypeIt 1D spectrum fits files

    Args:
        files (str) : The list of file names of your spec1d file
        objid (str) : The id (one per fits file) of the object you want to load. (default is the first object)
        order (int) : which order you want to load (default is None, loading all orders)
        extract (str) : 'OPT' or 'BOX'
        flux (bool) : default is True, loading fluxed spectra

    Returns:
        XSpectrum1D: spectrum_out
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


def load_multiext_fits(filename, ext):
    """
    Load data and primary header from a multi-extension FITS file

    Args:
        filename (:obj:`str`):
            Name of the file.
        ext (:obj:`str`, :obj:`int`, :obj:`list`):
            One or more file extensions with data to return.  The
            extension can be designated by its 0-indexed integer
            number or its name.

    Returns:
        tuple: Returns the image data from each provided extension.
        If return_header is true, the primary header is also
        returned.
    """
    # Format the input and set the tuple for an empty return
    _ext = ext if isinstance(ext, list) else [ext]
    n_ext = len(_ext)
    # Open the file
    hdu = fits.open(filename)
    head0 = hdu[0].header
    # Only one extension
    if n_ext == 1:
        data = hdu[_ext[0]].data.astype(np.float)
        return data, head0
    # Multiple extensions
    data = tuple([None if hdu[k].data is None else hdu[k].data.astype(np.float) for k in _ext])
    # Return
    return data+(head0,)

