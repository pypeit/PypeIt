"""
Module for the SpecObjs and SpecObj classes

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import os
import re

import numpy as np

from astropy.units import Quantity
from astropy.io import fits

from pypeit import msgs
from pypeit.core import save
from pypeit import specobj
from pypeit.io import initialize_header
from pypeit.spectrographs.util import load_spectrograph

from IPython import embed

class SpecObjs(object):
    """
    Object to hold a set of SpecObj objects

    Note that this class overloads:

        - ``__getitem__`` to allow one to pull an attribute or a portion
          of the SpecObjs list
        - ``__setattr__`` to force a custom assignment method
        - ``__getattr__`` to generate an array of attribute 'k' from the
          specobjs.

    Args:
        specobjs (ndarray or list, optional):  One or more SpecObj objects

    Attributes:
        summary (astropy.table.Table):
    """
    @classmethod
    def from_fitsfile(cls, fits_file):
        """
        Instantiate from a FITS file

        Also tag on the Header

        Args:
            fits_file (str):

        Returns:
            specobsj.SpecObjs

        """
        # HDUList
        hdul = fits.open(fits_file)
        nhdu = len(hdul)
        # Init
        slf = cls()
        # Add on the header
        slf.header = hdul[0].header
        # Loop on em
        for kk in range(1,nhdu):
            tbl = fits.connect.read_table_fits(hdul, hdu=kk)
            sobj = specobj.SpecObj.from_table(tbl)
            slf.add_sobj(sobj)

        # JFH I'm commenting this out below. I prefer to just directly write out attributes and reinstantiate them
        # from files. Doing things like this just leads to errors
        # since it is not in touch with the code that actually determined what these attributes should be.

        # PYPELINE specific
        #if slf[0].PYPELINE == 'Echelle':
        #    # Set ech_objid
        #    uni_frac = np.unique(slf.ECH_FRACPOS)
        #    for ii, ufrac in enumerate(uni_frac):
        #        idx = np.isclose(slf.ECH_FRACPOS, ufrac)
        #        slf[idx].ECH_OBJID = ii
        #    # Set ech_orderindx
        #    uni_order = np.unique(slf.ECH_ORDER)
        #    for ii, uorder in enumerate(uni_order):
        #        idx = slf.ECH_ORDER == uorder
        #        slf[idx].ech_orderindx = ii

        # Return
        return slf

    def __init__(self, specobjs=None, header=None):

        # Only one attribute is allowed for this Object -- specobjs
        if specobjs is None:
            self.specobjs = np.array([])
        else:
            if isinstance(specobjs, (list, np.ndarray)):
                specobjs = np.array(specobjs)
            self.specobjs = specobjs

        self.header = header if header is not None else None

        # Turn off attributes from here
        #   Anything else set will be on the individual specobj objects in the specobjs array
        self.__initialised = True

    def __setattr__(self, item, value):
        """
        Define a custom assignment method.

        Args:
            item (str):
                Item to set
            value (object):
                Value of the item

        """
        if not '_SpecObjs__initialised' in self.__dict__:  # this test allows attributes to be set in the __init__ method
            return dict.__setattr__(self, item, value)
        elif item in self.__dict__:  # any normal attributes are handled normally
            dict.__setattr__(self, item, value)
        else:
            # Special handling when the input is an array/list and the length matches that of the slice
            if isinstance(value, (list, np.ndarray)):
                if len(self.specobjs) == len(value):  # Assume these are to be paired up
                    for kk, specobj in enumerate(self.specobjs):
                        setattr(specobj, item, value[kk])
                    return
            #
            for specobj in self.specobjs:
                setattr(specobj, item, value)

    @property
    def nobj(self):
        """
        Return the number of SpecObj objects

        Returns:
            int

        """
        return len(self.specobjs)

    def unpack_object(self, ret_flam=False):
        """
        Utility function to unpack the sobjs for one object and
        return various numpy arrays describing the spectrum and meta
        data. The user needs to already have trimmed the Specobjs to
        the relevant indices for the object.

        Args:
           ret_flam (:obj:`bool`, optional):
              If True return the FLAM, otherwise return COUNTS.

        Returns:
            tuple: Returns the following where all numpy arrays
            returned have shape (nspec, norders) for Echelle data and
            (nspec,) for Multislit data.

                - wave (`numpy.ndarray`_): Wavelength grids
                - flux (`numpy.ndarray`_): Flambda or counts
                - flux_ivar (`numpy.ndarray`_): Inverse variance (of
                  Flambda or counts)
                - flux_gpm (`numpy.ndarray`_): Good pixel mask.
                  True=Good
                - meta_spec (dict:) Dictionary containing meta data.
                  The keys are defined by
                  spectrograph.header_cards_from_spec()
                - header (astropy.io.header object): header from
                  spec1d file
    """

        # Read in the spec1d file
        norddet = self.nobj
        if ret_flam:
            # TODO Should nspec be an attribute of specobj?
            nspec = self[0].OPT_FLAM.size
        else:
            nspec = self[0].OPT_COUNTS.size
        # Allocate arrays and unpack spectrum
        wave = np.zeros((nspec, norddet))
        flux = np.zeros((nspec, norddet))
        flux_ivar = np.zeros((nspec, norddet))
        flux_gpm = np.zeros((nspec, norddet), dtype=bool)
        detector = np.zeros(norddet, dtype=int)
        ech_orders = np.zeros(norddet, dtype=int)
        # TODO make the extraction that is desired OPT vs BOX and optional input variable.
        for iorddet in range(norddet):
            wave[:, iorddet] = self[iorddet].OPT_WAVE
            flux_gpm[:, iorddet] = self[iorddet].OPT_MASK
            detector[iorddet] = self[iorddet].DET
            if self[0].PYPELINE == 'Echelle':
                ech_orders[iorddet] = self[iorddet].ECH_ORDER
            if ret_flam:
                flux[:, iorddet] = self[iorddet].OPT_FLAM
                flux_ivar[:, iorddet] = self[iorddet].OPT_FLAM_IVAR
            else:
                flux[:, iorddet] = self[iorddet].OPT_COUNTS
                flux_ivar[:, iorddet] = self[iorddet].OPT_COUNTS_IVAR

        # Populate meta data
        # TODO Remove this hack is it needed? If PYP_SPEC is always written then it is not.
        try:
            spectrograph = load_spectrograph(self.header['PYP_SPEC'])
        except:
            # TODO JFH  This is a hack until a generic spectrograph is implemented.
            spectrograph = load_spectrograph('shane_kast_blue')

        meta_spec = {}
        core_keys = spectrograph.header_cards_for_spec()
        for key in core_keys:
            try:
                meta_spec[key.upper()] = self.header[key.upper()]
            except KeyError:
                msgs.warn('Core meta data is missing from the specobjs header ')
                pass
        # Add the pyp spec.
        # TODO JFH: Make this an atribute of the specobj by default.
        meta_spec['PYP_SPEC'] = self.header['PYP_SPEC']
        meta_spec['PYPELINE'] = self[0].PYPELINE
        meta_spec['DET'] = detector
        if self[0].PYPELINE == 'MultiSlit' and self.nobj == 1:
            meta_spec['ECH_ORDERS'] = None
            return wave.reshape(nspec), flux.reshape(nspec), flux_ivar.reshape(nspec), \
                   flux_gpm.reshape(nspec), meta_spec, self.header
        else:
            meta_spec['ECH_ORDERS'] = ech_orders
            return wave, flux, flux_ivar, flux_gpm, meta_spec, self.header






    def get_std(self, multi_spec_det=None):
        """
        Return the standard star from this Specobjs. For MultiSlit this
        will be a single specobj in SpecObjs container, for Echelle it
        will be the standard for all the orders.

        Args:
            multi_spec_det (list):
                If there are multiple detectors arranged in the spectral direction, return the sobjs for
                the standard on each detector.

        Returns:
            SpecObj or SpecObjs

        """
        # Is this MultiSlit or Echelle
        pypeline = (self.PYPELINE)[0]
        if 'MultiSlit' in pypeline:
            # Have to do a loop to extract the counts for all objects
            SNR = np.median(self.OPT_COUNTS*np.sqrt(self.OPT_COUNTS_IVAR), axis=1)
            # For multiple detectors grab the requested detectors
            if multi_spec_det is not None:
                sobjs_std = SpecObjs(header=self.header)
                # Now append the maximum S/N object on each detector
                for idet in multi_spec_det:
                    this_det = self.DET == idet
                    istd = SNR[this_det].argmax()
                    sobjs_std.add_sobj(self[this_det][istd])
            else: # For normal multislit take the brightest object
                istd = SNR.argmax()
                # Return
                sobjs_std = SpecObjs(specobjs=[self[istd]], header=self.header)
            sobjs_std.header = self.header
            return sobjs_std
        elif 'Echelle' in pypeline:
            uni_objid = np.unique(self.ECH_FRACPOS)  # A little risky using floats
            uni_order = np.unique(self.ECH_ORDER)
            nobj = len(uni_objid)
            norders = len(uni_order)
            # Build up S/N
            SNR = np.zeros((norders, nobj))
            for iobj in range(nobj):
                for iord in range(norders):
                    ind = (self.ECH_FRACPOS == uni_objid[iobj]) & (self.ECH_ORDER == uni_order[iord])
                    spec = self[ind]
                    SNR[iord, iobj] = np.median(spec[0].OPT_COUNTS*np.sqrt(spec[0].OPT_COUNTS_IVAR))
            # Maximize S/N
            SNR_all = np.sqrt(np.sum(SNR**2,axis=0))
            objid_std = uni_objid[SNR_all.argmax()]
            # Finish
            indx = self.ECH_FRACPOS == objid_std
            # Return
            sobjs_std = SpecObjs(specobjs=self[indx], header=self.header)
            sobjs_std.header = self.header
            return sobjs_std
        else:
            msgs.error('Unknown pypeline')

    def append_neg(self, sobjs_neg):
        """
        Append negative objects and change the sign of their objids for IR reductions

        Args:
            sobjs_neg (SpecObjs):

        """
        if sobjs_neg.nobj == 0:
            msgs.warn("No negative objects found...")
            return
        # Assign the sign and the objids
        sobjs_neg.sign = -1.0
        if sobjs_neg[0].PYPELINE == 'Echelle':
            sobjs_neg.ECH_OBJID = -sobjs_neg.ECH_OBJID
            sobjs_neg.OBJID = -sobjs_neg.OBJID
        elif sobjs_neg[0].PYPELINE == 'MultiSlit':
            sobjs_neg.OBJID = -sobjs_neg.OBJID
        else:
            msgs.error("Should not get here")
        self.add_sobj(sobjs_neg)

        # Sort objects according to their spatial location. Necessary for the extraction to properly work
        if self.nobj > 0:
            self.specobjs = self.specobjs[np.argsort(self.SPAT_PIXPOS)]

    def purge_neg(self):
        """
        Purge negative objects from specobjs for IR reductions

        """
        # Assign the sign and the objids
        if self.nobj > 0:
            if self[0].PYPELINE == 'Echelle':
                index = self.ECH_OBJID < 0
            elif self[0].PYPELINE == 'MultiSlit':
                index = self.OBJID < 0
            else:
                msgs.error("Should not get here")
            self.remove_sobj(index)


    def slitorder_indices(self, slitorder):
        """
        Return the set of indices matching the input slit/order
        """
        if self[0].PYPELINE == 'Echelle':
            indx = self.ECH_ORDERINDX == slitorder
        elif self[0].PYPELINE == 'MultiSlit':
            indx = self.SLITID == slitorder
        else:
            msgs.error("Should not get here")
        #
        return indx

    def name_indices(self, name):
        """
        Return the set of indices matching the input slit/order

        Parameters
        ----------
        name: str
           The name of the object

        Returns
        -------
        indx: array-like, shape (nobj,)
           Array of indices with the corresponding name.

        """
        if self[0].PYPELINE == 'Echelle':
            indx = self.ECH_NAME == name
        elif self[0].PYPELINE == 'MultiSlit':
            indx = self.NAME == name
        else:
            msgs.error("Should not get here")
        return indx


    def slitorder_objid_indices(self, slitorder, objid):
        """
        Return the set of indices matching the input slit/order and the input objid
        """
        if self[0].PYPELINE == 'Echelle':
            indx = (self.ECH_ORDERINDX == slitorder) & (self.ECH_OBJID == objid)
        elif self[0].PYPELINE == 'MultiSlit':
            indx = (self.SLITID == slitorder) & (self.OBJID == objid)
        else:
            msgs.error("Should not get here")
        #
        return indx

    def set_names(self):
        for sobj in self.specobjs:
            sobj.set_name()

    def add_sobj(self, sobj):
        """
        Add one or more SpecObj
        The summary table is rebuilt

        Args:
            sobj (SpecObj or list or ndarray):  On or more SpecObj objects

        Returns:


        """
        if isinstance(sobj, specobj.SpecObj):
            self.specobjs = np.append(self.specobjs, [sobj])
        elif isinstance(sobj, (np.ndarray,list)):
            self.specobjs = np.append(self.specobjs, sobj)
        elif isinstance(sobj, SpecObjs):
            for isobj in sobj:
                self.specobjs = np.append(self.specobjs, isobj)

    def remove_sobj(self, index):
        """
        Remove an object

        Args:
            index: int

        Returns:

        """
        msk = np.ones(self.specobjs.size, dtype=bool)
        msk[index] = False
        # Do it
        self.specobjs = self.specobjs[msk]

    def copy(self):
        """
        Generate a copy of self

        Returns:
            SpecObjs

        """
        sobj_copy = SpecObjs(header=self.header)
        for sobj in self.specobjs:
            sobj_copy.add_sobj(sobj.copy())
        return sobj_copy

    def grab_spec_arrays(self, obj_id, DET=None, ECH_ORDER=None, **kwargs):
        # In development
        pass


#    JFH This function is deprecated
#    def set_idx(self):
#        """
#        Set the idx in all the SpecObj
#        Update the summary Table
#
#        Returns:
#
#        """
#        for sobj in self.specobjs:
#            sobj.set_idx()

    def __getitem__(self, item):
        """
        Overload to allow one to pull an attribute or a portion of the
        SpecObjs list

        Args:
            item (:obj:`str`, :obj:`int`, :obj:`slice`)

        Returns:
            The selected items as either an object,
            :class:`pypeit.specobj.SpecObj`, or
            :class:`pypeit.specobjs.SpecObjs`, depending on the input
            item.
        """
        if isinstance(item, str):
            return self.__getattr__(item)
        elif isinstance(item, (int, np.integer)):
            return self.specobjs[item]   # TODO Is this using pointers or creating new data????
        elif (isinstance(item, slice) or  # Stolen from astropy.table
            isinstance(item, np.ndarray) or
            isinstance(item, list) or
            isinstance(item, tuple) and all(isinstance(x, np.ndarray) for x in item)):
            # here for the many ways to give a slice; a tuple of ndarray
            # is produced by np.where, as in t[np.where(t['a'] > 2)]
            # For all, a new table is constructed with slice of all columns
            return SpecObjs(specobjs=self.specobjs[item], header=self.header)

    def __getattr__(self, k):
        """
        Overloaded to generate an array of attribute 'k' from the
        :class:`pypeit.specobj.SpecObj` objects.

        First attempts to grab data from the Summary table, then the list
        """
        if len(self.specobjs) == 0:
            raise ValueError("Empty specobjs")
        try:
            lst = [getattr(specobj, k) for specobj in self.specobjs]
        except ValueError:
            raise ValueError("Attribute does not exist")
        # Recast as an array
        return lst_to_array(lst)

    # Printing
    def __repr__(self):
        txt = '<{:s}:'.format(self.__class__.__name__)
        if self.nobj == 0:
            txt += "Empty SpecObjs"
        else:
            txt += '\n'
            for sobj in self.specobjs:
                txt += '  {} \n'.format(sobj)
        txt += '>'
        return txt

    def __len__(self):
        return len(self.specobjs)

    def build_header(self, head_fitstbl, spectrograph):
        """
        Builds the spec1d file header from the fitstbl meta data and meta data from spectrograph
        Args:
            head_fitstbl (header):
               Header from fitstbl
            spectrograph (object):
               Spectrograph object

        Returns:

        """
        header = fits.PrimaryHDU().header
        try:
            header['MJD-OBS'] = head_fitstbl['mjd']  # recorded as 'mjd' in fitstbl
        except KeyError:
            header['MJD-OBS'] = head_fitstbl['MJD-OBS']
        core_keys = spectrograph.header_cards_for_spec()
        for key in core_keys:
            # Allow for fitstbl vs. header
            try:
                header[key.upper()] = head_fitstbl[key.upper()]
            except KeyError:
                header[key.upper()] = head_fitstbl[key]
        # Specify which pipeline created this file
        header['PYPELINE'] = spectrograph.pypeline
        header['PYP_SPEC'] = (spectrograph.spectrograph, 'PypeIt: Spectrograph name')
        # Observatory
        telescope = spectrograph.telescope
        header['LON-OBS'] = telescope['longitude']
        header['LAT-OBS'] = telescope['latitude']
        header['ALT-OBS'] = telescope['elevation']
        _ = initialize_header(header)
        return header



    def write_to_fits(self, header, outfile, overwrite=True, update_det=None):
        """
        Write the set of SpecObj objects to one multi-extension FITS file

        Args:
            outfile (str):
            header:
            overwrite (bool, optional):
            update_det (int or list, optional):
              If provided, do not clobber the existing file but only update
              the indicated detectors.  Useful for re-running on a subset of detectors

        """
        if os.path.isfile(outfile) and (not overwrite):
            msgs.warn("Outfile exists.  Set overwrite=True to clobber it")
            return

        # If the file exists and update_det is provided, use the existing header
        #   and load up all the other hdus so that we only over-write the ones
        #   we are updating
        if os.path.isfile(outfile) and (update_det is not None):
            hdus, prihdu = save.init_hdus(update_det, outfile)
        else:
            # Build up the Header
            prihdu = fits.PrimaryHDU()
            hdus = [prihdu]
            prihdu.header = header

        ext = len(hdus)-1
        # Loop on the SpecObj objects
        for sobj in self.specobjs:
            if sobj is None:
                continue
            ext += 1
            # Add header keyword
            keywd = 'EXT{:04d}'.format(ext)
            prihdu.header[keywd] = sobj.NAME

            # Table
            shdu = fits.table_to_hdu(sobj._data)
            shdu.name = sobj.NAME
            # Append
            hdus += [shdu]

        # A few more for the header
        prihdu.header['NSPEC'] = len(hdus) - 1
        #prihdu.header['NPIX'] = specObjs.trace_spat.shape[1]
        # Code versions
        #_ = initialize_header(prihdu.header)

        # Finish
        hdulist = fits.HDUList(hdus)
        hdulist.writeto(outfile, overwrite=overwrite)
        msgs.info("Wrote 1D spectra to {:s}".format(outfile))
        return


def lst_to_array(lst, mask=None):
    """
    Simple method to convert a list to an array

    Allows for a list of Quantity objects

    Args:
        lst : list
          Should be number or Quantities
        mask (ndarray of bool, optional):  Limit to a subset of the list.  True=good

    Returns:
        ndarray or Quantity array:  Converted list

    """
    if mask is None:
        mask = np.array([True]*len(lst))
    if isinstance(lst[0], Quantity):
        return Quantity(lst)[mask]
    else:
        return np.array(lst)[mask]
