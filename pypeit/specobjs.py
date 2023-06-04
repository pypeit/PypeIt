"""
Module for the SpecObjs and SpecObj classes

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import os
import re
from typing import List

from IPython import embed

import numpy as np

from astropy import units
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time

from pypeit import msgs
from pypeit import specobj
from pypeit import io
from pypeit.spectrographs.util import load_spectrograph
from pypeit.core import parse
from pypeit.images.detector_container import DetectorContainer
from pypeit.images.mosaic import Mosaic
from pypeit import slittrace
from pypeit import utils

class SpecObjs:
    """
    Object to hold a set of :class:`~pypeit.specobj.SpecObj` objects

    Note that this class overloads:

        - ``__getitem__`` to allow one to pull an attribute or a portion
          of the SpecObjs list
        - ``__setattr__`` to force a custom assignment method
        - ``__getattr__`` to generate an array of attribute 'k' from the
          specobjs.

    Args:
        specobjs (`numpy.ndarray`_, list, optional):
            One or more :class:`~pypeit.specobj.SpecObj`  objects

    Attributes:
        summary (astropy.table.Table):
    """
    version = '1.0.0'

    #TODO JFH This method only populates some of the underlying specobj attributes, for example RA and DEC are not
    # getting set. We should do our best to populate everything that got written out to the file. This is part of having
    # a rigid data model.
    @classmethod
    def from_fitsfile(cls, fits_file, det=None, chk_version=True):
        """
        Instantiate from a spec1d FITS file

        Also tag on the Header

        Args:
            fits_file (:obj:`str`):
                The name of the fits file to read.
            det (:obj:`str`, optional):
                The string identifier for a detector or mosaic used to select
                the loaded spectra.  If None, all spectra are loaded.
            chk_version (:obj:`bool`, optional):
                If False, allow a mismatch in datamodel to proceed

        Returns:
            :class:`~pypeit.specobjs.SpecObjs`: The loaded spectra from the
            provided fits file.
        """
        # HDUList
        hdul = io.fits_open(fits_file)
        # Init
        slf = cls()
        # Add on the header
        slf.header = hdul[0].header
        # Keep track of HDUList for closing later

        # Catch common error of trying to read a OneSpec file
        if 'DMODCLS' in hdul[1].header and hdul[1].header['DMODCLS'] == 'OneSpec':
            msgs.error('This is a OneSpec file.  You are treating it like a SpecObjs file.')

        detector_hdus = {}
        # Loop for Detectors first as we need to add these to the objects
        for hdu in hdul[1:]:
            if 'DETECTOR' not in hdu.name:
                continue
            if 'DMODCLS' not in hdu.header:
                msgs.error('HDUs with DETECTOR in the name must have DMODCLS in their header.')
            try:
                dmodcls = eval(hdu.header['DMODCLS'])
            except:
                msgs.error(f"Unknown detector type datamodel class: {hdu.header['DMODCLS']}")
            # NOTE: This requires that any "detector" datamodel class has a
            # from_hdu method, and the name of the HDU must have a known format
            # (e.g., 'DET01-DETECTOR').
            _det = hdu.name.split('-')[0]
            detector_hdus[_det] = dmodcls.from_hdu(hdu)

        # Now the objects
        for hdu in hdul[1:]:
            if 'DETECTOR' in hdu.name:
                continue
            sobj = specobj.SpecObj.from_hdu(hdu, chk_version=chk_version)
            # Restrict on det?
            if det is not None and sobj.DET != det:
                continue
            # Check for detector
            if sobj.DET in detector_hdus.keys():
                sobj.DETECTOR = detector_hdus[sobj.DET]
            # Append
            slf.add_sobj(sobj)
        # Return
        hdul.close()
        return slf


    def __init__(self, specobjs=None, header=None):

        # Only two attributes are allowed for this Object -- specobjs, header
        if specobjs is None:
            self.specobjs = np.array([])
        else:
            if isinstance(specobjs, (list, np.ndarray)):
                specobjs = np.array(specobjs)
            self.specobjs = specobjs

        self.header = header
        self.hdul = None

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

    def unpack_object(self, ret_flam=False, extract_type='OPT'):
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
                  spectrograph.parse_spec_header()
                - header (astropy.io.header object): header from
                  spec1d file
        """
        # Prep
        norddet = self.nobj
        flux_attr = 'FLAM' if ret_flam else 'COUNTS'
        flux_key = '{}_{}'.format(extract_type, flux_attr)
        wave_key = '{}_WAVE'.format(extract_type)
        # Test
        if getattr(self, flux_key)[0] is None:
            msgs.error("Flux not available for {}.  Try the other ".format(flux_key))
        #
        nspec = getattr(self, flux_key)[0].size
        # Allocate arrays and unpack spectrum
        wave = np.zeros((nspec, norddet))
        flux = np.zeros((nspec, norddet))
        flux_ivar = np.zeros((nspec, norddet))
        flux_gpm = np.zeros((nspec, norddet), dtype=bool)
        detector = [None]*norddet
        ech_orders = np.zeros(norddet, dtype=int)

        # TODO make the extraction that is desired OPT vs BOX an optional input variable.
        for iorddet in range(norddet):
            wave[:, iorddet] = getattr(self, wave_key)[iorddet]
            flux_gpm[:, iorddet] = getattr(self, '{}_MASK'.format(extract_type))[iorddet]
            detector[iorddet] = self[iorddet].DET
            if self[0].PYPELINE == 'Echelle':
                ech_orders[iorddet] = self[iorddet].ECH_ORDER
            flux[:, iorddet] = getattr(self, flux_key)[iorddet]
            flux_ivar[:, iorddet] = getattr(self, flux_key+'_IVAR')[iorddet] #OPT_FLAM_IVAR

        # Populate meta data
        spectrograph = load_spectrograph(self.header['PYP_SPEC'])

        meta_spec = spectrograph.parse_spec_header(self.header)
        # Add the pyp spec.
        # TODO JFH: Make this an atribute of the specobj by default.
        meta_spec['PYP_SPEC'] = self.header['PYP_SPEC']
        meta_spec['PYPELINE'] = self[0].PYPELINE
        meta_spec['DET'] = np.array(detector)
        meta_spec['DISPNAME'] = self.header['DISPNAME']
        # Return
        if self[0].PYPELINE in ['MultiSlit', 'IFU'] and self.nobj == 1:
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
                If there are multiple detectors arranged in the spectral
                direction, return the sobjs for the standard on each detector.

        Returns:
            SpecObj or SpecObjs or None

        """
        # Is this MultiSlit or Echelle
        pypeline = (self.PYPELINE)[0]
        if 'MultiSlit' in pypeline or 'IFU' in pypeline:
            # Have to do a loop to extract the counts for all objects
            if self.OPT_COUNTS[0] is not None:
                SNR = np.median(self.OPT_COUNTS * np.sqrt(self.OPT_COUNTS_IVAR), axis=1)
            elif self.BOX_COUNTS[0] is not None:
                SNR = np.median(self.BOX_COUNTS * np.sqrt(self.BOX_COUNTS_IVAR), axis=1)
            else:
                return None

            # For multiple detectors grab the requested detectors
            if multi_spec_det is not None:
                # TODO: This is a hack assuming the integers in multi_spec_det
                # are for *detectors*, not mosaics.
                if any([isinstance(d, int) for d in multi_spec_det]):
                    _multi_spec_det = [DetectorContainer.get_name(d) if isinstance(d, int) else d
                                            for d in multi_spec_det]
                else:
                    _multi_spec_det = multi_spec_det
                sobjs_std = SpecObjs(header=self.header)
                # Now append the maximum S/N object on each detector
                for idet in _multi_spec_det:
                    this_det = self.DET == idet
                    if not np.any(this_det):
                        unique_det = np.unique(self.DET)
                        msgs.error(f'No matches for {idet} in spec1d file.  Unique options found'
                                   f"are {', '.join(unique_det)}.  Check usage of multi_spec_det.")
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
                    # Grab SNR
                    if spec[0].OPT_COUNTS is not None:
                        SNR[iord, iobj] = np.median(spec[0].OPT_COUNTS*np.sqrt(spec[0].OPT_COUNTS_IVAR))
                    elif spec[0].BOX_COUNTS is not None:
                        SNR[iord, iobj] = np.median(spec[0].BOX_COUNTS * np.sqrt(spec[0].BOX_COUNTS_IVAR))
                    else:
                        return None
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
        elif sobjs_neg[0].PYPELINE == 'IFU':
            sobjs_neg.OBJID = -sobjs_neg.OBJID
        else:
            msgs.error("The '{0:s}' PYPELINE is not defined".format(self[0].PYPELINE))
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
            elif self[0].PYPELINE == 'IFU':
                index = self.OBJID < 0
            else:
                msgs.error("The '{0:s}' PYPELINE is not defined".format(self[0].PYPELINE))
            self.remove_sobj(index)


    def make_neg_pos(self):
        """
        Purge negative objects from specobjs for IR reductions

        """
        # Assign the sign and the objids
        if self.nobj > 0:
            if self[0].PYPELINE == 'Echelle':
                index = self.ECH_OBJID < 0
            elif self[0].PYPELINE == 'MultiSlit':
                index = self.OBJID < 0
            elif self[0].PYPELINE == 'IFU':
                index = self.OBJID < 0
            else:
                msgs.error("Should not get here")
            try:
                self[index].OPT_COUNTS *= -1
            except (TypeError,ValueError):
                pass
            try:
                self[index].BOX_COUNTS *= -1
            except (TypeError,ValueError):
                pass

    def slitorder_indices(self, slitorder):
        """
        Return the set of indices matching the input slit/order

        Args:
            slitorder (int):
        Returns:
            int:
        """
        if self[0].PYPELINE == 'Echelle':
            indx = self.ECH_ORDER == slitorder
        elif self[0].PYPELINE == 'MultiSlit':
            indx = self.SLITID == slitorder
        elif self[0].PYPELINE == 'IFU':
            indx = self.SLITID == slitorder
        else:
            msgs.error("The '{0:s}' PYPELINE is not defined".format(self[0].PYPELINE))
        #
        return indx

    def name_indices(self, name):
        """
        Return the set of indices matching the input slit/order

        Args:
            name (str): The name of the object

        Returns:
            `numpy.ndarray`_: Array of indices with the corresponding
            name. Shape is (nobj,).
        """
        if self[0].PYPELINE == 'Echelle':
            indx = self.ECH_NAME == name
        elif self[0].PYPELINE == 'MultiSlit':
            indx = self.NAME == name
        elif self[0].PYPELINE == 'IFU':
            indx = self.NAME == name
        else:
            msgs.error("The '{0:s}' PYPELINE is not defined".format(self[0].PYPELINE))
        return indx

    def slitorder_objid_indices(self, slitorder, objid, toler=5):
        """
        Return the set of indices matching the input slit/order and the input
        objid
        
        Args:
            slitorder (int):
                Order/Spatial pixel value for slit of interest.
            objid (int):
                ID value for object of interest.
            toler (int, optional):
                Tolerance for slit spatial pixel values used for slit
                identification. Default = 5

        Returns:
            :obj:`int`: Index value for input slit/order and object ID values
            for specobjs object.

        """

        if self[0].PYPELINE == 'Echelle':
            indx = (self.ECH_ORDER == slitorder) & (self.ECH_OBJID == objid)
        elif self[0].PYPELINE == 'MultiSlit':
            indx = (np.abs(self.SLITID - slitorder) <= toler) & (self.OBJID == objid)
        elif self[0].PYPELINE == 'IFU':
            indx = (self.SLITID == slitorder) & (self.OBJID == objid)
        else:
            msgs.error("The '{0:s}' PYPELINE is not defined".format(self[0].PYPELINE))
        #
        return indx

    def set_names(self):
        """
        Simple method to (re)set the names of all the SpecObj
        """
        for sobj in self.specobjs:
            sobj.set_name()

    def add_sobj(self, sobj):
        """
        Append one or more SpecObj to the existing set.

        Args:
            sobj (:class:`~pypeit.specobj.SpecObj`, :class:`~pypeit.specobjs.SpecObjs`, :obj:`list`, `numpy.ndarray`_):
                One or more SpecObj objects to append.  If a list or array, the
                elements of these must be instances of
                :class:`~pypeit.specobj.SpecObj`.
        """
        if isinstance(sobj, specobj.SpecObj):
            self.specobjs = np.append(self.specobjs, [sobj])
            return
        if isinstance(sobj, SpecObjs):
            # TODO: Possible to instead do np.append(self.specobjs, sobj.specobjs)?
            for isobj in sobj:
                self.specobjs = np.append(self.specobjs, isobj)
            return
        if not isinstance(sobj, (np.ndarray, list)):
            msgs.error(f'Unable to add {type(sobj)} objects to SpecObjs')
        if any([not isinstance(s, specobj.SpecObj) for s in sobj]):
            msgs.error('List or arrays of objects to add must all be of type SpecObj.')
        self.specobjs = np.append(self.specobjs, sobj)

    def remove_sobj(self, index):
        """
        Remove one or more SpecObj by index

        Args:
            index (int, `numpy.ndarray`_):
        """
        msk = np.ones(self.specobjs.size, dtype=bool)
        msk[index] = False
        # Do it
        self.specobjs = self.specobjs[msk]

    def ready_for_fluxing(self):
        # Fluxing
        required_header = ['EXPTIME', 'AIRMASS']  # These are sufficient to apply a sensitivity function
        required_header += ['DISPNAME', 'PYP_SPEC', 'RA', 'DEC']  # These are to generate one
        required_for_fluxing = ['_WAVE', '_COUNTS', '_IVAR']

        chk = True
        # Check header
        for key in required_header:
            chk &= key in self.header

        for sobj in self.specobjs:
            sub_box, sub_opt = True, True
            if sobj is not None:
                # Only need one of these but need them all
                for item in required_for_fluxing:
                    if hasattr(sobj, 'BOX'+item):
                        sub_box &= True
                    if hasattr(sobj, 'OPT'+item):
                        sub_opt &= True
                # chk
                chk &= (sub_box or sub_opt)
        return chk

    def copy(self):
        """
        Generate a copy of self

        Returns:
            :class:`SpecObjs`:

        """
        sobj_copy = SpecObjs(header=self.header)
        for sobj in self.specobjs:
            sobj_copy.add_sobj(sobj.copy())
        return sobj_copy

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

    def __getattr__(self, attr):
        """
        Overloaded to generate an array of attribute 'k' from the
        :class:`pypeit.specobj.SpecObj` objects.

        First attempts to grab data from the Summary table, then the list
        """
        if len(self.specobjs) == 0:
            raise ValueError('SpecObjs is empty!')
        if attr in self.__dict__:  # any normal attributes are handled normally
            return self.__dict__[attr]
        try:
            getattr(self.specobjs[0], attr)
        except NameError:
            raise NameError(f'{attr} is not an attribute of SpecObjs or SpecObj.')

        return lst_to_array([getattr(specobj, attr) for specobj in self.specobjs])

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

    @property
    def size(self):
        return self.specobjs.size

    @property
    def shape(self):
        return self.specobjs.shape

    def write_to_fits(self, subheader, outfile, overwrite=True, update_det=None,
                      slitspatnum=None, history=None, debug=False):
        """
        Write the set of SpecObj objects to one multi-extension FITS file

        Args:
            subheader (:obj:`dict`):
            outfile (str):
            overwrite (bool, optional):
            slitspatnum (:obj:`str` or :obj:`list`, optional):
              Restricted set of slits for reduction.
              If provided, do not clobber the existing file but only update
              the indicated slits.  Useful for re-running on a subset of slits
            update_det (int or list, optional):
              If provided, do not clobber the existing file but only update
              the indicated detectors.  Useful for re-running on a subset of detectors

        """
        if os.path.isfile(outfile) and not overwrite:
            msgs.warn(f'{outfile} exits. Set overwrite=True to overwrite it.')
            return

        # If the file exists and update_det (and slit_spat_num) is provided, use the existing header
        #   and load up all the other hdus so that we only over-write the ones
        #   we are updating
        if os.path.isfile(outfile) and (update_det is not None or slitspatnum is not None):
            _specobjs = SpecObjs.from_fitsfile(outfile)
            mask = np.ones(_specobjs.nobj, dtype=bool)
            # Update
            if slitspatnum is not None: # slitspatnum
                dets, spat_ids = parse.parse_slitspatnum(slitspatnum)
                for det, spat_id in zip(dets, spat_ids):
                    mask[(_specobjs.DET == det) & (_specobjs.SLITID == spat_id)] = False
            elif update_det is not None:
                # Pop out those with this detector (and slit if slit_spat_num is provided)
                for det in np.atleast_1d(update_det):
                    mask[_specobjs.DET == det] = False
            _specobjs = _specobjs[mask]
            # Add in the new
            # TODO: Is the loop necessary? add_sobj can take many SpecObj objects.
            for sobj in self.specobjs:
                _specobjs.add_sobj(sobj)
        else:
            _specobjs = self.specobjs

        # Build up the Header
        header = io.initialize_header()
        for key in subheader.keys():
            if key.upper() == 'HISTORY':
                if history is None:
                    for line in str(subheader[key.upper()]).split('\n'):
                        header[key.upper()] = line
            else:
                header[key.upper()] = subheader[key]

        # Init
        prihdu = fits.PrimaryHDU(header=header)
        hdus = [prihdu]

        # Add class info
        prihdu.header['DMODCLS'] = (self.__class__.__name__, 'Datamodel class')
        prihdu.header['DMODVER'] = (self.version, 'Datamodel version')

        # Add history
        if history is not None:
            history.write_to_header(prihdu.header)

        detector_hdus = {}
        nspec, ext = 0, 0
        # Loop on the SpecObj objects
        for sobj in _specobjs:
            if sobj is None:
                continue
            # HDUs
            if debug:
                raise NotImplementedError('Debugging for developers only.')
                #embed()
                #exit()
            shdul = sobj.to_hdu()
            if len(shdul) not in [1, 2]:
                msgs.error('CODING ERROR: SpecObj datamodel changed.  to_hdu should return 1 or 2 '
                           'HDUs.  If returned, the 2nd one should be the detector/mosaic.')
            if len(shdul) == 2:
                detector_hdus[sobj['DET']] = shdul[1]
                shdu = [shdul[0]]
            else:
                shdu = shdul

            if len(shdu) != 1 or not isinstance(shdu[0], fits.hdu.table.BinTableHDU):
                msgs.error('CODING ERROR: SpecObj datamodel changed.')

            # Name
            shdu[0].name = sobj.NAME
            # Extension
            keywd = 'EXT{:04d}'.format(ext)
            prihdu.header[keywd] = sobj.NAME
            ext += 1
            nspec += 1
            # Append
            hdus += shdu

        # Deal with Detectors
        for key, item in detector_hdus.items():
            prefix = item.header['name']
            # Name
            if prefix not in item.name:  # In case we are re-loading
                item.name = f'{prefix}-{item.name}'
            # Append
            hdus += [item]

        # A few more for the header
        prihdu.header['NSPEC'] = nspec

        # Finish
        hdulist = fits.HDUList(hdus)
        if debug:
            raise NotImplementedError('Debugging for developers only.')
             #embed()
             #exit()
        hdulist.writeto(outfile, overwrite=overwrite)
        msgs.info(f'Wrote 1D spectra to {outfile}')

    def write_info(self, outfile, pypeline):
        """
        Write a summary of items to an ASCII file

        Args:
            outfile (:obj:`str`):  Output filename
            pypeline (:obj:`str`): PypeIt pipeline mode
        """
        # TODO -- Deal with update_det
        # Lists for a Table
        slits, names, maskdef_id, objname, objra, objdec, spat_pixpos, spat_fracpos, boxsize, opt_fwhm, s2n = \
            [], [], [], [], [], [], [], [], [], [], []
        wave_rms = []
        maskdef_extract = []
        manual_extract = []
        # binspectral, binspatial = parse.parse_binning(binning)
        for specobj in self.specobjs:
            det = specobj.DET
            if specobj is None:
                continue
            # Detector items
            binspectral, binspatial = parse.parse_binning(specobj.DETECTOR.binning)
            platescale = specobj.DETECTOR.platescale
            # Append
            spat_pixpos.append(specobj.SPAT_PIXPOS)
            if pypeline == 'MultiSlit':
                spat_fracpos.append(specobj.SPAT_FRACPOS)
                slits.append(specobj.SLITID)
                names.append(specobj.NAME)
            elif pypeline == 'IFU':
                spat_fracpos.append(specobj.SPAT_FRACPOS)
                slits.append(specobj.SLITID)
                names.append(specobj.NAME)
            elif pypeline == 'Echelle':
                spat_fracpos.append(specobj.ECH_FRACPOS)
                slits.append(specobj.ECH_ORDER)
                names.append(specobj.ECH_NAME)
            # Wave RMS
            wave_rms.append(specobj.WAVE_RMS)
            # Boxcar width
            if specobj.BOX_RADIUS is not None:
                slit_pix = 2.0 * specobj.BOX_RADIUS
                # Convert to arcsec
                binspectral, binspatial = parse.parse_binning(specobj.DETECTOR.binning)
                #binspectral, binspatial = parse.parse_binning(binning)
                # JFH TODO This should be using the order_platescale for each order. Furthermore, not all detectors
                # have the same platescale, i.e. with GNIRS it is the same detector but a different camera hence a
                # different attribute. platescale should be a spectrograph attribute determined on the fly.
                # boxsize.append(slit_pix*binspatial*spectrograph.detector[specobj.DET-1]['platescale'])
                boxsize.append(slit_pix * binspatial * platescale)
            else:
                boxsize.append(0.)

            # Optimal profile (FWHM)
            # S2N -- default to boxcar
            if specobj.FWHMFIT is not None and specobj.OPT_COUNTS is not None:
                opt_fwhm.append(np.median(specobj.FWHMFIT) * binspatial * platescale)
            else:  # Optimal is not required to occur
                opt_fwhm.append(0.)
            # NOTE: Below requires that S2N not be None, otherwise the code will
            # fault.  If the code gets here and S2N is None, check that 1D
            # extractions have been performed.
            s2n.append(specobj.S2N)
            # Manual extraction?
            manual_extract.append(specobj.hand_extract_flag)
            # Slitmask info
            maskdef_id.append(specobj.MASKDEF_ID)
            objname.append(specobj.MASKDEF_OBJNAME)
            objra.append(specobj.RA)
            objdec.append(specobj.DEC)
            maskdef_extract.append(specobj.MASKDEF_EXTRACT)

        # Generate the table, if we have at least one source
        if len(names) > 0:
            obj_tbl = Table()
            if pypeline == 'MultiSlit':
                obj_tbl['slit'] = slits
                obj_tbl['slit'].format = 'd'
            elif pypeline == 'IFU':
                obj_tbl['slit'] = slits
                obj_tbl['slit'].format = 'd'
            elif pypeline == 'Echelle':
                obj_tbl['order'] = slits
                obj_tbl['order'].format = 'd'
            obj_tbl['name'] = names
            if not np.all(np.array(maskdef_id) == None):
                obj_tbl['maskdef_id'] = maskdef_id
            if not np.all(np.array(objname) == None):
                obj_tbl['objname'] = objname
            if not np.all(np.array(objra) == None):
                obj_tbl['objra'] = objra
                if None not in objra:
                    obj_tbl['objra'].format = '.5f'
                obj_tbl['objdec'] = objdec
                if None not in objdec:
                    obj_tbl['objdec'].format = '.5f'
            obj_tbl['spat_pixpos'] = spat_pixpos
            obj_tbl['spat_pixpos'].format = '.1f'
            obj_tbl['spat_fracpos'] = spat_fracpos
            obj_tbl['spat_fracpos'].format = '.3f'
            obj_tbl['box_width'] = boxsize
            obj_tbl['box_width'].format = '.2f'
            obj_tbl['box_width'].unit = units.arcsec
            obj_tbl['opt_fwhm'] = opt_fwhm
            obj_tbl['opt_fwhm'].format = '.3f'
            obj_tbl['opt_fwhm'].unit = units.arcsec
            obj_tbl['s2n'] = s2n
            obj_tbl['s2n'].format = '.2f'
            # is this a forced extraction at the expected position from slitmask design?
            if not np.all(np.array(maskdef_extract) == None):
                    obj_tbl['maskdef_extract'] = maskdef_extract
            # only if any manual extraction exists, print this
            if not np.all(np.array(manual_extract) == False):
                obj_tbl['manual_extract'] = manual_extract

            # Wavelengths
            if not np.all(np.array(wave_rms) == None):
                obj_tbl['wv_rms'] = wave_rms
                if None not in wave_rms:
                    obj_tbl['wv_rms'].format = '.3f'
            # Write
            obj_tbl.write(outfile,format='ascii.fixed_width', overwrite=True)

    def get_extraction_groups(self, model_full_slit=False) -> List[List[int]]:
        """
        Returns:
            List[List[int]]: A list of extraction groups, each of which is a list of integer
                object indices that should be extracted together by core.skysub.local_skysub_extract
        """
        nobj = len(self.specobjs)

        if model_full_slit:
            return [list(range(nobj))]

        # initialize adjacency matrix
        adj = np.full((nobj, nobj), dtype=bool, fill_value=False)
        ## build adjacency matrix
        # adj[i, j] is True iff objects i and j are touching each other
        for i in range(nobj):
            left_edge_i = self.specobjs[i].TRACE_SPAT - self.specobjs[i].maskwidth - 1
            righ_edge_i = self.specobjs[i].TRACE_SPAT + self.specobjs[i].maskwidth + 1
            for j in range(i + 1, nobj):
                left_edge_j = self.specobjs[j].TRACE_SPAT - self.specobjs[j].maskwidth - 1
                righ_edge_j = self.specobjs[j].TRACE_SPAT + self.specobjs[j].maskwidth + 1

                touch = (left_edge_j < righ_edge_i) & (left_edge_i < righ_edge_j)

                if touch.any():
                    adj[i, j] = True
                    adj[j, i] = True

        ## Find all connected components in the graph of objects.
        # One call to DFS will visit every object it can that is "connected" to
        # the starting object by the touching relation.
        visited = [False]*nobj
        groups = []
        while not all(visited):
            # pick a starting unvisited vertex
            v = visited.index(False)
            group = []
            # DFS starting at v. Afterwards, group contains every object that
            # is "connected" to v by the touching relation.
            utils.DFS(v, visited, group, adj)
            groups.append(group)

        return groups

#TODO Should this be a classmethod on specobjs??
def get_std_trace(detname, std_outfile, chk_version=True):
    """
     Returns the trace of the standard.

     Args:
         det (:obj:`int`, :obj:`tuple`):
             1-indexed detector(s) to process.
         std_outfile (:obj:`str`):
             Filename with the standard star spec1d file.  Can be None.
     Returns:
         `numpy.ndarray`_: Trace of the standard star on input detector.
         Will be None if ``std_outfile`` is None, or if the selected detector/mosaic is not available
         in the provided spec1d file.
     """

    sobjs = SpecObjs.from_fitsfile(std_outfile, chk_version=chk_version)
    pypeline = sobjs.PYPELINE
    # Does the detector match?
    # TODO: Instrument specific logic here could be implemented with the
    # parset. For example LRIS-B or LRIS-R we we would use the standard
    # from another detector.

    this_det = sobjs.DET == detname
    if np.any(this_det):
        sobjs_det = sobjs[this_det]
        sobjs_std = sobjs_det.get_std()
        # No standard extracted on this detector??
        if sobjs_std is None:
            return None
        std_trace = sobjs_std.TRACE_SPAT
        # flatten the array if this multislit
        if 'MultiSlit' in pypeline:
            std_trace = std_trace.flatten()
        elif 'Echelle' in pypeline:
            std_trace = std_trace.T
        else:
            msgs.error('Unrecognized pypeline')
    else:
        std_trace = None

    return std_trace

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
    if isinstance(lst[0], units.Quantity):
        return units.Quantity(lst)[mask]
    else:
        return np.array(lst)[mask]

