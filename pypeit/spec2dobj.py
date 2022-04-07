"""
Module for the Spec2DObj class

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import os
import inspect
import datetime

from copy import deepcopy

import numpy as np

from astropy.io import fits
from astropy.stats import mad_std
import astropy

from pypeit import msgs
from pypeit import io
from pypeit import datamodel
from pypeit import slittrace
from pypeit.core import parse 
from pypeit.images import imagebitmask
from pypeit.images.detector_container import DetectorContainer
from pypeit.images.mosaic import Mosaic

from IPython import embed

#def spec2d_hdu_prefix(det):
#    return 'DET{:02d}-'.format(det)


class Spec2DObj(datamodel.DataContainer):
    """Class to handle 2D spectral image outputs of PypeIt

    One generates one of these Objects for each detector in the exposure.

    See datamodel below and at :ref:`spec2dobj_datamodel`

    Args:

    Attributes:
        head0 (`astropy.fits.Header`):
            Primary header if instantiated from a FITS file

    """
    version = '1.0.4'

    # TODO 2d data model should be expanded to include:
    # waveimage  --  flexure and heliocentric corrections should be applied to the final waveimage and since this is unique to
    #                every exposure (i.e. it depneds on obstime, RA, DEC and the flexure incurred) it should be written out for
    #                each science frame.

    # Because we are including nested DataContainers, be careful not to
    # duplicate variable names!!
    datamodel = {'sciimg': dict(otype=np.ndarray, atype=np.floating,
                                descr='2D processed science image (float32)'),
                 'ivarraw': dict(otype=np.ndarray, atype=np.floating,
                                 descr='2D processed inverse variance image (float32)'),
                 'skymodel': dict(otype=np.ndarray, atype=np.floating,
                                  descr='2D sky model image (float32)'),
                 'objmodel': dict(otype=np.ndarray, atype=np.floating,
                                  descr='2D object model image (float32)'),
                 'ivarmodel': dict(otype=np.ndarray, atype=np.floating,
                                   descr='2D ivar model image (float32)'),
                 'tilts': dict(otype=np.ndarray, atype=np.floating,
                               descr='2D tilts image (float64)'),
                 'scaleimg': dict(otype=np.ndarray, atype=np.floating,
                                  descr='2D multiplicative scale image [or a single scalar as an array] that has been applied to '
                                        'the science image (float32)'),
                 'waveimg': dict(otype=np.ndarray, atype=np.floating,
                                 descr='2D wavelength image in vacuum (float64)'),
                 'bpmmask': dict(otype=np.ndarray, atype=np.integer,
                                 descr='2D bad-pixel mask for the image'),
                 'imgbitm': dict(otype=str, descr='List of BITMASK keys from ImageBitMask'),
                 'slits': dict(otype=slittrace.SlitTraceSet,
                               descr='SlitTraceSet defining the slits'),
                 'maskdef_designtab': dict(otype=astropy.table.Table,
                                           descr='Table with slitmask design and object info'),
                 'sci_spat_flexure': dict(otype=float,
                                          descr='Shift, in spatial pixels, between this image '
                                                'and SlitTrace'),
                 'sci_spec_flexure': dict(otype=astropy.table.Table,
                                          descr='Global shift of the spectrum to correct for spectral'
                                                'flexure (pixels). This is based on the sky spectrum at'
                                                'the center of each slit'),
                 'vel_type': dict(otype=str, descr='Type of reference frame correction (if any). '
                                                   'Options are listed in the routine: '
                                                   'WavelengthSolutionPar.valid_reference_frames() '
                                                   'Current list: observed, heliocentric, barycentric'),
                 'vel_corr': dict(otype=float,
                                  descr='Relativistic velocity correction for wavelengths'),
                 'med_chis': dict(otype=np.ndarray, atype=np.floating,
                               descr='Median of the chi image for each slit/order'),
                 'std_chis': dict(otype=np.ndarray, atype=np.floating,
                               descr='std of the chi image for each slit/order'),
                 'det': dict(otype=int, descr='Detector index'),
                 'detector': dict(otype=(DetectorContainer, Mosaic),
                                  descr='Detector or Mosaic metadata') }

    @classmethod
    def from_file(cls, file, detname, chk_version=True):
        """
        Overload :func:`pypeit.datamodel.DataContainer.from_file` to allow det
        input and to slurp the header

        Args:
            file (:obj:`str`):
                File name to read.
            detname (:obj:`str`):
                The string identifier for the detector or mosaic used to select
                the data that is read.
            chk_version (:obj:`bool`, optional):
                If False, allow a mismatch in datamodel to proceed

        Returns:
            :class:`~pypeit.spec2dobj.Spec2DObj`: 2D spectra object.

        """
        hdul = io.fits_open(file)
        # Quick check on det
        if not np.any([detname in hdu.name for hdu in hdul]):
            msgs.error(f'{detname} not available in any extension of {file}')
        slf = super().from_hdu(hdul, hdu_prefix=f'{detname}-', chk_version=chk_version)
        slf.head0 = hdul[0].header
        slf.chk_version = chk_version
        return slf

    def __init__(self, sciimg, ivarraw, skymodel, objmodel, ivarmodel,
                 scaleimg, waveimg, bpmmask, detector, sci_spat_flexure, sci_spec_flexure,
                 vel_type, vel_corr, slits, tilts, maskdef_designtab):
        # Slurp
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        _d = dict([(k,values[k]) for k in args[1:]])
        # Setup the DataContainer
        datamodel.DataContainer.__init__(self, d=_d)

    def _init_internals(self):
        self.process_steps = None
        self.head0 = None
        self.chk_version = None  # Mainly for viewing/using old versions

    def _validate(self):
        """
        Assert that the detector has been set

        Returns:

        """
        # Check the bitmask is current
        bitmask = imagebitmask.ImageBitMask()
        if self.imgbitm is None:
            self.imgbitm = ','.join(list(bitmask.keys()))
        else:
            # Validate
            if self.imgbitm != ','.join(list(bitmask.keys())) and self.chk_version:
                msgs.error("Input BITMASK keys differ from current data model!")

        # Check the detector/mosaic identifier has been provided (note this is a
        # property method)
        if self.detname is None:
            msgs.error('Detector/Mosaic string identifier must be set at instantiation.')

    def _bundle(self):
        """
        Over-write default _bundle() method to separate the DetectorContainer
        into its own HDU

        Returns:
            :obj:`list`: A list of dictionaries, each list element is
            written to its own fits extension. See the description
            above.
        """
        d = []
        # Rest of the datamodel
        for key in self.keys():
            # Skip Nones
            if self[key] is None:
                continue
            # Array?
            if self.datamodel[key]['otype'] == np.ndarray:
                tmp = {}
                if self.datamodel[key]['atype'] == np.floating and key not in ['waveimg', 'tilts']:
                    tmp[key] = self[key].astype(np.float32)
                else:
                    tmp[key] = self[key]
                d.append(tmp)
            # Detector
            elif key == 'detector':
                d.append(dict(detector=self.detector))
            # SlitTraceSet
            elif key == 'slits':
                d.append(dict(slits=self.slits))
            # maskdef_designtab
            elif key == 'maskdef_designtab':
                d.append(dict(maskdef_designtab=self.maskdef_designtab))
            # Spectral flexure
            elif key == 'sci_spec_flexure':
                d.append(dict(sci_spec_flexure=self.sci_spec_flexure))
            else: # Add to header of the primary image
                d[0][key] = self[key]
        # Return
        return d

    @property
    def detname(self):
        if self.detector is None:
            return None
        return self.detector.name

    @property
    def hdu_prefix(self):
        """
        Provides for a dynamic hdu_prefix based on our naming model.

        Returns:
            :obj:`str`: Detector/mosaic identifier

        """
        return f'{self.detname}-'

    def update_slits(self, spec2DObj):
        """
        Update the object at all good slits in the input object

        Args:
            spec2DObj (`Spec2DObj`):

        Returns:

        """
        # Quick checks
        if spec2DObj.detname != self.detname:
            msgs.error("Objects are not even the same detector!!")
        if not np.array_equal(spec2DObj.slits.spat_id, spec2DObj.slits.spat_id):
            msgs.error("SPAT_IDs are not in sync!")

        # Find the good ones on the input object
        bpm = spec2DObj.slits.mask.astype(bool)
        exc_reduce = np.invert(spec2DObj.slits.bitmask.flagged(
            spec2DObj.slits.mask, flag=spec2DObj.slits.bitmask.exclude_for_reducing))
        gpm = np.invert(bpm & exc_reduce)

        # Update slits.mask
        self.slits.mask[gpm] = spec2DObj.slits.mask[gpm]

        # Slitmask
        slitmask = spec2DObj.slits.slit_img(flexure=spec2DObj.sci_spat_flexure,
                                                 exclude_flag=spec2DObj.slits.bitmask.exclude_for_reducing)
        # Fill in the image
        for slit_idx, spat_id in enumerate(spec2DObj.slits.spat_id[gpm]):
            inmask = slitmask == spat_id
            # Get em all
            for imgname in ['sciimg','ivarraw','skymodel','objmodel','ivarmodel','waveimg','bpmmask']:
                self[imgname][inmask] = spec2DObj[imgname][inmask]

    def calc_chi_slit(self, slitidx:int, pad:int=None, remove_object:bool=True):
        """ Calculate a chi map and run some stats on it
        for a given slit/order

        Args:
            slitidx (int): Given slit/order
            pad (int, optional):  Ignore pixels within pad of edges. 
                Defaults to None.
            remove_object (bool, optional):  Remove object model (if 
                it exists)

        Returns:
            tuple: np.ndarray (chi image), median (float), std (float)
        """
        slit_select = self.slits.slit_img(pad=pad, slitidx=slitidx)
        skysub_img = self.sciimg - self.skymodel
        if remove_object and self.objmodel is not None:
            skysub_img -= self.objmodel

        # Chi
        chi = skysub_img * np.sqrt(self.ivarmodel) * (self.bpmmask == 0)
        chi_slit = chi * (slit_select == self.slits.spat_id[slitidx]) * (self.bpmmask == 0)

        # All bad?
        if np.all(chi_slit == 0):
            return None, 0., 0.
        
        # Stats
        median = np.median(chi_slit[chi_slit!=0]) 
        std = mad_std(chi_slit[chi_slit!=0])
        #
        return chi_slit, median, std

    def gen_qa(self):
        """ Generate QA for the slits/orders 

        Saved to the DataContainer
        """

        # Loop on slits to generate stats on chi^2
        med_chis = []
        std_chis = []
        for slitidx in range(self.slits.nslits):
            _, med, std = self.calc_chi_slit(slitidx)
            med_chis.append(med)
            std_chis.append(std)
        # Save
        self.med_chis = np.array(med_chis)
        self.std_chis = np.array(std_chis)
        return

class AllSpec2DObj:
    """
    Simple object to hold Spec2DObj objects
    and perform I/O

    Anything that goes into self['meta'] must be parseable into a FITS Header

    Restrict keys to be type int or 'meta'
    and items to be :class:`Spec2DObj`

    """
    hdr_prefix = 'ALLSPEC2D_'
    @classmethod
    def from_fits(cls, filename, chk_version=True):
        """

        Args:
            filename (:obj:`str`):
                Name of the file to read.
            chk_version (:obj:`bool`, optional):
                If True, demand the on-disk datamodel equals the current one.
                Passed to from_hdu() of DataContainer.

        Returns:
            :class:`~pypeit.spec2dobj.AllSpec2DObj`: The constructed object.
        """
        # Instantiate
        self = cls()
        # Open
        hdul = io.fits_open(filename)
        # Meta
        for key in hdul[0].header.keys():
            if key == self.hdr_prefix+'DETS':
                continue
            if self.hdr_prefix in key:
                meta_key = key.split(self.hdr_prefix)[-1].lower()
                self['meta'][meta_key] = hdul[0].header[key]
        # Detectors included
        detectors = hdul[0].header[self.hdr_prefix+'DETS']
        for detname in detectors.split(','):
            self[detname] = Spec2DObj.from_hdu(hdul, hdu_prefix=f'{detname}-',
                                               chk_version=chk_version)
        return self

    def __init__(self):
        # Init meta
        self['meta'] = {}

    # TODO -- Turn off attribute setting

    @property
    def detectors(self):
        """
        Return the list of detector/mosaic names, assuming they are the list of
        all keys except for meta.
        """
        dets = self.keys
        dets.remove('meta')
        dets.sort()
        return dets

    @property
    def keys(self):
        """
        Return the list of attributes
        """
        return list(self.__dict__.keys())

#    def __setitem__(self, item, value):
#        """
#        Over-load this to insist the item is either:
#          `meta`
#          1,2,3,...
#
#        Args:
#            item (:obj:`str` or :obj:`int`):
#            value (object or :class:`Spec2DObj`):
#                Should be FITS header write-able if going into `meta`
##        """
#        # Check item
#        if not isinstance(item, int) and item != 'meta':
#            raise KeyError('Key must be an integer, i.e. detector number or "meta"')
#        # Check value
#        if isinstance(item, int):
#            assert isinstance(value, Spec2DObj), 'Item must be a Spec2DObj'
#        self.__dict__[item] = value

    def __setitem__(self, item, value):
        # Check item
        if not isinstance(item, str):
            raise TypeError('Item key must be a string.')
        if item != 'meta':
            if not isinstance(value, Spec2DObj):
                raise KeyError('Any item not assigned to the meta dictionary must be a Spec2DObj.')
            if value.detname is not None and value.detname != item:
                msgs.warn(f'Mismatch between keyword used to define the Spec2DObj item ({item}) '
                          f'and the name of the detector/mosaic ({value.detname}).')
        self.__dict__[item] = value

    def __getitem__(self, item):
        """Get an item directly from the internal dict."""
        return self.__dict__[item]

    def build_primary_hdr(self, raw_header, spectrograph, master_key_dict=None, master_dir=None,
                          redux_path=None, subheader=None, history=None):
        """
        Build the primary header for a spec2d file

        Args:
            raw_header (`astropy.io.fits.Header`_):
                Header from the raw FITS file (i.e. original header)
            spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
                Spectrograph used to obtain the data.
            master_key_dict (:obj:`dict`, optional):
                Dictionary of master keys from :class:`~pypeit.calibrations.Calibrations`.
            master_dir (:obj:`str`):
                Path to the ``Masters`` folder
            redux_path (:obj:`str`, optional):
                Full path to the reduction output files.
            subheader (:obj:`dict`, optional):
                Generated by
                :func:`~pypeit.spectrographs.spectrograph.Spectrograph.subheader_for_spec`.

        Returns:
            `astropy.io.fits.Header`_: The primary header for the output fits file.
        """
        hdr = io.initialize_header(primary=True)

        hdukeys = ['BUNIT', 'COMMENT', '', 'BITPIX', 'NAXIS', 'NAXIS1', 'NAXIS2',
                   'HISTORY', 'EXTEND', 'DATASEC']
        for key in raw_header.keys():
            # Use new ones
            if key in hdukeys:
                continue
            # Update unused ones
            hdr[key] = raw_header[key]
        # History
        if history is not None:
            history.write_to_header(hdr)

        # Sub-header
        if subheader is not None:
            for key in subheader.keys():
                hdr[key.upper()] = subheader[key]

        # PYPEIT
        # TODO Should the spectrograph be written to the header?
        hdr['PIPELINE'] = str('PYPEIT')
        hdr['PYPELINE'] = spectrograph.pypeline
        hdr['PYP_SPEC'] = spectrograph.name
        hdr['DATE-RDX'] = str(datetime.date.today().strftime('%Y-%m-%d'))

        # MasterFrame info
        # TODO -- Should this be in the header of the individual HDUs ?
        if master_key_dict is not None:
            if 'bias' in master_key_dict.keys():
                hdr['BIASMKEY'] = master_key_dict['bias']
            if 'arc' in master_key_dict.keys():
                hdr['ARCMKEY'] = master_key_dict['arc']
            if 'trace' in master_key_dict.keys():
                hdr['TRACMKEY'] = master_key_dict['trace']
            if 'flat' in master_key_dict.keys():
                hdr['FLATMKEY'] = master_key_dict['flat']

        # Processing steps
        # TODO: Assumes processing steps for all detectors are the same...  Does
        # this matter?
        det = self.detectors[0]
        if self[det].process_steps is not None:
            hdr['PROCSTEP'] = (','.join(self[det].process_steps), 'Completed reduction steps')

        # Some paths
        if master_dir is not None:
            hdr['PYPMFDIR'] = str(master_dir)
        if redux_path is not None:
            hdr['PYPRDXP'] = redux_path
        # Sky sub mode
        if 'bkg_redux' in self['meta'] and self['meta']['bkg_redux']:
            hdr['SKYSUB'] = 'DIFF'
        else:
            hdr['SKYSUB'] = 'MODEL'
        # obj find mode
        if 'find_negative' in self['meta'] and self['meta']['find_negative']:
            hdr['FINDOBJ'] = 'POS_NEG'
        else:
            hdr['FINDOBJ'] = 'POS'
         
        return hdr

    def write_to_fits(self, outfile, pri_hdr=None, update_det=None, 
                      slitspatnum=None, overwrite=True):
        """
        Write the spec2d FITS file

        Args:
            outfile (:obj:`str`):
                Output filename
            pri_hdr (:class:`astropy.io.fits.Header`, optional):
                Header to be used in lieu of default
                Usually generated by :func:`pypeit,spec2dobj.AllSpec2DObj.build_primary_hdr`
            slitspatnum (:obj:`str` or :obj:`list`, optional):
              Restricted set of slits for reduction
              If provided, do not clobber the existing file but only update
              the indicated slits.  Useful for re-running on a subset of slits
            pri_hdr():
                Baseline primary header.  If None, initial primary header is
                empty.  Usually generated by
                :func:`pypeit,spec2dobj.AllSpec2DObj.build_primary_hdr`
            update_det (:obj:`list`, optional):
                If the output file already exists, this sets the list of
                detectors/mosaics to update with the data in this object.  If
                None, a new file is constructed from scratch, if overwrite is
                True.  Otherwise, the existing file is read and any detectors in
                that file but *not* in this one are added to this object.  I.e.,
                if ``update_det`` is not None, **this method can alter the
                object**.
            overwrite (:obj:`bool`, optional):
                If true and the output file already exists, overwrite it.  The
                combination of this and ``update_det`` may also alter this
                object based on the existing file.
        """
        if os.path.isfile(outfile):
            # Clobber?
            if not overwrite:
                msgs.warn("File {} exits.  Use -o to overwrite.".format(outfile))
                return
            # Load up the original
            _allspecobj = AllSpec2DObj.from_fits(outfile)
            # Replace the newly reduced detector?
            if update_det is not None:
                for det in _allspecobj.detectors:
                    if det in np.atleast_1d(update_det):
                        continue
                    else:
                        self[det] = _allspecobj[det]
            elif slitspatnum is not None: # Update specific slits!
                
                # Grab modified detectors and slits
                dets, spat_ids = parse.parse_slitspatnum(slitspatnum)

                # Loop on detectors to be fussed with
                for det in _allspecobj.detectors:
                    if det in dets:
                        # Check version 
                        if self[det].version != _allspecobj[det].version:
                            msgs.error("Original spec2D object has a different version.  Too risky to continue.  Rerun both")
                        # Generate the slit "mask"
                        slitmask = _allspecobj[det].slits.slit_img(
                            flexure=_allspecobj[det].sci_spat_flexure)
                        # Save the new one in a copy
                        new_Spec2DObj = deepcopy(self[det])
                        # Replace with the old
                        self[det] = _allspecobj[int(det)]
                        # Spat ids    
                        spats = spat_ids[dets==det]
                        for spat_id in spats:
                            # Find pixels to replace
                            replace_pix = slitmask == spat_id
                            # Fill em in
                            for key in new_Spec2DObj.datamodel.keys():
                                if new_Spec2DObj.datamodel[key]['otype'] == np.ndarray and (
                                    new_Spec2DObj[key].shape == slitmask.shape):
                                    self[det][key][replace_pix] = new_Spec2DObj[key][replace_pix]
                    else:
                        self[det] = _allspecobj[det]

        # Primary HDU for output
        prihdu = fits.PrimaryHDU()
        # Header
        if pri_hdr is not None:
            prihdu.header = pri_hdr

        # Add class name
        prihdu.header['PYP_CLS'] = self.__class__.__name__

        # Add meta to Primary Header
        for key in self['meta']:
            try:
                prihdu.header[self.hdr_prefix+key.upper()] = self['meta'][key]
            except:
                msgs.warn(f'Cannot add meta entry {key} to primary header!')
                continue
            if key.lower() != key:
                msgs.warn('Keywords in the meta dictionary are always read back in as lower case. '
                          f'Subsequent reads of {outfile} will have converted {key} to '
                          f'{key.lower()}!')

        # Loop on em (in order of detector)
        extnum = 1
        hdus = [prihdu]
        for det in self.detectors:
            hdul = self[det].to_hdu()
            # TODO -- Make adding EXT000X a default of DataContainer?
            # TODO: Why is this needed?
            for hdu in hdul:
                keywd = 'EXT{:04d}'.format(extnum)
                prihdu.header[keywd] = hdu.name
                extnum += 1
            # Add em in
            hdus += hdul

        # Detectors included
        prihdu.header[self.hdr_prefix+'DETS'] = ','.join(self.detectors)

        # Finish
        hdulist = fits.HDUList(hdus)
        hdulist.writeto(outfile, overwrite=overwrite)
        msgs.info("Wrote: {:s}".format(outfile))

    def __repr__(self):
        # Generate sets string
        txt = '<{:s}: '.format(self.__class__.__name__)
        txt += 'dets=('
        for det in self.detectors:
            txt += str(det)+','
        txt += ') >'
        return txt

