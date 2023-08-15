"""
Module for INT/IDS specific methods.

.. include:: ../include/links.rst
"""
import os

from IPython import embed

import numpy as np

from astropy.time import Time

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch, parse
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container


class INTIDSSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle INT/IDS specific code
    """
    ndet = 1
    telescope = telescopes.INTTelescopePar()
    url = 'https://www.ing.iac.es/astronomy/instruments/ids/'

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # Ignore PCA
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        # Bound the detector with slit edges if no edges are found
        par['calibrations']['slitedges']['bound_detector'] = True

        # Always correct for flexure, starting with default parameters
        par['flexure']['spec_method'] = 'boxcar'
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 0.001]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['pixelflatframe']['exprng'] = [0, None]
        par['calibrations']['traceframe']['exprng'] = [0, None]
        par['calibrations']['arcframe']['exprng'] = [None, 61]
        par['calibrations']['standardframe']['exprng'] = [1, 61]
        #
        par['scienceframe']['exprng'] = [61, None]
        return par

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='RA')
        self.meta['dec'] = dict(ext=0, card='DEC')
        self.meta['target'] = dict(ext=0, card='OBJECT')
        self.meta['decker'] = dict(ext=0, card='SLTWDSKY')
        self.meta['binning'] = dict(ext=0, card=None, compound=True)
        self.meta['mjd'] = dict(ext=0, card='MJD-OBS')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        self.meta['dispname'] = dict(ext=0, card='GRATNAME')
        self.meta['dispangle'] = dict(ext=0, card='GRATANGL', rtol=1e-3)  # TODO :: Need to check the tolerance is OK
        # Additional ones, generally for configuration determination or time
        self.meta['instrument'] = dict(ext=0, card='DETECTOR')
        self.meta['lampstat01'] = dict(ext=0, card=None, compound=True)
        self.meta['lampstat02'] = dict(ext=0, card=None, compound=True)

    def compound_meta(self, headarr, meta_key):
        """
        Methods to generate metadata requiring interpretation of the header
        data, instead of simply reading the value of a header card.

        Args:
            headarr (:obj:`list`):
                List of `astropy.io.fits.Header`_ objects.
            meta_key (:obj:`str`):
                Metadata keyword to construct.

        Returns:
            object: Metadata value read from the header(s).
        """
        if meta_key == 'binning':
            binspatial, binspec = parse.parse_binning(headarr[0]['CCDSUM'])
            binning = parse.binning2string(binspec, binspatial)
            return binning
        elif meta_key == 'lampstat01':  # Dome lamps
            if headarr[0]['COMPMPOS'].strip().lower() == 'in':
                if headarr[0]['AGARCLMP'].strip() == 'W':
                    return "on"
                else:
                    return "off"
            else:
                return "off"
        elif meta_key == 'lampstat02':  # Arc lamps
            if headarr[0]['COMPMPOS'].strip().lower() == 'in':
                tst = headarr[0]['AGARCLMP'].strip()
                if tst != "W" and tst != "":  # TODO :: Need to update this... currently this just checks that the lamp list is not empty and the tungsten lamp is not on
                    return "on"
                else:
                    return "off"
            else:
                return "off"
        msgs.error(f"Not ready for this compound meta: {meta_key}")

    def configuration_keys(self):
        """
        Return the metadata keys that define a unique instrument
        configuration.

        This list is used by :class:`~pypeit.metadata.PypeItMetaData` to
        identify the unique configurations among the list of frames read
        for a given reduction.

        Returns:
            :obj:`list`: List of keywords of data pulled from file headers
            and used to constuct the :class:`~pypeit.metadata.PypeItMetaData`
            object.
        """
        # decker is not included because arcs are often taken with a 0.5" slit
        return ['dispname', 'dispangle', 'instrument']

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.

        Args:
            ftype (:obj:`str`):
                Type of frame to check. Must be a valid frame type; see
                frame-type :ref:`frame_type_defs`.
            fitstbl (`astropy.table.Table`_):
                The table with the metadata for one or more frames to check.
            exprng (:obj:`list`, optional):
                Range in the allowed exposure time for a frame of type
                ``ftype``. See
                :func:`pypeit.core.framematch.check_frame_exptime`.

        Returns:
            `numpy.ndarray`_: Boolean array with the flags selecting the
            exposures in ``fitstbl`` that are ``ftype`` type frames.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype in ['science', 'standard']:
            return good_exp & self.lamps(fitstbl, 'off')
        if ftype == 'bias':
            return good_exp
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            # Flats and trace frames are typed together
            return good_exp & self.lamps(fitstbl, 'dome')
        if ftype in ['pinhole', 'dark']:
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            return good_exp & self.lamps(fitstbl, 'arcs')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)
  
    def lamps(self, fitstbl, status):
        """
        Check the lamp status.

        Args:
            fitstbl (`astropy.table.Table`_):
                The table with the fits header meta data.
            status (:obj:`str`):
                The status to check. Can be ``'off'``, ``'arcs'``, or
                ``'dome'``.

        Returns:
            `numpy.ndarray`_: A boolean array selecting fits files that meet
            the selected lamp status.

        Raises:
            ValueError:
                Raised if the status is not one of the valid options.
        """
        if status == 'off':
            # Check if all are off
            return np.all(np.array([ (fitstbl[k] == 'off') for k in fitstbl.keys() if 'lampstat' in k]), axis=0)
        if status == 'dome':
            # Check if any dome lamps are on
            dome_lamp_stat = ['lampstat01']
            return np.any(np.array([ fitstbl[k] == 'on' for k in fitstbl.keys() if k in dome_lamp_stat]), axis=0)
        if status == 'arcs':
            # Check if any arc lamps are on
            arc_lamp_stat = ['lampstat02']
            return np.any(np.array([ fitstbl[k] == 'on' for k in fitstbl.keys() if k in arc_lamp_stat]), axis=0)
        raise ValueError('No implementation for status = {0}'.format(status))


class INTIDSBlueSpectrograph(INTIDSSpectrograph):
    """
    Child to handle INT/IDS acquired with the EEV10 detector
    """

    name = 'int_ids_eev'
    camera = 'EEV10'
    supported = True
    header_name = 'ids_eev'

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        Args:
            det (:obj:`int`):
                1-indexed detector number.
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.  If not
                provided, frame-dependent parameters are set to a default.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """
        # Detector 1
        detector_dict = dict(
            binning='1,1' if hdu is None else self.get_meta_value(self.get_headarr(hdu), 'binning'),
            det=1,
            dataext=0,
            specaxis=1,
            specflip=False,
            spatflip=False,
            platescale=0.43,
            saturation=65535.,
            mincounts=-1e10,
            nonlinear=0.76,
            numamplifiers=2,
            gain=np.asarray([1.2, 1.2]),
            ronoise=np.asarray([3.7, 3.7]),
            xgap=0.,
            ygap=0.,
            ysize=1.,
            darkcurr=0.0,
            # These are rows, columns on the raw frame, 1-indexed
            datasec=np.asarray(['[:, 1:1024]', '[:, 1025:2048]']),
            oscansec=np.asarray(['[:, 2050:2080]', '[:, 2081:2111]']),
        )
        return detector_container.DetectorContainer(**detector_dict)

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # par['flexure']['spectrum'] = 'sky_kastb_600.fits'
        # 1D wavelength solution
        par['calibrations']['wavelengths']['sigdetect'] = 5.
        par['calibrations']['wavelengths']['rms_threshold'] = 0.20
        par['calibrations']['wavelengths']['lamps'] = ['CdI','HgI','HeI']

        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['n_first'] = 3
        par['calibrations']['wavelengths']['match_toler'] = 2.5

        # Set wave tilts order
        par['calibrations']['tilts']['spat_order'] = 3
        par['calibrations']['tilts']['spec_order'] = 5
        par['calibrations']['tilts']['maxdev_tracefit'] = 0.02
        par['calibrations']['tilts']['maxdev2d'] = 0.02

        return par

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        Args:
            scifile (:obj:`str`):
                File to use when determining the configuration and how
                to adjust the input parameters.
            inp_par (:class:`~pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`~pypeit.par.parset.ParSet`: The PypeIt parameter set
            adjusted for configuration specific parameter values.
        """
        par = super().config_specific_par(scifile, inp_par=inp_par)
        # TODO: Should we allow the user to override these?

        if self.get_meta_value(scifile, 'dispname') == '600/4310':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'blah.fits'
        elif self.get_meta_value(scifile, 'dispname') == '452/3306':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'blah.fits'
        else:
            msgs.error("NEED TO ADD YOUR GRISM HERE!")
        # Return
        return par

    def raw_header_cards(self):
        """
        Return additional raw header cards to be propagated in
        downstream output files for configuration identification.

        The list of raw data FITS keywords should be those used to populate
        the :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.configuration_keys`
        or are used in :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.config_specific_par`
        for a particular spectrograph, if different from the name of the
        PypeIt metadata keyword.

        This list is used by :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.subheader_for_spec`
        to include additional FITS keywords in downstream output files.

        Returns:
            :obj:`list`: List of keywords from the raw data files that should
            be propagated in output files.
        """
        return ['GRISM_N', 'BSPLIT_N']


class INTIDSRedSpectrograph(INTIDSSpectrograph):
    """
    Child to handle INT/IDS data acquired with the default (RED+2) detector
    """

    name = 'int_ids_red'
    camera = 'RED+2'
    supported = False
    header_name = 'int_ids_red'

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        .. warning::

            Many of the necessary detector parameters are read from the file
            header, meaning the ``hdu`` argument is effectively **required** for
            INT/IDS (RED+2).  The optional use of ``hdu`` is only viable for
            automatically generated documentation.

        Args:
            det (:obj:`int`):
                1-indexed detector number.
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.  If not
                provided, frame-dependent parameters are set to a default.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """
        # Binning
        # TODO: Could this be detector dependent??
        binning = '1,1' if hdu is None else self.get_meta_value(self.get_headarr(hdu), 'binning')

        # Detector 1
        detector_dict = dict(
            binning         = binning,
            det             = 1,
            dataext         = 0,
            specaxis        = 0,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.44,  # arcsec/pixel
            darkcurr        = 8.0,   # e-/pix/hour
            saturation      = 65535.,
            nonlinear       = 0.95,  # < 1% deviation from linear with <60k counts
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d([0.91]),
            ronoise         = np.atleast_1d([3.48]),
            datasec         = None,
            oscansec        = None
            )

        if hdu is None:
            return detector_container.DetectorContainer(**detector_dict)

        # TODO: I don't know how to handle the stuff below for the
        # auto-generated detector table...

        # Parse datasec, oscancsec from the header
        header = hdu[0].header
        naxis1 = header['NAXIS1']
        crval1u = header['CRVAL1U']
        nover = header['COVER']

        ndata = naxis1 - nover*2

        x1_0 = 1             # Amp 1
        x1_1 = 512 - crval1u
        x2_0 = max(x1_1+1,1)       # Amp 2
        x2_1 = ndata

        xo1_1 = x2_1+1
        xo1_2 = x2_1+nover
        xo2_1 = xo1_2+1
        xo2_2 = xo1_2+nover

        # Allow for reading only Amp 2!
        if x1_1 < 3:
            msgs.warn("Only Amp 2 data was written.  Ignoring Amp 1")
            detector_dict['numamplifiers'] = 1
            detector_dict['gain'] = np.atleast_1d(detector_dict['gain'][0])
            detector_dict['ronoise'] = np.atleast_1d(detector_dict['ronoise'][0])
            # These are rows, columns on the raw frame, 1-indexed
            datasec = ['[:,{}:{}]'.format(x2_0,x2_1)]
            oscansec = ['[:,{}:{}]'.format(xo2_1,xo2_2)]
        else:
            # These are rows, columns on the raw frame, 1-indexed
            datasec = ['[:,{}:{}]'.format(x1_0, x1_1), 
                    '[:,{}:{}]'.format(x2_0,x2_1)]
            oscansec = ['[:,{}:{}]'.format(xo1_1,xo1_2), 
                        '[:,{}:{}]'.format(xo2_1,xo2_2)]

        # Fill it up
        detector_dict['datasec'] = np.atleast_1d(datasec)
        detector_dict['oscansec'] = np.atleast_1d(oscansec)

        return detector_container.DetectorContainer(**detector_dict)

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # 1D wavelength solution
        par['calibrations']['wavelengths']['lamps'] = ['NeI','HgI','HeI','ArI']
        return par

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        .. todo::
            Document the changes made!

        Args:
            scifile (str):
                File to use when determining the configuration and how
                to adjust the input parameters.
            inp_par (:class:`pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`pypeit.par.parset.ParSet`: The PypeIt paramter set
            adjusted for configuration specific parameter values.
        """
        par = self.default_pypeit_par() if inp_par is None else inp_par
        # TODO: Should we allow the user to override these?

        par['calibrations']['wavelengths']['method'] = 'holy-grail'
        par['calibrations']['wavelengths']['lamps'] = ['NeI', 'ArI', 'ArII', 'CuI']
        # if self.get_meta_value(scifile, 'dispname') == '300/7500':
        #     # TODO -- Note in docs that a NoNe solution is available too
        #     par['calibrations']['wavelengths']['method'] = 'full_template'
        #     par['calibrations']['wavelengths']['reid_arxiv'] = 'blah.fits'
        #     # Add CdI
        #     par['calibrations']['wavelengths']['lamps'] = ['NeI', 'HgI', 'HeI', 'ArI', 'CdI']
        # elif self.get_meta_value(scifile, 'dispname') == '1200/5000':
        #     par['calibrations']['wavelengths']['method'] = 'full_template'
        #     par['calibrations']['wavelengths']['reid_arxiv'] = 'blah.fits'
        #     par['calibrations']['wavelengths']['lamps'] = ['NeI', 'HgI', 'HeI', 'ArI', 'CdI']
        # else:
        #     par['calibrations']['wavelengths']['use_instr_flag'] = True

        # Return
        return par
