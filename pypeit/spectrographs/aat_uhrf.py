"""
Module for Shane/Kast specific methods.

.. include:: ../include/links.rst
"""
import os

from IPython import embed

import numpy as np

from astropy.time import Time

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container
from pypeit import data



class AATUHRFSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle AAT/UHRF specific code
    """
    ndet = 1
    telescope = telescopes.AATTelescopePar()
    url = 'https://aat.anu.edu.au/science/instruments/decomissioned/uhrf/overview'
    ql_supported = False
    name = 'aat_uhrf'
    camera = 'UHRF'
    supported = True
    header_name = 'uhrf'

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
            darkcurr=0.0,  # e-/pixel/hour
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
        par['sensfunc']['IR']['telgridfile'] = 'TellPCA_3000_26000_R10000.fits'
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
        # dispname is arm specific (blue/red)
        self.meta['decker'] = dict(ext=0, card='SLIT_N')
        self.meta['binning'] = dict(ext=0, card=None, default='1,1')
        self.meta['mjd'] = dict(ext=0, card=None, compound=True)
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Additional ones, generally for configuration determination or time
        self.meta['dichroic'] = dict(ext=0, card='BSPLIT_N')
        self.meta['instrument'] = dict(ext=0, card='VERSION')
        lamp_names = [ '1', '2', '3', '4', '5',
                       'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K']
        for kk,lamp_name in enumerate(lamp_names):
            self.meta['lampstat{:02d}'.format(kk+1)] = dict(ext=0, card='LAMPSTA{0}'.format(lamp_name))

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
        if meta_key == 'mjd':
            time = headarr[0]['DATE']
            ttime = Time(time, format='isot')
            return ttime.mjd
        msgs.error("Not ready for this compound meta")

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
        return ['dispname', 'dichroic' ]

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
            return good_exp # & (fitstbl['target'] == 'Bias')
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            # Flats and trace frames are typed together
            return good_exp & self.lamps(fitstbl, 'dome') # & (fitstbl['target'] == 'Dome Flat')
        if ftype in ['pinhole', 'dark']:
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            return good_exp & self.lamps(fitstbl, 'arcs')#  & (fitstbl['target'] == 'Arcs')

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
            return np.all(np.array([ (fitstbl[k] == 'off') | (fitstbl[k] == 'None')
                                        for k in fitstbl.keys() if 'lampstat' in k]), axis=0)
        if status == 'arcs':
            # Check if any arc lamps are on
            arc_lamp_stat = [ 'lampstat{0:02d}'.format(i) for i in range(6,17) ]
            return np.any(np.array([ fitstbl[k] == 'on' for k in fitstbl.keys()
                                            if k in arc_lamp_stat]), axis=0)
        if status == 'dome':
            # Check if any dome lamps are on
            dome_lamp_stat = [ 'lampstat{0:02d}'.format(i) for i in range(1,6) ]
            return np.any(np.array([ fitstbl[k] == 'on' for k in fitstbl.keys()
                                            if k in dome_lamp_stat]), axis=0)
        raise ValueError('No implementation for status = {0}'.format(status))

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
            par['calibrations']['wavelengths']['reid_arxiv'] = 'shane_kast_blue_600.fits'
        elif self.get_meta_value(scifile, 'dispname') == '452/3306':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'shane_kast_blue_452.fits'
        elif self.get_meta_value(scifile, 'dispname') == '830/3460':  # NOT YET TESTED
            par['calibrations']['wavelengths']['reid_arxiv'] = 'shane_kast_blue_830.fits'
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
