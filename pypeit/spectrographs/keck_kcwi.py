"""
Implements KCWI-specific functions.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from IPython import embed

import numpy as np

from astropy import wcs, units
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import EarthLocation
from scipy.optimize import curve_fit

from pypeit import msgs
from pypeit import telescopes
from pypeit import utils
from pypeit import io
from pypeit.core import parse
from pypeit.core import procimg
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container


class KeckKCWIKCRMSpectrograph(spectrograph.Spectrograph):
    """
    Parent to handle Keck/KCWI+KCRM specific code

    .. todo::
        * Need to apply spectral flexure and heliocentric correction to waveimg -- done?
        * Copy fast_histogram code into PypeIt?
        * Re-write flexure code with datamodel + implement spectral flexure QA in find_objects.py
        * When making the datacube, add an option to apply a spectral flexure correction from a different frame?
        * Write some detailed docs about the corrections that can be used when making a datacube
        * Consider introducing a new method (par['flexure']['spec_method']) for IFU flexure corrections (see find-objects.py)
    """
    ndet = 1
    telescope = telescopes.KeckTelescopePar()
    pypeline = 'SlicerIFU'
    supported = True

    def __init__(self):
        super().__init__()

        # TODO :: Might need to change the tolerance of disperser angle in
        # pypeit setup (two BH2 nights where sufficiently different that this
        # was important).

        # TODO :: Might consider changing TelescopePar to use the astropy
        # EarthLocation. KBW: Fine with me!
        self.location = EarthLocation.of_site('Keck Observatory')

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}

        # Required (core)
        self.meta['ra'] = dict(ext=0, card=None, compound=True)
        self.meta['dec'] = dict(ext=0, card=None, compound=True)
        self.meta['target'] = dict(ext=0, card='TARGNAME')
        self.meta['decker'] = dict(ext=0, card='IFUNAM')
        self.meta['binning'] = dict(card=None, compound=True)

        self.meta['mjd'] = dict(ext=0, card='MJD')
        self.meta['exptime'] = dict(card=None, compound=True)
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        self.meta['posang'] = dict(card=None, compound=True)
        self.meta['ra_off'] = dict(ext=0, card='RAOFF')
        self.meta['dec_off'] = dict(ext=0, card='DECOFF')


        # Extras for config and frametyping
        self.meta['hatch'] = dict(ext=0, card='HATPOS')
        #        self.meta['idname'] = dict(ext=0, card='CALXPOS')
        self.meta['idname'] = dict(ext=0, card='IMTYPE')
        self.meta['calpos'] = dict(ext=0, card='CALMNAM')
        self.meta['slitwid'] = dict(card=None, compound=True)

        # Get atmospheric conditions (note, these are the conditions at the end of the exposure)
        self.meta['obstime'] = dict(card=None, compound=True, required=False)
        self.meta['pressure'] = dict(card=None, compound=True, required=False)
        self.meta['temperature'] = dict(card=None, compound=True, required=False)
        self.meta['humidity'] = dict(card=None, compound=True, required=False)
        self.meta['parangle'] = dict(card=None, compound=True, required=False)
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')

        # Lamps
        lamp_names = ['LMP0', 'LMP1', 'LMP2', 'LMP3']  # FeAr, ThAr, Aux, Continuum
        for kk, lamp_name in enumerate(lamp_names):
            self.meta['lampstat{:02d}'.format(kk + 1)] = dict(ext=0, card=lamp_name + 'STAT')
        for kk, lamp_name in enumerate(lamp_names):
            if lamp_name == 'LMP3':
                # There is no shutter on LMP3
                self.meta['lampshst{:02d}'.format(kk + 1)] = dict(ext=0, card=None, default=1)
                continue
            self.meta['lampshst{:02d}'.format(kk + 1)] = dict(ext=0, card=lamp_name + 'SHST')
        # Add in the dome lamp
        self.meta['lampstat{:02d}'.format(len(lamp_names) + 1)] = dict(ext=0, card='FLSPECTR')
        self.meta['lampshst{:02d}'.format(len(lamp_names) + 1)] = dict(ext=0, card=None, default=1)


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

        headarr = self.get_headarr(scifile)

        # Templates
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['lamps'] = ['FeI', 'ArI', 'ArII']
        if self.get_meta_value(headarr, 'dispname') == 'BH2':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_kcwi_BH2.fits'
        elif self.get_meta_value(headarr, 'dispname') == 'BM':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_kcwi_BM.fits'
        elif self.get_meta_value(headarr, 'dispname') == 'BL':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_kcwi_BL.fits'
        elif self.get_meta_value(headarr, 'dispname') == 'RL':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_kcrm_RL.fits'
        elif self.get_meta_value(headarr, 'dispname') == 'RM1':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_kcrm_RM1.fits'
        elif self.get_meta_value(headarr, 'dispname') == 'RM2':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_kcrm_RM2.fits'
        elif self.get_meta_value(headarr, 'dispname') == 'RH3':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_kcrm_RH3.fits'
        else:
            msgs.warn("Full template solution is unavailable")
            msgs.info("Adopting holy-grail algorithm - Check the wavelength solution!")
            par['calibrations']['wavelengths']['method'] = 'holy-grail'
        # FWHM
        # binning = parse.parse_binning(self.get_meta_value(headarr, 'binning'))
        # par['calibrations']['wavelengths']['fwhm'] = 6.0 / binning[1]

        # Return
        return par

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
        return ['dispname', 'decker', 'binning', 'cenwave']

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
            binspatial, binspec = parse.parse_binning(headarr[0]['BINNING'])
            binning = parse.binning2string(binspec, binspatial)
            return binning
        elif meta_key == 'exptime':
            try:
                return headarr[0]['ELAPTIME']
            except KeyError:
                return headarr[0]['TELAPSE']
        elif meta_key == 'slitwid':
            # Get the slice scale
            slicescale = 0.00037718  # Degrees per 'large slicer' slice
            ifunum = headarr[0]['IFUNUM']
            if ifunum == 2:
                slicescale /= 2.0
            elif ifunum == 3:
                slicescale /= 4.0
            return slicescale
        elif meta_key == 'ra' or meta_key == 'dec':
            try:
                if self.is_nasmask(headarr[0]):
                    hdrstr = 'RABASE' if meta_key == 'ra' else 'DECBASE'
                else:
                    hdrstr = 'RA' if meta_key == 'ra' else 'DEC'
            except KeyError:
                try:
                    hdrstr = 'TARGRA' if meta_key == 'ra' else 'TARGDEC'
                except KeyError:
                    msgs.error(f'Cannot determine the {meta_key} from the header')
            return headarr[0][hdrstr]
        elif meta_key == 'pressure':
            try:
                return headarr[0]['WXPRESS']  # Must be in astropy.units.mbar
            except KeyError:
                msgs.warn("Pressure is not in header")
                msgs.info("The default pressure will be assumed: 611 mbar")
                return 611.0
        elif meta_key == 'temperature':
            try:
                return headarr[0]['WXOUTTMP']  # Must be in astropy.units.deg_C
            except KeyError:
                msgs.warn("Temperature is not in header")
                msgs.info("The default temperature will be assumed: 1.5 deg C")
                return 1.5  # van Kooten & Izett, arXiv:2208.11794
        elif meta_key == 'humidity':
            try:
                # Humidity expressed as a percentage, not a fraction
                return headarr[0]['WXOUTHUM']
            except KeyError:
                msgs.warn("Humidity is not in header")
                msgs.info("The default relative humidity will be assumed: 20 %")
                return 20.0  # van Kooten & Izett, arXiv:2208.11794
        elif meta_key == 'parangle':
            try:
                # Parallactic angle expressed in radians
                return headarr[0]['PARANG'] * np.pi / 180.0
            except KeyError:
                msgs.error("Parallactic angle is not in header")
        elif meta_key == 'obstime':
            return Time(headarr[0]['DATE-END'])
        elif meta_key == 'posang':
            hdr = headarr[0]
            # Get rotator position
            if 'ROTPOSN' in hdr:
                rpos = hdr['ROTPOSN']
            else:
                rpos = 0.
            if 'ROTREFAN' in hdr:
                rref = hdr['ROTREFAN']
            else:
                rref = 0.
            # Get the offset and PA
            skypa = rpos + rref  # IFU position angle (degrees)
            return skypa
        else:
            msgs.error("Not ready for this compound meta")


    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 0.001]
        par['calibrations']['darkframe']['exprng'] = [0.01, None]

        # Set the number of alignments in the align frames
        par['calibrations']['alignment']['locations'] = [0.1, 0.3, 0.5, 0.7, 0.9]  # TODO:: Check this - is this accurate enough?

        # LACosmics parameters
        par['scienceframe']['process']['sigclip'] = 4.0
        par['scienceframe']['process']['objlim'] = 1.5
        # Illumination corrections
        par['scienceframe']['process']['use_illumflat'] = True  # illumflat is applied when building the relative scale image in reduce.py, so should be applied to scienceframe too.
        par['scienceframe']['process']['use_specillum'] = True  # apply relative spectral illumination
        par['scienceframe']['process']['spat_flexure_correct'] = False  # don't correct for spatial flexure - varying spatial illumination profile could throw this correction off. Also, there's no way to do astrometric correction if we can't correct for spatial flexure of the contbars frames
        par['scienceframe']['process']['use_biasimage'] = True  # Need to use bias frames for KCWI, because the bias level varies monotonically with spatial and spectral direction
        par['scienceframe']['process']['use_darkimage'] = False

        # Don't do 1D extraction for 3D data - it's meaningless because the DAR correction must be performed on the 3D data.
        par['reduce']['extraction']['skip_extraction'] = True  # Because extraction occurs before the DAR correction, don't extract

        # Make sure that this is reduced as a slit (as opposed to fiber) spectrograph
        par['reduce']['cube']['slit_spec'] = True
        par['reduce']['cube']['combine'] = False  # Make separate spec3d files from the input spec2d files

        # Sky subtraction parameters
        par['reduce']['skysub']['no_poly'] = False
        par['reduce']['skysub']['bspline_spacing'] = 0.6
        par['reduce']['skysub']['joint_fit'] = False

        # Don't correct flexure by default, but you should use slitcen,
        # because this is a slit-based IFU where no objects are extracted.
        par['flexure']['spec_method'] = 'skip'
        par['flexure']['spec_maxshift'] = 3  # Just in case someone switches on spectral flexure, this needs to be minimal

        # Flux calibration parameters
        par['sensfunc']['UVIS']['extinct_correct'] = False  # This must be False - the extinction correction is performed when making the datacube

        # If telluric is triggered
        par['sensfunc']['IR']['telgridfile'] = 'TellPCA_3000_26000_R15000.fits'

        return par

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard PypeIt file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
        return super().pypeit_file_keys() + ['posang', 'ra_off', 'dec_off', 'idname', 'calpos']

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
        if ftype == 'science':
            return good_exp & (fitstbl['idname'] == 'OBJECT') & (fitstbl['calpos'] == 'Sky') \
                    & self.lamps(fitstbl, 'off') & (fitstbl['hatch'] == 'Open')
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'BIAS')
        if ftype in ['pixelflat', 'scattlight']:  # Scattered light needs lots of counts - so set it to the pixelflat by default
            # Use internal lamp
            return good_exp & (fitstbl['idname'] == 'FLATLAMP') & (fitstbl['calpos'] == 'Mirror') \
                    & self.lamps(fitstbl, 'cont_noarc')
        if ftype in ['illumflat', 'trace']:
            # Use dome flats
            return good_exp & (fitstbl['idname'] == 'DOMEFLAT') & (fitstbl['calpos'] == 'Sky') \
                    & self.lamps(fitstbl, 'dome_noarc') & (fitstbl['hatch'] == 'Open')
        if ftype == 'dark':
            # Dark frames
            return good_exp & (fitstbl['idname'] == 'DARK') & self.lamps(fitstbl, 'off') \
                    & (fitstbl['hatch'] == 'Closed')
        if ftype == 'align':
            # Alignment frames
            # NOTE: Different from previous versions, this now only warns the user if everyth
            is_align = good_exp & (fitstbl['idname'] == 'CONTBARS') \
                        & (fitstbl['calpos'] == 'Mirror') & self.lamps(fitstbl, 'cont')
            if np.any(is_align & np.logical_not(self.lamps(fitstbl, 'cont_noarc'))):
                msgs.warn('Alignment frames have both the continuum and arc lamps on (although '
                          'arc-lamp shutter might be closed)!')
            return is_align
        if ftype == 'arc':
            # PypeIt is only setup to wavelength calibrate using the FeAr lamp.
            return good_exp & (fitstbl['idname'] == 'ARCLAMP') & (fitstbl['calpos'] == 'Mirror') \
                    & self.lamps(fitstbl, 'arcs')
        if ftype == 'tilt':
            # Check to see if ThAr frames are available. If so, use them. Otherwise, use the FeAr lamp.
            tmp = good_exp & (fitstbl['idname'] == 'ARCLAMP') & (fitstbl['calpos'] == 'Mirror') \
                    & self.lamps(fitstbl, 'tilts')
            if np.any(tmp):
                return tmp
            # No ThAr frames, so use the FeAr lamp
            return good_exp & (fitstbl['idname'] == 'ARCLAMP') & (fitstbl['calpos'] == 'Mirror') \
                   & self.lamps(fitstbl, 'arcs')
        if ftype == 'pinhole':
            # Don't type pinhole frames
            return np.zeros(len(fitstbl), dtype=bool)

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
            lampstat = np.array([np.isin(fitstbl[k], ['0', 'None', 'off'])
                                    for k in fitstbl.keys() if 'lampstat' in k])
            return np.all(lampstat, axis=0)  # Lamp has to be off
        if status in ['arcs', 'tilts']:
            # Check if any arc/tilt lamps are on
            if status == 'arcs':
                # Only FeAr is used for wavelength calibration
                arc_lamp_stat = ['lampstat{0:02d}'.format(1)]
                arc_lamp_shst = ['lampshst{0:02d}'.format(1)]
            else:
                # Only ThAr is used to calculate the tilt
                arc_lamp_stat = ['lampstat{0:02d}'.format(2)]
                arc_lamp_shst = ['lampshst{0:02d}'.format(2)]
            lamp_stat = np.array([fitstbl[k] == '1' for k in fitstbl.keys() if k in arc_lamp_stat])
            lamp_shst = np.array([fitstbl[k] == '1' for k in fitstbl.keys() if k in arc_lamp_shst])
            # Make sure the continuum frames are off
            dome_lamps = ['lampstat{0:02d}'.format(i) for i in range(4, 5)]
            dome_lamp_stat = np.array([fitstbl[k] == '0' for k in fitstbl.keys()
                                       if k in dome_lamps])
            return np.any(lamp_stat & lamp_shst & dome_lamp_stat, axis=0)  # i.e. lamp on and shutter open
        if status in ['cont_noarc', 'cont']:
            # Check if any internal lamps are on (Continuum) - Ignore lampstat03 (Aux) - not sure what this is used for
            cont_lamp_stat = ['lampstat{0:02d}'.format(4)]
            lamp_stat = np.array([fitstbl[k] == '1' for k in fitstbl.keys()
                                  if k in cont_lamp_stat])
            if status == 'cont_noarc':
                # Make sure arcs are off - it seems even with the shutter closed, the arcs
                arc_lamps = ['lampstat{0:02d}'.format(i) for i in range(1, 3)]
                arc_lamp_stat = np.array([fitstbl[k] == '0' for k in fitstbl.keys()
                                        if k in arc_lamps])
                lamp_stat = lamp_stat & arc_lamp_stat
            return np.any(lamp_stat, axis=0)  # i.e. lamp on
        if status in ['dome_noarc', 'dome']:
            # Check if any dome lamps are on (Continuum) - Ignore lampstat03 (Aux) - not sure what this is used for
            dome_lamp_stat = ['lampstat{0:02d}'.format(5)]
            lamp_stat = np.array([fitstbl[k] == 'on' for k in fitstbl.keys()
                                  if k in dome_lamp_stat])
            if status == 'dome_noarc':
                # Make sure arcs are off - it seems even with the shutter closed, the arcs
                arc_lamps = ['lampstat{0:02d}'.format(i) for i in range(1, 3)]
                arc_lamp_stat = np.array([fitstbl[k] == '0' for k in fitstbl.keys()
                                          if k in arc_lamps])
                lamp_stat = lamp_stat & arc_lamp_stat
            return np.any(lamp_stat, axis=0)  # i.e. lamp on
        raise ValueError('No implementation for status = {0}'.format(status))

    def get_lamps_status(self, headarr):
        """
        Return a string containing the information on the lamp status.

        Args:
            headarr (:obj:`list`):
                A list of 1 or more `astropy.io.fits.Header`_ objects.

        Returns:
            :obj:`str`: A string that uniquely represents the lamp status.
        """
        # Loop through all lamps and collect their status
        kk = 1
        lampstat = []
        while True:
            lampkey1 = 'lampstat{:02d}'.format(kk)
            if lampkey1 not in self.meta.keys():
                break
            ext1, card1 = self.meta[lampkey1]['ext'], self.meta[lampkey1]['card']
            lampkey2 = 'lampshst{:02d}'.format(kk)
            if self.meta[lampkey2]['card'] is None:
                lampstat += [str(headarr[ext1][card1])]
            else:
                ext2, card2 = self.meta[lampkey2]['ext'], self.meta[lampkey2]['card']
                lampstat += ["{0:s}-{1:s}".format(str(headarr[ext1][card1]), str(headarr[ext2][card2]))]
            kk += 1
        return "_".join(lampstat)


    def calc_pattern_freq(self, frame, rawdatasec_img, oscansec_img, hdu):
        """
        Calculate the pattern frequency using the overscan region that covers
        the overscan and data sections. Using a larger range allows the
        frequency to be pinned down with high accuracy.

        NOTE: The amplifiers are arranged as follows:

        |   (0,ny)  --------- (nx,ny)
        |           | 3 | 4 |
        |           ---------
        |           | 1 | 2 |
        |     (0,0) --------- (nx, 0)

        .. todo::

            PATTERN FREQUENCY ALGORITHM HAS NOT BEEN TESTED WHEN BINNING != 1x1

        Parameters
        ----------
        frame : `numpy.ndarray`_
            Raw data frame to be used to estimate the pattern frequency.
        rawdatasec_img : `numpy.ndarray`_
            Array the same shape as ``frame``, used as a mask to identify the
            data pixels (0 is no data, non-zero values indicate the amplifier
            number).
        oscansec_img : `numpy.ndarray`_
            Array the same shape as ``frame``, used as a mask to identify the
            overscan pixels (0 is no data, non-zero values indicate the
            amplifier number).
        hdu : `astropy.io.fits.HDUList`_
            Opened fits file.

        Returns
        -------
        patt_freqs : :obj:`list`
            List of pattern frequencies.
        """
        msgs.info("Calculating pattern noise frequency")

        # Make a copy of te original frame
        raw_img = frame.copy()

        # Get a unique list of the amplifiers
        unq_amps = np.sort(np.unique(oscansec_img[np.where(oscansec_img >= 1)]))
        num_amps = unq_amps.size

        # Loop through amplifiers and calculate the frequency
        patt_freqs = []
        for amp in unq_amps:
            # Grab the pixels where the amplifier has data
            pixs = np.where((rawdatasec_img == amp) | (oscansec_img == amp))
            rmin, rmax = np.min(pixs[1]), np.max(pixs[1])
            # Deal with the different locations of the overscan regions in 2- and 4- amp mode
            if num_amps == 2:
                cmin = 1+np.max(pixs[0])
                frame = raw_img[cmin:, rmin:rmax].astype(float)
            elif num_amps == 4:
                if amp in [1, 2]:
                    pixalt = np.where((rawdatasec_img == amp+2) | (oscansec_img == amp+2))
                    cmin = 1+np.max(pixs[0])
                    cmax = (np.min(pixalt[0]) + cmin)//2  # Average of the bottom of the top amp, and top of the bottom amp
                else:
                    pixalt = np.where((rawdatasec_img == amp-2) | (oscansec_img == amp-2))
                    cmax = 1+np.min(pixs[0])
                    cmin = (np.max(pixalt[0]) + cmax)//2
                frame = raw_img[cmin:cmax, rmin:rmax].astype(float)
            # Calculate the pattern frequency
            freq = procimg.pattern_frequency(frame)
            patt_freqs.append(freq)
            msgs.info("Pattern frequency of amplifier {0:d}/{1:d} = {2:f}".format(amp, num_amps, freq))

        # Return the list of pattern frequencies
        return patt_freqs

    def get_wcs(self, hdr, slits, platescale, wave0, dwv, spatial_scale=None):
        """
        Construct/Read a World-Coordinate System for a frame.

        Args:
            hdr (`astropy.io.fits.Header`_):
                The header of the raw frame. The information in this
                header will be extracted and returned as a WCS.
            slits (:class:`~pypeit.slittrace.SlitTraceSet`):
                Slit traces.
            platescale (:obj:`float`):
                The platescale of an unbinned pixel in arcsec/pixel (e.g.
                detector.platescale). See also 'spatial_scale'
            wave0 (:obj:`float`):
                The wavelength zeropoint.
            dwv (:obj:`float`):
                Change in wavelength per spectral pixel.
            spatial_scale (:obj:`float`, None, optional):
                The spatial scale (units=arcsec/pixel) of the WCS to be used.
                This variable is fixed, and is independent of the binning.
                If spatial_scale is set, it will be used for the spatial size
                of the WCS and the platescale will be ignored. If None, then
                the platescale will be used.

        Returns:
            `astropy.wcs.WCS`_: The world-coordinate system.
        """
        msgs.info(f"Generating {self.camera} WCS")
        # Get the x and y binning factors, and the typical slit length
        binspec, binspat = parse.parse_binning(self.get_meta_value([hdr], 'binning'))

        # Get the pixel and slice scales
        pxscl = platescale * binspat / 3600.0  # 3600 is to convert arcsec to degrees
        slscl = self.get_meta_value([hdr], 'slitwid')
        if spatial_scale is not None:
            if pxscl > spatial_scale / 3600.0:
                msgs.warn("Spatial scale requested ({0:f}'') is less than the pixel scale ({1:f}'')".format(spatial_scale, pxscl*3600.0))
            # Update the pixel scale
            pxscl = spatial_scale / 3600.0  # 3600 is to convert arcsec to degrees

        # Get the typical slit length (this changes by ~0.3% over all slits, so a constant is fine for now)
        slitlength = int(np.round(np.median(slits.get_slitlengths(initial=True, median=True))))

        # Get RA/DEC
        ra = self.compound_meta([hdr], 'ra')
        dec = self.compound_meta([hdr], 'dec')

        skypa = self.compound_meta([hdr], 'posang')
        rotoff = 0.0  # IFU-SKYPA offset (degrees)
        crota = np.radians(-(skypa + rotoff))

        # Calculate the fits coordinates
        cdelt1 = -slscl
        cdelt2 = pxscl

        # Calculate the CD Matrix
        cd11 = cdelt1 * np.cos(crota)                          # RA degrees per column
        cd12 = abs(cdelt2) * np.sign(cdelt1) * np.sin(crota)   # RA degrees per row
        cd21 = -abs(cdelt1) * np.sign(cdelt2) * np.sin(crota)  # DEC degress per column
        cd22 = cdelt2 * np.cos(crota)                          # DEC degrees per row
        # Get reference pixels (set these to the middle of the FOV)
        crpix1 = 24/2   # i.e. 24 slices/2
        crpix2 = slitlength / 2.
        crpix3 = 1.
        # Get the offset
        porg = hdr['PONAME']
        ifunum = hdr['IFUNUM']
        if 'IFU' in porg:
            # if ifunum == 1:  # Large slicer
            #     off1 = 1.0
            #     off2 = 4.0
            # elif ifunum == 2:  # Medium slicer
            #     off1 = 1.0
            #     off2 = 5.0
            # elif ifunum == 3:  # Small slicer
            #     off1 = 0.05
            #     off2 = 5.6
            # else:
            #     msgs.warn("Unknown IFU number: {0:d}".format(ifunum))
            off1 = 0.
            off2 = 0.
            off1 /= binspec
            off2 /= binspat
            crpix1 += off1
            crpix2 += off2

        # Create a new WCS object.
        w = wcs.WCS(naxis=3)
        w.wcs.equinox = hdr['EQUINOX']
        w.wcs.name = self.camera
        w.wcs.radesys = 'FK5'
        w.wcs.lonpole = 180.0  # Native longitude of the Celestial pole
        w.wcs.latpole = 0.0  # Native latitude of the Celestial pole
        # Insert the coordinate frame
        w.wcs.cname = ['RA', 'DEC', 'Wavelength']
        w.wcs.cunit = [units.degree, units.degree, units.Angstrom]
        w.wcs.ctype = ["RA---TAN", "DEC--TAN", "WAVE"]  # Note, WAVE is vacuum wavelength
        w.wcs.crval = [ra, dec, wave0]  # RA, DEC, and wavelength zeropoints
        w.wcs.crpix = [crpix1, crpix2, crpix3]  # RA, DEC, and wavelength reference pixels
        w.wcs.cd = np.array([[cd11, cd12, 0.0], [cd21, cd22, 0.0], [0.0, 0.0, dwv]])
        return w

    def get_datacube_bins(self, slitlength, minmax, num_wave):
        r"""
        Calculate the bin edges to be used when making a datacube.

        Args:
            slitlength (:obj:`int`):
                Length of the slit in pixels
            minmax (`numpy.ndarray`_):
                An array with the minimum and maximum pixel locations on each
                slit relative to the reference location (usually the centre
                of the slit). Shape must be :math:`(N_{\rm slits},2)`, and is
                typically the array returned by
                :func:`~pypeit.slittrace.SlitTraceSet.get_radec_image`.
            num_wave (:obj:`int`):
                Number of wavelength steps.  Given by::
                    int(round((wavemax-wavemin)/delta_wave))

        Args:
            :obj:`tuple`: Three 1D `numpy.ndarray`_ providing the bins to use
            when constructing a histogram of the spec2d files. The elements
            are :math:`(x,y,\lambda)`.
        """
        xbins = np.arange(1 + 24) - 24/2 - 0.5
        ybins = np.linspace(np.min(minmax[:, 0]), np.max(minmax[:, 1]), 1+slitlength) - 0.5
        spec_bins = np.arange(1+num_wave) - 0.5
        return xbins, ybins, spec_bins

    def bpm(self, filename, det, shape=None, msbias=None):
        """
        Generate a default bad-pixel mask for KCWI and KCRM.

        Even though they are both optional, either the precise shape for
        the image (``shape``) or an example file that can be read to get
        the shape (``filename`` using :func:`get_image_shape`) *must* be
        provided.

        Args:
            filename (:obj:`str` or None):
                An example file to use to get the image shape.
            det (:obj:`int`):
                1-indexed detector number to use when getting the image
                shape from the example file.
            shape (tuple, optional):
                Processed image shape
                Required if filename is None
                Ignored if filename is not None
            msbias (`numpy.ndarray`_, optional):
                Processed bias frame used to identify bad pixels. **This is
                ignored for KCWI.**

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """
        # Call the base-class method to generate the empty bpm; msbias is always set to None.
        bpm_img = super().bpm(filename, det, shape=shape, msbias=None)

        # Extract some header info
        head0 = fits.getheader(filename, ext=0)
        ampmode = head0['AMPMODE']
        binning = head0['BINNING']

        # Construct a list of the bad columns
        # KCWI --> AMPMODE = 'ALL', 'TBO', 'TUP'
        # KCRM --> AMPMODE = 'L2U2', 'L2U2L1U1'
        bc = None
        if ampmode == 'ALL':
            # TODO: There are several bad columns in this mode, but this is typically only used for arcs.
            #       It's the same set of bad columns seen in the TBO and TUP amplifier modes.
            if binning == '1,1':
                bc = [[3676, 3676, 2056, 2244]]
            elif binning == '2,2':
                bc = [[1838, 1838, 1028, 1121]]
        elif ampmode == 'TBO':
            if binning == '1,1':
                bc = [[2622, 2622,  619,  687],
                      [2739, 2739, 1748, 1860],
                      [3295, 3300, 2556, 2560],
                      [3675, 3676, 2243, 4111]]
            elif binning == '2,2':
                bc = [[1311, 1311,  310,  354],
                      [1369, 1369,  876,  947],
                      [1646, 1650, 1278, 1280],
                      [1838, 1838, 1122, 2055]]
        elif ampmode == 'TUP':
            if binning == '1,1':
#                bc = [[2622, 2622, 3492, 3528],
                bc = [[2622, 2622, 3492, 4111],   # Extending this BPM, as sometimes the bad column is larger than this.
                      [3295, 3300, 1550, 1555],
                      [3676, 3676, 1866, 4111]]
            elif binning == '2,2':
#                bc = [[1311, 1311, 1745, 1788],
                bc = [[1311, 1311, 1745, 2055],   # Extending this BPM, as sometimes the bad column is larger than this.
                      [1646, 1650,  775,  777],
                      [1838, 1838,  933, 2055]]
        elif ampmode == 'L2U2':
            if binning == '1,1':
                bc = [[649, 651, 0, 613]]  # This accounts for the spatflip - not sure if the 649-651 is too broad though...
            elif binning == '2,2':
                bc = [[325, 325, 0, 307]]  # This accounts for the spatflip
        elif ampmode == "L2U2L1U1":
            pass
            # Currently unchecked...
            # if binning == '1,1':
            #     bc = [[3460, 3460, 2064, 3520]]
            # elif binning == '2,2':
            #     bc = [[1838, 1838, 1028, 1121]]

        # Check if the bad columns haven't been set
        if bc is None:
            msgs.warn("KCRM bad pixel mask is not available for ampmode={0:s} binning={1:s}".format(ampmode, binning))
            bc = []

        # Apply these bad columns to the mask
        for bb in range(len(bc)):
            bpm_img[bc[bb][2]:bc[bb][3]+1, bc[bb][0]:bc[bb][1]+1] = 1

        return np.flipud(bpm_img)


class KeckKCWISpectrograph(KeckKCWIKCRMSpectrograph):
    """
    Child to handle Keck/KCWI specific code
    """
    name = 'keck_kcwi'
    camera = 'KCWI'
    url = 'https://www2.keck.hawaii.edu/inst/kcwi/'
    header_name = 'KCWI'
    comment = 'Supported setups: BL, BM, BH2; see :doc:`keck_kcwi`'

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        .. warning::

            Many of the necessary detector parameters are read from the file
            header, meaning the ``hdu`` argument is effectively **required** for
            KCWI.  The optional use of ``hdu`` is only viable for automatically
            generated documentation.

        Args:
            det (:obj:`int`):
                1-indexed detector number.
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """
        if hdu is None:
            binning = '2,2'
            specflip = None
            numamps = None
            gainarr = None
            ronarr = None
#            dsecarr = None
#            msgs.error("A required keyword argument (hdu) was not supplied")
        else:
            # Some properties of the image
            binning = self.compound_meta(self.get_headarr(hdu), "binning")
            numamps = hdu[0].header['NVIDINP']
            specflip = True if hdu[0].header['AMPID1'] == 2 else False
            gainmul, gainarr = hdu[0].header['GAINMUL'], np.zeros(numamps)
            ronarr = np.zeros(numamps)  # Set this to zero (determine the readout noise from the overscan regions)
#            dsecarr = np.array(['']*numamps)

            for ii in range(numamps):
                # Assign the gain for this amplifier
                gainarr[ii] = hdu[0].header["GAIN{0:1d}".format(ii + 1)]# * gainmul

        detector = dict(det             = det,
                        binning         = binning,
                        dataext         = 0,
                        specaxis        = 0,
                        specflip        = specflip,
                        spatflip        = False,
                        platescale      = 0.145728,  # arcsec/pixel
                        darkcurr        = 1.0,  # e-/hour/unbinned pixel
                        mincounts       = -1e10,
                        saturation      = 65535.,
                        nonlinear       = 0.95,       # For lack of a better number!
                        numamplifiers   = numamps,
                        gain            = gainarr,
                        ronoise         = ronarr,
                        # These are never used because the image reader sets these up using the file headers data.
                        # datasec         = dsecarr, #.copy(),     # <-- This is provided in the header
                        # oscansec        = dsecarr, #.copy(),     # <-- This is provided in the header
                        )
        # Return
        return detector_container.DetectorContainer(**detector)

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        super().init_meta()
        self.meta['dispname'] = dict(ext=0, card='BGRATNAM')
        self.meta['dispangle'] = dict(ext=0, card='BGRANGLE', rtol=0.01)
        self.meta['cenwave'] = dict(ext=0, card='BCWAVE', rtol=0.01)

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
        return ['BGRATNAM', 'IFUNAM', 'BGRANGLE']

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # Subtract the detector pattern from certain frames.
        # NOTE: The pattern subtraction is time-consuming, meaning we don't
        # perform it (by default) for the high S/N pixel flat images but we do
        # for everything else.
        par['calibrations']['biasframe']['process']['use_pattern'] = True
        par['calibrations']['darkframe']['process']['use_pattern'] = True
        par['calibrations']['pixelflatframe']['process']['use_pattern'] = False
        par['calibrations']['illumflatframe']['process']['use_pattern'] = True
        par['calibrations']['standardframe']['process']['use_pattern'] = True
        par['scienceframe']['process']['use_pattern'] = True
        # Subtract scattered light, but only for the pixel and illum flats,
        # as well as and science/standard star data.
        par['calibrations']['scattlight_pad'] = 6  # This is the unbinned number of pixels to pad
        par['calibrations']['pixelflatframe']['process']['subtract_scattlight'] = True
        par['calibrations']['illumflatframe']['process']['subtract_scattlight'] = True
        par['scienceframe']['process']['subtract_scattlight'] = True
        par['scienceframe']['process']['scattlight']['finecorr_method'] = 'median'
        par['scienceframe']['process']['scattlight']['finecorr_pad'] = 4  # This is the unbinned number of pixels to pad
        par['scienceframe']['process']['scattlight']['finecorr_order'] = 2
        # par['scienceframe']['process']['scattlight']['finecorr_mask'] = 12  # Mask the middle inter-slit region. It contains a strange scattered light feature that doesn't appear to affect any other inter-slit regions

        # Correct the illumflat for pixel-to-pixel sensitivity variations
        par['calibrations']['illumflatframe']['process']['use_pixelflat'] = True

        # Make sure the overscan is subtracted from the dark
        par['calibrations']['darkframe']['process']['use_overscan'] = True

        # Set the slit edge parameters
        par['calibrations']['slitedges']['fit_order'] = 4
        par['calibrations']['slitedges']['pad'] = 2  # Need to pad out the tilts for the astrometric transform when creating a datacube.
        par['calibrations']['slitedges']['edge_thresh'] = 5  # 5 works well with a range of setups tested by RJC (mostly 1x1 binning)

        # KCWI has non-uniform spectral resolution across the field-of-view
        par['calibrations']['wavelengths']['fwhm_spec_order'] = 1
        par['calibrations']['wavelengths']['fwhm_spat_order'] = 2

        # Alter the method used to combine pixel flats
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        par['calibrations']['flatfield']['spec_samp_coarse'] = 20.0
        par['calibrations']['flatfield']['spat_samp'] = 1.0  # This should give 1% accuracy in the spatial illumination correction for 2x2 binning, and <0.5% accuracy for 1x1 binning
        #par['calibrations']['flatfield']['tweak_slits'] = False  # Do not tweak the slit edges (we want to use the full slit)
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.0  # Make sure the full slit is used (i.e. when the illumination fraction is > 0.5)
        par['calibrations']['flatfield']['tweak_slits_maxfrac'] = 0.0  # Make sure the full slit is used (i.e. no padding)
        par['calibrations']['flatfield']['slit_trim'] = 3  # Trim the slit edges
        # Relative illumination correction
        par['calibrations']['flatfield']['slit_illum_relative'] = True  # Calculate the relative slit illumination
        par['calibrations']['flatfield']['slit_illum_ref_idx'] = 14  # The reference index - this should probably be the same for the science frame
        par['calibrations']['flatfield']['slit_illum_smooth_npix'] = 5  # Sufficiently small value so less structure in relative weights
        par['calibrations']['flatfield']['fit_2d_det_response'] = True  # Include the 2D detector response in the pixelflat.

        return par

    @staticmethod
    def is_nasmask(hdr):
        """
        Determine if a frame used nod-and-shuffle.

        Args:
            hdr (`astropy.io.fits.Header`_):
                The header of the raw frame.

        Returns:
            :obj:`bool`: True if NAS used.
        """
        return 'Mask' in hdr['BNASNAM']

    def get_rawimage(self, raw_file, det):
        """
        Read a raw KCWI data frame

        NOTE: The amplifiers are arranged as follows:

        |   (0,ny)  --------- (nx,ny)
        |           | 3 | 4 |
        |           ---------
        |           | 1 | 2 |
        |     (0,0) --------- (nx, 0)

        Parameters
        ----------
        raw_file : :obj:`str`
            File to read
        det : :obj:`int`
            1-indexed detector to read

        Returns
        -------
        detector_par : :class:`pypeit.images.detector_container.DetectorContainer`
            Detector metadata parameters.
        raw_img : `numpy.ndarray`_
            Raw image for this detector.
        hdu : `astropy.io.fits.HDUList`_
            Opened fits file
        exptime : :obj:`float`
            Exposure time read from the file header
        rawdatasec_img : `numpy.ndarray`_
            Data (Science) section of the detector as provided by setting the
            (1-indexed) number of the amplifier used to read each detector
            pixel. Pixels unassociated with any amplifier are set to 0.
        oscansec_img : `numpy.ndarray`_
            Overscan section of the detector as provided by setting the
            (1-indexed) number of the amplifier used to read each detector
            pixel. Pixels unassociated with any amplifier are set to 0.
        """
        fil = utils.find_single_file(f'{raw_file}*', required=True)

        # Read
        msgs.info(f'Reading KCWI file: {fil}')
        hdu = io.fits_open(fil)
        detpar = self.get_detector_par(det if det is not None else 1, hdu=hdu)
        head0 = hdu[0].header
        raw_img = hdu[detpar['dataext']].data.astype(float)

        # Some properties of the image
        numamps = head0['NVIDINP']
        # Exposure time (used by ProcessRawImage)
        headarr = self.get_headarr(hdu)
        exptime = self.get_meta_value(headarr, 'exptime')

        # Always assume normal FITS header formatting
        one_indexed = True
        include_last = True
        for section in ['DSEC', 'BSEC']:

            # Initialize the image (0 means no amplifier)
            pix_img = np.zeros(raw_img.shape, dtype=int)
            for aa, ampid in enumerate(1+np.arange(numamps)):
                # Get the data section
                sec = head0[section+"{0:1d}".format(ampid)]

                # Convert the data section from a string to a slice
                # TODO :: RJC - I think something has changed here... and the BPM is flipped (or not flipped) for different amp modes.
                # RJC - Note, KCWI records binned sections, so there's no need to pass binning in as an argument
                datasec = parse.sec2slice(sec, one_indexed=one_indexed, include_end=include_last, require_dim=2)
                # Flip the datasec
                datasec = datasec[::-1]

                # Assign the amplifier
                pix_img[datasec] = aa+1

            # Finish
            if section == 'DSEC':
                rawdatasec_img = pix_img.copy()
            elif section == 'BSEC':
                oscansec_img = pix_img.copy()

        # Return
        return detpar, raw_img, hdu, exptime, rawdatasec_img, oscansec_img

    def scattered_light_archive(self, binning, dispname):
        """ Archival model parameters for the scattered light. These are based on best fits to currently available data.

        For KCWI, the main contributor to the scattered light is referred to as the "narcissistic ghost"
        by Morrissey et al. (2018), ApJ, 864, 93. This scattered light is thought to be a reflection off the
        detector that travels back through the optical system. Some fraction gets sent back out to space, while
        the remained comes back through the optical system and a fuzzy version of this is re-imaged onto the
        detector. The current KCWI scattered light model is designed to account for these effects.

        Parameters
        ----------
        binning : :obj:`str`, optional
            Comma-separated binning along the spectral and spatial directions; e.g., ``2,1``
        dispname : :obj:`str`, optional
            Name of the disperser

        Returns
        -------
        x0 : `numpy.ndarray`_
            A 1D array containing the best-fitting model parameters
        bounds : :obj:`tuple`
            A tuple of two elements, containing two `numpy.ndarray`_ of the same length as x0. These
            two arrays contain the lower (first element of the tuple) and upper (second element of the tuple)
            bounds to consider on the scattered light model parameters.
        """
        # Grab the binning for convenience
        specbin, spatbin = parse.parse_binning(binning)

        # Get some starting parameters (these were determined by fitting spectra,
        # and should be close to the final fitted values to reduce computational time)
        # Note :: These values need to be originally based on data that uses 1x1 binning,
        # and are now scaled here according to the binning of the current data to be analysed.
        if dispname == 'BH2':
            # This solution had Cost: 1.0393e+08 and was based on a 1x1 dataset using pixelflat as the scattlight frame, and assuming pad=6
            x0 = np.array([67.15200530737414 / specbin, 157.1074288810557 / spatbin,  # Gaussian kernel widths
                           179.2999412601927 / specbin, 139.76705167365654 / spatbin,  # Lorentzian kernel widths
                           4.562250837489429 / specbin, 5.054084609129999 / spatbin,  # pixel offsets
                           0.9551212825380554, 0.9982713953567679,  # Zoom factor (spec, spat)
                           24.967114476694526,  # constant flux offset
                           0.09655699215523728,  # kernel angle
                           0.5524841628998713,  # Relative kernel scale (>1 means the kernel is more Gaussian, >0 but <1 makes the profile more lorentzian)
                           0.0035617904638365447, -0.004572414005692468, # Polynomial terms (coefficients of "spat" and "spat*spec")
                           0.10339051745377138, -0.011285432228341519, -0.007042406129643602])  # Polynomial terms (coefficients of spec**index)
        elif dispname == 'BM':
            # This solution had Cost: 4.8690e+07 and was based on a 1x1 dataset using pixelflat as the scattlight frame, and assuming pad=6
            x0 = np.array([57.52686698670778 / specbin, 44.22645738529251 / spatbin,  # Gaussian kernel widths
                           177.49996713255973 / specbin, 157.85206762558929 / spatbin,  # Lorentzian kernel widths
                           1.5547056520696672 / specbin, 5.1916115942048915 / spatbin,  # pixel offsets
                           0.9969235709861033, 0.9988876925628252,  # Zoom factor (spec, spat)
                           5.117391743273053,  # constant flux offset
                           0.1073126011932528,  # kernel angle
                           0.37915677855016305,  # Relative kernel scale (>1 means the kernel is more Gaussian, >0 but <1 makes the profile more lorentzian)
                           -0.0023288032996881753, 0.002167786497577728,  # Polynomial terms (coefficients of "spat" and "spat*spec")
                           0.08981952806519802, -0.07364263035160445, 0.04799106653657783])  # Polynomial terms (coefficients of spec**index)
        elif dispname == 'BL':
            # This solution had Cost: 7.0172e+06 and was based on a 2x2 dataset using pixelflat as the scattlight frame, and assuming pad=10
            x0 = np.array([54.843502304988725 / specbin, 71.36603219575882 / spatbin,  # Gaussian kernel widths
                           166.5990017834228 / specbin, 164.45188033168876 / spatbin,  # Lorentzian kernel widths
                           -5.759623374637964 / specbin, 5.01392929142184 / spatbin,  # pixel offsets
                           1.0017829374409521, 1.000312421855213,  # Zoom factor (spec, spat)
                           4.429458755393496,  # constant flux offset
                           -0.11853206385621386,  # kernel angle
                           0.4961668294341919,  # Relative kernel scale (>1 means the kernel is more Gaussian, >0 but <1 makes the profile more lorentzian)
                           -0.004790394657721825, 0.0032481886185675036,  # Polynomial terms (coefficients of "spat" and "spat*spec")
                           0.07823077510724392, -0.0644638013233617, 0.01819438897935518])  # Polynomial terms (coefficients of spec**index)
        else:
            msgs.warn(f"Initial scattered light model parameters have not been setup for grating {dispname}")
            x0 = np.array([54.843502304988725 / specbin, 71.36603219575882 / spatbin,  # Gaussian kernel widths
                           166.5990017834228 / specbin, 164.45188033168876 / spatbin,  # Lorentzian kernel widths
                           -5.759623374637964 / specbin, 5.01392929142184 / spatbin,  # pixel offsets
                           1.0017829374409521, 1.000312421855213,  # Zoom factor (spec, spat)
                           4.429458755393496,  # constant flux offset
                           -0.11853206385621386,  # kernel angle
                           0.4961668294341919,  # Relative kernel scale (>1 means the kernel is more Gaussian, >0 but <1 makes the profile more lorentzian)
                           -0.004790394657721825, 0.0032481886185675036,  # Polynomial terms (coefficients of "spat" and "spat*spec")
                           0.07823077510724392, -0.0644638013233617, 0.01819438897935518])  # Polynomial terms (coefficients of spec**index)

        # Now set the bounds of the fitted parameters
        bounds = ([# Lower bounds:
                      1, 1,  # Gaussian kernel widths
                      1, 1,  # Lorentzian kernel widths
                      -10 / specbin, -10 / spatbin,  # pixel offsets
                      0, 0,  # Zoom factor (spec, spat)
                      -1000, -(10 / 180) * np.pi, 0.0,  # constant flux offset, kernel angle, relative kernel scale
                      -10, -10, -10, -10, -10],  # Polynomial terms
                  # Upper bounds
                     [600 / specbin, 600 / spatbin,  # Gaussian kernel widths
                      600 / specbin, 600 / spatbin,  # Lorentzian kernel widths
                      10 / specbin, 10 / spatbin,  # pixel offsets
                      2, 2,  # Zoom factor (spec, spat)
                      1000.0, +(10 / 180) * np.pi, 1000.0,  # constant flux offset, kernel angle, relative kernel scale
                      10, 10, 10, 10, 10])  # Polynomial terms

        # Return the best-fitting archival parameters and the bounds
        return x0, bounds

    def fit_2d_det_response(self, det_resp, gpmask):
        r"""
        Perform a 2D model fit to the KCWI detector response.
        A few different setups were inspected (BH2 & BM with different
        grating angles), and a very similar response pattern was found for all
        setups, indicating that this structure is something to do with
        the detector. The starting parameters and functional form are
        assumed to be sufficient for all setups.

        Args:
            det_resp (`numpy.ndarray`_):
                An image of the flatfield structure.
            gpmask (`numpy.ndarray`_):
                Good pixel mask (True=good), the same shape as ff_struct.

        Returns:
            `numpy.ndarray`_: A model fit to the flatfield structure.
        """
        msgs.info("Performing a 2D fit to the detector response")

        # Define a 2D sine function, which is a good description of KCWI data
        def sinfunc2d(x, amp, scl, quad, phase, wavelength, angle):
            """
            2D Sine function
            """
            xx, yy = x
            angle *= np.pi / 180.0
            return 1 + (amp + xx*scl + xx*xx*quad) * np.sin(
                2 * np.pi * (xx * np.cos(angle) + yy * np.sin(angle)) / wavelength + phase)

        x = np.arange(det_resp.shape[0])
        y = np.arange(det_resp.shape[1])
        xx, yy = np.meshgrid(x, y, indexing='ij')
        # Prepare the starting parameters
        amp = 0.02  # Roughly a 2% effect
        scale, quad = 0.0, 0.0  # Assume the amplitude is constant over the detector
        wavelength = np.sqrt(det_resp.shape[0] ** 2 + det_resp.shape[1] ** 2) / 31.5  # 31-32 cycles of the pattern from corner to corner
        phase, angle = 0.0, -45.34  # No phase, and a decent guess at the angle
        p0 = [amp, scale, quad, phase, wavelength, angle]
        this_bpm = gpmask & (np.abs(det_resp-1) < 0.1)  # Only expect this to be a 5% effect
        popt, pcov = curve_fit(sinfunc2d, (xx[this_bpm], yy[this_bpm]), det_resp[this_bpm], p0=p0)
        return sinfunc2d((xx, yy), *popt)


class KeckKCRMSpectrograph(KeckKCWIKCRMSpectrograph):
    """
    Child to handle Keck/KCRM specific code

    """
    name = 'keck_kcrm'
    camera = 'KCRM'
    url = 'https://www2.keck.hawaii.edu/inst/kcwi/'  # TODO :: Need to update this website
    header_name = 'KCRM'
    comment = 'Supported setups: RL, RM1, RM2, RH3; see :doc:`keck_kcwi`'

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        .. warning::

            Many of the necessary detector parameters are read from the file
            header, meaning the ``hdu`` argument is effectively **required** for
            KCRM.  The optional use of ``hdu`` is only viable for automatically
            generated documentation.

        Args:
            det (:obj:`int`):
                1-indexed detector number.
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """
        if hdu is None:
            binning = '2,2'
            specflip = None
            numamps = None
            gainarr = None
            ronarr = None
#            dsecarr = None
#            msgs.error("A required keyword argument (hdu) was not supplied")
        else:
            # Some properties of the image
            binning = self.compound_meta(self.get_headarr(hdu), "binning")
            nampsxy = hdu[0].header['NAMPSXY'].split()
            numamps = int(nampsxy[0]) * int(nampsxy[1])
            amps = self.get_amplifiers(numamps)

            specflip = False
            gainarr = np.zeros(numamps)
            ronarr = np.zeros(numamps)  # Set this to zero (determine the readout noise from the overscan regions)
            for aa, amp in enumerate(amps):
                # Assign the gain for this amplifier
                gainarr[aa] = hdu[0].header["GAIN{0:1d}".format(amp)]

        detector = dict(det             = det,
                        binning         = binning,
                        dataext         = 0,
                        specaxis        = 0,
                        specflip        = specflip,
                        spatflip        = True,  # Due to the extra mirror, the slices are flipped relative to KCWI
                        platescale      = 0.145728,  # arcsec/pixel TODO :: Need to double check this
                        darkcurr        = None,  # e-/pixel/hour  TODO :: Need to check this.
                        mincounts       = -1e10,
                        saturation      = 65535.,
                        nonlinear       = 0.95,       # For lack of a better number!
                        numamplifiers   = numamps,
                        gain            = gainarr,
                        ronoise         = ronarr,
                        # These are never used because the image reader sets these up using the file headers data.
                        # datasec         = dsecarr, #.copy(),     # <-- This is provided in the header
                        # oscansec        = dsecarr, #.copy(),     # <-- This is provided in the header
                        )
        # Return
        return detector_container.DetectorContainer(**detector)

    def get_amplifiers(self, numamps):
        """
        Obtain a list of the amplifier ID numbers

        Args:
            numamps (:obj:`int`):
                Number of amplifiers used for readout

        Returns:
            :obj:`list`:
            A list (of length numamps) containing the ID number of the amplifiers used for readout
        """
        if numamps == 2:
            return [1, 3]
        elif numamps == 4:
            return [0, 1, 2, 3]
        else:
            msgs.error("PypeIt only supports 2 or 4 amplifier readout of KCRM data")

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        super().init_meta()
        self.meta['dispname'] = dict(ext=0, card='RGRATNAM')
        self.meta['dispangle'] = dict(ext=0, card='RGRANGLE', rtol=0.01)
        self.meta['cenwave'] = dict(ext=0, card='RCWAVE', rtol=0.01)

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
        return ['RGRATNAM', 'IFUNAM', 'RGRANGLE']

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # Subtract the detector pattern from certain frames.
        # NOTE: The pattern subtraction is time-consuming, meaning we don't
        # perform it (by default) for the high S/N pixel flat images but we do
        # for everything else.
        par['calibrations']['biasframe']['process']['use_pattern'] = False
        par['calibrations']['darkframe']['process']['use_pattern'] = False
        par['calibrations']['pixelflatframe']['process']['use_pattern'] = False
        par['calibrations']['illumflatframe']['process']['use_pattern'] = False
        par['calibrations']['standardframe']['process']['use_pattern'] = False
        par['scienceframe']['process']['use_pattern'] = False

        # Correct the illumflat for pixel-to-pixel sensitivity variations
        par['calibrations']['illumflatframe']['process']['use_pixelflat'] = True

        # Make sure the overscan is subtracted from the dark
        par['calibrations']['darkframe']['process']['use_overscan'] = True

        # Set the slit edge parameters
        par['calibrations']['slitedges']['fit_order'] = 4
        par['calibrations']['slitedges']['pad'] = 2  # Need to pad out the tilts for the astrometric transform when creating a datacube.
        par['calibrations']['slitedges']['edge_thresh'] = 5  # 5 works well with a range of setups tested by RJC (mostly 1x1 binning)

        # KCWI has non-uniform spectral resolution across the field-of-view
        par['calibrations']['wavelengths']['fwhm_spec_order'] = 1
        par['calibrations']['wavelengths']['fwhm_spat_order'] = 2

        # Alter the method used to combine pixel flats
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        par['calibrations']['flatfield']['spec_samp_coarse'] = 20.0
        #par['calibrations']['flatfield']['tweak_slits'] = False  # Do not tweak the slit edges (we want to use the full slit)
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.0  # Make sure the full slit is used (i.e. when the illumination fraction is > 0.5)
        par['calibrations']['flatfield']['tweak_slits_maxfrac'] = 0.0  # Make sure the full slit is used (i.e. no padding)
        par['calibrations']['flatfield']['slit_trim'] = 3  # Trim the slit edges
        # Relative illumination correction
        par['calibrations']['flatfield']['slit_illum_relative'] = True  # Calculate the relative slit illumination
        par['calibrations']['flatfield']['slit_illum_ref_idx'] = 14  # The reference index - this should probably be the same for the science frame
        par['calibrations']['flatfield']['slit_illum_smooth_npix'] = 5  # Sufficiently small value so less structure in relative weights
        par['calibrations']['flatfield']['fit_2d_det_response'] = True  # Include the 2D detector response in the pixelflat.

        # Sky subtraction parameters
        par['reduce']['skysub']['bspline_spacing'] = 0.4
        par['reduce']['skysub']['joint_fit'] = False

        return par

    @staticmethod
    def is_nasmask(hdr):
        """
        Determine if a frame used nod-and-shuffle.

        Args:
            hdr (`astropy.io.fits.Header`_):
                The header of the raw frame.

        Returns:
            :obj:`bool`: True if NAS used.
        """
        return 'Mask' in hdr['RNASNAM']

    def get_rawimage(self, raw_file, det):
        """
        Read a raw KCRM data frame

        Parameters
        ----------
        raw_file : :obj:`str`
            File to read
        det : :obj:`int`
            1-indexed detector to read

        Returns
        -------
        detector_par : :class:`pypeit.images.detector_container.DetectorContainer`
            Detector metadata parameters.
        raw_img : `numpy.ndarray`_
            Raw image for this detector.
        hdu : `astropy.io.fits.HDUList`_
            Opened fits file
        exptime : :obj:`float`
            Exposure time read from the file header
        rawdatasec_img : `numpy.ndarray`_
            Data (Science) section of the detector as provided by setting the
            (1-indexed) number of the amplifier used to read each detector
            pixel. Pixels unassociated with any amplifier are set to 0.
        oscansec_img : `numpy.ndarray`_
            Overscan section of the detector as provided by setting the
            (1-indexed) number of the amplifier used to read each detector
            pixel. Pixels unassociated with any amplifier are set to 0.
        """
        fil = utils.find_single_file(f'{raw_file}*', required=True)

        # Read
        msgs.info(f'Reading KCWI file: {fil}')
        hdu = io.fits_open(fil)
        detpar = self.get_detector_par(det if det is not None else 1, hdu=hdu)
        head0 = hdu[0].header
        raw_img = hdu[detpar['dataext']].data.astype(float)

        # Some properties of the image
        nampsxy = head0['NAMPSXY'].split()
        numamps = int(nampsxy[0]) * int(nampsxy[1])
        amps = self.get_amplifiers(numamps)
        # Exposure time (used by ProcessRawImage)
        headarr = self.get_headarr(hdu)
        exptime = self.get_meta_value(headarr, 'exptime')

        # get the x and y binning factors...
        #binning = self.get_meta_value(headarr, 'binning')

        # Always assume normal FITS header formatting
        one_indexed = True
        include_last = True
        for section in ['DSEC', 'BSEC']:

            # Initialize the image (0 means no amplifier)
            pix_img = np.zeros(raw_img.shape, dtype=int)
            for aa, ampid in enumerate(amps):
                # Get the data section
                sec = head0[section+"{0:1d}".format(ampid)]

                # Convert the data section from a string to a slice
                # TODO :: RJC - I think something has changed here... and the BPM is flipped (or not flipped) for different amp modes.
                # RJC - Note, KCWI records binned sections, so there's no need to pass binning in as an argument
                datasec = parse.sec2slice(sec, one_indexed=one_indexed, include_end=include_last, require_dim=2)
                # Flip the datasec
                datasec = datasec[::-1]

                # Assign the amplifier
                pix_img[datasec] = aa+1

            # Finish
            if section == 'DSEC':
                rawdatasec_img = pix_img.copy()
            elif section == 'BSEC':
                oscansec_img = pix_img.copy()

        # Return
        return detpar, raw_img, hdu, exptime, rawdatasec_img, oscansec_img
