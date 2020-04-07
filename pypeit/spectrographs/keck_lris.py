""" Module for LRIS specific codes
"""
import glob
import os
import numpy as np
from astropy.io import fits

from pkg_resources import resource_filename

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container

from IPython import embed

class KeckLRISSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Keck/LRIS specific code
    """
    ndet = 2

    def __init__(self):
        # Get it started
        super(KeckLRISSpectrograph, self).__init__()
        self.spectrograph = 'keck_lris_base'
        self.telescope = telescopes.KeckTelescopePar()

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for Keck LRISr reductions.
        """
        par = pypeitpar.PypeItPar()
        # Set wave tilts order
        par['calibrations']['slitedges']['edge_thresh'] = 15.
        par['calibrations']['slitedges']['fit_order'] = 3
        par['calibrations']['slitedges']['sync_center'] = 'gap'
        # TODO: I had to increase this from 1. to 2. to deal with
        # Keck_LRIS_red/multi_1200_9000_d680_1x2/ . May need a
        # different solution given that this is binned data and most of
        # the data in the dev suite is unbinned.
        # JXP -- Increased to 6 arcsec.  I don't know how 2 (or 1!) could have worked.
        par['calibrations']['slitedges']['minimum_slit_length'] = 6
        # 1D wavelengths
        par['calibrations']['wavelengths']['rms_threshold'] = 0.20  # Might be grism dependent
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['pixelflatframe']['exprng'] = [None, 30]    # This may be too low for LRISb
        par['calibrations']['traceframe']['exprng'] = [None, 30]

        # Flexure
        # Always correct for spectral flexure, starting with default parameters
        par['flexure']['spec_method'] = 'boxcar'
        # Always correct for spatial flexure on science images
        # TODO -- Decide whether to make the following defaults
        #   May not want to do them for LongSlit
        par['scienceframe']['process']['spat_flexure_correct'] = True
        par['calibrations']['standardframe']['process']['spat_flexure_correct'] = True

        par['scienceframe']['exprng'] = [29, None]
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

        # Ignore PCA if longslit
        #  This is a little risky as a user could put long into their maskname
        #  But they would then need to over-ride in their PypeIt file
        if scifile is None:
            msgs.error("You have not included a standard or science file in your PypeIt file to determine the configuration")
        if 'long' in self.get_meta_value(scifile, 'decker'):
            par['calibrations']['slitedges']['sync_predict'] = 'nearest'
            # This might only be required for det=2, but we'll see..
            if self.spectrograph == 'keck_lris_red':
                par['calibrations']['slitedges']['edge_thresh'] = 1000.

        return par

    def init_meta(self):
        """
        Generate the meta data dict
        Note that the children can add to this

        Returns:
            self.meta: dict (generated in place)

        """
        meta = {}
        # Required (core)
        meta['ra'] = dict(ext=0, card='RA')
        meta['dec'] = dict(ext=0, card='DEC')
        meta['target'] = dict(ext=0, card='TARGNAME')
        meta['decker'] = dict(ext=0, card='SLITNAME')
        meta['binning'] = dict(card=None, compound=True)

        meta['mjd'] = dict(ext=0, card='MJD-OBS')
        meta['exptime'] = dict(ext=0, card='ELAPTIME')
        meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        meta['dichroic'] = dict(ext=0, card='DICHNAME')
        meta['hatch'] = dict(ext=0, card='TRAPDOOR')
        # Red only, but grabbing here
        meta['dispangle'] = dict(ext=0, card='GRANGLE', rtol=1e-2)

        # Lamps
        lamp_names = ['MERCURY', 'NEON', 'ARGON', 'CADMIUM', 'ZINC', 'KRYPTON', 'XENON',
                      'FEARGON', 'DEUTERI', 'FLAMP1', 'FLAMP2', 'HALOGEN']
        for kk,lamp_name in enumerate(lamp_names):
            meta['lampstat{:02d}'.format(kk+1)] = dict(ext=0, card=lamp_name)
        # Ingest
        self.meta = meta

    def compound_meta(self, headarr, meta_key):
        if meta_key == 'binning':
#            return '1,1'
            binspatial, binspec = parse.parse_binning(headarr[0]['BINNING'])
            binning = parse.binning2string(binspec, binspatial)
            return binning
        else:
            msgs.error("Not ready for this compound meta")

    def configuration_keys(self):
        """
        Return the metadata keys that defines a unique instrument
        configuration.

        This list is used by :class:`pypeit.metadata.PypeItMetaData` to
        identify the unique configurations among the list of frames read
        for a given reduction.

        Returns:
            list: List of keywords of data pulled from file headers and
            used to constuct the :class:`pypeit.metadata.PypeItMetaData`
            object.
        """
        return ['dispname', 'dichroic', 'decker', 'binning']

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype == 'science':
            return good_exp & self.lamps(fitstbl, 'off') & (fitstbl['hatch'] == 'open')
        if ftype == 'bias':
            return good_exp & self.lamps(fitstbl, 'off') & (fitstbl['hatch'] == 'closed')
        if ftype in ['pixelflat', 'trace']:
            # Flats and trace frames are typed together
            return good_exp & self.lamps(fitstbl, 'dome') & (fitstbl['hatch'] == 'open')
        if ftype in ['pinhole', 'dark']:
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            return good_exp & self.lamps(fitstbl, 'arcs') & (fitstbl['hatch'] == 'closed')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)
  
    def lamps(self, fitstbl, status):
        """
        Check the lamp status.

        Args:
            fitstbl (:obj:`astropy.table.Table`):
                The table with the fits header meta data.
            status (:obj:`str`):
                The status to check.  Can be `off`, `arcs`, or `dome`.
        
        Returns:
            numpy.ndarray: A boolean array selecting fits files that
            meet the selected lamp status.

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
            arc_lamp_stat = [ 'lampstat{0:02d}'.format(i) for i in range(1,9) ]
            return np.any(np.array([ fitstbl[k] == 'on' for k in fitstbl.keys()
                                            if k in arc_lamp_stat]), axis=0)
        if status == 'dome':
            # Check if any dome lamps are on
            # Warning 9, 10 are FEARGON and DEUTERI
            dome_lamp_stat = [ 'lampstat{0:02d}'.format(i) for i in range(9,13) ]
            return np.any(np.array([ fitstbl[k] == 'on' for k in fitstbl.keys()
                                            if k in dome_lamp_stat]), axis=0)
        raise ValueError('No implementation for status = {0}'.format(status))

    def get_rawimage(self, raw_file, det):
        """
        Read a raw LRIS data frame (one or more detectors)
        Packed in a multi-extension HDU
        Based on readmhdufits.pro
        Parameters
        ----------
        raw_file : str
          Filename
        det (int or None):
          Detector number; Default = both
        Returns
        -------
        tuple
            See :func:`pypeit.spectrograph.spectrograph.get_rawimage`
        """
        # Check for file; allow for extra .gz, etc. suffix
        fil = glob.glob(raw_file + '*')
        if len(fil) != 1:
            msgs.error("Found {:d} files matching {:s}".format(len(fil)))

        # Read
        msgs.info("Reading LRIS file: {:s}".format(fil[0]))
        hdu = fits.open(fil[0])
        head0 = hdu[0].header

        # Get post, pre-pix values
        precol = head0['PRECOL']
        postpix = head0['POSTPIX']
        preline = head0['PRELINE']
        postline = head0['POSTLINE']

        # get the x and y binning factors...
        binning = head0['BINNING']
        xbin, ybin = [int(ibin) for ibin in binning.split(',')]

        # First read over the header info to determine the size of the output array...
        n_ext = len(hdu) - 1  # Number of extensions (usually 4)
        xcol = []
        xmax = 0
        ymax = 0
        xmin = 10000
        ymin = 10000
        for i in np.arange(1, n_ext + 1):
            theader = hdu[i].header
            detsec = theader['DETSEC']
            if detsec != '0':
                # parse the DETSEC keyword to determine the size of the array.
                x1, x2, y1, y2 = np.array(parse.load_sections(detsec, fmt_iraf=False)).flatten()

                # find the range of detector space occupied by the data
                # [xmin:xmax,ymin:ymax]
                xt = max(x2, x1)
                xmax = max(xt, xmax)
                yt = max(y2, y1)
                ymax = max(yt, ymax)

                # find the min size of the array
                xt = min(x1, x2)
                xmin = min(xmin, xt)
                yt = min(y1, y2)
                ymin = min(ymin, yt)
                # Save
                xcol.append(xt)

        # determine the output array size...
        nx = xmax - xmin + 1
        ny = ymax - ymin + 1

        # change size for binning...
        nx = nx // xbin
        ny = ny // ybin

        # Update PRECOL and POSTPIX
        precol = precol // xbin
        postpix = postpix // xbin

        # Deal with detectors
        if det in [1, 2]:
            nx = nx // 2
            n_ext = n_ext // 2
            det_idx = np.arange(n_ext, dtype=np.int) + (det - 1) * n_ext
        elif det is None:
            det_idx = np.arange(n_ext).astype(int)
        else:
            raise ValueError('Bad value for det')

        # change size for pre/postscan...
        nx += n_ext * (precol + postpix)
        ny += preline + postline

        # allocate output arrays...
        array = np.zeros((nx, ny))
        order = np.argsort(np.array(xcol))
        rawdatasec_img = np.zeros_like(array, dtype=int)
        oscansec_img = np.zeros_like(array, dtype=int)

        # insert extensions into master image...
        for amp, i in enumerate(order[det_idx]):

            # grab complete extension...
            data, predata, postdata, x1, y1 = lris_read_amp(hdu, i + 1)

            # insert predata...
            buf = predata.shape
            nxpre = buf[0]
            xs = amp * precol
            xe = xs + nxpre
            # predata (ignored)
            array[xs:xe, :] = predata

            # insert data...
            buf = data.shape
            nxdata = buf[0]
            xs = n_ext * precol + amp * nxdata  # (x1-xmin)/xbin
            xe = xs + nxdata
            array[xs:xe, :] = data
            rawdatasec_img[xs:xe, preline:ny-postline] = amp+1

            # ; insert postdata...
            buf = postdata.shape
            nxpost = buf[0]
            xs = nx - n_ext * postpix + amp * postpix
            xe = xs + nxpost
            array[xs:xe, :] = postdata
            oscansec_img[xs:xe, preline:ny-postline] = amp+1

        # Need the exposure time
        exptime = hdu[self.meta['exptime']['ext']].header[self.meta['exptime']['card']]
        # Return
        return self.get_detector_par(hdu, det if det is None else 1), \
                array.T, hdu, exptime, rawdatasec_img.T, oscansec_img.T

    def subheader_for_spec(self, row_fitstbl, raw_header):
        """
        See :func:`pypeit.spectrograph.spectrograph.Spectrograph.subheader_for_spec`
        for doc string

        Args:
            row_fitstbl (:class:`astropy.table.Row` or :class:`astropy.io.fits.Header`):
            raw_header (:class:`astropy.io.fits.Header`):

        Returns:
            :obj:`dict`: -- Used to generate a Header or table downstream

        """
        return super(KeckLRISSpectrograph, self).subheader_for_spec(
            row_fitstbl, raw_header, extra_header_cards=['GRANAME', 'GRISNAME', 'SLITNAME'])


class KeckLRISBSpectrograph(KeckLRISSpectrograph):
    """
    Child to handle Keck/LRISb specific code
    """
    def __init__(self):
        # Get it started
        super(KeckLRISBSpectrograph, self).__init__()
        self.spectrograph = 'keck_lris_blue'
        self.camera = 'LRISb'
        # Uses default timeunit
        # Uses default primary_hdrext
        self.sky_file = 'sky_LRISb_600.fits'

    def get_detector_par(self, hdu, det):
        # Binning
        binning = self.get_meta_value(self.get_headarr(hdu), 'binning')  # Could this be detector dependent??

        # Detector 1
        detector_dict1 = dict(
            binning         = binning,
            det             = 1,
            dataext         = 1,
            specaxis        = 0,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.135,
            darkcurr        = 0.0,
            saturation      = 65535.,
            nonlinear       = 0.86,
            mincounts       = -1e10,
            numamplifiers   = 2,
            gain            = np.atleast_1d([1.55, 1.56]),
            ronoise         = np.atleast_1d([3.9, 4.2]),
            )
        # Detector 2
        detector_dict2 = detector_dict1.copy()
        detector_dict2.update(dict(
            det=2,
            dataext=2,
            gain=np.atleast_1d([1.63, 1.70]),
            ronoise=np.atleast_1d([3.6, 3.6])
        ))

        # Instantiate
        detector_dicts = [detector_dict1, detector_dict2]
        detector = detector_container.DetectorContainer(**detector_dicts[det])

        # Deal with number of amps
        namps = hdu[0].header['NUMAMPS']
        # The website does not give values for single amp per detector so we take the mean
        #   of the values provided
        if namps == 2:
            detector.numamplifiers = 1
            detector.gain = np.atleast_1d(np.mean(detector.gain))
            detector.ronoise = np.atleast_1d(np.mean(detector.ronoise))
        elif namps == 4:
            pass
        else:
            msgs.error("Did not see this namps coming..")

        # Return
        return detector


    def default_pypeit_par(self):
        """
        Set default parameters for Keck LRISb reductions.
        """
        par = KeckLRISSpectrograph.default_pypeit_par()
        par['rdx']['spectrograph'] = 'keck_lris_blue'

        par['calibrations']['slitedges']['det_min_spec_length'] = 0.1
        par['calibrations']['slitedges']['fit_min_spec_length'] = 0.2

        # 1D wavelength solution -- Additional parameters are grism dependent
        par['calibrations']['wavelengths']['rms_threshold'] = 0.20  # Might be grism dependent..
        par['calibrations']['wavelengths']['sigdetect'] = 10.0

        par['calibrations']['wavelengths']['lamps'] = ['NeI', 'ArI', 'CdI', 'KrI', 'XeI', 'ZnI', 'HgI']
        #par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['n_first'] = 3
        par['calibrations']['wavelengths']['match_toler'] = 2.5
        par['calibrations']['wavelengths']['method'] = 'full_template'


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
        # Start with instrument wide
        par = super(KeckLRISBSpectrograph, self).config_specific_par(scifile, inp_par=inp_par)

        # Wavelength calibrations
        if self.get_meta_value(scifile, 'dispname') == '300/5000':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_lris_blue_300_d680.fits'
            par['flexure']['spectrum'] = os.path.join(resource_filename('pypeit', 'data/sky_spec/'),
                                                      'sky_LRISb_400.fits')
        elif self.get_meta_value(scifile, 'dispname') == '400/3400':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_lris_blue_400_d560.fits'
            par['flexure']['spectrum'] = os.path.join(resource_filename('pypeit', 'data/sky_spec/'),
                                                  'sky_LRISb_400.fits')
        elif self.get_meta_value(scifile, 'dispname') == '600/4000':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_lris_blue_600_d560.fits'
            par['flexure']['spectrum'] = os.path.join(resource_filename('pypeit', 'data/sky_spec/'),
                                                      'sky_LRISb_600.fits')
        elif self.get_meta_value(scifile, 'dispname') == '1200/3400':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_lris_blue_1200_d460.fits'
            par['flexure']['spectrum'] = os.path.join(resource_filename('pypeit', 'data/sky_spec/'),
                                                      'sky_LRISb_600.fits')

        # FWHM
        binning = parse.parse_binning(self.get_meta_value(scifile, 'binning'))
        par['calibrations']['wavelengths']['fwhm'] = 8.0 / binning[0]

        # Slit tracing
        # Reduce the slit parameters because the flux does not span the full detector
        #   It is primarily on the upper half of the detector (usually)
        if self.get_meta_value(scifile, 'dispname') == '300/5000':
            par['calibrations']['slitedges']['smash_range'] = [0.5, 1.]

        # Return
        return par

    def init_meta(self):
        """
        Meta data specific to Keck LRIS red

        Returns:

        """
        super(KeckLRISBSpectrograph, self).init_meta()
        # Add the name of the dispersing element
        self.meta['dispname'] = dict(ext=0, card='GRISNAME')

    def bpm(self, filename, det, shape=None, msbias=None):
        """ Generate a BPM

        Args:
            filename (str):
            det (int):
            msbias : numpy.ndarray, required if the user wishes to generate a BPM based on a master bias

        Returns:
            np.ndarray:

        """
        # Get the empty bpm: force is always True
        bpm_img = self.empty_bpm(filename, det, shape=shape)

        # Fill in bad pixels if a master bias frame is provided
        if msbias is not None:
            return self.bpm_frombias(msbias, det, bpm_img)

        # Only defined for det=1
        if det == 1:
            msgs.info("Using hard-coded BPM for det=1 on LRISb")
            bpm_img[:,:3] = 1

        return bpm_img


class KeckLRISRSpectrograph(KeckLRISSpectrograph):
    """
    Child to handle Keck/LRISr specific code
    """
    def __init__(self):
        # Get it started
        super(KeckLRISRSpectrograph, self).__init__()
        self.spectrograph = 'keck_lris_red'
        self.camera = 'LRISr'
#        # Allow for variable namps
#        if ifile is not None:
#            hdu = fits.open(ifile)
#            head0 = hdu[0].header
#            if head0['AMPPSIZE'] == '[1:1024,1:4096]':
#                namps = 2
#            elif head0['AMPPSIZE'] == '[1:2048,1:4096]':
#                namps = 1
#            else:
#                msgs.error("Bad amp size (windowed??)!!")
#        else: # Default
#            namps = 2
        #
    def get_detector_par(self, hdu, det):
        # Binning
        binning = self.get_meta_value(self.get_headarr(hdu), 'binning')  # Could this be detector dependent??

        # Detector 1
        detector_dict1 = dict(
            binning=binning,
            det=1,
            dataext=1,
            specaxis=0,
            specflip=False,
            spatflip=False,
            platescale=0.135,
            darkcurr=0.0,
            saturation=65535.,
            nonlinear=0.76,
            mincounts=-1e10,
            numamplifiers=2,
            gain=np.atleast_1d([1.255, 1.18]),
            ronoise=np.atleast_1d([4.64, 4.76]),
        )
        # Detector 2
        detector_dict2 = detector_dict1.copy()
        detector_dict2.update(dict(
            det=2,
            dataext=2,
            gain=np.atleast_1d([1.191, 1.162]),
            ronoise=np.atleast_1d([4.54, 4.62])
        ))

        # Instantiate
        detector_dicts = [detector_dict1, detector_dict2]
        detector = detector_container.DetectorContainer(**detector_dicts[det])

        # Deal with number of amps
        namps = hdu[0].header['NUMAMPS']
        # The website does not give values for single amp per detector so we take the mean
        #   of the values provided
        if namps == 2:
            detector.numamplifiers = 1
            # Long silt mode
            if hdu[0].header['AMPPSIZE'] == '[1:1024,1:4096]':
                idx = 0 if det==1 else 1  # Vid1 for det=1, Vid4 for det=2
                detector.gain = np.atleast_1d(detector.gain[idx])
                detector.ronoise = np.atleast_1d(detector.ronoise[idx])
            else:
                detector.gain = np.atleast_1d(np.mean(detector.gain))
                detector.ronoise = np.atleast_1d(np.mean(detector.ronoise))
        elif namps == 4:
            pass
        else:
            msgs.error("Did not see this namps coming..")

        # Return
        return detector

    def default_pypeit_par(self):
        """
        Set default parameters for Keck LRISr reductions.
        """
        par = KeckLRISSpectrograph.default_pypeit_par()
        par['rdx']['spectrograph'] = 'keck_lris_red'
        #
        par['calibrations']['slitedges']['edge_thresh'] = 20.

        # 1D wavelength solution
        par['calibrations']['wavelengths']['lamps'] = ['NeI', 'ArI', 'CdI', 'KrI', 'XeI', 'ZnI', 'HgI']
        #par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['sigdetect'] = 10.0
        # Tilts
        # These are the defaults
        par['calibrations']['tilts']['tracethresh'] = 25
        par['calibrations']['tilts']['spat_order'] = 4
        par['calibrations']['tilts']['spec_order'] = 7
        par['calibrations']['tilts']['maxdev2d'] = 1.0
        par['calibrations']['tilts']['maxdev_tracefit'] = 1.0
        par['calibrations']['tilts']['sigrej2d'] = 5.0

        #  Sky Subtraction
        par['reduce']['skysub']['bspline_spacing'] = 0.8

        # Defaults for anything other than 1,1 binning
        #  Rest config_specific_par below if binning is (1,1)
        par['scienceframe']['process']['sigclip'] = 5.
        par['scienceframe']['process']['objlim'] = 5.

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
        # Start with instrument wide
        par = super(KeckLRISRSpectrograph, self).config_specific_par(scifile, inp_par=inp_par)

        # Lacosmic CR settings
        #   Grab the defaults for LRISr
        binning = self.get_meta_value(scifile, 'binning')
        # Unbinned LRISr needs very aggressive LACosmics parameters for 1x1 binning
        if binning == '1,1':
            sigclip = 3.0
            objlim = 0.5
            par['scienceframe']['process']['sigclip'] = sigclip
            par['scienceframe']['process']['objlim'] = objlim

        # Wavelength calibrations
        if self.get_meta_value(scifile, 'dispname') == '400/8500':  # This is basically a reidentify
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_lris_red_400.fits'
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['sigdetect'] = 20.0
            par['calibrations']['wavelengths']['nsnippet'] = 1
        elif self.get_meta_value(scifile, 'dispname') == '600/5000':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_lris_red_600_5000.fits'
            par['calibrations']['wavelengths']['method'] = 'full_template'
        elif self.get_meta_value(scifile, 'dispname') == '600/7500':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_lris_red_600_7500.fits'
            par['calibrations']['wavelengths']['method'] = 'full_template'
        elif self.get_meta_value(scifile, 'dispname') == '1200/9000':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_lris_red_1200_9000.fits'
            par['calibrations']['wavelengths']['method'] = 'full_template'

        # FWHM
        binning = parse.parse_binning(self.get_meta_value(scifile, 'binning'))
        par['calibrations']['wavelengths']['fwhm'] = 8.0 / binning[0]


        # Return
        return par


    def init_meta(self):
        """
        Meta data specific to Keck LRIS red

        Returns:

        """
        super(KeckLRISRSpectrograph, self).init_meta()
        # Add the name of the dispersing element
        self.meta['dispname'] = dict(ext=0, card='GRANAME')

    def configuration_keys(self):
        """
        Return the metadata keys that defines a unique instrument
        configuration.

        This list is used by :class:`pypeit.metadata.PypeItMetaData` to
        identify the unique configurations among the list of frames read
        for a given reduction.

        Returns:

            list: List of keywords of data pulled from meta
        """
        cfg_keys = super(KeckLRISRSpectrograph, self).configuration_keys()
        # Add grating tilt
        return cfg_keys+['dispangle']

    def bpm(self, filename, det, shape=None, msbias=None):
        """ Generate a BPM

        Args:
            filename (str):
            det (int):
            msbias : numpy.ndarray, required if the user wishes to generate a BPM based on a master bias

        Returns:
            np.ndarray

        """
        # Get the empty bpm: force is always True
        bpm_img = self.empty_bpm(filename, det, shape=shape)

        # Fill in bad pixels if a master bias frame is provided
        if msbias is not None:
            return self.bpm_frombias(msbias, det, bpm_img)

        # Only defined for det=2
        if det == 2:
            msgs.info("Using hard-coded BPM for det=2 on LRISr")

            # Get the binning
            hdu = fits.open(filename)
            binning = hdu[0].header['BINNING']
            hdu.close()

            # Apply the mask
            xbin = int(binning.split(',')[0])
            badc = 16//xbin
            bpm_img[:,:badc] = 1

            # Mask the end too (this is risky as an edge may appear)
            #  But there is often weird behavior at the ends of these detectors
            bpm_img[:,-10:] = 1

        return bpm_img


#class KeckLRISRLSpectrograph(KeckLRISRSpectrograph):
#    """
#    Child to handle Keck/LRISr in Long-slit readout mode (Vid1, Vid4)
#    """
#    def __init__(self):
#        # Get it started
#        super(KeckLRISRSpectrograph, self).__init__()
#        self.spectrograph = 'keck_lris_red_longonly'
#        self.camera = 'LRISr'
#        self.detector = [
#                # Detector 1
#                pypeitpar.DetectorPar(
#                            dataext         =1,
#                            specaxis        =0,
#                            specflip        = False,
#                            xgap            =0.,
#                            ygap            =0.,
#                            ysize           =1.,
#                            platescale      =0.135,
#                            darkcurr        =0.0,
#                            saturation      =65535.,  # Gain applied
#                            nonlinear       =0.86,          # Modified by JXP to go higher
#                            numamplifiers   =1,
#                            gain            =[1.255],
#                            ronoise         =[4.64],
#                            datasec         = ['',''],
#                            oscansec        = ['',''],
#                            suffix          ='_01red'
#                            ),
#                #Detector 2
#                pypeitpar.DetectorPar(
#                            dataext         =2,
#                            specaxis        =0,
#                            specflip        = False,
#                            xgap            =0.,
#                            ygap            =0.,
#                            ysize           =1.,
#                            platescale      =0.135,
#                            darkcurr        =0.,
#                            saturation      =65535.,  # Gain applied
#                            nonlinear       =0.86,
#                            numamplifiers   =1,
#                            gain            =[1.162],
#                            ronoise         =[4.62],
#                            datasec         = ['',''],
#                            oscansec        = ['',''],
#                            suffix          ='_02red'
#                            )]
#        self.numhead = 3
#        # Uses default timeunit


def lris_read_amp(inp, ext):
    """
    Read one amplifier of an LRIS multi-extension FITS image

    Args:
        inp (str, astropy.io.fits.HDUList):
            filename or HDUList
        ext (int):
            Extension index

    Returns:
        tuple:
            data
            predata
            postdata
            x1
            y1

    """
    # Parse input
    if isinstance(inp, str):
        hdu = fits.open(inp)
    else:
        hdu = inp
    n_ext = len(hdu) - 1  # Number of extensions (usually 4)

    # Get the pre and post pix values
    # for LRIS red POSTLINE = 20, POSTPIX = 80, PRELINE = 0, PRECOL = 12
    head0 = hdu[0].header
    precol = head0['precol']
    postpix = head0['postpix']

    # Deal with binning
    binning = head0['BINNING']
    xbin, ybin = [int(ibin) for ibin in binning.split(',')]
    precol = precol//xbin
    postpix = postpix//xbin

    # get entire extension...
    temp = hdu[ext].data.transpose() # Silly Python nrow,ncol formatting
    tsize = temp.shape
    nxt = tsize[0]

    # parse the DETSEC keyword to determine the size of the array.
    header = hdu[ext].header
    detsec = header['DETSEC']
    x1, x2, y1, y2 = np.array(parse.load_sections(detsec, fmt_iraf=False)).flatten()

    # parse the DATASEC keyword to determine the size of the science region (unbinned)
    datasec = header['DATASEC']
    xdata1, xdata2, ydata1, ydata2 = np.array(parse.load_sections(datasec, fmt_iraf=False)).flatten()

    # grab the components...
    predata = temp[0:precol, :]
    # datasec appears to have the x value for the keywords that are zero
    # based. This is only true in the image header extensions
    # not true in the main header.  They also appear inconsistent between
    # LRISr and LRISb!
    #data     = temp[xdata1-1:xdata2-1,*]
    #data = temp[xdata1:xdata2+1, :]
    if (xdata1-1) != precol:
        msgs.error("Something wrong in LRIS datasec or precol")
    xshape = 1024 // xbin * (4//n_ext)  # Allow for single amp
    if (xshape+precol+postpix) != temp.shape[0]:
        msgs.warn("Unexpected size for LRIS detector.  We expect you did some windowing...")
        xshape = temp.shape[0] - precol - postpix
    data = temp[precol:precol+xshape,:]
    postdata = temp[nxt-postpix:nxt, :]

    # flip in X as needed...
    if x1 > x2:
        xt = x2
        x2 = x1
        x1 = xt
        data = np.flipud(data)

    # flip in Y as needed...
    if y1 > y2:
        yt = y2
        y2 = y1
        y1 = yt
        data = np.fliplr(data)
        predata = np.fliplr(predata)
        postdata = np.fliplr(postdata)

    return data, predata, postdata, x1, y1


def convert_lowredux_pixelflat(infil, outfil):
    """ Convert LowRedux pixelflat to PYPIT format
    Returns
    -------

    """
    # Read
    hdu = fits.open(infil)
    data = hdu[0].data

    #
    prihdu = fits.PrimaryHDU()
    hdus = [prihdu]
    prihdu.header['FRAMETYP'] = 'pixelflat'

    # Detector 1
    img1 = data[:,:data.shape[1]//2]
    hdu = fits.ImageHDU(img1)
    hdu.name = 'DET1'
    prihdu.header['EXT0001'] = 'DET1-pixelflat'
    hdus.append(hdu)

    # Detector 2
    img2 = data[:,data.shape[1]//2:]
    hdu = fits.ImageHDU(img2)
    hdu.name = 'DET2'
    prihdu.header['EXT0002'] = 'DET2-pixelflat'
    hdus.append(hdu)

    # Finish
    hdulist = fits.HDUList(hdus)
    hdulist.writeto(outfil, clobber=True)
    print('Wrote {:s}'.format(outfil))

