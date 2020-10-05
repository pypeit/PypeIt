""" Implements KCWI-specific functions.
"""

import glob
import numpy as np
from IPython import embed

from astropy import wcs, units
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import parse
from pypeit.core import procimg
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container


class KeckKCWISpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Keck/KCWI specific code

    TODO :: Need to apply spectral flexure and heliocentric correction to waveimg
    """
    ndet = 1

    def __init__(self):
        # Get it started
        # TODO :: Might need to change the tolerance of disperser angle in pypeit setup (two BH2 nights where sufficiently different that this was important).
        super(KeckKCWISpectrograph, self).__init__()
        self.spectrograph = 'keck_kcwi'
        self.telescope = telescopes.KeckTelescopePar()
        self.camera = 'KCWI'
        self.location = EarthLocation.of_site('Keck Observatory')  # TODO :: Might consider changing TelescopePar to use the astropy EarthLocation
        # Uses default timeunit
        # Uses default primary_hdrext
        # self.sky_file ?

        # Don't instantiate these until they're needed
        self.grating = None
        self.optical_model = None
        self.detector_map = None

    @property
    def pypeline(self):
        return 'IFU'

    def get_detector_par(self, hdu, det):
        """
        Return a DectectorContainer for the current image

        Args:
            hdu (`astropy.io.fits.HDUList`):
                HDUList of the image of interest.
                Ought to be the raw file, or else..
            det (int):

        Returns:
            :class:`pypeit.images.detector_container.DetectorContainer`:

        """
        # Some properties of the image
        head0 = hdu[0].header
        binning = self.compound_meta(self.get_headarr(hdu), "binning")
        numamps = head0['NVIDINP']
        specflip = True if head0['AMPID1'] == 2 else False
        gainmul, gainarr = head0['GAINMUL'], np.zeros(numamps)
        ronarr = np.zeros(numamps)  # Set this to zero (determine the readout noise from the overscan regions)
        dsecarr = np.array(['']*numamps)

        for ii in range(numamps):
            # Assign the gain for this amplifier
            gainarr[ii] = head0["GAIN{0:1d}".format(ii + 1)]# * gainmul

        detector = dict(det             = det,
                        binning         = binning,
                        dataext         = 0,
                        specaxis        = 0,
                        specflip        = specflip,
                        spatflip        = False,
                        platescale      = 0.145728,  # arcsec/pixel
                        darkcurr        = None,  # <-- TODO : Need to set this
                        mincounts       = -1e10,
                        saturation      = 65535.,
                        nonlinear       = 0.95,       # For lack of a better number!
                        numamplifiers   = numamps,
                        gain            = gainarr,
                        ronoise         = ronarr,
                        datasec         = dsecarr.copy(),     # <-- This is provided in the header
                        oscansec        = dsecarr.copy(),     # <-- This is provided in the header
                        )
        # Return
        return detector_container.DetectorContainer(**detector)

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

        headarr = self.get_headarr(scifile)

        # Templates
        if self.get_meta_value(headarr, 'dispname') == 'BH2':
            par['calibrations']['wavelengths']['method'] = 'full_template'  # 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_kcwi_BH2_4200.fits'
            par['calibrations']['wavelengths']['lamps'] = ['FeI', 'ArI', 'ArII']
        elif self.get_meta_value(headarr, 'dispname') == 'BM':
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_kcwi_BM.fits'
            par['calibrations']['wavelengths']['lamps'] = ['FeI', 'ArI', 'ArII']

        # FWHM
        # binning = parse.parse_binning(self.get_meta_value(headarr, 'binning'))
        # par['calibrations']['wavelengths']['fwhm'] = 6.0 / binning[1]

        # Return
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
        meta['ra'] = dict(ext=0, card=None, compound=True)
        meta['dec'] = dict(ext=0, card=None, compound=True)
        meta['target'] = dict(ext=0, card='TARGNAME')
        meta['dispname'] = dict(ext=0, card='BGRATNAM')
        meta['decker'] = dict(ext=0, card='IFUNAM')
        meta['binning'] = dict(card=None, compound=True)

        meta['mjd'] = dict(ext=0, card='MJD')
        meta['exptime'] = dict(ext=0, card='ELAPTIME')
        meta['airmass'] = dict(ext=0, card='AIRMASS')

        # Extras for config and frametyping
        meta['hatch'] = dict(ext=0, card='HATNUM')
        meta['idname'] = dict(ext=0, card='CALXPOS')
        meta['dispangle'] = dict(ext=0, card='BGRANGLE', rtol=0.01)
        meta['slitwid'] = dict(card=None, compound=True)

        # Get atmospheric conditions (note, these are the conditions at the end of the exposure)
        meta['obstime'] = dict(card=None, compound=True, required=False)
        meta['pressure'] = dict(card=None, compound=True, required=False)
        meta['temperature'] = dict(card=None, compound=True, required=False)
        meta['humidity'] = dict(card=None, compound=True, required=False)

        # Lamps
        lamp_names = ['LMP0', 'LMP1', 'LMP2', 'LMP3']  # FeAr, ThAr, Aux, Continuum
        for kk, lamp_name in enumerate(lamp_names):
            meta['lampstat{:02d}'.format(kk + 1)] = dict(ext=0, card=lamp_name+'STAT')
        for kk, lamp_name in enumerate(lamp_names):
            if lamp_name == 'LMP3':
                # There is no shutter on LMP3
                meta['lampshst{:02d}'.format(kk + 1)] = dict(ext=0, card=None, default=1)
                continue
            meta['lampshst{:02d}'.format(kk + 1)] = dict(ext=0, card=lamp_name+'SHST')
        # Ingest
        self.meta = meta

    def default_pypeit_par(self):
        """
        Set default parameters for Keck KCWI reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'keck_kcwi'

        # Subtract the detector pattern from certain frames
        par['calibrations']['biasframe']['process']['use_pattern'] = True
        par['calibrations']['darkframe']['process']['use_pattern'] = True
        par['calibrations']['pixelflatframe']['process']['use_pattern'] = False
        par['calibrations']['illumflatframe']['process']['use_pattern'] = True
        par['calibrations']['standardframe']['process']['use_pattern'] = True
        par['scienceframe']['process']['use_pattern'] = True

        # Make sure the overscan is subtracted from the dark
        par['calibrations']['darkframe']['process']['use_overscan'] = True

        # Set the slit edge parameters
        par['calibrations']['slitedges']['fit_order'] = 4

        # Always correct for flexure, starting with default parameters
        # slitcen must be used, because this is a slit-based IFU where
        # no objects are extracted.
        par['flexure']['spec_method'] = 'slitcen'

        # Alter the method used to combine pixel flats
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        par['calibrations']['pixelflatframe']['process']['sig_lohi'] = [10., 10.]
        par['calibrations']['flatfield']['spec_samp_coarse'] = 20.0
        #par['calibrations']['flatfield']['tweak_slits'] = False  # Do not tweak the slit edges (we want to use the full slit)
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.0  # Make sure the full slit is used (i.e. when the illumination fraction is > 0.5)
        par['calibrations']['flatfield']['tweak_slits_maxfrac'] = 0.0  # Make sure the full slit is used (i.e. no padding)
        par['calibrations']['flatfield']['slit_trim'] = 3  # Trim the slit edges
        par['calibrations']['flatfield']['slit_illum_relative'] = True  # Calculate the relative slit illumination

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 0.01]
        par['calibrations']['darkframe']['exprng'] = [0.01, None]
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['pixelflatframe']['exprng'] = [None, 30]
        par['calibrations']['traceframe']['exprng'] = [None, 30]
        par['scienceframe']['exprng'] = [30, None]

        # Set the number of alignments in the align frames
        par['calibrations']['alignment']['locations'] = [0.1, 0.3, 0.5, 0.7, 0.9]  # TODO:: Check this!!

        # LACosmics parameters
        par['scienceframe']['process']['sigclip'] = 4.0
        par['scienceframe']['process']['objlim'] = 1.5
        par['scienceframe']['process']['use_illumflat'] = True  # illumflat is applied when building the relative scale image in reduce.py, so should be applied to scienceframe too.
        par['scienceframe']['process']['use_specillum'] = True  # apply relative spectral illumination
        par['scienceframe']['process']['spat_flexure_correct'] = True  # correct for spatial flexure
        par['scienceframe']['process']['use_biasimage'] = False
        par['scienceframe']['process']['use_darkimage'] = False

        # Don't do optimal extraction for 3D data.
        par['reduce']['extraction']['skip_optimal'] = True

        # Make sure that this is reduced as a slit (as opposed to fiber) spectrograph
        par['reduce']['cube']['slit_spec'] = True

        # Sky subtraction parameters
        par['reduce']['skysub']['no_poly'] = True
        par['reduce']['skysub']['bspline_spacing'] = 0.6
        par['reduce']['skysub']['joint_fit'] = True

        return par

    def compound_meta(self, headarr, meta_key):
        if meta_key == 'binning':
            binspatial, binspec = parse.parse_binning(headarr[0]['BINNING'])
            binning = parse.binning2string(binspec, binspatial)
            return binning
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
                    hdrstr = ''
            return headarr[0][hdrstr]
        elif meta_key == 'pressure':
            return headarr[0]['WXPRESS'] * 0.001 * units.bar
        elif meta_key == 'temperature':
            return headarr[0]['WXOUTTMP'] * units.deg_C
        elif meta_key == 'humidity':
            return headarr[0]['WXOUTHUM'] / 100.0
        elif meta_key == 'obstime':
            return Time(headarr[0]['DATE-END'])
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
        return ['dispname', 'decker', 'binning', 'dispangle']

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype == 'science':
            return good_exp & self.lamps(fitstbl, 'off') & (fitstbl['hatch'] == '1')  #hatch=1,0=open,closed
        if ftype == 'bias':
            return good_exp & self.lamps(fitstbl, 'off') & (fitstbl['hatch'] == '0')
        if ftype in ['pixelflat', 'illumflat', 'trace']:
            # Flats and trace frames are typed together
            return good_exp & self.lamps(fitstbl, 'dome_noarc') & (fitstbl['hatch'] == '0') & (fitstbl['idname'] == '6')
        if ftype in ['dark']:
            # Dark frames
            return good_exp & self.lamps(fitstbl, 'off') & (fitstbl['hatch'] == '0')
        if ftype in ['align']:
            # Alignment frames
            return good_exp & self.lamps(fitstbl, 'dome') & (fitstbl['hatch'] == '0') & (fitstbl['idname'] == '4')
        if ftype in ['arc', 'tilt']:
            return good_exp & self.lamps(fitstbl, 'arcs') & (fitstbl['hatch'] == '0')
        if ftype in ['pinhole']:
            # Don't type pinhole frames
            return np.zeros(len(fitstbl), dtype=bool)

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
            lampstat = np.array([(fitstbl[k] == '0') | (fitstbl[k] == 'None')
                                    for k in fitstbl.keys() if 'lampstat' in k])
            lampshst = np.array([(fitstbl[k] == '0') | (fitstbl[k] == 'None')
                                    for k in fitstbl.keys() if 'lampshst' in k])
            return np.all(lampstat, axis=0)  # Lamp has to be off
            # return np.all(lampstat | lampshst, axis=0)  # i.e. either the shutter is closed or the lamp is off
        if status == 'arcs':
            # Check if any arc lamps are on (FeAr | ThAr)
            arc_lamp_stat = ['lampstat{0:02d}'.format(i) for i in range(1, 3)]
            arc_lamp_shst = ['lampshst{0:02d}'.format(i) for i in range(1, 3)]
            lamp_stat = np.array([fitstbl[k] == '1' for k in fitstbl.keys()
                                  if k in arc_lamp_stat])
            lamp_shst = np.array([fitstbl[k] == '1' for k in fitstbl.keys()
                                  if k in arc_lamp_shst])
            # Make sure the continuum frames are off
            dome_lamps = ['lampstat{0:02d}'.format(i) for i in range(4, 5)]
            dome_lamp_stat = np.array([fitstbl[k] == '0' for k in fitstbl.keys()
                                       if k in dome_lamps])
            return np.any(lamp_stat & lamp_shst & dome_lamp_stat, axis=0)  # i.e. lamp on and shutter open
        if status in ['dome_noarc', 'dome']:
            # Check if any dome lamps are on (Continuum) - Ignore lampstat03 (Aux) - not sure what this is used for
            dome_lamp_stat = ['lampstat{0:02d}'.format(i) for i in range(4, 5)]
            lamp_stat = np.array([fitstbl[k] == '1' for k in fitstbl.keys()
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
        Return a string containing the information on the lamp status

        Args:
            headarr (list of fits headers):
              list of headers

        Returns:
            str: A string that uniquely represents the lamp status
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
        raw_file : str
            Filename
        det (int or None):
            Detector number

        Returns
        -------
        array : ndarray
            Combined image
        hdu : HDUList
            Opened fits file.
        sections : list
            List of datasec, oscansec, ampsec sections. datasec,
            oscansec needs to be for an *unbinned* image as per
            standard convention
        """
        # Check for file; allow for extra .gz, etc. suffix
        fil = glob.glob(raw_file + '*')
        if len(fil) != 1:
            msgs.error("Found {:d} files matching {:s}".format(len(fil), raw_file))

        # Read
        msgs.info("Reading KCWI file: {:s}".format(fil[0]))
        hdu = fits.open(fil[0])
        detpar = self.get_detector_par(hdu, det if det is None else 1)
        head0 = hdu[0].header
        raw_img = hdu[detpar['dataext']].data.astype(float)

        # Some properties of the image
        numamps = head0['NVIDINP']
        # Exposure time (used by ProcessRawImage)
        headarr = self.get_headarr(hdu)
        exptime = self.get_meta_value(headarr, 'exptime')

        # get the x and y binning factors...
        binning = self.get_meta_value(headarr, 'binning')

        # Always assume normal FITS header formatting
        one_indexed = True
        include_last = True
        for section in ['DSEC', 'BSEC']:

            # Initialize the image (0 means no amplifier)
            pix_img = np.zeros(raw_img.shape, dtype=int)
            for i in range(numamps):
                # Get the data section
                sec = head0[section+"{0:1d}".format(i+1)]

                # Convert the data section from a string to a slice
                # TODO :: I fear something has changed here... and the BPM is flipped (ot not flipped) for different amp modes.
                datasec = parse.sec2slice(sec, one_indexed=one_indexed,
                                          include_end=include_last, require_dim=2,
                                          binning=binning)
                # Flip the datasec
                datasec = datasec[::-1]

                # Assign the amplifier
                pix_img[datasec] = i+1

            # Finish
            if section == 'DSEC':
                rawdatasec_img = pix_img.copy()
            elif section == 'BSEC':
                oscansec_img = pix_img.copy()

        # Calculate the pattern frequency
        hdu = self.calc_pattern_freq(raw_img, rawdatasec_img, oscansec_img, hdu)

        # Return
        return detpar, raw_img, hdu, exptime, rawdatasec_img, oscansec_img

    def calc_pattern_freq(self, frame, rawdatasec_img, oscansec_img, hdu):
        """Calculate the pattern frequency using the overscan region that
        covers the overscan and data sections. Using a larger range allows
        the frequency to be pinned down with high accuracy.

        NOTE: The amplifiers are arranged as follows:

        |   (0,ny)  --------- (nx,ny)
        |           | 3 | 4 |
        |           ---------
        |           | 1 | 2 |
        |     (0,0) --------- (nx, 0)

        TODO :: PATTERN FREQUENCY ALGORITHM HAS NOT BEEN TESTED WHEN BINNING != 1x1

        Parameters
        ----------
        frame : ndarray
            Raw data frame to be used to estimate the pattern frequency
        rawdatasec_img : ndarray
            array the same shape as frame, used as a mask to identify the
            data pixels (0 is no data, non-zero values indicate the amplifier number)
        oscansec_img : ndarray
            array the same shape as frame, used as a mask to identify the
            overscan pixels (0 is no data, non-zero values indicate the amplifier number)
        hdu : HDUList
            Opened fits file.

        Returns
        -------
        hdu : HDUList
            The input HDUList, with header updated to include the frequency of each amplifier
        """
        msgs.info("Calculating pattern noise frequency")

        # Make a copy of te original frame
        raw_img = frame.copy()

        # Get a unique list of the amplifiers
        unq_amps = np.sort(np.unique(oscansec_img[np.where(oscansec_img >= 1)]))
        num_amps = unq_amps.size

        # Loop through amplifiers and calculate the frequency
        for amp in unq_amps:
            # Grab the pixels where the amplifier has data
            pixs = np.where((rawdatasec_img == amp) | (oscansec_img == amp))
            rmin, rmax = np.min(pixs[1]), np.max(pixs[1])
            # Deal with the different locations of the overscan regions in 2- and 4- amp mode
            if num_amps == 2:
                cmin = 1+np.max(pixs[0])
                frame = raw_img[cmin:, rmin:rmax].astype(np.float64)
            elif num_amps == 4:
                if amp in [1, 2]:
                    pixalt = np.where((rawdatasec_img == amp+2) | (oscansec_img == amp+2))
                    cmin = 1+np.max(pixs[0])
                    cmax = (np.min(pixalt[0]) + cmin)//2  # Average of the bottom of the top amp, and top of the bottom amp
                else:
                    pixalt = np.where((rawdatasec_img == amp-2) | (oscansec_img == amp-2))
                    cmax = 1+np.min(pixs[0])
                    cmin = (np.max(pixalt[0]) + cmax)//2
                frame = raw_img[cmin:cmax, rmin:rmax].astype(np.float64)
            # Calculate the pattern frequency
            freq = procimg.pattern_frequency(frame)
            msgs.info("Pattern frequency of amplifier {0:d}/{1:d} = {2:f}".format(amp, num_amps, freq))
            # Add the frequency to the zeroth header
            hdu[0].header['PYPFRQ{0:02d}'.format(amp)] = freq

        # Return the updated HDU
        return hdu

    def bpm(self, filename, det, shape=None, msbias=None):
        """
        Override parent bpm function with BPM specific to KCWI.

        Parameters
        ----------
        det : int, REQUIRED
        **null_kwargs:
            Captured and never used

        Returns
        -------
        bpix : ndarray
            0 = ok; 1 = Mask

        """
        bpm_img = self.empty_bpm(filename, det, shape=shape)

        # Fill in bad pixels if a master bias frame is provided
        # if msbias is not None:
        #     return self.bpm_frombias(msbias, det, bpm_img)

        # Extract some header info
        #msgs.info("Reading AMPMODE and BINNING from KCWI file: {:s}".format(filename))
        head0 = fits.getheader(filename, ext=0)
        ampmode = head0['AMPMODE']
        binning = head0['BINNING']

        # Construct a list of the bad columns
        # Note: These were taken from v1.1.0 (REL) Date: 2018/06/11 of KDERP
        #       KDERP store values and in the code (stage1) subtract 1 from the badcol data files.
        #       Instead of this, I have already pre-subtracted the values in the following arrays.
        bc = None
        if ampmode == 'ALL':
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
        if ampmode == 'TUP':
            if binning == '1,1':
                bc = [[2622, 2622, 3492, 3528],
                      [3295, 3300, 1550, 1555],
                      [3676, 3676, 1866, 4111]]
            elif binning == '2,2':
                bc = [[1311, 1311, 1745, 1788],
                      [1646, 1650,  775,  777],
                      [1838, 1838,  933, 2055]]
        if bc is None:
            msgs.warn("Bad pixel mask is not available for ampmode={0:s} binning={1:s}".format(ampmode, binning))
            bc = []

        # Apply these bad columns to the mask
        for bb in range(len(bc)):
            bpm_img[bc[bb][2]:bc[bb][3]+1, bc[bb][0]:bc[bb][1]+1] = 1

        return np.flipud(bpm_img)

    @staticmethod
    def is_nasmask(hdr):
        """Determine if a frame used NAS

        Parameters
        ----------
        hdr : fits header
            The header of the raw frame.

        Returns
        -------
        bool : True if NAS used
        """
        if 'Mask' in hdr['BNASNAM']:
            return True
        else:
            return False

    def get_wcs(self, hdr, slits, platescale, wave0, dwv):
        """Get the WCS for a frame

        Parameters
        ----------
        detector : DetectorContainer
            A DetectorContainer object containing parameters of the detector
        hdr : fits header
            The header of the raw frame. The information in this
            header will be extracted and returned as a WCS.
        slits : :class:`pypeit.slittrace.SlitTraceSet`
            Master slit edges
        platescale : float
            platescale of an unbinned pixel in arcsec/pixel (e.g. detector.platescale)
        wave0 : float
            wavelength zeropoint
        dwv : float
            delta wavelength per spectral pixel

        Returns
        -------
        astropy.wcs : An astropy WCS object.
        """
        msgs.info("Calculating the WCS")
        # Get the x and y binning factors, and the typical slit length
        binspec, binspat = parse.parse_binning(self.get_meta_value([hdr], 'binning'))

        # Get the pixel and slice scales
        pxscl = platescale * binspat / 3600.0  # Need to convert arcsec to degrees
        slscl = self.get_meta_value([hdr], 'slitwid')

        # Get the typical slit length (this changes by ~0.3% over all slits, so a constant is fine for now)
        slitlength = int(np.round(np.median(slits.get_slitlengths(initial=True, median=True))))

        # Get RA/DEC
        raval = self.compound_meta([hdr], 'ra')
        decval = self.compound_meta([hdr], 'dec')

        # Create a coordinate
        coord = SkyCoord(raval, decval, unit=(units.deg, units.deg))

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
        rotoff = 0.0  # IFU-SKYPA offset (degrees)
        skypa = rpos + rref  # IFU position angle (degrees)
        crota = np.radians(-(skypa + rotoff))

        # Calculate the fits coordinates
        cdelt1 = -slscl
        cdelt2 = pxscl
        if coord is None:
            ra = 0.
            dec = 0.
            crota = 1
        else:
            ra = coord.ra.degree
            dec = coord.dec.degree
        # Calculate the CD Matrix
        cd11 = cdelt1 * np.cos(crota)                          # RA degrees per column
        cd12 = abs(cdelt2) * np.sign(cdelt1) * np.sin(crota)   # RA degrees per row
        cd21 = -abs(cdelt1) * np.sign(cdelt2) * np.sin(crota)  # DEC degress per column
        cd22 = cdelt2 * np.cos(crota)                          # DEC degrees per row
        # Get reference pixels (set these to the middle of the FOV)
        crpix1 = 12.   # i.e. 24 slices/2
        crpix2 = slitlength / 2.
        crpix3 = 1.
        # Get the offset
        porg = hdr['PONAME']
        ifunum = hdr['IFUNUM']
        if 'IFU' in porg:
            if ifunum == 1:  # Large slicer
                off1 = 1.0
                off2 = 4.0
            elif ifunum == 2:  # Medium slicer
                off1 = 1.0
                off2 = 5.0
            elif ifunum == 3:  # Small slicer
                off1 = 0.05
                off2 = 5.6
            else:
                msgs.warn("Unknown IFU number: {0:d}".format(ifunum))
                off1 = 0.
                off2 = 0.
            off1 /= binspec
            off2 /= binspat
            crpix1 += off1
            crpix2 += off2

        # Create a new WCS object.
        msgs.info("Generating KCWI WCS")
        w = wcs.WCS(naxis=3)
        w.wcs.equinox = hdr['EQUINOX']
        w.wcs.name = 'KCWI'
        w.wcs.radesys = 'FK5'
        # Insert the coordinate frame
        w.wcs.cname = ['KCWI RA', 'KCWI DEC', 'KCWI Wavelength']
        w.wcs.cunit = [units.degree, units.degree, units.Angstrom]
        w.wcs.ctype = ["RA---TAN", "DEC--TAN", "AWAV"]
        w.wcs.crval = [ra, dec, wave0]  # RA, DEC, and wavelength zeropoints
        w.wcs.crpix = [crpix1, crpix2, crpix3]  # RA, DEC, and wavelength reference pixels
        w.wcs.cd = np.array([[cd11, cd12, 0.0], [cd21, cd22, 0.0], [0.0, 0.0, dwv]])
        w.wcs.lonpole = 180.0  # Native longitude of the Celestial pole
        w.wcs.latpole = 0.0  # Native latitude of the Celestial pole

        return w

    def get_datacube_bins(self, slitlength, minmax, num_wave):
        """Calculate the bin edges to be used when making a datacube

        Parameters
        ----------
        slitlength : int
            Length of the slit in pixels
        minmax : `numpy.ndarray`_
            An array of size (nslits, 2), listing the minimum and maximum pixel
            locations on each slit relative to the reference location (usually
            the centre of the slit). This array is returned by the function
            `slittrace.SlitTraceSet.get_radec_image`_
        num_wave : int
            Number of wavelength steps = int(round((wavemax-wavemin)/delta_wave))

        Returns
        -------
        tuple : Three 1D numpy.ndarray providing the bins to use when constructing a histogram
                of the spec2d files. The elements are (x, y, lambda).
        """
        xbins = np.arange(1 + 24) - 12.0 - 0.5
        ybins = np.linspace(np.min(minmax[:, 0]), np.max(minmax[:, 1]), 1+slitlength) - 0.5
        spec_bins = np.arange(1+num_wave) - 0.5
        return xbins, ybins, spec_bins
