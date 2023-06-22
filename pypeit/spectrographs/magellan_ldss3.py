""" Module for Magellan/LDSS3 specific codes
"""
import glob
import numpy as np

from pypeit import msgs
from pypeit import telescopes
from pypeit import io
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit.core import parse
from pypeit.images import detector_container

from IPython import embed

class MagellanLDSS3Spectrograph(spectrograph.Spectrograph):
    """
    Child to handle Magellan/LDSS3 specific code
    """
    ndet = 1
    name = 'magellan_ldss3'
    telescope = telescopes.MMTTelescopePar()
    camera = 'LDSS3'
    supported = True
    comment = 'Longslit only (so far)'

    def get_detector_par(self, hdu, det):
        """
        Return metadata for the selected detector.

        Args:
            hdu (`astropy.io.fits.HDUList`_):
                The open fits file with the raw image of interest.
            det (:obj:`int`):
                1-indexed detector number.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """

        # Binning
        binning = self.get_meta_value(self.get_headarr(hdu), 'binning')

        detector_dict = dict(
                            binning         = binning,
                            det             = 1,
                            dataext         = 1,
                            specaxis        = 0,
                            specflip        = False,
                            spatflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.189,
                            darkcurr        = 25.0,
                            saturation      = 205000.,
                            nonlinear       = 0.85,
                            mincounts       = -1e10,
                            numamplifiers   = 2,
                            gain            = np.atleast_1d([1.65,1.47]),
                            ronoise         = np.atleast_1d([4.67,5.06]),
                            )

        # Instantiate
        return detector_container.DetectorContainer(**detector_dict)

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the ``PypeIt``-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.

        """
        meta = {}
        # Required (core)
        meta['ra'] = dict(ext=0, card='RA')
        meta['dec'] = dict(ext=0, card='DEC')
        meta['target'] = dict(ext=0, card='OBJECT')
        meta['decker'] = dict(ext=0, card='APERTURE')
        meta['binning'] = dict(ext=0, card='BINNING', compound=True)

        meta['mjd'] = dict(ext=0, card='JD')
        meta['exptime'] = dict(ext=0, card='EXPTIME')
        meta['airmass'] = dict(ext=0, card='AIRMASS')

        meta['dispname'] = dict(ext=0, card='GRISM')
        meta['idname'] = dict(ext=0, card='EXPTYPE')
        meta['amp'] = dict(ext=0, card='OPAMP') # used for distinguish different amplifiers

        # Ingest
        self.meta = meta


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
            binspatial, binspec = parse.parse_binning(headarr[1]['BINNING'])
            binning = parse.binning2string(binspec, binspatial)
            return binning

    def default_pypeit_par(self):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of ``PypeIt`` methods.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'magellan_ldss3'
        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.5
        par['calibrations']['wavelengths']['sigdetect'] = 10.
        par['calibrations']['wavelengths']['fwhm']= 5.0
        par['calibrations']['wavelengths']['method'] = 'holy-grail'

        # Tilt and slit parameters
        par['calibrations']['tilts']['tracethresh'] =  20.0
        par['calibrations']['tilts']['spat_order'] = 6
        par['calibrations']['tilts']['spec_order'] = 6

        # edges
        par['calibrations']['slitedges']['edge_thresh'] = 300.

        # Processing steps
        turn_off = dict(use_biasimage=False, use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)

        # Extraction
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['extraction']['sn_gauss'] = 4.0
        ## Do not perform global sky subtraction for standard stars
        par['reduce']['skysub']['global_sky_std']  = False

        # Flexure
        par['flexure']['spec_method'] = 'boxcar'

        # cosmic ray rejection parameters for science frames
        par['scienceframe']['process']['sigclip'] = 5.0
        par['scienceframe']['process']['objlim'] = 2.0

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [10, 500]
        par['calibrations']['arcframe']['exprng'] = [0, 10]
        par['calibrations']['tiltframe']['exprng'] = [0, 10]
        par['calibrations']['darkframe']['exprng'] = [0, None]
        par['scienceframe']['exprng'] = [200, None]

        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 7
        par['sensfunc']['IR']['telgridfile'] = 'TelFit_LasCampanas_3100_26100_R20000.fits'

        return par

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the ``PypeIt`` parameters to hard-wired values used for
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
        par = self.default_pypeit_par() if inp_par is None else inp_par

        headarr = self.get_headarr(scifile)

        # Turn PCA off for long slits
#        if ('center' in self.get_meta_value(headarr, 'decker')) or \
#                ('Center' in self.get_meta_value(headarr, 'decker')) or \
#                ('red' in self.get_meta_value(headarr, 'decker')) or \
#                ('Red' in self.get_meta_value(headarr, 'decker')) or \
#                ('blue' in self.get_meta_value(headarr, 'decker')) or \
#                ('Blue' in self.get_meta_value(headarr, 'decker')) :
        if 'longslit' in self.get_meta_value(headarr, 'decker'):
            par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        # Turn on the use of mask design
        else:
            pass
            #par['calibrations']['slitedges']['use_maskdesign'] = False
            # Since we use the slitmask info to find the alignment boxes, I don't need `minimum_slit_length_sci`
            #par['calibrations']['slitedges']['minimum_slit_length_sci'] = None
            # Sometime the added missing slits at the edge of the detector are to small to be useful.
            #par['calibrations']['slitedges']['minimum_slit_length'] = 2.
            # Since we use the slitmask info to add and remove traces, 'minimum_slit_gap' may undo the matching effort.
            #par['calibrations']['slitedges']['minimum_slit_gap'] = 0.

        # Templates
        if self.get_meta_value(headarr, 'dispname') == 'VPH-Red':
            #par['calibrations']['wavelengths']['method'] = 'full_template'
            #par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_deimos_600.fits'
            par['calibrations']['wavelengths']['lamps'] = ['OH_MODS']
        else:
            #par['calibrations']['wavelengths']['method'] = 'full_template'
            #par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_deimos_830G.fits'
            par['calibrations']['wavelengths']['lamps'] = ['HeI', 'NeI', 'ArI', 'ArII']

        # FWHM
        binning = parse.parse_binning(self.get_meta_value(headarr, 'binning'))
        par['calibrations']['wavelengths']['fwhm'] = 6.0 / binning[1]

        # Return
        return par

    def bpm(self, filename, det, shape=None, msbias=None):
        """
        Generate a default bad-pixel mask.

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
                Master bias frame used to identify bad pixels

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """
        # Get the empty bpm: force is always True
        try:
            bpm_img = self.empty_bpm(filename, det, shape=shape)
        except:
            embed(header='259 of magellan_ldss3')

        # Fill in bad pixels if a master bias frame is provided
        if msbias is not None:
            return self.bpm_frombias(msbias, det, bpm_img)

        return bpm_img

    def configuration_keys(self):
        return ['dispname']

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
            return good_exp & (fitstbl['idname'] == 'Object') #& (fitstbl['amp'] == '1')
        if ftype == 'standard':
            return good_exp & (fitstbl['idname'] == 'Object') #& (fitstbl['amp'] == '1')
        if ftype in ['arc', 'tilt']:
            arc1 = good_exp & (fitstbl['idname'] == 'Object')
            #arc2 = (fitstbl['dispname'] == 'VPH-Red') & (fitstbl['idname'] == 'Object') &  (fitstbl['exptime'] >600)
            arc = (fitstbl['dispname'] == 'VPH-Red') & (fitstbl['idname'] == 'Object') &  (fitstbl['exptime'] >600)
            #return (arc1 | arc2) #& (fitstbl['amp'] == '1')
            return arc1
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            return good_exp & (fitstbl['idname'] == 'Flat') #& (fitstbl['amp'] == '1')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)


    def get_rawimage(self, raw_file, det):
        """
        Read raw images and generate a few other bits and pieces
        that are key for image processing.

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

        raw_file2 = raw_file.replace('c1.fits','c2.fits')

        # Check for file; allow for extra .gz, etc. suffix
        fil = glob.glob(raw_file + '*')
        if len(fil) != 1:
            msgs.error("Found {:d} files matching {:s}".format(len(fil)))

        # Check for file; allow for extra .gz, etc. suffix
        fil2 = glob.glob(raw_file2 + '*')
        if len(fil2) != 1:
            msgs.warn("Found {:d} files matching {:s}".format(len(fil2), raw_file))
            msgs.warn("Proceeding without the second amplifier")
            have_file2 = False
        else:
            have_file2 = True


        # detector par
        hdu = io.fits_open(fil[0])
        detector_par = self.get_detector_par(hdu, det if det is None else 1)

        data1, overscan1, datasec1, biassec1, x1_1, x2_1, nxb1 = ldss3_read_amp(fil[0])
        if have_file2:
            data2, overscan2, datasec2, biassec2, x1_2, x2_2, nxb2 = ldss3_read_amp(fil2[0])
        else:
            nxb2 = nxb1
            x2_2 = x2_1

        nx, ny = x2_1 + x2_2 + nxb1 + nxb2, data1.shape[1]

        # allocate output array...
        array = np.zeros((nx, ny))
        rawdatasec_img = np.zeros_like(array, dtype=int)
        oscansec_img = np.zeros_like(array, dtype=int)

        ## For amplifier 1
        array[:nxb1,:] = overscan1
        array[nxb1:x2_1+nxb1,:] = data1
        rawdatasec_img[nxb1:x2_1+nxb1, :] = 1
        oscansec_img[:nxb1,:] = 1 # exclude the first pixel since it always has problem

        ## For amplifier 2
        if have_file2:
            array[x2_1+nxb1+x2_2:x2_1+nxb1+x2_2+nxb2,:] = np.flipud(overscan2)
            array[x2_1+nxb1:x2_1+nxb1+x2_2,:] = np.flipud(data2)
        rawdatasec_img[x2_1+nxb1:x2_1+nxb1+x2_2, :] = 2
        oscansec_img[x2_1+nxb1+x2_2:x2_1+nxb1+x2_2+nxb2,:] = 2 # exclude the first pixel since it always has problem

        # Transpose now (helps with debuggin)
        array = array.T
        rawdatasec_img = rawdatasec_img.T
        oscansec_img = oscansec_img.T

        # Need the exposure time
        exptime = hdu[self.meta['exptime']['ext']].header[self.meta['exptime']['card']]
        # Return, transposing array back to orient the overscan properly
        return detector_par,array, hdu, exptime, rawdatasec_img, oscansec_img

def ldss3_read_amp(fil:str):
    """ Read a single amp of LDSS3 data

    Args:
        fil (str): filename

    Returns:
        tuple: data, overscan, datasec, biassec, x1, x2, nxb
    """

    msgs.info("Reading LDSS3 file: {:s}".format(fil))
    hdu = io.fits_open(fil)
    head1 = hdu[0].header

    # ToDo: Need to check binned data
    # get the x and y binning factors...
    binning = head1['BINNING']
    xbin, ybin = [int(ibin) for ibin in binning.split('x')]

    # First read over the header info to determine the size of the output array...
    datasec = head1['DATASEC']
    x1, x2, y1, y2 = np.array(parse.load_sections(datasec, fmt_iraf=False)).flatten()
    biassec = head1['BIASSEC']
    b1, b2, b3, b4 = np.array(parse.load_sections(biassec, fmt_iraf=False)).flatten()
    nxb = b2 - b1 + 1

    # determine the output array size...
    nx = (x2 - x1 + 1) + nxb
    ny = y2 - y1 + 1

    # allocate output array...
    array = hdu[0].data.T[:, :ny] * 1.0
    data = array[:nx-nxb,:]
    overscan = array[nx-nxb:,:]

    return data, overscan, datasec, biassec, x1, x2, nxb

