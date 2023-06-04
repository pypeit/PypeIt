"""
Module for LRIS specific methods.

.. include:: ../include/links.rst
"""
import glob
import os

from IPython import embed

import numpy as np

from astropy.io import fits
from astropy import time
from astropy.coordinates import SkyCoord 
from astropy import units

import linetools.utils

from pypeit import msgs
from pypeit import telescopes
from pypeit import io
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.spectrographs import slitmask
from pypeit.images import detector_container
from pypeit import data


class KeckLRISSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Keck/LRIS specific code
    """
    ndet = 2
    telescope = telescopes.KeckTelescopePar()
    url = 'https://www2.keck.hawaii.edu/inst/lris/'

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # Set wave tilts order
        par['calibrations']['slitedges']['edge_thresh'] = 15.
        par['calibrations']['slitedges']['fit_order'] = 3
        par['calibrations']['slitedges']['sync_center'] = 'gap'
        # TODO: I had to increase this from 1. to 2. to deal with
        # Keck_LRIS_red/multi_1200_9000_d680_1x2/ . May need a
        # different solution given that this is binned data and most of
        # the data in the dev suite is unbinned.
        # JXP -- Increased to 6 arcsec.  I don't know how 2 (or 1!) could have worked.
        par['calibrations']['slitedges']['minimum_slit_length_sci'] = 6
        # Remove slits that are too short
        par['calibrations']['slitedges']['minimum_slit_length'] = 4.
        # 1D wavelengths
        par['calibrations']['wavelengths']['rms_threshold'] = 0.20  # Might be grism dependent
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['pixelflatframe']['exprng'] = [None, 60]
        par['calibrations']['traceframe']['exprng'] = [None, 60]
        par['calibrations']['standardframe']['exprng'] = [None, 30]

        # Flexure
        # Always correct for spectral flexure, starting with default parameters
        par['flexure']['spec_method'] = 'boxcar'
        # Always correct for spatial flexure on science images
        # TODO -- Decide whether to make the following defaults
        #   May not want to do them for LongSlit
        par['scienceframe']['process']['spat_flexure_correct'] = True
        par['calibrations']['standardframe']['process']['spat_flexure_correct'] = True

        par['scienceframe']['exprng'] = [60, None]


        # If telluric is triggered
        par['sensfunc']['IR']['telgridfile'] = 'TelFit_MaunaKea_3100_26100_R20000.fits'

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

        # Ignore PCA if longslit
        #  This is a little risky as a user could put long into their maskname
        #  But they would then need to over-ride in their PypeIt file
        if scifile is None:
            msgs.error("You have not included a standard or science file in your PypeIt file to determine the configuration")
        if 'long' in self.get_meta_value(scifile, 'decker'):
            par['calibrations']['slitedges']['sync_predict'] = 'nearest'
            # This might only be required for det=2, but we'll see..
            # TODO: Why is this here and not in KeckLRISRSpectrograph???
            if self.name == 'keck_lris_red':
                par['calibrations']['slitedges']['edge_thresh'] = 1000.

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
        self.meta['target'] = dict(ext=0, card='TARGNAME')
        self.meta['decker'] = dict(ext=0, card='SLITNAME')
        self.meta['binning'] = dict(card=None, compound=True)
        # 
        self.meta['mjd'] = dict(ext=0, card='MJD-OBS')
        self.meta['exptime'] = dict(ext=0, card='ELAPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['dichroic'] = dict(ext=0, card='DICHNAME')
        self.meta['hatch'] = dict(ext=0, card='TRAPDOOR')
        # Red only, but grabbing here
        self.meta['dispangle'] = dict(ext=0, card='GRANGLE', rtol=1e-2)
        self.meta['cenwave'] = dict(ext=0, card='WAVELEN', rtol=2.0)
        self.meta['frameno'] = dict(ext=0, card='FRAMENO')
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')

        # Extras for pypeit file
        if self.name == 'keck_lris_red_mark4':
            self.meta['amp'] = dict(ext=0, card='TAPLINES')
        else:
            self.meta['amp'] = dict(ext=0, card='NUMAMPS')

        # Lamps -- Have varied in time..
        for kk in range(12): # This needs to match the length of LAMPS below
            self.meta['lampstat{:02d}'.format(kk+1)] = dict(card=None, compound=True)

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
        elif 'lampstat' in meta_key:
            idx = int(meta_key[-2:])
            try:
                curr_date = time.Time(self.get_meta_value(headarr, 'mjd'), format='mjd')
            except:
                msgs.warn('No mjd in header. You either have bad headers, '
                          'or incorrectly specified the wrong spectrograph, '
                          'or are reading in other files from your directory.  '
                          'Using 2022-01-01 as the date for parsing lamp info from headers')
                curr_date =  time.Time("2022-01-01", format='isot')
            # Modern -- Assuming the change occurred with the new red detector
            t_newlamp = time.Time("2014-02-15", format='isot')  # LAMPS changed in Header
            if curr_date > t_newlamp:
                lamp_names = ['MERCURY', 'NEON', 'ARGON', 'CADMIUM', 'ZINC', 'KRYPTON', 'XENON',
                              'FEARGON', 'DEUTERI', 'FLAMP1', 'FLAMP2', 'HALOGEN']
                return headarr[0][lamp_names[idx-1]]  # Use this index is offset by 1
            else:  # Original lamps
                plamps = headarr[0]['LAMPS'].split(',')
                # https: // www2.keck.hawaii.edu / inst / lris / instrument_key_list.html
                old_lamp_names = ['MERCURY', 'NEON', 'ARGON', 'CADMIUM', 'ZINC', 'HALOGEN']
                if idx <= 5: # Arcs
                    return ('off' if plamps[idx - 1] == '0' else 'on')
                elif idx == 10:  # Current FLAMP1
                    return headarr[0]['FLIMAGIN'].strip()
                elif idx == 11:  # Current FLAMP2
                    return headarr[0]['FLSPECTR'].strip()
                elif idx == 12:  # Current Halogen slot
                    return ('off' if plamps[len(old_lamp_names)-1] == '0' else 'on')
                else:  # Lamp didn't exist.  Set to None
                    return 'None'
        else:
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
        return super().configuration_keys() + ['binning']

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard PypeIt file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
        return super().pypeit_file_keys() + ['frameno']

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
            return good_exp & self.lamps(fitstbl, 'off') & (fitstbl['hatch'] == 'open')
        if ftype == 'standard':
            return good_exp & self.lamps(fitstbl, 'off') & (fitstbl['hatch'] == 'open')
        if ftype == 'bias':
            return good_exp & self.lamps(fitstbl, 'off') & (fitstbl['hatch'] == 'closed')
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            # Allow for dome or internal
            good_dome = self.lamps(fitstbl, 'dome') & (fitstbl['hatch'] == 'open')
            good_internal = self.lamps(fitstbl, 'halogen') & (fitstbl['hatch'] == 'closed')
            # Flats and trace frames are typed together
            return good_exp & (good_dome + good_internal)
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
        elif status == 'arcs':
            # Check if any arc lamps are on
            arc_lamp_stat = [ 'lampstat{0:02d}'.format(i) for i in range(1,9) ]
            return np.any(np.array([ fitstbl[k] == 'on' for k in fitstbl.keys()
                                            if k in arc_lamp_stat]), axis=0)
        elif status == 'dome':
            # Check if any dome lamps are on
            # Warning 9, 10 are FEARGON and DEUTERI
            dome_lamp_stat = [ 'lampstat{0:02d}'.format(i) for i in range(9,13) ]
            return np.any(np.array([ fitstbl[k] == 'on' for k in fitstbl.keys()
                                            if k in dome_lamp_stat]), axis=0)
        elif status == 'halogen':
            halo_lamp_stat = ['lampstat12']
            return np.any(np.array([ fitstbl[k] == 'on' for k in fitstbl.keys()
                                            if k in halo_lamp_stat]), axis=0)
        else:
            msgs.error(f"Bad status option! {status}")

        raise ValueError('No implementation for status = {0}'.format(status))

    def get_rawimage(self, raw_file, det):
        """
        Read raw images and generate a few other bits and pieces
        that are key for image processing.

        Based on readmhdufits.pro

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
        # Check for file; allow for extra .gz, etc. suffix
        fil = glob.glob(raw_file + '*')
        if len(fil) != 1:
            msgs.error("Found {:d} files matching {:s}".format(len(fil)))

        # Read
        msgs.info("Reading LRIS file: {:s}".format(fil[0]))
        hdu = io.fits_open(fil[0])
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
        extensions = []
        for kk, ihdu in enumerate(hdu):
            if 'VidInp' in ihdu.name:
                extensions.append(kk)
        n_ext = len(extensions)
        xcol = []
        xmax = 0
        ymax = 0
        xmin = 10000
        ymin = 10000

        for i in extensions:
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
            det_idx = np.arange(n_ext, dtype=int) + (det - 1) * n_ext
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

        # insert extensions into calibration image...
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
        return self.get_detector_par(det if det is not None else 1, hdu=hdu), \
                array.T, hdu, exptime, rawdatasec_img.T, oscansec_img.T

    def subheader_for_spec(self, row_fitstbl, raw_header, extra_header_cards=None,
                           allow_missing=False):
        """
        Generate a dict that will be added to the Header of spectra files
        generated by PypeIt (e.g. :class:`~pypeit.specobjs.SpecObjs`).

        Args:
            row_fitstbl (dict-like):
                Typically an `astropy.table.Row`_ or
                `astropy.io.fits.Header`_ with keys defined by
                :func:`~pypeit.core.meta.define_core_meta`.
            raw_header (`astropy.io.fits.Header`_):
                Header that defines the instrument and detector, meaning that
                the header must contain the ``INSTRUME``, ``DETECTOR``,
                ``GRANAME``, ``GRISNAME``, and ``SLITNAME`` header cards. If
                provided, this must also contain the header cards provided by
                ``extra_header_cards``.
            extra_header_cards (:obj:`list`, optional):
                Additional header cards from ``raw_header`` to include in the
                output dictionary. Can be an empty list or None.
            allow_missing (:obj:`bool`, optional):
                Ignore any keywords returned by
                :func:`~pypeit.core.meta.define_core_meta` are not present in
                ``row_fitstbl``. Otherwise, raise ``PypeItError``.

        Returns:
            :obj:`dict`: Dictionary with data to include an output fits
            header file or table downstream.
        """
        _extra_header_cards = ['GRANAME', 'GRISNAME', 'SLITNAME']
        if extra_header_cards is not None:
            _extra_header_cards += extra_header_cards
        return super().subheader_for_spec(row_fitstbl, raw_header,
                                          extra_header_cards=_extra_header_cards,
                                          allow_missing=allow_missing)

    def get_slitmask(self, filename:str):
        """
        Parse the slitmask data from a LRIS file into :attr:`slitmask`, a
        :class:`~pypeit.spectrographs.slitmask.SlitMask` object.

        Args:
            filename (:obj:`str`):
                Name of the file to read.

        Returns:
            :class:`~pypeit.spectrographs.slitmask.SlitMask`: The slitmask
            data read from the file. The returned object is the same as
            :attr:`slitmask`.
        """
        self.slitmask = slitmask.load_keck_deimoslris(filename, self.name)
        return self.slitmask

    def get_maskdef_slitedges(self, ccdnum=None, filename=None, debug=None,
                              trc_path=None, binning=None):
        """
        Provides the slit edges positions predicted by the slitmask design using
        the mask coordinates already converted from mm to pixels by the method
        `mask_to_pixel_coordinates`.

        If not already instantiated, the :attr:`slitmask`, :attr:`amap`,
        and :attr:`bmap` attributes are instantiated.  If so, a file must be provided.

        Args:
            ccdnum (:obj:`int`):
                Detector number
            filename (:obj:`str`):
                The filename to use to (re)instantiate the :attr:`slitmask` and :attr:`grating`.
                Default is None, i.e., to use previously instantiated attributes.
            debug (:obj:`bool`, optional):
                Run in debug mode.

        Returns:
            :obj:`tuple`: Three `numpy.ndarray`_ and a :class:`~pypeit.spectrographs.slitmask.SlitMask`.
            Two arrays are the predictions of the slit edges from the slitmask design and
            one contains the indices to order the slits from left to right in the PypeIt orientation

        """
        # Re-initiate slitmask
        if filename is not None:
            self.get_slitmask(filename)
        else:
            msgs.error('The name of a science file should be provided for the mask info')

        if self.slitmask is None:
            msgs.error('Unable to read slitmask design info. Provide a file.')

        platescale = self.get_detector_par(det=1)['platescale']


        hdu = fits.open(filename)
        binning = self.get_meta_value(self.get_headarr(hdu), 'binning')
        _, bin_spat = parse.parse_binning(binning)

        # Slit center
        slit_coords = SkyCoord(ra=self.slitmask.onsky[:,0], 
                               dec=self.slitmask.onsky[:,1], unit='deg')
        mask_coord = SkyCoord(ra=self.slitmask.mask_radec[0],
                              dec=self.slitmask.mask_radec[1], unit='deg')

        # build an array of values containing the bottom (right) edge of the slits
        # starting edge
        left_edges = []
        #for islit in x_order:
        for islit in range(self.slitmask.nslits):
            sep = mask_coord.separation(slit_coords[islit])
            PA = mask_coord.position_angle(slit_coords[islit])
            #
            alpha = sep.to('arcsec') * np.cos(PA-self.slitmask.posx_pa*units.deg)
            #delta = sep.to('arcsec') * np.sin(PA-self.slitmask.posx_pa*units.deg)
            dx_pix = (alpha.value-self.slitmask.onsky[islit,2]/2.) / (platescale*bin_spat)
            # target is the slit number
            left_edges.append(np.round(dx_pix))
        left_edges = np.array(left_edges, dtype=int)

        # Build up the right edges
        right_edges = left_edges + np.round(
            self.slitmask.onsky[:,2]/(platescale*bin_spat)).astype(int)

        # Center of slit
        centers = (left_edges + right_edges)/2.

        # Trim down by detector
        # TODO -- Deal with Mark4
        max_spat = 2048//bin_spat
        if ccdnum == 1:
            if self.name == 'keck_lris_red':
                good = centers < 0.
                xstart = max_spat + 160//bin_spat  # The 160 is for the chip gap
            elif self.name == 'keck_lris_blue':
                good = centers < 0.
                xstart = max_spat + 30//bin_spat  
            else:
                msgs.error(f'Not ready to use slitmasks for {self.name}.  Develop it!')
        else:
            if self.name in ['keck_lris_red', 'keck_lris_blue']:
                good = centers >= 0.
                xstart = -48//bin_spat
            else:             
                msgs.error(f'Not ready to use slitmasks for {self.name}.  Develop it!')
        left_edges += xstart
        right_edges += xstart
        left_edges[~good] = -1

        # Toss any left edges off the right-side of the detector
        keep = left_edges < max_spat
        left_edges[~keep] = -1

        right_edges[left_edges == -1] = -1
        # Deal with right edge off the detector
        for islit in range(self.slitmask.nslits):
            if left_edges[islit] != -1 and right_edges[islit] > max_spat:
                right_edges[islit] = max_spat
        # Now the left
        for islit in range(self.slitmask.nslits):
            if right_edges[islit] != -1 and left_edges[islit] < 0:
                left_edges[islit] = 0

        # Order from left to right
        # Highest x is leftmost on DET=1
        x_order = np.argsort(self.slitmask.corners[:,1,0])
        sortindx = x_order[::-1]

        # Return
        return left_edges.astype(float), right_edges.astype(float), sortindx, self.slitmask



class KeckLRISBSpectrograph(KeckLRISSpectrograph):
    """
    Child to handle Keck/LRISb specific code
    """

    name = 'keck_lris_blue'
    camera = 'LRISb'
    header_name = 'LRISBLUE'
    supported = True
    comment = 'Blue camera; see :doc:`lris`'
    
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
        # Binning
        # TODO: Could this be detector dependent?
        binning = '1,1' if hdu is None else self.get_meta_value(self.get_headarr(hdu), 'binning')

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
        detector = detector_container.DetectorContainer(**detector_dicts[det-1])

        if hdu is None:
            return detector

        # TODO: I don't know how to handle the hack below for the auto-generated
        # detector table...

        # Deal with number of amps
        namps = hdu[0].header['NUMAMPS']
        # The website does not give values for single amp per detector so we take the mean
        #   of the values provided
        if namps == 2 or (namps==4 and len(hdu)==3):  # Longslit readout mode is the latter.  This is a hack..
            detector.numamplifiers = 1
            detector.gain = np.atleast_1d(np.mean(detector.gain))
            detector.ronoise = np.atleast_1d(np.mean(detector.ronoise))
        elif namps == 4:
            pass
        else:
            msgs.error("Did not see this namps coming..")

        # Return
        return detector

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

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

        # Allow for longer exposure times on blue side (especially if using the Dome lamps)
        par['calibrations']['pixelflatframe']['exprng'] = [None, 300]
        par['calibrations']['traceframe']['exprng'] = [None, 300]

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
        # Start with instrument wide
        par = super().config_specific_par(scifile, inp_par=inp_par)

        # Wavelength calibrations
        if self.get_meta_value(scifile, 'dispname') == '300/5000':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_lris_blue_300_d680.fits'
            par['flexure']['spectrum'] = 'sky_LRISb_400.fits'
        elif self.get_meta_value(scifile, 'dispname') == '400/3400':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_lris_blue_400_d560.fits'
            par['flexure']['spectrum'] = 'sky_LRISb_400.fits'
        elif self.get_meta_value(scifile, 'dispname') == '600/4000':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_lris_blue_600_d560.fits'
            par['flexure']['spectrum'] = 'sky_LRISb_600.fits'
        elif self.get_meta_value(scifile, 'dispname') == '1200/3400':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_lris_blue_1200_d460.fits'
            par['flexure']['spectrum'] = 'sky_LRISb_600.fits'

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
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        super().init_meta()
        # Add the name of the dispersing element
        self.meta['dispname'] = dict(ext=0, card='GRISNAME')

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
        return ['GRISNAME', 'DICHNAME', 'SLITNAME']

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
                Processed bias frame used to identify bad pixels

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """
        # Call the base-class method to generate the empty bpm
        bpm_img = super().bpm(filename, det, shape=shape, msbias=msbias)

        # Only defined for det=1
        if det == 1:
            msgs.info("Using hard-coded BPM for det=1 on LRISb")
            bpm_img[:,:3] = 1

        return bpm_img


class KeckLRISBOrigSpectrograph(KeckLRISBSpectrograph):
    """
    Child to handle the LRISb detector packed prior to 01 JUL 2009
    """
    ndet = 2
    name = 'keck_lris_blue_orig'
    camera = 'LRISb'
    supported = True    # TODO: Is this true?
    comment = 'Original detector; replaced in 20??; see :doc:`lris`'

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        super().init_meta()
        # Remove the lamps
        keys = list(self.meta.keys())
        for key in keys:
            if 'lampstat' in key:
                self.meta.pop(key)

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
        detpar = super().get_detector_par(det, hdu=hdu)
        detpar['specflip'] = True
        return detpar

    def get_rawimage(self, raw_file, det):
        """
        Read raw images and generate a few other bits and pieces
        that are key for image processing.

        Over-ride standard get_rawimage() for LRISb to deal
        with the original approach

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
        # Image info
        image, hdul, elaptime, rawdatasec_img, oscansec_img = get_orig_rawimage(raw_file)
        # Cut down
        if np.max(rawdatasec_img) != 4:
            msgs.error("Deal with not 2 AMP mode!!")
        if det == 1:
            bad_amp = rawdatasec_img > 2
            rawdatasec_img[bad_amp] = 0
            bad_amp = oscansec_img > 2
            oscansec_img[bad_amp] = 0
        elif det == 2:
            # Kludge this to be 1 and 2's
            for timage in [rawdatasec_img, oscansec_img]:
                # Zero out
                bad_amp = timage <= 2
                timage[bad_amp] = 0
                # Offset
                good_amp = timage > 2
                timage[good_amp] -= 2
        else:
            msgs.error("Should not be here in keck_lris!")

        # Detector
        detector_par = self.get_detector_par(det-1, hdu=hdul)

#        # Flip the spectral axis
#        detector_par['specflip'] = True

        # Return
        return detector_par, image, hdul, elaptime, rawdatasec_img, oscansec_img


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
                Processed bias frame used to identify bad pixels. **This is
                always ignored.**

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """
        # Call the base-class method to generate the empty bpm
        return super().bpm(filename, det, shape=shape, msbias=None)


class KeckLRISRSpectrograph(KeckLRISSpectrograph):
    """
    Child to handle Keck/LRISr specific code
    """
    name = 'keck_lris_red'
    camera = 'LRISr'
    header_name = 'LRIS'
    supported = True
    ql_supported = True
    comment = 'Red camera;  LBNL detector, 2kx4k; see :doc:`lris`'
    
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
        # Binning
        # TODO: Could this be detector dependent??
        binning = '1,1' if hdu is None else self.get_meta_value(self.get_headarr(hdu), 'binning')

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

        if hdu is not None:
            # Allow for post COVID detector issues
            t2020_1 = time.Time("2020-06-30", format='isot')  # First run
            t2020_2 = time.Time("2020-07-29", format='isot')  # Second run
            # Check for the new detector (Mark4) upgrade
            t2021_upgrade = time.Time("2021-04-15", format='isot') 
            date = time.Time(self.get_meta_value(self.get_headarr(hdu), 'mjd'), 
                             format='mjd')

            if date < t2020_1:
                pass
            elif date < t2020_2: # This is for the June 30 2020 run
                msgs.warn("We are using LRISr gain/RN values based on WMKO estimates.")
                detector_dict1['gain'] = np.atleast_1d([37.6])
                detector_dict2['gain'] = np.atleast_1d([1.26])
                detector_dict1['ronoise'] = np.atleast_1d([99.])
                detector_dict2['ronoise'] = np.atleast_1d([5.2])
            elif date > t2021_upgrade: 
                # Note:  We are unlikely to trip this.  Other things probably failed first
                msgs.error("This is the new detector.  Use keck_lris_red_mark4")
            else: # This is the 2020 July 29 run
                msgs.warn("We are using LRISr gain/RN values based on WMKO estimates.")
                detector_dict1['gain'] = np.atleast_1d([1.45])
                detector_dict2['gain'] = np.atleast_1d([1.25])
                detector_dict1['ronoise'] = np.atleast_1d([4.47])
                detector_dict2['ronoise'] = np.atleast_1d([4.75])

        # Instantiate
        detector_dicts = [detector_dict1, detector_dict2]
        detector = detector_container.DetectorContainer(**detector_dicts[det-1])

        if hdu is None:
            return detector

        # Deal with number of amps
        namps = hdu[0].header['NUMAMPS']
        # The website does not give values for single amp per detector so we take the mean
        #   of the values provided
        if namps == 2 or (namps==4 and len(hdu)==3):  # Longslit readout mode is the latter.  This is a hack..
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



    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

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

        # Sensitivity function defaults
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 9

        return par

    # NOTE: This function is used by the dev-suite
    def get_ql_calib_dir(self, file):
        """
        Returns calibration file directory for quicklook reductions.

        Args:
            file (str):
              Image file

        Returns:
            :obj:`str`: Quicklook calibrations directory

        """
        lris_grating = self.get_meta_value(file, 'dispname')
        lris_dichroic = self.get_meta_value(file, 'dichroic')
        setup_path = lris_grating.replace('/','_') + '_d' + lris_dichroic
        return os.path.join(self.name, setup_path)

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
        # Start with instrument wide
        par = super().config_specific_par(scifile, inp_par=inp_par)

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
            if self.name == 'keck_lris_red_mark4':
                par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_lris_red_mark4_R400.fits'
            else:
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
        elif self.get_meta_value(scifile, 'dispname') == '600/10000':  # d680
            par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_lris_red_600_10000.fits'
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
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        super().init_meta()
        # Add the name of the dispersing element
        self.meta['dispname'] = dict(ext=0, card='GRANAME')

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
        return super().configuration_keys() + ['dispangle', 'cenwave', 'amp', 'binning']

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
        return ['GRANAME', 'DICHNAME', 'SLITNAME', 'GRANGLE', 'WAVELEN', 'TAPLINES', 'NUMAMPS']

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
                Processed bias frame used to identify bad pixels.

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """
        # Call the base-class method to generate the empty bpm
        bpm_img = super().bpm(filename, det, shape=shape, msbias=msbias)

        # Only defined for det=2
        if det == 2:
            msgs.info("Using hard-coded BPM for det=2 on LRISr")

            # Get the binning
            hdu = io.fits_open(filename)
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

class KeckLRISRMark4Spectrograph(KeckLRISRSpectrograph):
    """
    Child to handle the new Mark4 detector
    """
    ndet = 1
    name = 'keck_lris_red_mark4'
    supported = True
    comment = 'New Mark4 detector, circa Spring 2021; Supported setups = R400'

    def init_meta(self):
        super().init_meta()
        # Over-ride a pair
        self.meta['mjd'] = dict(ext=0, card='MJD')
        self.meta['exptime'] = dict(ext=0, card='TELAPSE')

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
        # Binning
        binning = '1,1' if hdu is None else self.get_meta_value(self.get_headarr(hdu), 'binning')

        # Detector 1
        detector_dict1 = dict(
            binning=binning,
            det=1,
            dataext=0,
            specaxis=0,
            specflip=True,  
            spatflip=False,
            platescale=0.123,  # From the web page
            darkcurr=0.0,
            saturation=65535.,
            nonlinear=0.76,
            mincounts=-1e10,
            numamplifiers=2,  # These are defaults but can modify below
            gain=np.atleast_1d([1.61, 1.67*0.959]), # Assumes AMPMODE=HSPLIT,VUP;  Corrected by JXP using 2x1 binned flats
            ronoise=np.atleast_1d([3.65, 3.52]),
        )

        if hdu is None:
            return detector_container.DetectorContainer(**detector_dict1)

        # Date of Mark4 installation
        t2021_upgrade = time.Time("2021-04-15", format='isot') 
        # TODO -- Update with the date we transitioned to the correct ones
        t_gdhead = time.Time("2029-01-01", format='isot')
        date = time.Time(hdu[0].header['MJD'], format='mjd')

        if date < t2021_upgrade:
            msgs.error("This is not the Mark4 detector.  Use a different keck_lris_red spectrograph")

        # Deal with the intermediate headers
        if date < t_gdhead:
            amp_mode = hdu[0].header['AMPMODE']
            msgs.info("AMPMODE = {:s}".format(amp_mode))
            # Load up translation dict
            ampmode_translate_file = (
                data.Paths.data / 'spectrographs' / 'keck_lris_red_mark4' / 'dict_for_ampmode.json'
            )
            # Force any possible pathlib.Path object to string before `loadjson`
            ampmode_translate_dict = linetools.utils.loadjson(str(ampmode_translate_file))
            # Load up the corrected header
            swap_binning = f"{binning[-1]},{binning[0]}"  # LRIS convention is oppopsite ours
            header_file = (
                data.Paths.data /
                'spectrographs' /
                'keck_lris_red_mark4' /
                f'header{ampmode_translate_dict[amp_mode]}_{swap_binning.replace(",","_")}.fits'
            )
            correct_header = fits.getheader(header_file)
        else:
            correct_header = hdu[0].header

        # Deal with number of amps
        detector_dict1['numamplifiers'] = correct_header['TAPLINES']

        # The website does not give values for single amp per detector so we take the mean
        #   of the values provided
        if detector_dict1['numamplifiers'] == 2: 
            pass
        elif detector_dict1['numamplifiers'] == 4:
            # From the web page on 2021-10-04 (L1, L2, U1, U2)
            #  Corrected by JXP and SS using chk_lris_mark4_gain.py in the DevSuite
            detector_dict1['gain'] = np.atleast_1d([1.710, 
                                                    1.64*1.0245,  # L2
                                                    1.61*1.0185,  # U1
                                                    1.67*1.0052]) # U2
            detector_dict1['ronoise'] = np.atleast_1d([3.64, 3.45, 3.65, 3.52])
        else:
            msgs.error("Did not see this namps coming..")

        detector_dict1['datasec'] = []
        detector_dict1['oscansec'] = []

        # Parse which AMPS were used
        used_amps = []
        for amp in range(4):
            if f'AMPNM{amp}' in correct_header.keys():
                used_amps.append(amp)
        # Check
        assert detector_dict1['numamplifiers'] == len(used_amps)

        # Reverse engenieering to translate LRIS DSEC, BSEC
        #  into ones friendly for PypeIt...
        binspec = int(binning[0])
        binspatial = int(binning[-1])
        
        for iamp in used_amps:
            # These are column, row
            dsecs = correct_header[f'DSEC{iamp}'].split(',')
            d_rows = [int(item) for item in dsecs[1][:-1].split(':')]
            d_cols = [int(item) for item in dsecs[0][1:].split(':')]
            bsecs = correct_header[f'BSEC{iamp}'].split(',')
            o_rows = [int(item) for item in bsecs[1][:-1].split(':')]
            o_cols = [int(item) for item in bsecs[0][1:].split(':')]

            # Deal with binning (heaven help me!!)
            d_rows = [str(item*binspec) if item != 1 else str(item) for item in d_rows]
            o_rows = [str(item*binspec) if item != 1 else str(item) for item in o_rows]
            d_cols = [str(item*binspatial) if item != 1 else str(item) for item in d_cols]
            o_cols = [str(item*binspatial) if item != 1 else str(item) for item in o_cols]

            # These are now row, column
            #  And they need to be native!!  i.e. no binning accounted for..
            detector_dict1['datasec'] += [f"[{':'.join(d_rows)},{':'.join(d_cols)}]"]
            detector_dict1['oscansec'] += [f"[{':'.join(o_rows)},{':'.join(o_cols)}]"]

        detector_dict1['datasec'] = np.array(detector_dict1['datasec'])
        detector_dict1['oscansec'] = np.array(detector_dict1['oscansec'])

        # Instantiate
        detector = detector_container.DetectorContainer(**detector_dict1)

        #print('Dict1:', detector_dict1)
        #print('Binning:', binning)

        # Return
        return detector

    def get_rawimage(self, raw_file, det):
        """
        Read raw images and generate a few other bits and pieces
        that are key for image processing.

        Over-ride standard get_rawimage() for LRIS

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
        # Note:  There is no way we know to super super super
        return spectrograph.Spectrograph.get_rawimage(self, raw_file, det)

class KeckLRISROrigSpectrograph(KeckLRISRSpectrograph):
    """
    Child to handle the original LRISr detector (pre 01 JUL 2009)
    """
    ndet = 1
    name = 'keck_lris_red_orig'
    camera = 'LRISr'
    supported = True
    comment = 'Original detector; replaced in 2009; see :doc:`lris`'

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
        par['calibrations']['wavelengths']['lamps'] = ['NeI', 'ArI', 'KrI', 'XeI', 'HgI']

        return par

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
        # Binning
        # TODO: Could this be detector dependent??
        binning = '1,1' if hdu is None else self.get_meta_value(self.get_headarr(hdu), 'binning')

        # Detector 1
        detector_dict1 = dict(
            binning=binning,
            det=1,
            dataext=1,
            specaxis=1,
            specflip=False,
            spatflip=False,
            platescale=0.21,  # TO BE UPDATED!!
            darkcurr=0.0,
            saturation=65535.,
            nonlinear=0.76,
            mincounts=-1e10,
            numamplifiers=2,
            gain=np.atleast_1d([1.98, 2.17]),    # TAKEN FROM LOWREDUX
            ronoise=np.atleast_1d([6.1, 6.3]),   # TAKEN FROM LOWREDUX
        )
        # Instantiate
        detector = detector_container.DetectorContainer(**detector_dict1)

        # Deal with number of amps
        if hdu is not None and hdu[0].header['NUMAMPS'] != 2:
            msgs.error("Did not see this namps coming..")

        # Return
        return detector

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        super().init_meta()
        # Remove the lamps
        keys = list(self.meta.keys())
        for key in keys:
            if 'lampstat' in key:
                self.meta.pop(key)

    def get_rawimage(self, raw_file, det):
        """
        Read raw images and generate a few other bits and pieces
        that are key for image processing.

        Over-ride standard get_rawimage() for LRIS

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
        # Image info
        image, hdul, elaptime, rawdatasec_img, oscansec_img = get_orig_rawimage(raw_file)
        # Detector
        detector_par = self.get_detector_par(det, hdu=hdul)
        # Return
        return detector_par, image, hdul, elaptime, rawdatasec_img, oscansec_img

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
                Processed bias frame used to identify bad pixels. **This is
                always ignored.**

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """
        # Call the base-class method to generate the empty bpm
        return super().bpm(filename, det, shape=shape, msbias=None)


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
        hdu = io.fits_open(inp)
    else:
        hdu = inp
    # Count the number of extensions
    n_ext = np.sum(['VidInp' in h.name for h in hdu])

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
    hdu = io.fits_open(infil)
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


def get_orig_rawimage(raw_file, debug=False):
    """
    Read a raw, original LRIS data frame.

    Ported from LOWREDUX long_oscan.pro lris_oscan()

    Parameters
    ----------
    raw_file : :obj:`str`
        Filename
    debug : :obj:`bool`, optional
        Run in debug mode (doesn't do anything)

    Returns
    -------
    raw_img : `numpy.ndarray`_
        Raw image for this detector.
    hdu : `astropy.io.fits.HDUList`_
        Opened fits file
    exptime : :obj:`float`
        Exposure time read from the file header
    rawdatasec_img : `numpy.ndarray`_
        Data (Science) section of the detector as provided by setting the
        (1-indexed) number of the amplifier used to read each detector pixel.
        Pixels unassociated with any amplifier are set to 0.
    oscansec_img : `numpy.ndarray`_
        Overscan section of the detector as provided by setting the
        (1-indexed) number of the amplifier used to read each detector pixel.
        Pixels unassociated with any amplifier are set to 0.
    """
    # Open
    hdul = io.fits_open(raw_file)
    head0 = hdul[0].header
    # TODO -- Check date here and error/warn if not after the upgrade
    image = hdul[0].data.astype(float)

    # Get post, pre-pix values
    postpix = head0['POSTPIX']
    prepix = head0['PREPIX']
    post_buffer1 = 4
    post_buffer2 = 8
    namps = head0['NUMAMPS']

    # get the x and y binning factors...
    binning = head0['BINNING']
    xbin, ybin = [int(ibin) for ibin in binning.split(',')]

    rawdatasec_img = np.zeros_like(image, dtype=int)
    oscansec_img = np.zeros_like(image, dtype=int)

    datacol = namps * (prepix // xbin) + np.arange(namps) * 1024 // xbin
    postcol = datacol[namps - 1] + (1024 + post_buffer1) // xbin
    for iamp in range(namps): #= 0, namps - 1L
        biascols = np.arange((postpix - post_buffer2) // xbin) + (
                iamp * postpix) // xbin + postcol
        oscansec_img[:, biascols] = iamp+1
        imagecols = np.arange(1024 // xbin) + iamp * 1024 // xbin
        rawdatasec_img[:,imagecols + namps*(prepix // xbin)] = iamp+1
    return image, hdul, float(head0['ELAPTIME']), \
           rawdatasec_img, oscansec_img


