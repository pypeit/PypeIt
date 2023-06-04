"""
Module for Keck/MOSFIRE specific methods.

.. include:: ../include/links.rst
"""
import copy
import os
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch, meta
from pypeit import utils
from pypeit import io
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container
from scipy import special
from pypeit.spectrographs.slitmask import SlitMask

from pypeit.utils import index_of_x_eq_y

from IPython import embed

class KeckMOSFIRESpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Keck/MOSFIRE specific code
    """
    ndet = 1
    name = 'keck_mosfire'
    telescope = telescopes.KeckTelescopePar()
    camera = 'MOSFIRE'
    url = 'https://www2.keck.hawaii.edu/inst/mosfire/home.html'
    header_name = 'MOSFIRE'
    supported = True
    comment = 'Gratings tested: Y, J, J2, H, K; see :doc:`mosfire`'

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
            binning         = '1,1',
            det             = 1,
            dataext         = 0,
            specaxis        = 1,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.1798,
            darkcurr        = 0.8,
            saturation      = 1e9, # ADU, this is hacked for now
            nonlinear       = 1.00,  # docs say linear to 90,000 but our flats are usually higher
            numamplifiers   = 1,
            mincounts       = -1e10,
            gain            = np.atleast_1d(2.15),  # Taken from MOSFIRE detector webpage
            ronoise         = np.atleast_1d(5.8), # This is for 16 non-destructuve reads, the default readout mode
            datasec         = np.atleast_1d('[5:2044,5:2044]'),
            #oscansec        = np.atleast_1d('[:,:]')
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

        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.30 #0.20  # Might be grating dependent..
        par['calibrations']['wavelengths']['sigdetect']=5.0
        par['calibrations']['wavelengths']['fwhm']= 5.0
        par['calibrations']['wavelengths']['n_final']= 4
        par['calibrations']['wavelengths']['lamps'] = ['OH_NIRES']
        #par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['method'] = 'holy-grail'
        # Reidentification parameters
        #par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_nires.fits'
        par['calibrations']['slitedges']['edge_thresh'] = 50.
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Flats
        # Do not illumination correct. We should also not be flat fielding given the bars.
        # TODO Implement imaging flats for MOSFIRE. Do test with/without illumination flats.
        # Turn of illumflat
        turn_off = dict(use_biasimage=False, use_overscan=False, use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)

        # Extraction
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['extraction']['sn_gauss'] = 4.0

        # Flexure
        par['flexure']['spec_method'] = 'skip'

        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] ='nothing'

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, 20]
        par['calibrations']['arcframe']['exprng'] = [1, None]
        par['calibrations']['darkframe']['exprng'] = [1, None]
        par['scienceframe']['exprng'] = [20, None]

        # Sensitivity function parameters
        par['sensfunc']['extrap_blu'] = 0.0  # Y-band contaminated by higher order so don't extrap much
        par['sensfunc']['extrap_red'] = 0.0
        par['fluxcalib']['extrap_sens'] = True
        par['sensfunc']['extrap_red'] = 0.0
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 13
        par['sensfunc']['IR']['maxiter'] = 2
        par['sensfunc']['IR']['telgridfile'] = 'TelFit_MaunaKea_3100_26100_R20000.fits'
        return par

    # NOTE: This function is used by the dev-suite
    def get_ql_calib_dir(self, file):
        """
        Returns calibrations file directory for quicklook reductions.

        Args:
            file (str):
              Image file

        Returns:
            :obj:`str`: Quicklook calibrations directory

        """
        mosfire_filter = self.get_meta_value(file, 'filter1')
        return os.path.join(self.name, mosfire_filter)

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
        decker = self.get_meta_value(headarr, 'decker')

        if 'LONGSLIT' in decker:
            # turn PCA off
            par['calibrations']['slitedges']['sync_predict'] = 'nearest'
            # if "x" is not in the maskname, the maskname does not include the number of CSU
            # used for the longslit and the length of the longslit cannot be determined
            if ('LONGSLIT-46x' not in decker) and ('x' in decker):
                # find the spat pixel positions where the longslit starts and ends
                pix_start, pix_end = self.find_longslit_pos(scifile)
                # exclude the random slits outside the longslit from slit tracing
                par['calibrations']['slitedges']['exclude_regions'] = ['1:0:{}'.format(pix_start),
                                                                       '1:{}:2040'.format(pix_end)]
                par['calibrations']['slitedges']['det_buffer'] = 0
                # artificially add left and right edges
                par['calibrations']['slitedges']['bound_detector'] = True
            # set offsets for coadd2d
            par['coadd2d']['offsets'] = 'header'

        # Turn on the use of mask design
        else:
            par['calibrations']['slitedges']['use_maskdesign'] = True
            # use dither info in the header as the default offset
            par['reduce']['slitmask']['use_dither_offset'] = True
            # Assign RA, DEC, OBJNAME to detected objects
            par['reduce']['slitmask']['assign_obj'] = True
            # force extraction of undetected objects
            par['reduce']['slitmask']['extract_missing_objs'] = True
            if 'long2pos' in decker:
                # exclude the random slits outside the long2pos from slit tracing
                pix_start, pix_end = self._long2pos_pos()
                par['calibrations']['slitedges']['exclude_regions'] = ['1:0:{}'.format(pix_start),
                                                                       '1:{}:2040'.format(pix_end)]
                # assume that the main target is always detected, i.e., skipping force extraction
                par['reduce']['slitmask']['extract_missing_objs'] = False
            # set offsets for coadd2d
            par['coadd2d']['offsets'] = 'maskdef_offsets'

        # wavelength calibration
        supported_filters = ['Y', 'J', 'J2', 'H', 'K']
        filter = self.get_meta_value(headarr, 'filter1')
        # using OH lines
        if 'long2pos_specphot' not in decker and filter in supported_filters:
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['fwhm_fromlines'] = True
            par['calibrations']['wavelengths']['sigdetect'] = 10.
            # templates
            if filter == 'Y':
                par['calibrations']['wavelengths']['lamps'] = ['OH_MOSFIRE_Y']
                par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_mosfire_OH_Y.fits'
            elif filter == 'J':
                par['calibrations']['wavelengths']['lamps'] = ['OH_MOSFIRE_J']
                par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_mosfire_OH_J.fits'
            elif filter == 'J2':
                par['calibrations']['wavelengths']['lamps'] = ['OH_MOSFIRE_J']
                par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_mosfire_OH_J2.fits'
            elif filter == 'H':
                par['calibrations']['wavelengths']['lamps'] = ['OH_MOSFIRE_H']
                par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_mosfire_OH_H.fits'
            elif filter == 'K':
                par['calibrations']['wavelengths']['lamps'] = ['OH_MOSFIRE_K']
                par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_mosfire_OH_K.fits'

        # using arc lines (we use this as default only for long2pos_specphot mask)
        elif 'long2pos_specphot' in decker and filter in supported_filters:
            par['calibrations']['wavelengths']['lamps'] = ['Ar_IR_MOSFIRE', 'Ne_IR_MOSFIRE']
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['fwhm_fromlines'] = True
            # templates
            if filter == 'Y':
                par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_mosfire_arcs_Y.fits'
            elif filter == 'J':
                par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_mosfire_arcs_J.fits'
            elif filter == 'J2':
                par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_mosfire_arcs_J2.fits'
            elif filter == 'H':
                par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_mosfire_arcs_H.fits'
            elif filter == 'K':
                par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_mosfire_arcs_K.fits'

        # Return
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
        self.meta['decker'] = dict(ext=0, card='MASKNAME')
        self.meta['binning'] = dict(ext=0, card=None, default='1,1')

        self.meta['mjd'] = dict(ext=0, card='MJD-OBS')
        self.meta['exptime'] = dict(ext=0, card='TRUITIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='OBSMODE')
        self.meta['idname'] = dict(card=None, compound=True)
        self.meta['frameno'] = dict(ext=0, card='FRAMENUM')
        self.meta['object'] = dict(ext=0, card='OBJECT')
        # The following 3 metas (decker_secondary, slitwid, slitlength) are introduced
        # only to reduce data (LONGSLIT and long2pos) with calibrations taken with
        # a different decker ('MASKNAME')
        # decker_secondary is different than decker (MASKNAME) only for 'LONGSLIT' masks
        self.meta['decker_secondary'] = dict(card=None, compound=True)
        # slit width, defined only for 'LONGSLIT' masks
        self.meta['slitwid'] = dict(card=None, compound=True, rtol=0.1)
        # slit length in numbers of CSU, defined only for only for 'LONGSLIT' masks
        self.meta['slitlength'] = dict(card=None, compound=True, rtol=0.1)
        # Filter
        self.meta['filter1'] = dict(ext=0, card='FILTER')
        # Lamps on/off or Ar/Ne
        self.meta['lampstat01'] = dict(card=None, compound=True)

        # Dithering
        self.meta['dithpat'] = dict(ext=0, card='PATTERN')
        self.meta['dithpos'] = dict(ext=0, card='FRAMEID')
        self.meta['dithoff'] = dict(card=None, compound=True)
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')

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
        if meta_key == 'decker_secondary':
            # decker_secondary is different than decker (MASKNAME) only for 'LONGSLIT' masks
            maskname = headarr[0].get('MASKNAME')
            if 'LONGSLIT' in maskname:
                return maskname.split('(')[0].split('-')[0]
            else:
                return maskname
        if meta_key == 'slitlength':
            # slitlength is defined only for 'LONGSLIT' masks since this info is generally
            # included in the slitmask name (MASKNAME) of 'LONGSLIT' masks and
            # it's useful to associate science frames to calibrations taken with different MASKNAME
            maskname = headarr[0].get('MASKNAME')
            if 'LONGSLIT' in maskname and 'x' in maskname:
                return maskname.split('(')[0].split('x')[0].split('-')[1]
            else:
                return None

        if meta_key == 'slitwid':
            # slitwid is defined only for 'LONGSLIT' masks since this info is generally
            # included in the slitmask name (MASKNAME) of 'LONGSLIT' masks and
            # it's useful to associate science frames to calibrations taken with different MASKNAME
            maskname = headarr[0].get('MASKNAME')
            if 'LONGSLIT' in maskname and 'x' in maskname:
                return maskname.split('(')[0].split('x')[1]
            else:
                return None

        if meta_key == 'idname':
            FLATSPEC = headarr[0].get('FLATSPEC')
            PWSTATA7 = headarr[0].get('PWSTATA7')
            PWSTATA8 = headarr[0].get('PWSTATA8')
            if FLATSPEC == 0 and PWSTATA7 == 0 and PWSTATA8 == 0:
                if 'Flat' in headarr[0].get('OBJECT'):
                    return 'flatlampoff'
                else:
                    return 'object'
            elif FLATSPEC == 0 and headarr[0].get('FILTER') == 'Dark':
                return 'dark'
            elif FLATSPEC == 1:
                return 'flatlamp'
            elif PWSTATA7 == 1 or PWSTATA8 == 1:
                return 'arclamp'
            else:
                msgs.warn('Header keyword FLATSPEC, PWSTATA7, or PWSTATA8 may not exist')
                return 'unknown'
        if meta_key == 'lampstat01':
            if headarr[0].get('PWSTATA7') == 1 or headarr[0].get('PWSTATA8') == 1:
                lamps = []
                if headarr[0].get('PWSTATA7') == 1:
                    lamps.append(headarr[0].get('PWLOCA7')[:2])
                if headarr[0].get('PWSTATA8') == 1:
                    lamps.append(headarr[0].get('PWLOCA8')[:2])
                return ','.join(lamps)
            elif headarr[0].get('FLATSPEC') == 1 or headarr[0].get('FLSPECTR') == 'on':
                return 'on'
            else:
                return 'off'

        if meta_key == 'dithoff':
            if headarr[0].get('YOFFSET') is not None:
                return headarr[0].get('YOFFSET')
            else:
                return 0.0
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
        return ['decker_secondary', 'slitlength', 'slitwid', 'dispname', 'filter1']

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
        return ['MASKNAME', 'OBSMODE', 'FILTER']

    def modify_config(self, row, cfg):
        """

        Modify the configuration dictionary for a given frame. This method is
        used in :func:`~pypeit.metadata.PypeItMetaData.set_configurations` to
        modify in place the configuration requirement to assign a specific frame
        to the current setup.

        This is needed for the reduction of 'LONGSLIT' and 'long2pos' data,
        which often use calibrations taken with a different decker (MASKNAME).

            - For the 'LONGSLIT' masks, when we are assigning a configuration to
              a calibration file that was taken with the longest slit available
              (46 CSUs), since these calibrations are generally used for the
              reduction of science frames with shorter slits, we remove the
              configuration requirement on the slit length for the current file.

            - For the 'long2pos' masks, when we are assigning a configuration to
              a calibration file that was taken with the 'long2pos' mask, since
              these calibrations are generally used for the reduction of science
              frames taken with 'long2pos_specphot' masks, we modify the
              configuration requirement on the decker_secondary for the current
              file.

        Args:
            row (`astropy.table.Row`_):
                The table row with the metadata for one frame.
            cfg (:obj:`dict`):
                Dictionary with metadata associated to a specific configuration.

        Returns:
            :obj:`dict`: modified dictionary with metadata associated to a
            specific configuration.
        """
        if row['decker'] is not None \
                and cfg['decker_secondary'] is not None \
                and 'LONGSLIT' in row['decker'] \
                and 'LONGSLIT' in cfg['decker_secondary'] \
                and 'science' not in row['frametype'] \
                and 'standard' not in row['frametype'] \
                and row['slitlength'] == 46.:
            cfg2 = copy.deepcopy(cfg)
            cfg2.pop('slitlength')
            return cfg2
        if row['decker'] is not None \
                and cfg['decker_secondary'] is not None \
                and 'long2pos' in row['decker'] \
                and 'long2pos' in cfg['decker_secondary'] \
                and 'science' not in row['frametype'] \
                and 'standard' not in row['frametype']:
            cfg2 = copy.deepcopy(cfg)
            cfg2['decker_secondary'] = 'long2pos'
            return cfg2
        return cfg

    def get_comb_group(self, fitstbl):
        """
        Automatically assign combination groups and background images by parsing
        known dither patterns.

        This method is used in
        :func:`~pypeit.metadata.PypeItMetaData.set_combination_groups`, and
        directly modifies the ``comb_id`` and ``bkg_id`` columns in the provided
        table.

        Specifically here, this method parses the dither pattern of the
        science/standard frames in a given calibration group and assigns to each
        of them a comb_id and a bkg_id. The known dither patterns are: "Slit
        Nod", "Mask Nod", "ABA'B'", "ABAB", "ABBA", "long2pos_specphot", and
        "Stare". Note that the frames in the same dither positions (A positions
        or B positions) of each "ABAB" or "ABBA" sequence are 2D coadded
        (without optimal weighting) before the background subtraction, while for
        the other dither patterns, the frames in the same dither positions are
        not coadded.

        For "long2pos_specphot" masks, the ``comb_id`` and a ``bkg_id`` are
        assigned such that one of the two frames with spectra taken using the
        narrower slit is used as the background frame and subtracted from the
        frame with spectra taken using the wider slit.

        Args:
            fitstbl(`astropy.table.Table`_):
                The table with the metadata for all the frames.

        Returns:
            `astropy.table.Table`_: modified fitstbl.
        """
        #TODO incorporate parse_dither_pattern() here.

        # find index of fitstbl that contains science and standard frames
        # where science
        sci_idx = np.array(['science' in _tab for _tab in fitstbl['frametype']])
        # where standard
        std_idx = np.array(['standard' in _tab for _tab in fitstbl['frametype']])

        sci_std_idx = [sci_idx, std_idx]
        # loop over the science and standard frames
        for idx in sci_std_idx:
            setups = np.unique(fitstbl[idx]['setup'])
            # loop over the setups
            for setup in setups:
                in_cfg = idx & np.array([setup in _set for _set in fitstbl['setup']])
                if len(fitstbl[in_cfg]) == 1:
                    continue
                # how many dither patterns are used for the selected science/standard frames?
                uniq_dithpats = np.unique(fitstbl[in_cfg]['dithpat'])
                # loop through the dither patterns
                for dpat in uniq_dithpats:
                    if dpat == 'none':
                        continue
                    # where this dpat
                    dpat_idx = in_cfg & (fitstbl['dithpat'] == dpat) & \
                               np.array([dithpos in ["A", "B", "A'", "B'"] for dithpos in fitstbl['dithpos']])

                    # compute comb_id
                    if len(fitstbl[dpat_idx]) > 1:
                        # get default combid and bkgid
                        combid = np.copy(fitstbl['comb_id'][dpat_idx].data)
                        bkgid = np.copy(fitstbl['bkg_id'][dpat_idx].data)
                        dpos = fitstbl[dpat_idx]['dithpos']

                        if "long2pos_specphot" in fitstbl[dpat_idx]['decker']:
                            doff = fitstbl[dpat_idx]['dithoff']
                            # find the starting index of the BAA sequence
                            dpos_idx = np.where((dpos == "A") & (np.roll(dpos, -1) == "B") &
                                                (np.roll(dpos, -2) == "A"))[0]
                            for i in dpos_idx:
                                # make sure that that dither offsets are correct
                                if i < len(dpos)-2 and doff[i] == 0. and abs(doff[i+1]) > 0. and doff[i+1] == -doff[i+2]:
                                    bkgid[i] = combid[i+1]
                                    bkgid[i+1] = combid[i+2]
                                    bkgid[i+2] = combid[i+1]

                        elif "long2pos" in fitstbl[dpat_idx]['decker']:
                            # find the starting index of the BA sequence
                            dpos_idx = np.where((dpos == "B") & (np.roll(dpos, -1) == "A"))[0]
                            for i in dpos_idx:
                                # exclude when np.roll counts the 1st element of dpos to be in a
                                # sequence with the last element
                                if i < len(dpos) - 1:
                                    bkgid[i] = combid[i + 1]
                                    bkgid[i + 1] = combid[i]

                        elif dpat in ["Slit Nod", "Mask Nod"]:
                            # find the starting index of the AB sequence
                            dpos_idx = np.where((dpos == "A") & (np.roll(dpos, -1) == "B"))[0]
                            for i in dpos_idx:
                                # exclude when np.roll counts the 1st element of dpos to be in a
                                # sequence with the last element
                                if i < len(dpos)-1:
                                    bkgid[i] = combid[i+1]
                                    bkgid[i+1] = combid[i]

                        elif dpat == "ABA'B'":
                            # find the starting index of the ABA'B' sequence
                            dpos_idx = np.where((dpos == "A") & (np.roll(dpos, -1) == "B") &
                                                (np.roll(dpos, -2) == "A'") & (np.roll(dpos, -3) == "B'"))[0]
                            for i in dpos_idx:
                                if i < len(dpos) - 3:
                                    bkgid[i] = combid[i+1]
                                    bkgid[i+1] = combid[i]
                                    bkgid[i+2] = bkgid[i+3]
                                    bkgid[i+3] = bkgid[i+2]

                        elif dpat == "ABAB":
                            # find the starting index of the ABAB sequence
                            dpos_idx = np.where((dpos == "A") & (np.roll(dpos, -1) == "B") &
                                                (np.roll(dpos, -2) == "A") & (np.roll(dpos, -3) == "B"))[0]
                            for i in dpos_idx:
                                if i < len(dpos) - 3:
                                    bkgid[i] = combid[i+1]
                                    bkgid[i+1] = combid[i]
                                    combid[i+2] = combid[i]
                                    bkgid[i+2] = bkgid[i]
                                    combid[i+3] = combid[i+1]
                                    bkgid[i+3] = bkgid[i+1]

                        elif dpat == "ABBA":
                            # find the starting index of the ABBA sequence
                            dpos_idx = np.where((dpos == "A") & (np.roll(dpos, -1) == "B") &
                                                (np.roll(dpos, -2) == "B") & (np.roll(dpos, -3) == "A"))[0]
                            for i in dpos_idx:
                                if i < len(dpos) - 3:
                                    bkgid[i] = combid[i+1]
                                    bkgid[i+1] = combid[i]
                                    combid[i+2] = combid[i+1]
                                    bkgid[i+2] = bkgid[i+1]
                                    combid[i+3] = combid[i]
                                    bkgid[i+3] = bkgid[i]

                        # if dpat is "Stare" try to find a sequence using dpos
                        elif dpat == "Stare":
                            # find the starting index of a possible ABBA sequence
                            dpos_idx = np.where((dpos == "A") & (np.roll(dpos, -1) == "B") &
                                                (np.roll(dpos, -2) == "B") & (np.roll(dpos, -3) == "A"))[0]
                            if dpos_idx.size > 0:
                                for i in dpos_idx:
                                    if i < len(dpos) - 3:
                                        bkgid[i] = combid[i+1]
                                        bkgid[i+1] = combid[i]
                                        combid[i+2] = combid[i+1]
                                        bkgid[i+2] = bkgid[i+1]
                                        combid[i+3] = combid[i]
                                        bkgid[i+3] = bkgid[i]
                            # find the starting index of a possible ABA'B' sequence
                            dpos_idx = np.where((dpos == "A") & (np.roll(dpos, -1) == "B") &
                                                (np.roll(dpos, -2) == "A'") & (np.roll(dpos, -3) == "B'"))[0]
                            if dpos_idx.size > 0:
                                for i in dpos_idx:
                                    if i < len(dpos) - 3:
                                        bkgid[i] = combid[i+1]
                                        bkgid[i+1] = combid[i]
                                        bkgid[i+2] = bkgid[i+3]
                                        bkgid[i+3] = bkgid[i+2]
                            # find the starting index of a possible AB sequence
                            dpos_idx = np.where((dpos == "A") & (np.roll(dpos, -1) == "B"))[0]
                            if dpos_idx.size > 0:
                                for i in dpos_idx:
                                    # exclude when np.roll counts the 1st element of dpos to be in a
                                    # sequence with the last element
                                    if i < len(dpos)-1:
                                        bkgid[i] = combid[i+1]
                                        bkgid[i+1] = combid[i]

                        # assign bkgid for files that deviate from general a sequence
                        for i in range(len(fitstbl[dpat_idx])):
                            # if A frame doesn't have bkgid assigned
                            if bkgid[i] == -1 and \
                            (fitstbl[dpat_idx]['dithpos'][i] == "A" or fitstbl[dpat_idx]['dithpos'][i] == "A'"):
                                # find closest (in mjd) B frame to subtract from this A
                                if fitstbl[dpat_idx]['dithpos'][i] == "A":
                                    pos_idx = fitstbl[dpat_idx]['dithpos'] == "B"
                                elif fitstbl[dpat_idx]['dithpos'][i] == "A'":
                                    pos_idx = fitstbl[dpat_idx]['dithpos'] == "B'"
                                if np.any(pos_idx):
                                    close_idx = np.argmin(np.absolute(fitstbl[dpat_idx][pos_idx]['mjd'] - fitstbl[dpat_idx]['mjd'][i]))
                                    bkgid[i] = combid[pos_idx][close_idx]
                            # if B frame doesn't have bkgid assigned
                            if bkgid[i] == -1 and \
                            (fitstbl[dpat_idx]['dithpos'][i] == "B" or fitstbl[dpat_idx]['dithpos'][i] == "B'"):
                                # find closest (in mjd) A frame to subtract from this B
                                if fitstbl[dpat_idx]['dithpos'][i] == "B":
                                    pos_idx = np.where(fitstbl[dpat_idx]['dithpos'] == "A")[0]
                                elif fitstbl[dpat_idx]['dithpos'][i] == "B'":
                                    pos_idx = np.where(fitstbl[dpat_idx]['dithpos'] == "A'")[0]
                                if np.any(pos_idx):
                                    close_idx = np.argmin(np.absolute(fitstbl[dpat_idx][pos_idx]['mjd'] - fitstbl[dpat_idx]['mjd'][i]))
                                    bkgid[i] = combid[pos_idx][close_idx]
                        fitstbl['bkg_id'][dpat_idx] = bkgid
                        fitstbl['comb_id'][dpat_idx] = combid

        return fitstbl

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard PypeIt file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
#        pypeit_keys = super().pypeit_file_keys()
#        # TODO: Why are these added here? See
#        # pypeit.metadata.PypeItMetaData.set_pypeit_cols
#        pypeit_keys += [calib', 'comb_id', 'bkg_id']
#        return pypeit_keys
        pypeit_keys = super().pypeit_file_keys()
        pypeit_keys.remove('decker_secondary')
        pypeit_keys.remove('slitwid')
        pypeit_keys.remove('slitlength')
        return pypeit_keys + ['lampstat01', 'dithpat', 'dithpos', 'dithoff', 'frameno']

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
            return good_exp & (fitstbl['idname'] == 'object')
        if ftype in ['bias', 'dark']:
            return good_exp & (fitstbl['lampstat01'] == 'off') & (fitstbl['idname'] == 'dark')
        if ftype in ['lampoffflats']:
            return good_exp & (fitstbl['lampstat01'] == 'off') & (fitstbl['idname'] == 'flatlampoff')
        if ftype in ['illumflat', 'pixelflat', 'trace']:
            # Flats and trace frames are typed together
            return good_exp & (fitstbl['lampstat01'] == 'on') & (fitstbl['idname'] == 'flatlamp')
        if ftype == 'pinhole':
            # Don't type pinhole frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            # TODO: This is a kludge.  Allow science frames to also be
            # classified as arcs
            is_arc = fitstbl['idname'] == 'arclamp'
            is_obj = (fitstbl['lampstat01'] == 'off') & (fitstbl['idname'] == 'object') & ('long2pos_specphot' not in fitstbl['decker'])
            return good_exp & (is_arc | is_obj)
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    # TODO: Is this supposed to be deprecated in favor of get_comb_group?
    def parse_dither_pattern(self, file_list, ext=None):
        """
        Parse headers from a file list to determine the dither pattern.

        Parameters
        ----------
        file_list (list of strings):
            List of files for which dither pattern is desired
        ext (int, optional):
            Extension containing the relevant header for these files. Default=None. If None, code uses
            self.primary_hdrext

        Returns
        -------
        dither_pattern, dither_id, offset_arcsec

        dither_pattern (str `numpy.ndarray`_):
            Array of dither pattern names
        dither_id (str `numpy.ndarray`_):
            Array of dither pattern IDs
        offset_arc (float `numpy.ndarray`_):
            Array of dither pattern offsets
        """
        nfiles = len(file_list)
        offset_arcsec = np.zeros(nfiles)
        dither_pattern = []
        dither_id = []
        for ifile, file in enumerate(file_list):
            hdr = fits.getheader(file, self.primary_hdrext if ext is None else ext)
            dither_pattern.append(hdr['PATTERN'])
            dither_id.append(hdr['FRAMEID'])
            offset_arcsec[ifile] = hdr['YOFFSET']
        return np.array(dither_pattern), np.array(dither_id), np.array(offset_arcsec)

    def tweak_standard(self, wave_in, counts_in, counts_ivar_in, gpm_in, meta_table, debug=False):
        """

        This routine is for performing instrument/disperser specific tweaks to standard stars so that sensitivity
        function fits will be well behaved. For example, masking second order light. For instruments that don't
        require such tweaks it will just return the inputs, but for isntruments that do this function is overloaded
        with a method that performs the tweaks.

        Parameters
        ----------
        wave_in: `numpy.ndarray`_
            Input standard star wavelengths (:obj:`float`, ``shape = (nspec,)``)
        counts_in: `numpy.ndarray`_
            Input standard star counts (:obj:`float`, ``shape = (nspec,)``)
        counts_ivar_in: `numpy.ndarray`_
            Input inverse variance of standard star counts (:obj:`float`, ``shape = (nspec,)``)
        gpm_in: `numpy.ndarray`_
            Input good pixel mask for standard (:obj:`bool`, ``shape = (nspec,)``)
        meta_table: :obj:`dict`
            Table containing meta data that is slupred from the :class:`~pypeit.specobjs.SpecObjs`
            object.  See :meth:`~pypeit.specobjs.SpecObjs.unpack_object` for the
            contents of this table.

        Returns
        -------
        wave_out: `numpy.ndarray`_
            Output standard star wavelengths (:obj:`float`, ``shape = (nspec,)``)
        counts_out: `numpy.ndarray`_
            Output standard star counts (:obj:`float`, ``shape = (nspec,)``)
        counts_ivar_out: `numpy.ndarray`_
            Output inverse variance of standard star counts (:obj:`float`, ``shape = (nspec,)``)
        gpm_out: `numpy.ndarray`_
            Output good pixel mask for standard (:obj:`bool`, ``shape = (nspec,)``)
        """

        # Could check the wavelenghts here to do something more robust to header/meta data issues
        if 'Y-spectroscopy' in meta_table['DISPNAME']:
            #wave_out = np.copy(wave_in)
            #counts_out = np.copy(counts_in)
            #counts_ivar_out = np.copy(counts_ivar_in)
            #gpm_out = np.copy(gpm_in)
            # The blue edge and red edge of the detector are contaiminated by higher order light. These are masked
            # by hand.
            #second_order_region= (wave_in < 9520.0) | (wave_in > 11256.0)
            #gpm_out = gpm_in & np.logical_not(second_order_region)

            # Use a sigmoid function to apodize the spectrum smoothly in the regions where it is bad.
            #dlam = 10.0 # width that determines how shaprly  apodization occurs
            #sigmoid_blue_arg = (wave_in - wave_blue)/dlam
            #sigmoid_red_arg = (wave_red - wave_in)/dlam

            #sigmoid_apodize = special.expit(sigmoid_blue_arg) * special.expit(sigmoid_red_arg)
            #counts = counts_in*sigmoid_apodize
            # No we apodize only the flux. Since there are more counts in the in the unapodized spectrum, there is also
            # more variance in the original counts_ivar_in, and so in this way the S/N ratio is naturally reduced
            # in this region. There is not an obvious way to tweak the error vector here, and we don't to just mask
            # since then the polynomial fits go crazy at the boundaries. This is a reasonable compromise to not mask the
            # counts_ivar_in. The flux is bogus in the regions we are apodizing, so it does not really matter what we do,
            # it is better than operating on the original bogus flux.

            # Inflat the errors in the apodized region so that they don't inform the fits much
            #apo_pix = (sigmoid_apodize < 0.95)
            #mean_counts, med_counts, sigma_counts = sigma_clipped_stats(counts_in[np.logical_not(apo_pix)], sigma=3.0)
            #sigma_apo = med_counts/20.0 # roughly S/N ratio 20 in the apodized region
            #counts_ivar[apo_pix] = 1.0/sigma_apo**2
            #sigma_apo = np.sqrt(np.abs(counts[apo_pix]))
            #counts_ivar[apo_pix] = utils.inverse(sigma_apo**2)
            #counts_ivar[apo_pix] = utils.clip_ivar(counts[apo_pix], counts_ivar_in[apo_pix], 10.0, mask=gpm_in[apo_pix])
            wave_blue = 9520.0  # blue wavelength below which there is contamination
            wave_red = 11256.0  # red wavelength above which the spectrum is containated

        elif 'J2-spectroscopy' in meta_table['DISPNAME']:
            wave_blue = 11170.0  # blue wavelength below which there is contamination
            wave_red = 12600.0  # red wavelength above which the spectrum is containated

        else:
            # keep everything the same
            wave_blue = -np.inf
            wave_red = np.inf

        second_order_region= (wave_in < wave_blue) | (wave_in > wave_red)
        wave = wave_in.copy()
        counts = counts_in.copy()
        gpm = gpm_in.copy()
        counts_ivar = counts_ivar_in.copy()
        wave[second_order_region] = 0.0
        counts[second_order_region] = 0.0
        counts_ivar[second_order_region] = 0.0
        # By setting the wavelengths to zero, we guarantee that the sensitvity function will only be computed
        # over the valid wavelength region. While we could mask, this would still produce a wave_min and wave_max
        # for the zeropoint that includes the bad regions, and the polynomial fits will extrapolate crazily there
        gpm[second_order_region] = False
        #if debug:
        #    from matplotlib import pyplot as plt
        #    counts_sigma = np.sqrt(utils.inverse(counts_ivar_in))
        #    plt.plot(wave_in, counts, color='red', alpha=0.7, label='apodized flux')
        #    plt.plot(wave_in, counts_in, color='black', alpha=0.7, label='flux')
        #    plt.plot(wave_in, counts_sigma, color='blue', alpha=0.7, label='flux')
        #    plt.axvline(wave_blue, color='blue')
        #    plt.axvline(wave_red, color='red')
        #    plt.legend()
        #    plt.show()
        return wave, counts, counts_ivar, gpm

    def list_detectors(self, mosaic=False):
        """
        List the *names* of the detectors in this spectrograph.

        This is primarily used :func:`~pypeit.slittrace.average_maskdef_offset`
        to measure the mean offset between the measured and expected slit
        locations.

        Detectors separated along the dispersion direction should be ordered
        along the first axis of the returned array.  For example, Keck/DEIMOS
        returns:
        
        .. code-block:: python
        
            dets = np.array([['DET01', 'DET02', 'DET03', 'DET04'],
                             ['DET05', 'DET06', 'DET07', 'DET08']])

        such that all the bluest detectors are in ``dets[0]``, and the slits
        found in detectors 1 and 5 are just from the blue and red counterparts
        of the same slit.

        Args:
            mosaic (:obj:`bool`, optional):
                Is this a mosaic reduction?
                It is used to determine how to list the detector, i.e., 'DET' or 'MSC'.

        Returns:
            `numpy.ndarray`_: The list of detectors in a `numpy.ndarray`_.  If
            the array is 2D, there are detectors separated along the dispersion
            axis.
        """
        return np.array([detector_container.DetectorContainer.get_name(i+1) 
                            for i in range(self.ndet)])

    def get_slitmask(self, filename):
        """
        Parse the slitmask data from a MOSFIRE file into :attr:`slitmask`, a
        :class:`~pypeit.spectrographs.slitmask.SlitMask` object.

        This can be used for multi-object slitmask, but it it's not good
        for "LONGSLIT" nor "long2pos". Both "LONGSLIT" and "long2pos" have emtpy/incomplete
        binTable where the slitmask data are stored.


        Args:
            filename (:obj:`str`):
                Name of the file to read.

        Returns:
            :class:`~pypeit.spectrographs.slitmask.SlitMask`: The slitmask
            data read from the file. The returned object is the same as
            :attr:`slitmask`.
        """
        # Open the file
        hdu = io.fits_open(filename)

        # load slitmask info
        targs = hdu['Target_List'].data
        ssl = hdu['Science_Slit_List'].data
        msl = hdu['Mechanical_Slit_List'].data
        # some needed cleanup
        ssl = ssl[ssl['Slit_Number'] != ' ']
        msl = msl[msl['Slit_Number'] != ' ']

        # Book keeping: Count and check that the # of objects in the SSL matches that of the MSL
        # and if we recover the total number of CSUs
        numslits = np.zeros(len(ssl))
        for i in range(len(ssl)):
            slit = ssl[i]
            numslits[i] = np.where(slit['Target_Name'] == msl['Target_in_Slit'])[0].size

        if (numslits.sum() != self._CSUnumslits()) and ('LONGSLIT' not in self.get_meta_value(filename, 'decker')) \
                and ('long2pos' not in self.get_meta_value(filename, 'decker')):
            msgs.error('The number of allocated CSU slits does not match the number of possible slits. '
                       'Slitmask design matching not possible. Turn parameter `use_maskdesign` off')

        targ_dist_center = np.array(ssl['Target_to_center_of_slit_distance'], dtype=float)

        slit_centers = np.array(ssl['Slit_length'], dtype=float) / 2.
        # Some # adjustment for long2pos
        if 'long2pos' in self.get_meta_value(filename, 'decker'):
            slit_centers_wrong = np.array(ssl['Slit_length'], dtype=float) / 2.
            # correct slit length
            ssl['Slit_length'] = self._long2pos_slits_length()
            slit_centers = np.array(ssl['Slit_length'], dtype=float) / 2.
            centers_diff = slit_centers - slit_centers_wrong
            targ_dist_center[0] -= centers_diff[0]
            targ_dist_center[2] += centers_diff[2]

        # Projected distance (in arcsec) of the object from the left and right (top and bot) edges of the slit
        topdist = np.round(slit_centers + targ_dist_center, 3)
        botdist = np.round(slit_centers - targ_dist_center, 3)

        # Find the index to map the objects in the Science Slit List and the Target list
        indx = index_of_x_eq_y(targs['Target_Name'], ssl['Target_Name'])
        targs_mtch = targs[indx]
        obj_ra = targs_mtch['RA_Hours']+' '+targs_mtch['RA_Minutes']+' '+targs_mtch['RA_Seconds']
        obj_dec = targs_mtch['Dec_Degrees']+' '+targs_mtch['Dec_Minutes']+' '+targs_mtch['Dec_Seconds']
        obj_ra, obj_dec = meta.convert_radec(obj_ra, obj_dec)
        objname = [item.strip() for item in ssl['Target_Name']]
        #   - Pull out the slit ID, object ID, name, object coordinates, top and bottom distance
        objects = np.array([np.array(ssl['Slit_Number'], dtype=int),
                           np.zeros(ssl['Slit_Number'].size, dtype=int),   # no object ID
                           obj_ra,
                           obj_dec,
                           objname,
                           np.array(targs_mtch['Magnitude'], dtype=float),
                           ['None']*ssl['Slit_Number'].size,       # no magnitude band
                           topdist,
                           botdist]).T

        # PA corresponding to positive x on detector (spatial)
        posx_pa = hdu[0].header['SKYPA2']
        if posx_pa < 0.:
            posx_pa += 360.

        slit_ra = ssl['Slit_RA_Hours']+' '+ssl['Slit_RA_Minutes']+' '+ssl['Slit_RA_Seconds']
        slit_dec = ssl['Slit_Dec_Degrees']+' '+ssl['Slit_Dec_Minutes']+' '+ssl['Slit_Dec_Seconds']
        slit_ra, slit_dec = meta.convert_radec(slit_ra, slit_dec)

        # Instantiate the slit mask object and return it
        self.slitmask = SlitMask(np.array([np.zeros(ssl['Slit_Number'].size),   # mosfire maskdef has not slit corners
                                           np.zeros(ssl['Slit_Number'].size),
                                           np.zeros(ssl['Slit_Number'].size),
                                           np.zeros(ssl['Slit_Number'].size),
                                           np.zeros(ssl['Slit_Number'].size),
                                           np.zeros(ssl['Slit_Number'].size),
                                           np.zeros(ssl['Slit_Number'].size),
                                           np.zeros(ssl['Slit_Number'].size)]).T.reshape(-1,4,2),
                                 slitid=np.array(ssl['Slit_Number'], dtype=int),
                                 align=ssl['Target_Name'] == 'posB',
                                 science=ssl['Target_Name'] != 'posB',
                                 onsky=np.array([slit_ra,
                                                 slit_dec,
                                                 np.array(ssl['Slit_length'], dtype=float),
                                                 np.array(ssl['Slit_width'], dtype=float),
                                                 np.array([round(hdu[0].header['SKYPA3'],2)]*ssl['Slit_Number'].size)]).T,
                                 objects=objects,
                                 posx_pa=posx_pa)
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
            msgs.error('The name of a science file should be provided')

        if self.slitmask is None:
            msgs.error('Unable to read slitmask design info. Provide a file.')

        platescale = self.get_detector_par(det=1)['platescale']
        slit_gap = self._slit_gap(platescale)

        # build an array of values containing the bottom (right) edge of the slits
        # starting edge
        edge = self._starting_edge(filename)
        bot_edges = np.array([edge], dtype=int)
        for i in range(self.slitmask.nslits - 1):
            # target is the slit number
            edge -= (self.slitmask.onsky[:,2][i]/platescale + slit_gap)
            bot_edges = np.append(bot_edges, np.round(edge))
        if bot_edges[-1] < 0:
            bot_edges[-1] = 1

        # build an array of values containing the top (left) edge of the slits
        top_edges = bot_edges - np.round(self.slitmask.onsky[:,2]/platescale)
        if top_edges[-1] < 0:
            top_edges[-1] = 1

        # Sort slits from left to right
        sortindx = np.argsort(self.slitmask.slitid[::-1])

        # This print a QA table with info on the slits sorted from left to right.
        if not debug:
            num = 0
            msgs.info('Expected slits')
            msgs.info('*' * 18)
            msgs.info('{0:^6s} {1:^12s}'.format('N.', 'Slit_Number'))
            msgs.info('{0:^6s} {1:^12s}'.format('-' * 5, '-' * 13))
            for i in range(sortindx.shape[0]):
                msgs.info('{0:^6d} {1:^12d}'.format(num, self.slitmask.slitid[sortindx][i]))
                num += 1
            msgs.info('*' * 18)

        # If instead we run this method in debug mode, we print more info
        if debug:
            num = 0
            msgs.info('Expected slits')
            msgs.info('*' * 92)
            msgs.info('{0:^5s} {1:^10s} {2:^12s} {3:^12s} {4:^16s} {5:^16s}'.format('N.', 'Slit_Number',
                                                                                    'slitLen(arcsec)',
                                                                                    'slitWid(arcsec)',
                                                                                    'top_edges(pix)',
                                                                                    'bot_edges(pix)'))
            msgs.info('{0:^5s} {1:^10s} {2:^12s} {3:^12s} {4:^16s} {5:^14s}'.format('-' * 4, '-' * 13, '-' * 11,
                                                                                    '-' * 11, '-' * 18, '-' * 15))
            for i in range(sortindx.size):
                msgs.info('{0:^5d}{1:^14d} {2:^9.3f} {3:^12.3f}    {4:^16.2f} {5:^14.2f}'.format(num,
                            self.slitmask.slitid[sortindx][i], self.slitmask.onsky[:,2][sortindx][i],
                            self.slitmask.onsky[:,3][sortindx][i], top_edges[sortindx][i], bot_edges[sortindx][i]))
                num += 1
            msgs.info('*' * 92)

        return top_edges, bot_edges, sortindx, self.slitmask

    @staticmethod
    def _CSUnumslits():
        """
        Returns:
            :obj:int: Number of CSUs always used by MOSFIRE in the slitmask

        """
        return 46

    @staticmethod
    def _slit_gap(platescale):
        """
        Args:
            platescale (:obj:`float`): platescale for the current detector

        Returns:
            :obj:float: Gap between each slit. The unit is pixels if platescale is provided otherwise it is arcsec.

        """

        return round(0.97 / platescale) if platescale is not None else 0.97

    @staticmethod
    def _CSUlength(platescale):
        """
        Args:
            platescale (:obj:`float`): platescale for the current detector

        Returns:
            :obj:float: Nominal length of each CSU. The unit is pixels if platescale is provided otherwise it is arcsec.

        """

        return 7.01/platescale if platescale is not None else 7.01

    @staticmethod
    def _starting_edge(scifile):
        """
        Provides the slit edge from where to start when extracting the prediction from the slitmask design

        Args:
            scifile (:obj:`str`):
                The filename of the science frame.

        Returns:
            :obj:int: Pixel position of the starting edge
        """
        hdu = io.fits_open(scifile)
        decker = hdu[0].header['MASKNAME']

        return 2034 if 'long2pos' not in decker else 1188

    @staticmethod
    def _long2pos_slits_length():
        """

        Returns:
            :obj:`tuple`: Three float numbers indicating the length in arcsec of the three slits in the long2pos mask

        """
        return '22.970', '7.010', '22.970'

    @staticmethod
    def _long2pos_pos():
        """

        Returns:
            :obj:`tuple`: Two integer number indicating the x position of the
            beginning and the end of the three slits forming the long2pos mask

        """
        return 880, 1190

    @staticmethod
    def find_longslit_pos(scifile):
        """
        Given a MOSFIRE science raw file, find the position of the slit
        in the LONGSLIT slitmask

        Args:
            scifile: (:obj:`str`):
                Name of the science file to read.

        Returns:
            :obj:`tuple`: Two integer number indicating the x position of the
            beginning and the end of the slit.

        """
        # Read some values from header
        hdu = io.fits_open(scifile)
        decker = hdu[0].header['MASKNAME']
        platescale = hdu[0].header['PSCALE']

        slit_gap = KeckMOSFIRESpectrograph._slit_gap(platescale)  # pixels
        CSUlength = KeckMOSFIRESpectrograph._CSUlength(platescale)  # pixels

        # Number of CSU used to make this longslit is recorded in the MASKNAME
        CSUnum = int(decker.split("x")[0].split('-')[1])
        slit_length = CSUnum * CSUlength + (CSUnum-1)*slit_gap
        if CSUnum % 2 == 0:
            pix_start = hdu[0].header['CRPIX2'] - (slit_length/2. + (CSUlength+slit_gap)/2. + 1)
            pix_end = hdu[0].header['CRPIX2'] + (slit_length/2. - (CSUlength+slit_gap)/2. + 1)
        else:
            pix_start = hdu[0].header['CRPIX2'] - (slit_length/2. + 1)
            pix_end = hdu[0].header['CRPIX2'] + (slit_length/2. + 1)

        return int(round(pix_start)), int(round(pix_end))
