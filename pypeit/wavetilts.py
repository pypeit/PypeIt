"""
Module for guiding Arc/Sky line tracing

.. _numpy.ndarray: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html

"""
import os
import copy
import inspect

import numpy as np

from pypeit import msgs
from pypeit import masterframe
from pypeit import ginga
from pypeit.core import arc
from pypeit.core import tracewave, pixels
from pypeit.core import save
from pypeit.core import load

from IPython import embed


class WaveTilts(masterframe.MasterFrame):
    """
    Class to guide slit/order tracing

    Args:
        msarc (ndarray): Arc image
        tslits_dict (dict or None): dict from TraceSlits class (e.g. slitpix)
        spectrograph (:obj:`pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph object
        par (:class:`pypeit.par.pypeitpar.WaveTiltsPar` or None):
            The parameters used to fuss with the tilts
        wavepar (:class:`pypeit.par.pypeitpar.WaveSolutionPar` or None):
            The parameters used for the wavelength solution
        det (int): Detector index
        master_key (:obj:`str`, optional):
            The string identifier for the instrument configuration.  See
            :class:`pypeit.masterframe.MasterFrame`.
        master_dir (:obj:`str`, optional):
            Path to master frames.
        reuse_masters (:obj:`bool`, optional):
            Load master files from disk, if possible.
        qa_path (:obj:`str`, optional):
            Directory for QA output.
        msbpm (`numpy.ndarray`_, optional):
            Bad pixel mask.  If not provided, a dummy array with no
            masking is generated.


    Attributes:
        tilts_dict (dict):
            Holds the tilts data
        steps : list
        mask : ndarray, bool
          True = Ignore this slit
        all_trcdict : list of dict
          All trace dict's
        tilts : ndarray
          Tilts for a single slit/order
        all_ttilts : list of tuples
          Tuple of tilts ndarray's
        final_tilts : ndarray
          Final tilts image
        gpm (np.ndarray):
            Good pixel mask
            Eventually, we might attach this to self.msarc although that would then
            require that we write it to disk with self.msarc.image
    """
    # Frametype is a class attribute
    master_type = 'Tilts'

    @classmethod
    def from_master_file(cls, master_file):
        """

        Args:
            master_file (str):

        Returns:
            waveimage.WaveImage:

        """
        # Spectrograph
        spectrograph, extras = masterframe.items_from_master_file(master_file)
        head0 = extras[0]
        # Master info
        master_dir = head0['MSTRDIR']
        master_key = head0['MSTRKEY']
        # Instantiate
        slf = cls(None, None, spectrograph, None, None, master_dir=master_dir, master_key=master_key,
                  reuse_masters=True)
        # Load
        slf.tilts_dict = slf.load(ifile=master_file)
        # Return
        return slf

    # TODO This needs to be modified to take an inmask
    def __init__(self, msarc, tslits_dict, spectrograph, par, wavepar, det=1, master_key=None,
                 master_dir=None, reuse_masters=False, qa_path=None, msbpm=None):

        self.spectrograph = spectrograph
        self.par = par
        self.wavepar = wavepar

        # MasterFrame
        masterframe.MasterFrame.__init__(self, self.master_type, master_dir=master_dir,
                                         master_key=master_key, reuse_masters=reuse_masters)

        self.msarc = msarc
        self.tslits_dict = tslits_dict
        self.msbpm = msbpm
        self.det = det
        self.qa_path = qa_path

        # --------------------------------------------------------------
        # TODO: Build another base class that does these things for both
        # WaveTilts and WaveCalib?

        # Get the non-linear count level
        self.nonlinear_counts = 1e10 if self.spectrograph is None \
                                    else self.spectrograph.nonlinear_counts(det=self.det)

        # Set the slitmask and slit boundary related attributes that the
        # code needs for execution. This also deals with arcimages that
        # have a different binning then the trace images used to defined
        # the slits
        if self.tslits_dict is not None and self.msarc is not None:
            self.slitmask_science = pixels.tslits2mask(self.tslits_dict)
            gpm = (self.msbpm == 0) if self.msbpm is not None \
                                        else np.ones_like(self.slitmask_science, dtype=bool)
            self.shape_science = self.slitmask_science.shape
            self.shape_arc = self.msarc.image.shape
            self.nslits = self.tslits_dict['slit_left'].shape[1]
            self.slit_left = arc.resize_slits2arc(self.shape_arc, self.shape_science, self.tslits_dict['slit_left'])
            self.slit_righ = arc.resize_slits2arc(self.shape_arc, self.shape_science, self.tslits_dict['slit_righ'])
            self.slitcen   = arc.resize_slits2arc(self.shape_arc, self.shape_science, self.tslits_dict['slitcen'])
            self.slitmask  = arc.resize_mask2arc(self.shape_arc, self.slitmask_science)
            self.gpm = (arc.resize_mask2arc(self.shape_arc, gpm)) & (self.msarc.image < self.nonlinear_counts)
        else:
            self.slitmask_science = None
            self.shape_science = None
            self.shape_arc = None
            self.nslits = 0
            self.slit_left = None
            self.slit_righ = None
            self.slitcen = None
            self.slitmask = None
            self.gpm = None
        # --------------------------------------------------------------

        # Key Internals
        self.mask = None
        self.all_trace_dict = [None]*self.nslits
        self.tilts = None
        # 2D fits are stored as a dictionary rather than list because we will jsonify the dict
        self.all_fit_dict = [None]*self.nslits
        self.steps = []
        # Main outputs
        self.final_tilts = None
        self.fit_dict = None
        self.trace_dict = None

    def extract_arcs(self):
        """
        Extract the arcs down each slit/order

        Wrapper to arc.get_censpec()

        Args:

        Returns:
            np.ndarray, np.ndarray:  Extracted arcs

        """
        arccen, arc_maskslit = arc.get_censpec(self.slitcen, self.slitmask, self.msarc.image, gpm=self.gpm)
            #, nonlinear_counts=self.nonlinear_counts)
        # Step
        self.steps.append(inspect.stack()[0][3])
        return arccen, arc_maskslit

    def find_lines(self, arcspec, slit_cen, slit, debug=False):
        """
        Find the lines for tracing

        Wrapper to tracewave.tilts_find_lines()

        Args:
            arcspec:
            slit_cen:
            slit (int):
            debug:

        Returns:
            ndarray, ndarray:  Spectral, spatial positions of lines to trace

        """
        if self.par['idsonly'] is not None:
            # Put in some hook here for getting the lines out of the
            # wave calib for i.e. LRIS ghosts.
            only_these_lines = None
            pass
        else:
            only_these_lines = None

        tracethresh = self._parse_param(self.par, 'tracethresh', slit)
        lines_spec, lines_spat \
                = tracewave.tilts_find_lines(arcspec, slit_cen, tracethresh=tracethresh,
                                             sig_neigh=self.par['sig_neigh'],
                                             nfwhm_neigh=self.par['nfwhm_neigh'],
                                             only_these_lines=only_these_lines,
                                             fwhm=self.wavepar['fwhm'],
                                             nonlinear_counts=self.nonlinear_counts,
                                             debug_peaks=False, debug_lines=debug)

        self.steps.append(inspect.stack()[0][3])
        return lines_spec, lines_spat



    def fit_tilts(self, trc_tilt_dict, thismask, slit_cen, spat_order, spec_order, slit,
                  show_QA=False, doqa=True, debug=False):
        """
        Fit the tilts

        Args:
            trc_tilt_dict (dict): Contains information from tilt tracing
            slit_cen (ndarray): (nspec,) Central trace for this slit
            spat_order (int): Order of the 2d polynomial fit for the spatial direction
            spec_order (int): Order of the 2d polytnomial fit for the spectral direction
            slit (int): integer index for the slit in question

        Optional Args:
            show_QA: bool, default = False
                show the QA instead of writing it out to the outfile
            doqa: bool, default = True
                Construct the QA plot
            debug: bool, default = False
                Show additional plots useful for debugging.

        Returns:
           (tilts, coeffs)
            tilts: ndarray (nspec, nspat)
               tilts image
            coeff: ndarray (spat_order + 1, spec_order+1)
               Array containing the coefficients for the 2d legendre polynomial fit
        """
        # Now perform a fit to the tilts
        tilt_fit_dict, trc_tilt_dict_out \
                = tracewave.fit_tilts(trc_tilt_dict, thismask, slit_cen, spat_order=spat_order,
                                      spec_order=spec_order,maxdev=self.par['maxdev2d'],
                                      sigrej=self.par['sigrej2d'], func2d=self.par['func2d'],
                                      doqa=doqa, master_key=self.master_key, slit=slit,
                                      show_QA=show_QA, out_dir=self.qa_path, debug=debug)

        # Evaluate the fit
        #tilts = tracewave.fit2tilts((tilt_fit_dict['nspec'], tilt_fit_dict['nspat']), slit_cen,
        #                            tilt_fit_dict['coeff2'], tilt_fit_dict['func'])

        # Populate the fit dict, and update the all_trace_dict
        self.all_fit_dict[slit] = copy.deepcopy(tilt_fit_dict)
        self.all_trace_dict[slit] = copy.deepcopy(trc_tilt_dict_out)

        self.steps.append(inspect.stack()[0][3])
        return tilt_fit_dict['coeff2']

    def trace_tilts(self, arcimg, lines_spec, lines_spat, thismask, slit_cen):
        """
        Trace the tilts

        Args:

            arcimg (`numpy.ndarray`_):
                Arc image.  Shape is (nspec, nspat).
            lines_spec (`numpy.ndarray`_):
                Array containing the spectral pixel location of each
                line found for this slit.  Shape is (nlines,).
            lines_spat (`numpy.ndarray`_):
               Array containing the spatial pixel location of each line,
               which is the slitcen evaluate at the spectral position
               position of the line stored in lines_spec. Shape is
               (nlines,).
            thismask (`numpy.ndarray`_):
               Image indicating which pixels lie on the slit in
               equation. True = on the slit. False = not on slit.  Shape
               is (nspec, nspat) with dtype=bool.
            slit_cen (:obj:`int`):
                Integer index indicating the slit in question.

        Returns:
            dict: Dictionary containing information on the traced tilts required to fit the filts.

        """
        trace_dict = tracewave.trace_tilts(arcimg, lines_spec, lines_spat, thismask, slit_cen,
                                           inmask=self.gpm, fwhm=self.wavepar['fwhm'],
                                           spat_order=self.par['spat_order'],
                                           maxdev_tracefit=self.par['maxdev_tracefit'],
                                           sigrej_trace=self.par['sigrej_trace'])

        # Return
        self.steps.append(inspect.stack()[0][3])
        return trace_dict


    def run(self, maskslits=None, doqa=True, debug=False, show=False):
        """
        Main driver for tracing arc lines

        Code flow::
            1.  Extract an arc spectrum down the center of each slit/order
            2.  Loop on slits/orders
                i. Trace and fit the arc lines (This is done twice, once
                   with trace_crude as the tracing crutch, then again
                   with a PCA model fit as the crutch).
                ii. Repeat trace.
                iii.  2D Fit to the offset from slitcen
                iv. Save

        Keyword Args:
            maskslits (`numpy.ndarray`_, optional):
                Boolean array to ignore slits.
            doqa (bool):
            debug (bool):
            show (bool):

        Returns:
            dict, ndarray:  Tilts dict and maskslits array
        """

        if maskslits is None:
            maskslits = np.zeros(self.nslits, dtype=bool)

        # Extract the arc spectra for all slits
        self.arccen, self.arc_maskslit = self.extract_arcs()#self.slitcen, self.slitmask, self.inmask)

        # maskslit
        self.mask = maskslits & (self.arc_maskslit==1)
        gdslits = np.where(np.invert(self.mask))[0]

        # Final tilts image
        self.final_tilts = np.zeros(self.shape_science,dtype=float)
        max_spat_dim = (np.asarray(self.par['spat_order']) + 1).max()
        max_spec_dim = (np.asarray(self.par['spec_order']) + 1).max()
        self.coeffs = np.zeros((max_spec_dim, max_spat_dim,self.nslits))
        self.spat_order = np.zeros(self.nslits, dtype=int)
        self.spec_order = np.zeros(self.nslits, dtype=int)

        # TODO sort out show methods for debugging
        #if show:
        #    viewer,ch = ginga.show_image(self.msarc*(self.slitmask > -1),chname='tilts')

        # Loop on all slits
        for slit in gdslits:
            msgs.info('Computing tilts for slit {:d}/{:d}'.format(slit,self.nslits-1))
            # Identify lines for tracing tilts
            msgs.info('Finding lines for tilt analysis')
            self.lines_spec, self.lines_spat = self.find_lines(self.arccen[:,slit], self.slitcen[:,slit], slit, debug=debug)
            if self.lines_spec is None:
                self.mask[slit] = True
                maskslits[slit] = True
                continue

            thismask = self.slitmask == slit
            # Trace
            msgs.info('Trace the tilts')
            self.trace_dict = self.trace_tilts(self.msarc.image, self.lines_spec, self.lines_spat,
                                               thismask, self.slitcen[:,slit])
            #if show:
            #    ginga.show_tilts(viewer, ch, self.trace_dict)

            self.spat_order[slit] = self._parse_param(self.par, 'spat_order', slit)
            self.spec_order[slit] = self._parse_param(self.par, 'spec_order', slit)
            # 2D model of the tilts, includes construction of QA
            coeff_out = self.fit_tilts(self.trace_dict, thismask, self.slitcen[:,slit],
                                       self.spat_order[slit], self.spec_order[slit], slit,
                                       doqa=doqa, show_QA=show, debug=show)
            self.coeffs[0:self.spec_order[slit]+1, 0:self.spat_order[slit]+1 , slit] = coeff_out

            # Tilts are created with the size of the original slitmask,
            # which corresonds to the same binning as the science
            # images, trace images, and pixelflats etc.
            self.tilts = tracewave.fit2tilts(self.slitmask_science.shape, coeff_out,
                                             self.par['func2d'])
            # Save to final image
            thismask_science = self.slitmask_science == slit
            self.final_tilts[thismask_science] = self.tilts[thismask_science]

        self.tilts_dict = {'tilts':self.final_tilts, 'coeffs':self.coeffs, 'slitcen':self.slitcen,
                           'func2d':self.par['func2d'], 'nslit':self.nslits,
                           'spat_order':self.spat_order, 'spec_order':self.spec_order}
        return self.tilts_dict, maskslits

    def save(self, outfile=None, overwrite=True):
        """
        Save the wavelength tilts data to a master frame

        Args:
            outfile (:obj:`str`, optional):
                Name for the output file.  Defaults to
                :attr:`master_file_path`.
            overwrite (:obj:`bool`, optional):
                Overwrite any existing file.
        """
        _outfile = self.master_file_path if outfile is None else outfile
        # Check if it exists
        if os.path.exists(_outfile) and not overwrite:
            msgs.warn('Master file exists: {0}'.format(_outfile) + msgs.newline()
                      + 'Set overwrite=True to overwrite it.')
            return

        # Log
        msgs.info('Saving master frame to {0}'.format(_outfile))

        # Build the header
        hdr = self.build_master_header(steps=self.steps)
        #   - Set the master frame type
        hdr['FRAMETYP'] = (self.master_type, 'PypeIt: Master calibration frame type')
        #   - Tilts metadata
        hdr['FUNC2D'] = self.tilts_dict['func2d']
        hdr['NSLIT'] = self.tilts_dict['nslit']

        # Write the fits file
        data = [self.tilts_dict['tilts'], self.tilts_dict['coeffs'], self.tilts_dict['slitcen'],
                self.tilts_dict['spat_order'], self.tilts_dict['spec_order']]
        extnames = ['TILTS', 'COEFFS', 'SLITCEN', 'SPAT_ORDER', 'SPEC_ORDER']
        save.write_fits(hdr, data, _outfile, extnames=extnames)
        #fits.HDUList([fits.PrimaryHDU(header=hdr),
        #              fits.ImageHDU(data=self.tilts_dict['tilts'], name='TILTS'),
        #              fits.ImageHDU(data=self.tilts_dict['coeffs'], name='COEFFS'),
        #              fits.ImageHDU(data=self.tilts_dict['slitcen'], name='SLITCEN'),
        #              fits.ImageHDU(data=self.tilts_dict['spat_order'], name='SPAT_ORDER'),
        #              fits.ImageHDU(data=self.tilts_dict['spec_order'], name='SPEC_ORDER')
        #             ]).writeto(_outfile, overwrite=True)

    def load(self, ifile=None):
        """
        Load the tilts data.

        This is largely a wrapper for :func:`pypeit.wavetilts.WaveTilts.load_from_file`.

        Args:
            ifile (:obj:`str`, optional):
                Name of the master frame file.  Defaults to
                :attr:`master_file_path`.
            return_header (:obj:`bool`, optional):
                Return the header.

        Returns:
            dict: Returns the tilts dictionary.  If nothing is
            loaded, either because :attr:`reuse_masters` is `False` or
            the file does not exist, everything is returned as None (one
            per expected return object).
        """
        # Check on whether to reuse and whether the file exists
        master_file = self.chk_load_master(ifile)
        if master_file is None:
            return
        msgs.info('Loading Master frame: {0}'.format(master_file))
        # Load
        extnames = ['TILTS', 'COEFFS', 'SLITCEN', 'SPAT_ORDER', 'SPEC_ORDER']
        *data, head0 = load.load_multiext_fits(master_file, extnames)

        # Fill the dict
        self.tilts_dict = {}
        keys = ['func2d', 'nslit']
        for k in keys:
            self.tilts_dict[k] = head0[k.upper()]
        # Data
        for ii,ext in enumerate(extnames):
            self.tilts_dict[ext.lower()] = data[ii]
        # Return
        return self.tilts_dict

    @staticmethod
    def load_from_file(filename, return_header=False):
        """
        Load the tilts data, without the benefit of the rest of the
        class.

        Args:
            ifile (:obj:`str`, optional):
                Name of the master frame file.  Defaults to
                :attr:`master_file_path`.
            return_header (:obj:`bool`, optional):
                Return the header, which will include the TraceImage
                metadata if available.

        Returns:
            Returns the tilts dictionary.  If return_header is
            true, the primary header is also returned.  If nothing is
            loaded, either because :attr:`reuse_masters` is `False` or
            the file does not exist, everything is returned as None (one
            per expected return object).
        """
        embed(header='471')

        # Check it exists
        if not os.path.isfile(filename):
            msgs.error('File does not exist: {0}'.format(filename))

    def _parse_param(self, par, key, slit):
        """
        Grab a parameter for a given slit

        Args:
            par (ParSet):
            key (str):
            slit (int):

        Returns:
            object:  Value of the parameter

        """
        param_in = par[key]
        if isinstance(param_in, (float, int)):
            param = param_in
        elif isinstance(param_in, (list, np.ndarray)):
            param = param_in[slit]
        else:
            raise ValueError('Invalid input for parameter {:s}'.format(key))
        return param

    def show(self, attr, slit=None, display='ginga', cname=None):
        """
        Display an image or spectrum in TraceSlits

        Parameters
        ----------
        attr : str
          'fweight'  -- Show the msarc image and the tilts traced by fweight
          'model'    -- Show the msarc image and the poylynomial model fits to the individual arc lines that
                        were traced by fweight.
          'arcmodel -- This illustrates the global final 2-d model fit to the indivdiaul models of each traced fweight arc line
                       tilts evaluated at the location of the specific arclines that wered use for the fit.
          'final_tilts' -- Show the final 2-d tilt model for all the slits that were fit.
        slit : int, optional
                    -- The slit to plot. This needs to be an integer between 1 and nslit
        display : str (optional)
          'ginga' -- Display to an RC Ginga
        """

        viewer, ch = ginga.show_image(self.arcimg*(self.slitmask == slit), chname='Tilts')
        ginga.show_tilts(viewer, ch, self.trace_dict,
                         sedges=(self.tslits_dict['slit_left'][:,slit],
                         self.tslits_dict['slit_righ'][:,slit]), points=True, clear_canvas=True)

    def __repr__(self):
        # Generate sets string
        txt = '<{:s}: '.format(self.__class__.__name__)
        if len(self.steps) > 0:
            txt+= ' steps: ['
            for step in self.steps:
                txt += '{:s}, '.format(step)
            txt = txt[:-2]+']'  # Trim the trailing comma
        txt += '>'
        return txt

