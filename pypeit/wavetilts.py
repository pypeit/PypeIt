"""
Module for guiding Arc/Sky line tracing

.. _numpy.ndarray: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html

"""
import os
import copy
import inspect

import numpy as np

from astropy import stats

from pypeit import msgs
from pypeit import masterframe
from pypeit import ginga
from pypeit import utils
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
        Instantiate from a master_file

        Args:
            master_file (str):

        Returns:
            wavetilts.WaveTilts:
                With tilts_dict loaded up

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
        arccen, arccen_bpm, arc_maskslit = arc.get_censpec(self.slitcen, self.slitmask,
                                                           self.msarc.image, gpm=self.gpm)
            #, nonlinear_counts=self.nonlinear_counts)
        # Step
        self.steps.append(inspect.stack()[0][3])
        return arccen, arccen_bpm, arc_maskslit

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


    def model_arc_continuum(self):
        """
        Model the continuum of the arc image.
        """
        # TODO: Instead check that extract arcs has been run using the
        # "steps" attribute?
        if self.arccen is None:
            # Extract the arc spectra for all slits
            self.arccen, self.arccen_bpm, self.arc_maskslit = self.extract_arcs()

        fit_order = 3
        fit_function = 'legendre'
        upper_rej = 1.5
        lower_rej = 3.

        # Fit the continuum of the extracted arc spectra for each slit
        # TODO: use iter_continuum
        nspec, nslits = self.arccen.shape
        spec = np.arange(nspec, dtype=float)
        arc_continuum = np.zeros(self.arccen.shape, dtype=float)
        arc_fitmask = np.zeros(self.arccen.shape, dtype=bool)
        for i in range(nslits):
            arc_fitmask[:,i], coeff \
                    = utils.robust_polyfit_djs(spec, self.arccen[:,i], fit_order,
                                               function=fit_function, minx=spec[0], maxx=spec[-1],
                                               inmask=np.invert(self.arccen_bpm[:,i]),
                                               upper=upper_rej, lower=lower_rej, use_mad=True,
                                               sticky=True)
            arc_continuum[:,i] = utils.func_val(coeff, spec, fit_function, minx=spec[0],
                                                maxx=spec[-1])

        # For each slit, rescale the continuum to the spectrum at a
        # given spatial position along the slit/order. The
        # implementation below may be too simplistic.
        nspat = self.slitmask.shape[1]
        cont_image = np.zeros(self.msarc.image.shape, dtype=float)
        ii = np.tile(np.arange(nspec), (nspat,1)).T
        jj = np.tile(np.arange(nspat), (nspec,1))
        # TODO: Can probably do this without the for loop
        for i in range(nslits):
            # Find the pixels in this slit
            indx = self.slitmask == i

            # Get a single width for the slit to simplify the calculation
            width = np.sum(indx, axis=1)
            width = int(np.amax(width[np.invert(self.arccen_bpm[:,i])]))

            # Get the spatial indices for spectral pixels in the
            # spatial dimension that follow the curvature of the slit
            # center.  TODO: May need to be more sophisticated.
            _jj = (self.slitcen[:,i,None] + np.arange(width)[None,:] - width//2).astype(int)

            # Set the index of masked pixels or those off the detector
            # to -1 so that they don't cause the image indexing to
            # fault and can be selected for masking
            _jj[(_jj < 0) | (_jj >= nspat) | self.arccen_bpm[:,i,None]] = -1

            # Pull out the slit pixels into a square array and mask
            # pixels off of the slit
            aligned_spec = np.tile(np.arange(nspec), (width,1)).T
            aligned_flux = np.ma.MaskedArray(self.msarc.image[aligned_spec, _jj], mask=_jj==-1)

            # Use a sigma-clipped median to determine the scaling of
            # the continuum fit to the central extracted spectrum to
            # match the spectrum at each spatial position.
            # TODO: Instead of determining the scale factor directly,
            # use the slit illumination profile?
            cont_renorm = stats.sigma_clipped_stats(aligned_flux/arc_continuum[:,i,None],
                                                    sigma_lower=lower_rej, sigma_upper=upper_rej,
                                                    axis=0)[1]

            # Fill the image with the continuum for this slit
            indx = np.invert(aligned_flux.mask)
            cont_image[aligned_spec[indx], _jj[indx]] \
                    = (arc_continuum[:,i,None] * cont_renorm[None,:])[indx]

        # Remove continuum measurements that are off any slits (because
        # of the fixed width assumption)
        cont_image[self.slitmask == -1] = 0.
        return cont_image

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
        self.arccen, self.arccen_bpm, self.arc_maskslit = self.extract_arcs()

        # Subtract arc continuum
        sub_continuum = True
        cont_image = self.model_arc_continuum() if sub_continuum \
                        else np.zeros(self.shape_science, dtype=float)

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

        # TODO: Leave for now.  Used for debugging
#        debug = True
#        show = True
#        import pdb
#        from matplotlib import pyplot

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
            self.trace_dict = self.trace_tilts(self.msarc.image - cont_image,
                                               self.lines_spec, self.lines_spat,
                                               thismask, self.slitcen[:,slit])
            #if show:
            #    ginga.show_tilts(viewer, ch, self.trace_dict)

            self.spat_order[slit] = self._parse_param(self.par, 'spat_order', slit)
            self.spec_order[slit] = self._parse_param(self.par, 'spec_order', slit)
            # 2D model of the tilts, includes construction of QA
            # NOTE: This also fills in self.all_fit_dict and self.all_trace_dict
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

            # TODO: Leave this commented code block for now. It plots
            # the tilt traced positions and the fit against the
            # continuum-subtracted image. This was used for testing.
#            pyplot.imshow(self.msarc.image - cont_image, origin='lower', interpolation='nearest',
#                          aspect='auto', vmin=0, vmax=1e3)
#            spat = self.all_trace_dict[slit]['tilts_spat']
#            spec = self.all_trace_dict[slit]['tilts']
#            spec_fit = self.all_trace_dict[slit]['tilts_fit']
#            in_fit = self.all_trace_dict[slit]['tot_mask']
#            not_fit = np.invert(in_fit)
#            fit_rej = in_fit & np.invert(self.all_trace_dict[slit]['fit_mask'])
#            fit_keep = in_fit & self.all_trace_dict[slit]['fit_mask']
#            pyplot.scatter(spat[not_fit], spec[not_fit], color='C1', marker='.', s=50, lw=0)
#            pyplot.scatter(spat[fit_rej], spec[fit_rej], color='C3', marker='.', s=50, lw=0)
#            pyplot.scatter(spat[fit_keep], spec[fit_keep], color='k', marker='.', s=50, lw=0)
#            with_fit = np.invert(np.all(np.invert(fit_keep), axis=0))
#            pyplot.plot(spat[:,with_fit], spec_fit[:,with_fit], color='k')
#            pyplot.show()

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

