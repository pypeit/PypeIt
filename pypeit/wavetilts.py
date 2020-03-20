"""
Module for guiding Arc/Sky line tracing

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import os
import copy
import inspect

import numpy as np
from matplotlib import pyplot as plt

from astropy import stats, visualization

from pypeit import msgs
from pypeit import datamodel
from pypeit import ginga
from pypeit.core import arc
from pypeit.core import tracewave

from IPython import embed


class WaveTilts(datamodel.DataContainer):
    """
    Simple DataContainer for the output from BuildWaveTilts

    All of the items in the datamodel are required for instantiation,
      although they can be None (but shouldn't be)

    """
    version = '1.0.0'

    # I/O
    output_to_disk = None  # This writes all items that are not None
    hdu_prefix = None      # None required for this DataContainer

    # MasterFrame fun
    master_type = 'Tilts'
    file_format = 'fits'

    datamodel = {
        'tilts':  dict(otype=np.ndarray, atype=np.floating, desc='Image of the tilts (nspec, nspat)'),
        'coeffs': dict(otype=np.ndarray, atype=np.floating, desc='2D coefficents for the fit.' \
                       'One set per slit/order (3D array).'),
        'slitcen': dict(otype=np.ndarray, atype=np.floating,
                        desc='Location of the slit center.  2D array (nspec, nslit)'),
        'nslit': dict(otype=int, desc='Number of slits'),
        'spat_order': dict(otype=np.ndarray, atype=np.integer,
                           desc='Order for spatial fit'),
        'spec_order': dict(otype=np.ndarray, atype=np.integer,
                           desc='Order for spectral fit'),
        'func2d': dict(otype=str,desc='Function used for the 2D fit'),
    }
    def __init__(self, tilts, coeffs, slitcen, nslit, spat_order, spec_order, func2d):

        # Parse
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        d = dict([(k,values[k]) for k in args[1:]])
        # Setup the DataContainer
        datamodel.DataContainer.__init__(self, d=d)

    def fit2tiltimg(self, slitmask, flexure=0.):
        """
        Generate a tilt image from the fit parameters

        Mainly to allow for flexure

        Args:
            slitmask (`np.ndarray`_):
            flexure (float, optional):
                Spatial shift of the tilt image onto the desired frame (typically a science image)

        Returns:
            `np.ndarray`_:  New tilt image

        """
        # TODO -- worry about skipped/broken slits (i.e. maskslits)
        final_tilts = np.zeros_like(self.tilts)
        for slit in range(self.nslit):
            coeff_out = self.coeffs[:self.spec_order[slit]+1,:self.spat_order[slit]+1,slit]
            _tilts = tracewave.fit2tilts(final_tilts.shape, coeff_out, self.func2d, spat_shift=-1*flexure)
            #
            thismask_science = slitmask == slit
            final_tilts[thismask_science] = _tilts[thismask_science]
        # Return
        return final_tilts



class BuildWaveTilts(object):
    """
    Class to guide slit/order tracing

    Args:
        mstilt (:class:`pypeit.tilitimage.TiltImage`): Tilt image
        slits (:class:`pypeit.slittrace.SlitTraceSet`, None):
            Slit edges
        spectrograph (:obj:`pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph object
        par (:class:`pypeit.par.pypeitpar.WaveTiltsPar` or None):
            The parameters used to fuss with the tilts
        wavepar (:class:`pypeit.par.pypeitpar.WaveSolutionPar` or None):
            The parameters used for the wavelength solution
        det (int): Detector index
        qa_path (:obj:`str`, optional):
            Directory for QA output.
        msbpm (`numpy.ndarray`_, optional):
            Bad pixel mask.  If not provided, a dummy array with no
            masking is generated.
        master_key (:obj:`str`, optional):  For naming QA only


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
            Eventually, we might attach this to self.mstilt although that would then
            require that we write it to disk with self.mstilt.image
    """

#    @classmethod
#    def from_master_file(cls, master_file):
#        """
#        Instantiate from a master_file
#
#        Args:
#            master_file (str):
#
#        Returns:
#            wavetilts.WaveTilts:
#                With tilts_dict loaded up
#
#        """
#        # Spectrograph
#        spectrograph, extras = masterframe.items_from_master_file(master_file)
#        #head0 = extras[0]
#        # Master info
#        #master_dir = head0['MSTRDIR']
#        #master_key = head0['MSTRKEY']
#        # Instantiate
#        slf = cls(None, None, spectrograph, None, None)#, master_dir=master_dir, master_key=master_key, reuse_masters=True)
#        # Load
#        slf.tilts_dict = slf.load(ifile=master_file)
#        # Return
#        return slf

    # TODO This needs to be modified to take an inmask
    def __init__(self, mstilt, slits, spectrograph, par, wavepar, det=1, qa_path=None, msbpm=None,
                 master_key=None):

        # TODO: Perform type checking
        self.spectrograph = spectrograph
        self.par = par
        self.wavepar = wavepar

        self.mstilt = mstilt
        self.slits = slits
        self.msbpm = msbpm
        self.det = det
        self.qa_path = qa_path
        self.master_key = master_key

        # --------------------------------------------------------------
        # TODO: Build another base class that does these things for both
        # WaveTilts and WaveCalib?

        # Get the non-linear count level
        self.nonlinear_counts = 1e10 if self.spectrograph is None \
            else self.spectrograph.nonlinear_counts(self.mstilt.detector)

        # Set the slitmask and slit boundary related attributes that the
        # code needs for execution. This also deals with arcimages that
        # have a different binning then the trace images used to defined
        # the slits
        if self.slits is not None and self.mstilt is not None:
            # NOTE: This uses the interneral definition of `pad`
            # TODO -- This needs to deal with flexure
            self.slitmask_science = self.slits.slit_img()
            gpm = (self.msbpm == 0) if self.msbpm is not None \
                                        else np.ones_like(self.slitmask_science, dtype=bool)
            self.shape_science = self.slitmask_science.shape
            self.shape_arc = self.mstilt.image.shape
            self.nslits = self.slits.nslits
            # Grab left, right and resize as need be
            left, right = self.slits.select_edges()
            self.slit_left = arc.resize_slits2arc(self.shape_arc, self.shape_science, left)
            self.slit_righ = arc.resize_slits2arc(self.shape_arc, self.shape_science, right)
            self.slitcen = arc.resize_slits2arc(self.shape_arc, self.shape_science, (left+right)/2)
            self.slitmask  = arc.resize_mask2arc(self.shape_arc, self.slitmask_science)
            self.gpm = (arc.resize_mask2arc(self.shape_arc, gpm)) & (self.mstilt.image < self.nonlinear_counts)
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
                                                           self.mstilt.image, gpm=self.gpm)
            #, nonlinear_counts=self.nonlinear_counts)
        # Step
        self.steps.append(inspect.stack()[0][3])
        return arccen, arccen_bpm, arc_maskslit

    def find_lines(self, arcspec, slit_cen, slit, bpm=None, debug=False):
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
        # TODO: Implement this!
        only_these_lines = None
        if self.par['idsonly']:
            # Put in some hook here for getting the lines out of the
            # wave calib for i.e. LRIS ghosts.
            raise NotImplementedError('Select lines with IDs for tracing not yet implemented.')

        tracethresh = self._parse_param(self.par, 'tracethresh', slit)
        lines_spec, lines_spat, good \
                = tracewave.tilts_find_lines(arcspec, slit_cen, tracethresh=tracethresh,
                                             sig_neigh=self.par['sig_neigh'],
                                             nfwhm_neigh=self.par['nfwhm_neigh'],
                                             only_these_lines=only_these_lines,
                                             fwhm=self.wavepar['fwhm'],
                                             nonlinear_counts=self.nonlinear_counts,
                                             bpm=bpm, debug_peaks=False, debug_lines=debug)

        if debug:
            mean, median, stddev = stats.sigma_clipped_stats(self.mstilt.image, sigma=3.)
#            vmin, vmax = visualization.ZScaleInterval().get_limits(self.mstilt.image)
            vmin = median - 2*stddev
            vmax = median + 2*stddev
            plt.imshow(self.mstilt.image, origin='lower', interpolation='nearest', aspect='auto',
                       vmin=vmin, vmax=vmax)
            plt.scatter(lines_spat[good], lines_spec[good], marker='x', color='k', lw=2, s=50)
            plt.scatter(lines_spat[np.invert(good)], lines_spec[np.invert(good)], marker='x', color='C3', lw=2, s=50)
            plt.show()

        self.steps.append(inspect.stack()[0][3])
        return (None, None) if lines_spec is None else (lines_spec[good], lines_spat[good])



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
#        tilt_fit_dict, trc_tilt_dict_out \
        self.all_fit_dict[slit], self.all_trace_dict[slit] \
                = tracewave.fit_tilts(trc_tilt_dict, thismask, slit_cen, spat_order=spat_order,
                                      spec_order=spec_order,maxdev=self.par['maxdev2d'],
                                      sigrej=self.par['sigrej2d'], func2d=self.par['func2d'],
                                      doqa=doqa, master_key=self.master_key, slit=slit,
                                      minmax_extrap=self.par['minmax_extrap'],
                                      show_QA=show_QA, out_dir=self.qa_path, debug=debug)

        # Evaluate the fit
        #tilts = tracewave.fit2tilts((tilt_fit_dict['nspec'], tilt_fit_dict['nspat']), slit_cen,
        #                            tilt_fit_dict['coeff2'], tilt_fit_dict['func'])

        # Populate the fit dict, and update the all_trace_dict
#        self.all_fit_dict[slit] = copy.deepcopy(tilt_fit_dict)
#        self.all_trace_dict[slit] = copy.deepcopy(trc_tilt_dict_out)

        self.steps.append(inspect.stack()[0][3])
#        return tilt_fit_dict['coeff2']
        return self.all_fit_dict[slit]['coeff2']

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
                                           #, debug_pca=True, show_tracefits=True)

        # Return
        self.steps.append(inspect.stack()[0][3])
        return trace_dict

    def model_arc_continuum(self, debug=False):
        """
        Model the continuum of the arc image.

        The method uses the arc spectra extracted using
        :attr:`extract_arcs` and fits a characteristic low-order
        continuum for each slit/order using
        :func:`pypeit.util.robust_polyfit_djs` and the parameters
        `cont_function`, `cont_order`, and `cont_rej` from
        :attr:`par`. The characteristic continuum is then rescaled to
        match the continuum at each spatial position in the
        slit/order.
        
        .. note::
            The approach used here may be too simplistic (in the
            robustness of the continuum fit and then how the
            continuum is rescaled and projected for each slit/order).
            Tests should be performed to make sure that the approach
            is good enough to trace the centroid of the arc/sky lines
            without biasing the centroids due to a non-zero continuum
            level.

        Args:
            debug (:obj:`bool`, optional):
                Run the method in debug mode.

        Returns:
            numpy.ndarray: Returns a 2D image with the same shape as
            :attr:`mstilt` with the model continuum.
        """
        # TODO: Instead check that extract arcs has been run using the
        # "steps" attribute?
        if self.arccen is None:
            # Extract the arc spectra for all slits
            self.arccen, self.arccen_bpm, self.arc_maskslit = self.extract_arcs()

        # TODO: Should make this operation part of WaveTiltsPar ...
        # Parse the upper and lower sigma rejection thresholds; used
        # when rescaling continuum from center spectrum.
        lower_rej, upper_rej = self.par['cont_rej'] if hasattr(self.par['cont_rej'], '__len__') \
                                    else np.repeat(self.par['cont_rej'], 2)

        # Fit the continuum of the extracted arc spectra for each slit
        nspec, nslits = self.arccen.shape
        spec = np.arange(nspec, dtype=float)
        arc_continuum = np.zeros(self.arccen.shape, dtype=float)
        arc_fitmask = np.zeros(self.arccen.shape, dtype=bool)
        for i in range(nslits):
            # TODO: What to do with the following iter_continuum parameters?:
            #       sigthresh, sigrej, niter_cont, cont_samp, cont_frac_fwhm
            arc_continuum[:,i], arc_fitmask[:,i] \
                    = arc.iter_continuum(self.arccen[:,i], inmask=np.invert(self.arccen_bpm[:,i]),
                                         fwhm=self.wavepar['fwhm'])
            # TODO: Original version.  Please leave it for now.
#            arc_fitmask[:,i], coeff \
#                    = utils.robust_polyfit_djs(spec, self.arccen[:,i], self.par['cont_order'],
#                                               function=self.par['cont_function'],
#                                               minx=spec[0], maxx=spec[-1],
#                                               upper=upper_rej, lower=lower_rej, use_mad=True,
#                                               sticky=True)
#            arc_continuum[:,i] = utils.func_val(coeff, spec, self.par['cont_function'],
#                                                minx=spec[0], maxx=spec[-1])

            if debug:
                plt.plot(spec, self.arccen[:,i], color='k', label='Arc')
                plt.scatter(spec, np.ma.MaskedArray(self.arccen[:,i], mask=arc_fitmask[:,i]),
                            marker='x', s=10, color='C1', label='Ignored')
                plt.plot(arc_continuum[:,i], color='C3', label='Cont. Fit')
                plt.xlabel('Spectral pixel')
                plt.ylabel('Counts')
                plt.legend()
                plt.show()

        # For each slit, rescale the continuum to the spectrum at a
        # given spatial position along the slit/order. This
        # implementation may be too simplistic in how it treats the
        # spatial axis.
        nspat = self.slitmask.shape[1]
        cont_image = np.zeros(self.mstilt.image.shape, dtype=float)
        # TODO: Can probably do this without the for loop but this
        # still may be faster.
        for i in range(nslits):
            # Find the pixels in this slit
            indx = self.slitmask == i

            # Set a single width for the slit to simplify the
            # calculation
            width = np.sum(indx, axis=1)
            width = int(np.amax(width[np.invert(self.arccen_bpm[:,i])]))

            # Get the spatial indices for spectral pixels in the
            # spatial dimension that follow the curvature of the slit
            # center.
            # TODO: May need to be more sophisticated.
            _spat = (self.slitcen[:,i,None] + np.arange(width)[None,:] - width//2).astype(int)

            # Set the index of masked pixels or those off the detector
            # to -1 so that they don't cause the image indexing to
            # fault and can be selected for masking
            _spat[(_spat < 0) | (_spat >= nspat) | self.arccen_bpm[:,i,None]] = -1

            # Pull out the slit pixels into a square array and mask
            # pixels off of the slit
            aligned_spec = np.tile(np.arange(nspec), (width,1)).T
            aligned_flux = np.ma.MaskedArray(self.mstilt.image[aligned_spec, _spat],
                                             mask=_spat==-1)

            # Use a sigma-clipped median to determine the scaling of
            # the continuum fit to the central extracted spectrum to
            # match the spectrum at each spatial position.
            # TODO: Instead of determining the scale factor directly,
            # use the slit-illumination profile?
            cont_renorm = stats.sigma_clipped_stats(aligned_flux/arc_continuum[:,i,None],
                                                    sigma_lower=lower_rej, sigma_upper=upper_rej,
                                                    axis=0)[1]

            # Fill the image with the continuum for this slit
            indx = np.invert(aligned_flux.mask)
            cont_image[aligned_spec[indx], _spat[indx]] \
                    = (arc_continuum[:,i,None] * cont_renorm[None,:])[indx]

        # Remove continuum measurements that are off any slits (because
        # of the fixed width assumption)
        cont_image[self.slitmask == -1] = 0.
        return cont_image

    def run(self, maskslits=None, doqa=True, debug=False, show=False):
        """
        Main driver for tracing arc lines

        Code flow:
            1. Extract an arc spectrum down the center of each slit/order
            2. Loop on slits/orders
                i. Trace and fit the arc lines (This is done twice, once
                   with trace_crude as the tracing crutch, then again
                   with a PCA model fit as the crutch).
                ii. Repeat trace.
                iii.  2D Fit to the offset from slitcen
                iv. Save

        Args:
            maskslits (`numpy.ndarray`_, optional):
                Boolean array to ignore slits.
            doqa (bool):
            debug (bool):
            show (bool):

        Returns:
            tuple:
                :class:`WaveTilts`_
                `np.ndarray`_

        """

        if maskslits is None:
            maskslits = np.zeros(self.nslits, dtype=bool)

        # Extract the arc spectra for all slits
        self.arccen, self.arccen_bpm, self.arc_maskslit = self.extract_arcs()

        # TODO: Leave for now.  Used for debugging
#        self.par['rm_continuum'] = True
#        debug = True
#        show = True

        # Subtract arc continuum
        _mstilt = self.mstilt.image.copy()
        if self.par['rm_continuum']:
            continuum = self.model_arc_continuum(debug=debug)
            _mstilt -= continuum
            if debug:
                # TODO: Put this into a function
                vmin, vmax = visualization.ZScaleInterval().get_limits(_mstilt)
                w,h = plt.figaspect(1)
                fig = plt.figure(figsize=(3*w,h))
                ax = fig.add_axes([0.15/3, 0.1, 0.8/3, 0.8])
                ax.imshow(self.mstilt.image, origin='lower', interpolation='nearest',
                          aspect='auto', vmin=vmin, vmax=vmax)
                ax.set_title('MasterArc')
                ax = fig.add_axes([1.15/3, 0.1, 0.8/3, 0.8])
                ax.imshow(continuum, origin='lower', interpolation='nearest',
                          aspect='auto', vmin=vmin, vmax=vmax)
                ax.set_title('Continuum')
                ax = fig.add_axes([2.15/3, 0.1, 0.8/3, 0.8])
                ax.imshow(_mstilt, origin='lower', interpolation='nearest',
                          aspect='auto', vmin=vmin, vmax=vmax)
                ax.set_title('MasterArc - Continuum')
                plt.show()

        # maskslit
        self.mask = np.any([maskslits, self.arc_maskslit==1], axis=0)
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
        #    viewer,ch = ginga.show_image(self.mstilt*(self.slitmask > -1),chname='tilts')

        # Loop on all slits
        for slit in gdslits:
            msgs.info('Computing tilts for slit {0}/{1}'.format(slit, self.nslits-1))
            # Identify lines for tracing tilts
            msgs.info('Finding lines for tilt analysis')
            self.lines_spec, self.lines_spat \
                    = self.find_lines(self.arccen[:,slit], self.slitcen[:,slit], slit,
                                      bpm=self.arccen_bpm[:,slit], debug=False) #debug)
            if self.lines_spec is None:
                self.mask[slit] = True
                maskslits[slit] = True
                continue

            thismask = self.slitmask == slit

            # Performs the initial tracing of the line centroids as a
            # function of spatial position resulting in 1D traces for
            # each line.
            msgs.info('Trace the tilts')
            self.trace_dict = self.trace_tilts(_mstilt, self.lines_spec, self.lines_spat,
                                               thismask, self.slitcen[:,slit])

            # TODO: Show the traces before running the 2D fit

            #if show:
            #    ginga.show_tilts(viewer, ch, self.trace_dict)

            self.spat_order[slit] = self._parse_param(self.par, 'spat_order', slit)
            self.spec_order[slit] = self._parse_param(self.par, 'spec_order', slit)
            # 2D model of the tilts, includes construction of QA
            # NOTE: This also fills in self.all_fit_dict and self.all_trace_dict
            coeff_out = self.fit_tilts(self.trace_dict, thismask, self.slitcen[:,slit],
                                       self.spat_order[slit], self.spec_order[slit], slit,
                                       doqa=doqa, show_QA=show, debug=show)
            self.coeffs[:self.spec_order[slit]+1,:self.spat_order[slit]+1,slit] = coeff_out

            # Tilts are created with the size of the original slitmask,
            # which corresonds to the same binning as the science
            # images, trace images, and pixelflats etc.
            self.tilts = tracewave.fit2tilts(self.slitmask_science.shape, coeff_out,
                                             self.par['func2d'])
            # Save to final image
            thismask_science = self.slitmask_science == slit
            self.final_tilts[thismask_science] = self.tilts[thismask_science]

        if debug:
            # TODO: Add this to the show method?
            vmin, vmax = visualization.ZScaleInterval().get_limits(_mstilt)
            plt.imshow(_mstilt, origin='lower', interpolation='nearest', aspect='auto',
                       vmin=vmin, vmax=vmax)
            for slit in gdslits:
                spat = self.all_trace_dict[slit]['tilts_spat']
                spec = self.all_trace_dict[slit]['tilts']
                spec_fit = self.all_trace_dict[slit]['tilts_fit']
                in_fit = self.all_trace_dict[slit]['tot_mask']
                not_fit = np.invert(in_fit) & (spec > 0)
                fit_rej = in_fit & np.invert(self.all_trace_dict[slit]['fit_mask'])
                fit_keep = in_fit & self.all_trace_dict[slit]['fit_mask']
                plt.scatter(spat[not_fit], spec[not_fit], color='C1', marker='.', s=30, lw=0)
                plt.scatter(spat[fit_rej], spec[fit_rej], color='C3', marker='.', s=30, lw=0)
                plt.scatter(spat[fit_keep], spec[fit_keep], color='k', marker='.', s=30, lw=0)
                with_fit = np.invert(np.all(np.invert(fit_keep), axis=0))
                for t in range(in_fit.shape[1]):
                    if not with_fit[t]:
                        continue
                    l, r = np.nonzero(in_fit[:,t])[0][[0,-1]]
                    plt.plot(spat[l:r+1,t], spec_fit[l:r+1,t], color='k')
            plt.show()

        # Build and return DataContainer
        tilts_dict = {'tilts':self.final_tilts, 'coeffs':self.coeffs, 'slitcen':self.slitcen,
                      'func2d':self.par['func2d'], 'nslit':self.nslits,
                      'spat_order':self.spat_order, 'spec_order':self.spec_order}
        return WaveTilts(**tilts_dict), maskslits

#    def save(self, outfile, master_dir, overwrite=True):
#        """
#        Save the wavelength tilts data to a master frame
#
#        Args:
#            outfile (:obj:`str`):
#                Name for the output file.  Defaults to
#                :attr:`master_file_path`.
#            overwrite (:obj:`bool`, optional):
#                Overwrite any existing file.
#        """
#        #_outfile = self.master_file_path if outfile is None else outfile
#        _outfile = outfile #is None else outfile
#        # Check if it exists
#        if os.path.exists(_outfile) and not overwrite:
#            msgs.warn('Master file exists: {0}'.format(_outfile) + msgs.newline()
#                      + 'Set overwrite=True to overwrite it.')
#            return
#
#        # Log
#        msgs.info('Saving master frame to {0}'.format(_outfile))
#
#        hdr = masterframe.build_master_header(self, self.master_key,
#                                              master_dir, self.spectrograph.spectrograph,
#                                              steps=self.steps)
#
#        # Build the header
#        #hdr = self.build_master_header(steps=self.steps)
#        #   - Set the master frame type
#        hdr['FRAMETYP'] = (self.master_type, 'PypeIt: Master calibration frame type')
#        #   - Tilts metadata
#        hdr['FUNC2D'] = self.tilts_dict['func2d']
#        hdr['NSLIT'] = self.tilts_dict['nslit']
#
#        # Write the fits file
#        data = [self.tilts_dict['tilts'], self.tilts_dict['coeffs'], self.tilts_dict['slitcen'],
#                self.tilts_dict['spat_order'], self.tilts_dict['spec_order']]
#        extnames = ['TILTS', 'COEFFS', 'SLITCEN', 'SPAT_ORDER', 'SPEC_ORDER']
#        embed(header='643 of wavetilts')
#        save.write_fits(hdr, data, _outfile, extnames=extnames)

#    def load(self, ifile):
#        """
#        Load the tilts data.
#
#        This is largely a wrapper for :func:`pypeit.wavetilts.WaveTilts.load_from_file`.
#
#        Args:
#            ifile (:obj:`str`):
#                Name of the master frame file.  Defaults to
#                :attr:`master_file_path`.
#            return_header (:obj:`bool`, optional):
#                Return the header.
#
#        Returns:
#            dict: Returns the tilts dictionary.  If nothing is
#            loaded, either because :attr:`reuse_masters` is `False` or
#            the file does not exist, everything is returned as None (one
#            per expected return object).
#        """
#        # Check on whether to reuse and whether the file exists
#        #master_file = self.chk_load_master(ifile)
#        #if master_file is None:
#        #    return
#        msgs.info('Loading Master frame: {0}'.format(ifile))
#        # Load
#        extnames = ['TILTS', 'COEFFS', 'SLITCEN', 'SPAT_ORDER', 'SPEC_ORDER']
#        *data, head0 = load.load_multiext_fits(ifile, extnames)
#
#        # Fill the dict
#        self.tilts_dict = {}
#        keys = ['func2d', 'nslit']
#        for k in keys:
#            self.tilts_dict[k] = head0[k.upper()]
#        # Data
#        for ii,ext in enumerate(extnames):
#            self.tilts_dict[ext.lower()] = data[ii]
#        # Return
#        return self.tilts_dict

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
            Options are:
                - ``'fweight'``: Show the mstilt image and the tilts
                  traced by fweight
                - ``'model'``: Show the mstilt image and the poylynomial
                  model fits to the individual arc lines that were
                  traced by fweight.
                - ``'arcmodel'``: This illustrates the global final 2-d
                  model fit to the indivdiaul models of each traced
                  fweight arc line tilts evaluated at the location of
                  the specific arclines that wered use for the fit.
                - ``'final_tilts'``: Show the final 2-d tilt model for
                  all the slits that were fit.
        slit : int, optional
            The slit to plot. This needs to be an integer between 1 and nslit
        display : str (optional)
            Use 'ginga' to display to an RC Ginga.

        """

        viewer, ch = ginga.show_image(self.arcimg*(self.slitmask == slit), chname='Tilts')
        ginga.show_tilts(viewer, ch, self.trace_dict,
                         sedges=(self.slits.left[:,slit], self.slits.right[:,slit]), points=True,
                         clear_canvas=True)

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

