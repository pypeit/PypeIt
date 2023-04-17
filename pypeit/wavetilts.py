"""
Module for guiding Arc/Sky line tracing

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import os
import copy
import inspect

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

from astropy import stats, visualization
from astropy.stats import sigma_clipped_stats
from astropy import table

from pypeit import msgs, datamodel
from pypeit import utils
from pypeit import slittrace
from pypeit import wavecalib
from pypeit import masterframe
from pypeit.display import display
from pypeit.core import arc
from pypeit.core import tracewave
from pypeit.images import buildimage

from IPython import embed


class WaveTilts(datamodel.DataContainer):
    """
    Simple DataContainer for the output from BuildWaveTilts

    All of the items in the datamodel are required for instantiation, although
    they can be None (but shouldn't be)

    The datamodel attributes are:

    .. include:: ../include/class_datamodel_wavetilts.rst

    When written to an output-file HDU, all `numpy.ndarray`_ elements are
    bundled into an `astropy.io.fits.BinTableHDU`_, and the other elements are
    written as header keywords.  Any datamodel elements that are None are *not*
    included in the output.

    """
    version = '1.2.0'

    # MasterFrame fun
    master_type = 'Tilts'
    master_file_format = 'fits'

    datamodel = {'coeffs': dict(otype=np.ndarray, atype=np.floating,
                                descr='2D coefficents for the fit on the initial slits.  One '
                                      'set per slit/order (3D array).'),
                 'bpmtilts': dict(otype=np.ndarray, atype=np.integer,
                                  descr='Bad pixel mask for tilt solutions. Keys are taken from '
                                        'SlitTraceSetBitmask'),
                 'nslit': dict(otype=int,
                               descr='Total number of slits.  This can include masked slits'),
                 'spat_id': dict(otype=np.ndarray, atype=np.integer, descr='Slit spat_id '),
                 'spat_order': dict(otype=np.ndarray, atype=np.integer,
                                    descr='Order for spatial fit (nslit)'),
                 'spec_order': dict(otype=np.ndarray, atype=np.integer,
                                    descr='Order for spectral fit (nslit)'),
                 'func2d': dict(otype=str, descr='Function used for the 2D fit'),
                 'PYP_SPEC': dict(otype=str, descr='PypeIt spectrograph name'),
                 'slits_calib_key': dict(otype=str, descr='SlitTraceSet calibration key. '
                                                          'This helps to find the Slits calibration '
                                                          'file when running pypeit_chk_tilts()'),
                 'spat_flexure': dict(otype=float, descr='Flexure shift from the input TiltImage'),
                 'tilt_traces': dict(otype=table.Table, descr='Table with the positions of the '
                                                              'traced and fitted tilts for all the slits. '
                                                              'see :func:`make_tbl_tilt_traces` for more details. ')}

    def __init__(self, coeffs, nslit, spat_id, spat_order, spec_order, func2d, bpmtilts=None,
                 spat_flexure=None, PYP_SPEC=None, slits_calib_key=None, tilt_traces=None):

        # Parse
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        d = dict([(k,values[k]) for k in args[1:]])
        # Setup the DataContainer
        datamodel.DataContainer.__init__(self, d=d)

    def _init_internals(self):
        self.filename = None
        # Master stuff
        self.master_key = None
        self.master_dir = None

    def _bundle(self):
        """
        Bundle the data in preparation for writing to a fits file.

        See :func:`pypeit.datamodel.DataContainer._bundle`. Data is
        always written to a 'TILTS' extension.
        """
        _d = super(WaveTilts, self)._bundle(ext='TILTS')

        # similarly to SlitTraceSet
        # save the table
        tab_detached = _d[0]['TILTS']['tilt_traces']
        # remove `tab_detached` from the dict
        _d[0]['TILTS'].pop('tilt_traces')
        # add the table as a separate extension
        return [_d[0], {'tilt_traces': tab_detached}]

    def is_synced(self, slits):
        """
        Confirm the slits in WaveTilts are aligned to that in SlitTraceSet

        Barfs if not

        Args:
            slits (:class:`pypeit.slittrace.SlitTraceSet`):

        """
        if not np.array_equal(self.spat_id, slits.spat_id):
            msgs.error("Your tilt solutions are out of sync with your slits.  Remove Masters and start from scratch")

    def fit2tiltimg(self, slitmask, flexure=None):
        """
        Generate a tilt image from the fit parameters

        Mainly to allow for flexure

        Args:
            slitmask (`numpy.ndarray`_):
            flexure (float, optional):
                Spatial shift of the tilt image onto the desired frame
                (typically a science image)

        Returns:
            `numpy.ndarray`_:  New tilt image

        """
        _flexure = 0. if flexure is None else flexure

        final_tilts = np.zeros_like(slitmask).astype(float)
        gdslit_spat = np.unique(slitmask[slitmask >= 0]).astype(int)
        # Loop
        for slit_spat in gdslit_spat:
            slit_idx = self.spatid_to_zero(slit_spat)
            # Calculate
            coeff_out = self.coeffs[:self.spec_order[slit_idx]+1,:self.spat_order[slit_idx]+1,slit_idx]
            _tilts = tracewave.fit2tilts(final_tilts.shape, coeff_out, self.func2d, spat_shift=-1*_flexure)
            # Fill
            thismask_science = slitmask == slit_spat
            final_tilts[thismask_science] = _tilts[thismask_science]
        # Return
        return final_tilts

    def spatid_to_zero(self, spat_id):
        """
        Convert slit spat_id to zero-based
        Mainly for coeffs

        Args:
            spat_id (int):

        Returns:
            int:

        """
        mtch = self.spat_id == spat_id
        return np.where(mtch)[0][0]

    def show(self, tilt_img=None, slits=None, waveimg=None, wcs_match=True, in_ginga=True, show_traces=False):
        """
        """
        # Warning: This parses the filename, not the Header!
        master_key, master_dir = masterframe.grab_key_mdir(self.filename, from_filename=True)
        # get slits if not provided
        if slits is None:
            slit_masterframe_name = masterframe.construct_file_name(slittrace.SlitTraceSet, self.slits_calib_key,
                                                                    master_dir=master_dir)
            if os.path.exists(slit_masterframe_name):
                slits = slittrace.SlitTraceSet.from_file(slit_masterframe_name)
            else:
                msgs.warn('Could not load slits to show with tilts image.')
        # slits can still be None if Slits file does not exist
        if slits is not None:
            # get slitmask, used for generating waveimg and to mask tiltimg
            slitmask = slits.slit_img(initial=True, flexure=self.spat_flexure)
        # get waveimg if not provided and if showing in ginga
        if waveimg is None and in_ginga:
            wave_masterframe_name = masterframe.construct_file_name(wavecalib.WaveCalib, master_key,
                                                                    master_dir=master_dir)
            if os.path.exists(wave_masterframe_name) and slits is not None:
                wv_calib = wavecalib.WaveCalib.from_file(wave_masterframe_name)
                tilts = self.fit2tiltimg(slitmask, flexure=self.spat_flexure)
                waveimg = wv_calib.build_waveimg(tilts, slits, spat_flexure=self.spat_flexure)
            else:
                msgs.warn('Could not load Wave image to show with tilts image.')
        # get tiltimg if not provided
        tiltimg_masterframe_name = None
        if tilt_img is None:
            tiltimg_masterframe_name = masterframe.construct_file_name(buildimage.TiltImage, master_key,
                                                                       master_dir=master_dir)
        else:
            tiltimg_masterframe_name = tilt_img
        if tiltimg_masterframe_name is not None and os.path.exists(tiltimg_masterframe_name):
            tilt_img_dict = buildimage.TiltImage.from_file(tiltimg_masterframe_name)
            tilt_img = tilt_img_dict.image * (slitmask > -1) if slits is not None else tilt_img_dict.image
        else:
            msgs.err('Tiltimg not found.')

        zmax = stats.sigma_clip(tilt_img, sigma=10, return_bounds=True)[2]
        zmin = stats.sigma_clip(tilt_img, sigma=5, return_bounds=True)[1] * 2
        cut = (zmin, zmax)
        # show in ginga
        if in_ginga:
            # connect to ginga
            display.connect_to_ginga(raise_err=True, allow_new=True)
            # display image
            viewer, ch = display.show_image(tilt_img, chname='Tilts', cuts=cut,
                                            wcs_match=wcs_match, waveimg=waveimg, clear=True)
            # show tilts
            display.show_tilts(viewer, ch, self.tilt_traces, points=show_traces, nspec=tilt_img.shape[0])
            # shoe slit edges
            if slits is not None:
                left, right, mask = slits.select_edges()
                gpm = mask == 0
                display.show_slits(viewer, ch, left[:, gpm], right[:, gpm], slit_ids=slits.spat_id[gpm])

        # show in matplotlib
        else:
            show_tilts_mpl(tilt_img, self.tilt_traces, slits=slits, show_traces=show_traces, cut=cut)


class BuildWaveTilts:
    """
    Class to guide arc/sky tracing

    Args:
        mstilt (:class:`pypeit.images.buildimage.TiltImage`): Tilt image
        slits (:class:`pypeit.slittrace.SlitTraceSet`):
            Slit edges
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph object
        par (:class:`pypeit.par.pypeitpar.WaveTiltsPar` or None):
            The parameters used to fuss with the tilts
        wavepar (:class:`pypeit.par.pypeitpar.WaveSolutionPar` or None):
            The parameters used for the wavelength solution
        det (int): Detector index
        qa_path (:obj:`str`, optional):
            Directory for QA output.
        master_key (:obj:`str`, optional):  For naming QA only
        spat_flexure (float, optional):
            If input, the slitmask and slit edges are shifted prior
            to tilt analysis.


    Attributes:
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
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
        gpm (`numpy.ndarray`_):
            Good pixel mask
            Eventually, we might attach this to self.mstilt although that would then
            require that we write it to disk with self.mstilt.image
    """

    # TODO This needs to be modified to take an inmask
    def __init__(self, mstilt, slits, spectrograph, par, wavepar, det=1, qa_path=None,
                 master_key=None, spat_flexure=None):

        # TODO: Perform type checking
        self.spectrograph = spectrograph
        self.par = par
        self.wavepar = wavepar

        self.mstilt = mstilt
        self.slits = slits
        self.det = det
        self.qa_path = qa_path
        self.master_key = master_key
        self.spat_flexure = spat_flexure

        # Get the non-linear count level
        # TODO: This is currently hacked to deal with Mosaics
        try:
            self.nonlinear_counts = self.mstilt.detector.nonlinear_counts()
        except:
            self.nonlinear_counts = 1e10

#        # Get the non-linear count level
#        self.nonlinear_counts = 1e10 if self.spectrograph is None \
#            else self.mstilt.detector.nonlinear_counts()
##            else self.spectrograph.nonlinear_counts(self.mstilt.detector)

        # Set the slitmask and slit boundary related attributes that the
        # code needs for execution. This also deals with arcimages that
        # have a different binning then the trace images used to defined
        # the slits

        # TODO -- Tidy this up into one or two methods?
        # Load up all slits
        # TODO -- Discuss further with JFH
        all_left, all_right, mask = self.slits.select_edges(initial=True, flexure=self.spat_flexure)  # Grabs all, initial slits
        # self.tilt_bpm = np.invert(mask == 0)
        # At this point of the reduction the only bitmask flags that may have been generated are 'USERIGNORE',
        # 'SHORTSLIT', 'BOXSLIT' and 'BADWVCALIB'. Here we use only 'USERIGNORE' and 'SHORTSLIT' to create the bpm mask
        self.tilt_bpm = self.slits.bitmask.flagged(mask, flag=['SHORTSLIT', 'USERIGNORE'])
        self.tilt_bpm_init = self.tilt_bpm.copy()
        # Slitmask
        # TODO -- Discuss further with JFH
        self.slitmask_science = self.slits.slit_img(initial=True, flexure=self.spat_flexure, exclude_flag=['BOXSLIT'])  # All unmasked slits
        # Resize
        # TODO: Should this be the bpm or *any* flag?
        gpm = self.mstilt.select_flag(flag='BPM', invert=True) if self.mstilt is not None \
                    else np.ones_like(self.slitmask_science, dtype=bool)
        self.shape_science = self.slitmask_science.shape
        self.shape_tilt = self.mstilt.image.shape
        self.slitcen = arc.resize_slits2arc(self.shape_tilt, self.shape_science, (all_left+all_right)/2)
        self.slitmask = arc.resize_mask2arc(self.shape_tilt, self.slitmask_science)
        self.gpm = (arc.resize_mask2arc(self.shape_tilt, gpm)) & (self.mstilt.image < self.nonlinear_counts)
    # --------------------------------------------------------------

        # Key Internals
        self.mask = None
        self.all_trace_dict = [None]*self.slits.nslits
        self.tilts = None
        # 2D fits are stored as a dictionary rather than list because we will jsonify the dict
        self.all_fit_dict = [None]*self.slits.nslits
        self.steps = []
        # Main outputs
        self.final_tilts = None
        self.fit_dict = None
        self.trace_dict = None

    def extract_arcs(self):
        """
        Extract the arcs down each slit/order

        Wrapper to arc.get_censpec()

        Returns:
            :obj:`tuple`: Extracted arcs in two `numpy.ndarray`_ objects
        """
        arccen, arccen_bpm, arc_maskslit = arc.get_censpec(self.slitcen, self.slitmask,
                                                           self.mstilt.image, gpm=self.gpm,
                                                           slit_bpm=self.tilt_bpm)
            #, nonlinear_counts=self.nonlinear_counts)
        # Step
        self.steps.append(inspect.stack()[0][3])

        # Update the mask
        self.tilt_bpm |= arc_maskslit

        return arccen, arccen_bpm

    def find_lines(self, arcspec, slit_cen, slit_idx, bpm=None, debug=False):
        """
        Find the lines for tracing

        Wrapper to tracewave.tilts_find_lines()

        Args:
            arcspec:
            slit_cen:
            slit_idx (int):
                Slit index, zero-based
            bpm (`numpy.ndarray`_, optional):
            debug (bool, optional):

        Returns:
            tuple:  2 objectcs
                - `numpy.ndarray`_ or None:  Spectral positions of lines to trace
                - `numpy.ndarray`_ or None:  Spatial positions of lines to trace

        """
        # TODO: Implement this!
        only_these_lines = None
        if self.par['idsonly']:
            # Put in some hook here for getting the lines out of the
            # wave calib for i.e. LRIS ghosts.
            raise NotImplementedError('Select lines with IDs for tracing not yet implemented.')

        # TODO -- This should be order not slit!
        tracethresh = self._parse_param(self.par, 'tracethresh', slit_idx)
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



    def fit_tilts(self, trc_tilt_dict, thismask, slit_cen, spat_order, spec_order, slit_idx,
                  show_QA=False, doqa=True):
        """
        Fit the tilts

        all_fit_dict and all_trace_dict are filled in place

        Args:
            trc_tilt_dict (dict): Contains information from tilt tracing
            slit_cen (ndarray): (nspec,) Central trace for this slit
            spat_order (int): Order of the 2d polynomial fit for the spatial direction
            spec_order (int): Order of the 2d polytnomial fit for the spectral direction
            slit_idx (int): zero-based, integer index for the slit in question

        Optional Args:
            show_QA: bool, default = False
                show the QA instead of writing it out to the outfile
            doqa: bool, default = True
                Construct the QA plot

        Returns:
            `numpy.ndarray`_: coeff: ndarray (spat_order + 1, spec_order+1)
               Array containing the coefficients for the 2d legendre polynomial fit
        """
        # Index
        self.all_fit_dict[slit_idx], self.all_trace_dict[slit_idx] \
                = tracewave.fit_tilts(trc_tilt_dict, thismask, slit_cen, spat_order=spat_order,
                                      spec_order=spec_order,maxdev=self.par['maxdev2d'],
                                      sigrej=self.par['sigrej2d'], func2d=self.par['func2d'],
                                      doqa=doqa, master_key=self.master_key,
                                      slitord_id=self.slits.slitord_id[slit_idx],
                                      minmax_extrap=self.par['minmax_extrap'],
                                      show_QA=show_QA, out_dir=self.qa_path)

        self.steps.append(inspect.stack()[0][3])
        return self.all_fit_dict[slit_idx]['coeff2']

    def trace_tilts(self, arcimg, lines_spec, lines_spat, thismask, slit_cen,
                    debug_pca=False, show_tracefits=False):
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
                                           sigrej_trace=self.par['sigrej_trace'],
                                           debug_pca=debug_pca, show_tracefits=show_tracefits)

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
            if self.tilt_bpm[i]:
                continue
            # TODO: What to do with the following iter_continuum parameters?:
            #       sigthresh, sigrej, niter_cont, cont_samp, cont_frac_fwhm
            arc_continuum[:,i], arc_fitmask[:,i] \
                    = arc.iter_continuum(self.arccen[:,i], gpm=np.invert(self.arccen_bpm[:,i]),
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
            # Masked?
            if self.tilt_bpm[i]:
                continue
            # Find the pixels in this slit
            indx = self.slitmask == self.slits.spat_id[i]

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

    def run(self, doqa=True, debug=False, show=False):
        """
        Main driver for tracing arc lines

        Code flow:

            #. Extract an arc spectrum down the center of each slit/order
            #. Loop on slits/orders
                #. Trace and fit the arc lines (This is done twice, once
                   with trace_crude as the tracing crutch, then again
                   with a PCA model fit as the crutch).
                #. Repeat trace.
                #.  2D Fit to the offset from slitcen
                #. Save

        Args:
            doqa (bool):
            debug (bool):
            show (bool):

        Returns:
            :class:`WaveTilts`:

        """
        # Extract the arc spectra for all slits
        self.arccen, self.arccen_bpm = self.extract_arcs()

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

        # Final tilts image
        self.final_tilts = np.zeros(self.shape_science,dtype=float)
        max_spat_dim = (np.asarray(self.par['spat_order']) + 1).max()
        max_spec_dim = (np.asarray(self.par['spec_order']) + 1).max()
        self.coeffs = np.zeros((max_spec_dim, max_spat_dim,self.slits.nslits))
        self.spat_order = np.zeros(self.slits.nslits, dtype=int)
        self.spec_order = np.zeros(self.slits.nslits, dtype=int)

        # Loop on all slits
        for slit_idx, slit_spat in enumerate(self.slits.spat_id):
            if self.tilt_bpm[slit_idx]:
                msgs.info('Skipping bad slit {0}/{1}'.format(slit_idx, self.slits.nslits))
                self.slits.mask[slit_idx] = self.slits.bitmask.turn_on(self.slits.mask[slit_idx], 'BADTILTCALIB')
                continue
            #msgs.info('Computing tilts for slit {0}/{1}'.format(slit, self.slits.nslits-1))
            msgs.info('Computing tilts for slit {0}/{1}'.format(slit_idx, self.slits.nslits))
            # Identify lines for tracing tilts
            msgs.info('Finding lines for tilt analysis')
            self.lines_spec, self.lines_spat \
                    = self.find_lines(self.arccen[:,slit_idx], self.slitcen[:,slit_idx],
                                      slit_idx,
                                      bpm=self.arccen_bpm[:,slit_idx], debug=debug)

            if self.lines_spec is None:
                msgs.warn('Did not recover any lines for slit/order = {:d}'.format(self.slits.slitord_id[slit_idx]) +
                          '. This slit/order will not reduced!')
                self.slits.mask[slit_idx] = self.slits.bitmask.turn_on(self.slits.mask[slit_idx], 'BADTILTCALIB')
                continue

            thismask = self.slitmask == slit_spat

            # Performs the initial tracing of the line centroids as a
            # function of spatial position resulting in 1D traces for
            # each line.
            msgs.info('Trace the tilts')
            self.trace_dict = self.trace_tilts(_mstilt, self.lines_spec, self.lines_spat,
                                               thismask, self.slitcen[:, slit_idx])
            # IF there are < 2 usable arc lines for tilt tracing, PCA fit does not work and the reduction crushes
            # TODO investigate why some slits have <2 usable arc lines
            if np.sum(self.trace_dict['use_tilt']) < 2:
                msgs.warn('Less than 2 usable arc lines for slit/order = {:d}'.format(self.slits.slitord_id[slit_idx]) +
                          '. This slit/order will not reduced!')
                self.slits.mask[slit_idx] = self.slits.bitmask.turn_on(self.slits.mask[slit_idx], 'BADTILTCALIB')
                continue
            # Check the spectral coverage of the usable arc lines for tilts. If the coverage is small,
            # it will affect the wavelength range in waveimg (i.e, self.wv_calib.build_waveimg()) and
            # crash the reduction later on.
            # Here we mask slits that computed the tilts with arc lines coverage <10%
            use_tilt_spec_cov = (self.trace_dict['tilts_spec'][:, self.trace_dict['use_tilt']].max() -
                                 self.trace_dict['tilts_spec'][:, self.trace_dict['use_tilt']].min()) / self.arccen.shape[0]
            if use_tilt_spec_cov < 0.1:
                msgs.warn(f'The spectral coverage of the usable arc lines is {use_tilt_spec_cov:.3f} (less than 10%).' +
                          ' This slit/order will not be reduced!')
                self.slits.mask[slit_idx] = self.slits.bitmask.turn_on(self.slits.mask[slit_idx], 'BADTILTCALIB')
                continue

            # TODO: Show the traces before running the 2D fit

            self.spat_order[slit_idx] = self._parse_param(self.par, 'spat_order', slit_idx)
            self.spec_order[slit_idx] = self._parse_param(self.par, 'spec_order', slit_idx)
            # 2D model of the tilts, includes construction of QA
            # NOTE: This also fills in self.all_fit_dict and self.all_trace_dict
            coeff_out = self.fit_tilts(self.trace_dict, thismask, self.slitcen[:,slit_idx],
                                       self.spat_order[slit_idx], self.spec_order[slit_idx],
                                       slit_idx,
                                       doqa=doqa, show_QA=show)
            self.coeffs[:self.spec_order[slit_idx]+1,:self.spat_order[slit_idx]+1,slit_idx] = coeff_out

            # TODO: Need a way to assess the success of fit_tilts and
            # flag the slit if it fails

            # Tilts are created with the size of the original slitmask,
            # which corresonds to the same binning as the science
            # images, trace images, and pixelflats etc.
            self.tilts = tracewave.fit2tilts(self.slitmask_science.shape, coeff_out,
                                             self.par['func2d'])
            # Save to final image
            thismask_science = self.slitmask_science == slit_spat
            self.final_tilts[thismask_science] = self.tilts[thismask_science]

        if show:
            viewer, ch = display.show_image(self.mstilt.image * (self.slitmask > -1), chname='tilts')
            display.show_tilts(viewer, ch, self.make_tbl_tilt_traces())

        if debug:
            show_tilts_mpl(self.mstilt.image*(self.slitmask > -1), self.make_tbl_tilt_traces())

        # Record the Mask
        bpmtilts = np.zeros_like(self.slits.mask, dtype=self.slits.bitmask.minimum_dtype())
        for flag in ['BADTILTCALIB']:
            bpm = self.slits.bitmask.flagged(self.slits.mask, flag)
            if np.any(bpm):
                bpmtilts[bpm] = self.slits.bitmask.turn_on(bpmtilts[bpm], flag)

        # Build and return DataContainer
        tilts_dict = {'coeffs':self.coeffs,
                      'func2d':self.par['func2d'], 'nslit':self.slits.nslits,
                      'spat_order':self.spat_order, 'spec_order':self.spec_order,
                      'spat_id':self.slits.spat_id, 'bpmtilts': bpmtilts,
                      'spat_flexure': self.spat_flexure, 'PYP_SPEC': self.spectrograph.name,
                      'slits_calib_key': self.slits.master_key, 'tilt_traces': self.make_tbl_tilt_traces()}
        return WaveTilts(**tilts_dict)

    def make_tbl_tilt_traces(self):
        """
        Make a Table of the  traced and fitted tilts from self.all_trace_dict
        to save into the WaveTilts object

        Returns:
            `astropy.table.Table`_: Table including the traced and fitted tilts for each slit.
            Columns are:

                - slit_ids: Slit IDs for which the tilts were traced and fit
                - goodpix_spat: Good spatial pixels of the traced tilts
                - goodpix_tilt: Good spectral pixels  of the traced tilts
                - goodpix_slitid: Slit IDs for each goodpix. This is needed to associate goodpix to each slit
                - badpix_spat: Masked spatial pixels of the traced tilts
                - badpix_tilt: Masked spectral pixels of the traced tilts
                - badpix_slitid: Slit IDs for each badpix. This is needed to associate badpix to each slit
                - good2dfit_spat: Good spatial pixels of the 2D fit of the tilts
                - good2dfit_tilt: Good spectral pixels of the 2D fit of the tilts
                - good2dfit_slitid: Slit IDs for each good2dfit. This is needed to associate good2dfit to each slit
                - bad2dfit_spat: Rejected spatial pixels of the 2D fit of the tilts
                - bad2dfit_tilt: Rejected spectral pixels of the 2D fit of the tilts
                - bad2dfit_slitid: Slit IDs for each bad2dfit. This is needed to associate bad2dfit to each slit

        """

        if self.all_trace_dict is None:
            msgs.error('No tilts have been traced and fit yet. Run the run() method first.')

        # slit_ids
        slit_ids = np.array([])

        # good&bad traces
        goodpix_spat = np.array([])
        goodpix_tilt = np.array([])
        goodpix_slitid = np.array([])
        badpix_spat = np.array([])
        badpix_tilt = np.array([])
        badpix_slitid = np.array([])

        # good&bad fits
        good2dfit_spat = np.array([])
        good2dfit_tilt = np.array([])
        good2dfit_slitid = np.array([])
        bad2dfit_spat = np.array([])
        bad2dfit_tilt = np.array([])
        bad2dfit_slitid = np.array([])

        # fill the arrays with the traced and fitted tilts
        # each trc is a slit
        for i, trc in enumerate(self.all_trace_dict):
            if trc is not None:
                slit_ids = np.append(slit_ids, self.slits.spat_id[i])
                gpix = trc['tot_mask']
                goodpix_spat = np.append(goodpix_spat, trc['tilts_spat'].flatten()[gpix.flatten()])
                goodpix_tilt = np.append(goodpix_tilt, trc['tilts'].flatten()[gpix.flatten()])
                goodpix_slitid = np.append(goodpix_slitid, np.repeat(self.slits.spat_id[i], sum(gpix.flatten())))
                bpix = np.invert(gpix) & (trc['tilts'] > 0)
                badpix_spat = np.append(badpix_spat, trc['tilts_spat'].flatten()[bpix.flatten()])
                badpix_tilt = np.append(badpix_tilt, trc['tilts'].flatten()[bpix.flatten()])
                badpix_slitid = np.append(badpix_slitid, np.repeat(self.slits.spat_id[i], sum(bpix.flatten())))
                gfit = gpix & trc['fit_mask']
                good2dfit_spat = np.append(good2dfit_spat, trc['tilts_spat'].flatten()[gfit.flatten()])
                good2dfit_tilt = np.append(good2dfit_tilt, trc['tilt_2dfit'].flatten()[gfit.flatten()])
                good2dfit_slitid = np.append(good2dfit_slitid, np.repeat(self.slits.spat_id[i], sum(gfit.flatten())))
                bfit = gpix & np.invert(trc['fit_mask'])
                bad2dfit_spat = np.append(bad2dfit_spat, trc['tilts_spat'].flatten()[bfit.flatten()])
                bad2dfit_tilt = np.append(bad2dfit_tilt, trc['tilt_2dfit'].flatten()[bfit.flatten()])
                bad2dfit_slitid = np.append(bad2dfit_slitid, np.repeat(self.slits.spat_id[i], sum(bfit.flatten())))

        # fill the table
        tbl_tilt_traces = table.Table()
        # slit ids
        tbl_tilt_traces['slit_ids'] = np.expand_dims(slit_ids, axis=0)
        # traced tilts
        tbl_tilt_traces['goodpix_spat'] = np.expand_dims(goodpix_spat, axis=0)
        tbl_tilt_traces['goodpix_tilt'] = np.expand_dims(goodpix_tilt, axis=0)
        tbl_tilt_traces['goodpix_slitid'] = np.expand_dims(goodpix_slitid, axis=0)
        tbl_tilt_traces['badpix_spat'] = np.expand_dims(badpix_spat, axis=0)
        tbl_tilt_traces['badpix_tilt'] = np.expand_dims(badpix_tilt, axis=0)
        tbl_tilt_traces['badpix_slitid'] = np.expand_dims(badpix_slitid, axis=0)
        # fitted tilts
        tbl_tilt_traces['good2dfit_spat'] = np.expand_dims(good2dfit_spat, axis=0)
        tbl_tilt_traces['good2dfit_tilt'] = np.expand_dims(good2dfit_tilt, axis=0)
        tbl_tilt_traces['good2dfit_slitid'] = np.expand_dims(good2dfit_slitid, axis=0)
        tbl_tilt_traces['bad2dfit_spat'] = np.expand_dims(bad2dfit_spat, axis=0)
        tbl_tilt_traces['bad2dfit_tilt'] = np.expand_dims(bad2dfit_tilt, axis=0)
        tbl_tilt_traces['bad2dfit_slitid'] = np.expand_dims(bad2dfit_slitid, axis=0)

        return tbl_tilt_traces

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


def show_tilts_mpl(tilt_img, tilt_traces, slits=None, show_traces=False, cut=None):
    """Show the TiltImage with the traced and 2D fitted tilts overlaid

    Args:
        tilt_img (`np.ndarray`_):
            TiltImage
        tilt_traces (`astropy.table.Table`_):
            Table containing the traced and fitted tilts.
            See :func:`make_tbl_tilt_traces` for information on the table columns.
        slits (`slittrace.SlitTraceSet`_, optional):
            SlitTraceSet object
        show_traces (bool, optional):
            Show the traced tilts
        cut (tuple, optional):
            Lower and upper levels cut for the image display
    """

    if cut is None:
        cut = utils.growth_lim(tilt_img, 0.98, fac=1)

    w, h = plt.figaspect(1)
    fig = plt.figure(figsize=(1.5 * w, 1.5 * h))

    plt.imshow(tilt_img, origin='lower', interpolation='nearest', aspect='auto', cmap='gray',
               vmin=cut[0], vmax=cut[1], zorder=0)

    # traced tilts
    if show_traces:
        if tilt_traces['goodpix_tilt'][0].size > 0:
            plt.scatter(tilt_traces['goodpix_spat'][0], tilt_traces['goodpix_tilt'][0], color='cyan',
                        marker='s', s=10, lw=0, label='Good pixel', zorder=1)
    # masked tilts
    if tilt_traces['badpix_tilt'][0].size > 0:
        plt.scatter(tilt_traces['badpix_spat'], tilt_traces['badpix_tilt'], color='red',
                    marker='s', s=10, lw=0, label='Masked pixel', zorder=2)
    # 2D fit tilts
    # loop over each slit/line so that we can plot lines instead of points
    for islit in tilt_traces['slit_ids'][0]:
        this_slit = tilt_traces['good2dfit_slitid'][0] == islit
        if np.any(this_slit):
            good2dfit_spat = tilt_traces['good2dfit_spat'][0][this_slit]
            good2dfit_tilt = tilt_traces['good2dfit_tilt'][0][this_slit]
            # try to split the array into each line. This is for improving speed
            # sort tilts
            tsort = np.argsort(good2dfit_tilt)
            # find values that are close together
            close = np.diff(good2dfit_tilt[tsort]) <= 1.
            # divide into slices, each including a line
            lines = utils.contiguous_true(close)
            # loop over each line, this allows to use type='path' and therefore a faster plotting
            for iline in lines:
                x = good2dfit_spat[tsort][iline]
                y = good2dfit_tilt[tsort][iline]
                plt.plot(x, y, color='blue', lw=1, zorder=3)

    # rejected point in 2D fit
    if tilt_traces['bad2dfit_tilt'][0].size > 0:
        plt.scatter(tilt_traces['bad2dfit_spat'][0], tilt_traces['bad2dfit_tilt'][0], color='yellow',
                    marker='s', s=5, lw=0, label='Rejected in fit', zorder=4)

    if slits is not None:
        pstep = 50
        left, right, mask = slits.select_edges()
        gpm = mask == 0
        for i in range(left.shape[1]):
            if gpm[i]:
                spec = np.arange(left[:, i].size)
                plt.plot(left[::pstep, i], spec[::pstep], color='green', lw=2, zorder=5)
                plt.plot(right[::pstep, i], spec[::pstep], color='magenta', lw=2, zorder=5)
                x_spatid = left[spec.size//2, i] + 0.10*(right[spec.size//2, i] - left[spec.size//2, i])
                y_spatid = spec[spec.size//2]
                plt.text(x_spatid, y_spatid, str(slits.spat_id[i]), color='aquamarine', fontsize=10, alpha=1,
                         weight='bold', rotation='vertical', zorder=5)

    # add legend for 2D fitted tilts
    legend_element = [Line2D([0], [0], color='blue', ls='-', lw=1, label='Good tilt fit')]
    handles, labels = plt.gca().get_legend_handles_labels()
    handles.append(legend_element[0])

    lcolors = ['cyan', 'red', 'yellow', 'blue'] if show_traces else ['red', 'yellow', 'blue']
    plt.legend(handles=handles, labelcolor=lcolors)
    plt.ylabel('Spectral pixel index')
    plt.xlabel('Spatial pixel index')
    plt.show()

