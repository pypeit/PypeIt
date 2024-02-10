"""
Module for guiding Arc/Sky line tracing

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import os
import copy
import inspect

from IPython import embed
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

from astropy import stats, visualization
from astropy import table

from pypeit import msgs, datamodel, utils
from pypeit import calibframe
from pypeit import slittrace, wavecalib
from pypeit.display import display
from pypeit.core import arc
from pypeit.core import tracewave
from pypeit.images import buildimage


class WaveTilts(calibframe.CalibFrame):
    """
    Calibration frame containing the wavelength tilt calibration.

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

    # Calibration frame attributes
    calib_type = 'Tilts'
    calib_file_format = 'fits'

    # NOTE:
    #   - Internals are identical to the base class
    #   - Datamodel already contains CalibFrame base elements, so no need to
    #     include it here.

    datamodel = {'PYP_SPEC': dict(otype=str, descr='PypeIt spectrograph name'),
                 'coeffs': dict(otype=np.ndarray, atype=np.floating,
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
                 'spat_flexure': dict(otype=float, descr='Flexure shift from the input TiltImage'),
                 'slits_filename': dict(otype=str, descr='Path to SlitTraceSet file. This helps to '
                                                         'find the Slits calibration file when running '
                                                         'pypeit_chk_tilts()'),
                 'tiltimg_filename': dict(otype=str, descr='Path to Tiltimg file. This helps to '
                                                          'find Tiltimg file when running pypeit_chk_tilts()'),
                 'tilt_traces': dict(otype=table.Table, descr='Table with the positions of the '
                                                              'traced and fitted tilts for all the slits. '
                                                              'see :func:`~pypeit.wavetilts.BuildWaveTilts.make_tbl_tilt_traces` for more details. ')
                 }

    def __init__(self, coeffs, nslit, spat_id, spat_order, spec_order, func2d, bpmtilts=None,
                 spat_flexure=None, PYP_SPEC=None, slits_filename=None, tiltimg_filename=None,
                 tilt_traces=None):

        # Parse
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        d = dict([(k,values[k]) for k in args[1:]])
        # Setup the DataContainer
        datamodel.DataContainer.__init__(self, d=d)

    def _bundle(self):
        """
        Bundle the data in preparation for writing to a fits file.

        See :func:`~pypeit.datamodel.DataContainer._bundle`. Data is
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
            slits (:class:`~pypeit.slittrace.SlitTraceSet`):

        """
        if not np.array_equal(self.spat_id, slits.spat_id):
            msgs.error('Your tilt solutions are out of sync with your slits.  Remove calibrations '
                       'and restart from scratch.')

    def fit2tiltimg(self, slitmask, flexure=None):
        """
        Generate a tilt image from the fit parameters

        Mainly to allow for flexure

        Args:
            slitmask (`numpy.ndarray`_):
                ??
            flexure (float, optional):
                Spatial shift of the tilt image onto the desired frame
                (typically a science image)

        Returns:
            `numpy.ndarray`_:  New tilt image

        """
        msgs.info("Generating a tilts image from the fit parameters")

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

    def show(self, waveimg=None, wcs_match=True, in_ginga=True, show_traces=False,
             chk_version=True):
        """
        Show in ginga or mpl Tiltimg with the tilts traced and fitted overlaid

        Args:
            waveimg (`numpy.ndarray`_, optional):
                Image with the wavelength solution.
            wcs_match (bool, optional):
                If True, use this image as a reference image for the WCS and match all
                image in other channels to it.
            in_ginga (bool, optional):
                If True, show the image in ginga. Otherwise, use matplotlib.
            show_traces (bool, optional):
                If True, show the traces of the tilts on the image.
            chk_version (:obj:`bool`, optional):
                When reading in existing files written by PypeIt, perform strict
                version checking to ensure a valid file.  If False, the code
                will try to keep going, but this may lead to faults and quiet
                failures.  User beware!
        """
        # get tilt_img_dict
        if (Path(self.calib_dir).resolve() / self.tiltimg_filename).exists():
            cal_file = Path(self.calib_dir).resolve() / self.tiltimg_filename
            tilt_img_dict = buildimage.TiltImage.from_file(cal_file, chk_version=chk_version)
        else:
            msgs.error(f'Tilt image {str((Path(self.calib_dir).resolve() / self.tiltimg_filename))} NOT FOUND.')

        # get slits
        slitmask = None
        if (Path(self.calib_dir).resolve() / self.slits_filename).exists():
            cal_file = Path(self.calib_dir).resolve() / self.slits_filename
            slits = slittrace.SlitTraceSet.from_file(cal_file, chk_version=chk_version)
            _slitmask = slits.slit_img(initial=True, flexure=self.spat_flexure)
            _left, _right, _mask = slits.select_edges(flexure=self.spat_flexure)
            gpm = _mask == 0
            # resize
            slitmask = arc.resize_mask2arc(tilt_img_dict.image.shape, _slitmask)
            left = arc.resize_slits2arc(tilt_img_dict.image.shape, _slitmask.shape, _left)
            right = arc.resize_slits2arc(tilt_img_dict.image.shape, _slitmask.shape, _right)
        else:
            slits = None
            msgs.warn('Could not load slits to show with tilts image.')

        # get waveimg
        same_size = (slits.nspec, slits.nspat) == tilt_img_dict.image.shape
        if waveimg is None and slits is not None and same_size and in_ginga:
            wv_calib_name = wavecalib.WaveCalib.construct_file_name(self.calib_key, calib_dir=self.calib_dir)
            if Path(wv_calib_name).resolve().exists():
                wv_calib = wavecalib.WaveCalib.from_file(wv_calib_name, chk_version=chk_version)
                tilts = self.fit2tiltimg(slitmask, flexure=self.spat_flexure)
                waveimg = wv_calib.build_waveimg(tilts, slits, spat_flexure=self.spat_flexure)
            else:
                msgs.warn('Could not load Wave image to show with tilts image.')

        # Show
        # tilt image
        tilt_img = tilt_img_dict.image * (slitmask > -1) if slitmask is not None else tilt_img_dict.image
        # set cuts
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
            # show slit edges
            if slits is not None:
                display.show_slits(viewer, ch, left[:, gpm], right[:, gpm], slit_ids=slits.spat_id[gpm])

        # show in matplotlib
        else:
            if slits is not None:
                show_tilts_mpl(tilt_img, self.tilt_traces, show_traces=show_traces, left_edges=left[:, gpm],
                               right_edges=right[:, gpm], slit_ids=slits.spat_id[gpm], cut=cut)
            else:
                show_tilts_mpl(tilt_img, self.tilt_traces, show_traces=show_traces, cut=cut)


class BuildWaveTilts:
    """
    Class to guide arc/sky tracing

    Args:
        mstilt (:class:`~pypeit.images.buildimage.TiltImage`):
            Tilt image.  QA file naming inherits the calibration key
            (``calib_key``) from this object.
        slits (:class:`~pypeit.slittrace.SlitTraceSet`):
            Slit edges
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph object
        par (:class:`~pypeit.par.pypeitpar.WaveTiltsPar` or None):
            The parameters used to fuss with the tilts
        wavepar (:class:`~pypeit.par.pypeitpar.WavelengthSolutionPar` or None):
            The parameters used for the wavelength solution
        det (int): Detector index
        qa_path (:obj:`str`, optional):
            Directory for QA output.
        spat_flexure (float, optional):
            If input, the slitmask and slit edges are shifted prior
            to tilt analysis.


    Attributes:
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            ??
        steps (list):
            ??
        mask (`numpy.ndarray`_):
            boolean array; True = Ignore this slit
        all_trcdict (list):
            All trace dict's
        tilts (`numpy.ndarray`_):
            Tilts for a single slit/order
        all_tilts (list): 
            Tuple of tilts `numpy.ndarray`_ objects
        final_tilts (`numpy.ndarray`_):
            Final tilts image
        gpm (`numpy.ndarray`_):
            Good pixel mask.  Eventually, we might attach this to self.mstilt
            although that would then require that we write it to disk with
            self.mstilt.image
    """

    # TODO This needs to be modified to take an inmask
    def __init__(self, mstilt, slits, spectrograph, par, wavepar, det=1, qa_path=None,
                 spat_flexure=None):

        # TODO: Perform type checking
        self.spectrograph = spectrograph
        self.par = par
        self.wavepar = wavepar

        self.mstilt = mstilt
        self.slits = slits
        self.det = det
        self.qa_path = qa_path
        self.spat_flexure = spat_flexure

        # Get the non-linear count level
        # TODO: This is currently hacked to deal with Mosaics
        try:
            self.nonlinear_counts = self.mstilt.detector.nonlinear_counts()
        except:
            self.nonlinear_counts = 1e10

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
            arcspec ():
                ??
            slit_cen ():
                ??
            slit_idx (int):
                Slit index, zero-based
            bpm (`numpy.ndarray`_, optional):
                ??
            debug (bool, optional):
                ??

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
            slit_cen (`numpy.ndarray`_): (nspec,) Central trace for this slit
            spat_order (int): Order of the 2d polynomial fit for the spatial direction
            spec_order (int): Order of the 2d polytnomial fit for the spectral direction
            slit_idx (int): zero-based, integer index for the slit in question
            show_QA (bool, optional):
                show the QA instead of writing it out to the outfile
            doqa (bool, optional):
                Construct the QA plot

        Returns:
            `numpy.ndarray`_: Array containing the coefficients for the 2d
            legendre polynomial fit.  Shape is (spat_order + 1, spec_order+1).
        """
        # Index
        self.all_fit_dict[slit_idx], self.all_trace_dict[slit_idx] \
                = tracewave.fit_tilts(trc_tilt_dict, thismask, slit_cen, spat_order=spat_order,
                                      spec_order=spec_order,maxdev=self.par['maxdev2d'],
                                      sigrej=self.par['sigrej2d'], func2d=self.par['func2d'],
                                      doqa=doqa, calib_key=self.mstilt.calib_key,
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
            dict: Dictionary containing information on the traced tilts required
            to fit the filts.

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
        :func:`~pypeit.core.fitting.robust_fit` and the parameters
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
            `numpy.ndarray`_: Returns a 2D image with the same shape as
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
        #self.par['rm_continuum'] = True
        #debug = True
        #show = True

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
                ax.set_title('Arc')
                ax = fig.add_axes([1.15/3, 0.1, 0.8/3, 0.8])
                ax.imshow(continuum, origin='lower', interpolation='nearest',
                          aspect='auto', vmin=vmin, vmax=vmax)
                ax.set_title('Continuum')
                ax = fig.add_axes([2.15/3, 0.1, 0.8/3, 0.8])
                ax.imshow(_mstilt, origin='lower', interpolation='nearest',
                          aspect='auto', vmin=vmin, vmax=vmax)
                ax.set_title('Arc - Continuum')
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
                msgs.info(f'Skipping bad slit/order {self.slits.slitord_id[slit_idx]} ({slit_idx+1}/{self.slits.nslits})')
                self.slits.mask[slit_idx] = self.slits.bitmask.turn_on(self.slits.mask[slit_idx], 'BADTILTCALIB')
                continue
            msgs.info(f'Computing tilts for slit/order {self.slits.slitord_id[slit_idx]} ({slit_idx+1}/{self.slits.nslits})')
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
            bpm = self.slits.bitmask.flagged(self.slits.mask, flag=flag)
            if np.any(bpm):
                bpmtilts[bpm] = self.slits.bitmask.turn_on(bpmtilts[bpm], flag)

        # grab slits file_name
        slits_filename = self.slits.construct_file_name(self.slits.calib_key)
        # grab titimg file_name
        tiltimg_filename = self.mstilt.construct_file_name(self.mstilt.calib_key)

        # Build and return the calibration frame
        tilts = WaveTilts(self.coeffs, self.slits.nslits, self.slits.spat_id, self.spat_order,
                          self.spec_order, self.par['func2d'], bpmtilts=bpmtilts,
                          spat_flexure=self.spat_flexure, PYP_SPEC=self.spectrograph.name,
                          slits_filename=slits_filename, tiltimg_filename=tiltimg_filename,
                          tilt_traces=self.make_tbl_tilt_traces())
        # Inherit the calibration frame naming from self.mstilt
        # TODO: Should throw an error here if these calibration frame naming
        # elements are not defined by self.mstilts...
        tilts.copy_calib_internals(self.mstilt)
        return tilts

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
                - goodpix_lid: Line IDs for each goodpix. This is needed to associate goodpix to each line
                - badpix_spat: Masked spatial pixels of the traced tilts
                - badpix_tilt: Masked spectral pixels of the traced tilts
                - badpix_lid: Line IDs for each badpix. This is needed to associate badpix to each line
                - good2dfit_spat: Good spatial pixels of the 2D fit of the tilts
                - good2dfit_tilt: Good spectral pixels of the 2D fit of the tilts
                - good2dfit_lid: Line IDs for each good2dfit. This is needed to associate good2dfit to each line
                - bad2dfit_spat: Rejected spatial pixels of the 2D fit of the tilts
                - bad2dfit_tilt: Rejected spectral pixels of the 2D fit of the tilts
                - bad2dfit_lid: Line IDs for each bad2dfit. This is needed to associate bad2dfit to each line

        """

        if self.all_trace_dict is None:
            msgs.error('No tilts have been traced and fit yet. Run the run() method first.')

        # slit_ids
        slit_ids = np.array([])

        # good&bad traces
        goodpix_lid = np.array([])
        goodpix_spat = np.array([])
        goodpix_tilt = np.array([])

        badpix_lid = np.array([])
        badpix_spat = np.array([])
        badpix_tilt = np.array([])

        # good&bad fits
        good2dfit_lid = np.array([])
        good2dfit_spat = np.array([])
        good2dfit_tilt = np.array([])

        bad2dfit_lid = np.array([])
        bad2dfit_spat = np.array([])
        bad2dfit_tilt = np.array([])

        # fill the arrays with the traced and fitted tilts
        # each trc is a slit
        for i, trc in enumerate(self.all_trace_dict):
            if trc is not None and self.slits.mask[i] == 0:
                # slit ids
                slit_ids = np.append(slit_ids, self.slits.spat_id[i])
                # good pixels
                gpix = trc['tot_mask']
                # bad pixels
                bpix = np.invert(gpix) & (trc['tilts'] > 0)
                # good 2d fit
                gfit = gpix & trc['fit_mask']
                # bad 2d fit
                bfit = gpix & np.invert(trc['fit_mask'])

                for l in range(trc['tilts_spat'].shape[1]):
                    # good pixels
                    goodpix_lid = np.append(goodpix_lid, np.repeat(f'{self.slits.spat_id[i]}_{l + 1}', sum(gpix[:, l])))
                    goodpix_spat = np.append(goodpix_spat, trc['tilts_spat'][:,l][gpix[:,l]])
                    goodpix_tilt = np.append(goodpix_tilt, trc['tilts'][:,l][gpix[:,l]])
                    # bad pixels
                    badpix_lid = np.append(badpix_lid, np.repeat(f'{self.slits.spat_id[i]}_{l + 1}', sum(bpix[:, l])))
                    badpix_spat = np.append(badpix_spat, trc['tilts_spat'][:,l][bpix[:,l]])
                    badpix_tilt = np.append(badpix_tilt, trc['tilts'][:,l][bpix[:,l]])
                    # good 2d fit
                    good2dfit_lid = np.append(good2dfit_lid, np.repeat(f'{self.slits.spat_id[i]}_{l + 1}', sum(gfit[:, l])))
                    good2dfit_spat = np.append(good2dfit_spat, trc['tilts_spat'][:,l][gfit[:,l]])
                    good2dfit_tilt = np.append(good2dfit_tilt, trc['tilt_2dfit'][:,l][gfit[:,l]])
                    # bad 2d fit
                    bad2dfit_lid = np.append(bad2dfit_lid, np.repeat(f'{self.slits.spat_id[i]}_{l + 1}', sum(bfit[:, l])))
                    bad2dfit_spat = np.append(bad2dfit_spat, trc['tilts_spat'][:,l][bfit[:,l]])
                    bad2dfit_tilt = np.append(bad2dfit_tilt, trc['tilt_2dfit'][:,l][bfit[:,l]])

        # fill the table
        tbl_tilt_traces = table.Table()
        trc_arrays = [slit_ids, goodpix_spat, goodpix_tilt, badpix_lid, badpix_spat, badpix_tilt,
                      good2dfit_lid, good2dfit_spat, good2dfit_tilt, bad2dfit_lid, bad2dfit_spat, bad2dfit_tilt]
        tbl_keys = ['slit_ids', 'goodpix_spat', 'goodpix_tilt', 'badpix_lid', 'badpix_spat', 'badpix_tilt',
                    'good2dfit_lid', 'good2dfit_spat', 'good2dfit_tilt', 'bad2dfit_lid', 'bad2dfit_spat', 'bad2dfit_tilt']
        for i, arr in enumerate(trc_arrays):
            if arr.size > 0:
                tbl_tilt_traces[tbl_keys[i]] = np.expand_dims(arr, axis=0)

        if len(tbl_tilt_traces) == 0:
            msgs.warn('No traced and fitted tilts have been found.')
            return None

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


def show_tilts_mpl(tilt_img, tilt_traces, show_traces=False, left_edges=None,
                   right_edges=None, slit_ids=None, cut=None):
    """Show the TiltImage with the traced and 2D fitted tilts overlaid

    Args:
        tilt_img (`numpy.ndarray`_):
            TiltImage
        tilt_traces (`astropy.table.Table`_):
            Table containing the traced and fitted tilts.  See
            :func:`~pypeit.wavetilts.BuiltWaveTilts.make_tbl_tilt_traces` for
            information on the table columns.
        show_traces (bool, optional):
            Show the traced tilts
        left_edges (`numpy.ndarray`_, optional):
            Left edges of the slits
        right_edges (`numpy.ndarray`_, optional):
            Right edges of the slits
        slit_ids (`numpy.ndarray`_, optional):
            Slit IDs
        cut (tuple, optional):
            Lower and upper levels cut for the image display
    """

    if tilt_traces is None:
        return msgs.error('No tilts have been traced or fitted')

    if cut is None:
        cut = utils.growth_lim(tilt_img, 0.98, fac=1)

    w, h = plt.figaspect(1)
    fig = plt.figure(figsize=(1.5 * w, 1.5 * h))

    plt.imshow(tilt_img, origin='lower', interpolation='nearest', aspect='auto', cmap='gray',
               vmin=cut[0], vmax=cut[1], zorder=0)

    # traced tilts
    if show_traces:
        if 'goodpix_tilt' in tilt_traces.keys() and tilt_traces['goodpix_tilt'][0].size > 0:
            plt.scatter(tilt_traces['goodpix_spat'][0], tilt_traces['goodpix_tilt'][0], color='cyan',
                        edgecolors=None, marker='s', s=10, alpha=0.3, lw=0, label='Good pixel', zorder=1)
    # masked tilts
    if 'badpix_tilt' in tilt_traces.keys() and tilt_traces['badpix_tilt'][0].size > 0:
        plt.scatter(tilt_traces['badpix_spat'], tilt_traces['badpix_tilt'], color='red',
                    edgecolors=None, marker='s', s=10, alpha=0.3, lw=0, label='Masked pixel', zorder=2)
    # 2D fit tilts
    # loop over each line so that we can plot lines instead of points
    if 'good2dfit_lid' in tilt_traces.keys():
        for iline in np.unique(tilt_traces['good2dfit_lid'][0]):
            this_line = tilt_traces['good2dfit_lid'][0] == iline
            if np.any(this_line):
                good2dfit_spat = tilt_traces['good2dfit_spat'][0][this_line]
                good2dfit_tilt = tilt_traces['good2dfit_tilt'][0][this_line]
                plt.plot(good2dfit_spat, good2dfit_tilt, color='blue', lw=1, alpha=0.4, zorder=3)

    # rejected point in 2D fit
    if 'bad2dfit_tilt' in tilt_traces.keys() and tilt_traces['bad2dfit_tilt'][0].size > 0:
        plt.scatter(tilt_traces['bad2dfit_spat'][0], tilt_traces['bad2dfit_tilt'][0], color='yellow',
                    edgecolors=None, marker='s', s=5, alpha=0.3, lw=0, label='Rejected in fit', zorder=4)

    if left_edges is not None and right_edges is not None:
        pstep = 50
        for i in range(left_edges.shape[1]):
            spec = np.arange(left_edges[:, i].size)
            plt.plot(left_edges[::pstep, i], spec[::pstep], color='green', lw=2, zorder=5)
            plt.plot(right_edges[::pstep, i], spec[::pstep], color='magenta', lw=2, zorder=5)
            x_spatid = left_edges[spec.size//2, i] + 0.10*(right_edges[spec.size//2, i] - left_edges[spec.size//2, i])
            y_spatid = spec[spec.size//2]
            if slit_ids is not None:
                plt.text(x_spatid, y_spatid, str(slit_ids[i]), color='aquamarine', fontsize=10, alpha=1,
                         weight='bold', rotation='vertical', zorder=5)

    # add legend for 2D fitted tilts
    legend_element = [Line2D([0], [0], color='blue', ls='-', lw=1, label='Good tilt fit')]
    handles, labels = plt.gca().get_legend_handles_labels()
    handles.append(legend_element[0])

    lcolors = ['cyan', 'red', 'yellow', 'blue'] if show_traces else ['red', 'yellow', 'blue']
    legend = plt.legend(handles=handles, labelcolor=lcolors, loc=3, markerscale=2, frameon=False)
    for h in legend.legendHandles:
        h.set_alpha(1)
    
    plt.ylabel('Spectral pixel index')
    plt.xlabel('Spatial pixel index')
    plt.show()

