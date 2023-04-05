"""
Module for generating an Alignment image to map constant spatial locations

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import inspect
import numpy as np
from IPython import embed
from scipy.interpolate import interp1d, RegularGridInterpolator

from pypeit.display import display
from pypeit.core import findobj_skymask
from pypeit import datamodel
from pypeit import calibframe
from pypeit import msgs


class Alignments(calibframe.CalibFrame):
    """
    Calibration frame holding result of slit alignment processing.

    All of the items in the datamodel are required for instantiation, although
    they can be None (but shouldn't be)

    The datamodel attributes are:

    .. include:: ../include/class_datamodel_alignments.rst
    """
    version = '1.1.0'

    # Calibration frame attributes
    calib_type = 'Alignment'
    calib_file_format = 'fits'

    # Datamodel already includes PYP_SPEC, so no need to combine it with the
    # CalibFrame base datamodel.
    datamodel = {'PYP_SPEC': dict(otype=str, descr='PypeIt spectrograph name'),
                 'alignframe': dict(otype=np.ndarray, atype=np.floating,
                                    descr='Processed, combined alignment frames'),
                 'nspec': dict(otype=int, descr='The number of spectral elements'),
                 'nalign': dict(otype=int, descr='Number of alignment traces in each slit'),
                 'nslits': dict(otype=int, descr='The number of slits'),
                 'traces': dict(otype=np.ndarray, atype=np.floating,
                                descr='Traces of the alignment frame'),
                 'spat_id': dict(otype=np.ndarray, atype=np.integer, descr='Slit spat_id ')}

    def __init__(self, alignframe=None, nspec=None, nalign=None, nslits=None,
                 traces=None, PYP_SPEC=None, spat_id=None):
        # Parse
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        d = dict([(k,values[k]) for k in args[1:]])
        # Setup the DataContainer
        datamodel.DataContainer.__init__(self, d=d)

    def _validate(self):
        # TBC - need to check that all alignment traces have been correctly traced
        pass

    # NOTE: If you make changes to how this object is bundled into the output
    # datamodel, make sure you update the documentation in
    # doc/calibrations/align.rst!
    def _bundle(self):
        """
        Override the base class method simply to set the HDU extension name.
        """
        return super()._bundle(ext='ALIGN')

    def is_synced(self, slits):
        """
        Confirm the slits in the alignment are the same as that in SlitTraceSet

        Barfs if not

        Args:
            slits (:class:`pypeit.slittrace.SlitTraceSet`):

        """
        if not np.array_equal(self.spat_id, slits.spat_id):
            msgs.error('Your alignment solutions are out of sync with your slits.  Remove '
                       'Calibrations and restart from scratch.')

    def show(self, slits=None):
        """
        Simple wrapper for :func:`show_alignment`.

        Args:
            slits (:class:`pypeit.slittrace.SlitTraceSet`, optional):
                Slit properties, including traces.
        """
        # Show
        show_alignment(self.alignframe, align_traces=self.traces, slits=slits)


class TraceAlignment:
    """
    Class to guide the determination of the alignment traces

    Args:
        rawalignimg (:class:`~pypeit.images.pypeitimage.PypeItImage`):
            Align image, created by the AlignFrame class. Can be
            None.
        slits (:class:`~pypeit.slittrace.SlitTraceSet`):
            Slit edge traces.  Can be None.
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            The `Spectrograph` instance that sets the instrument used
            to take the observations. Can be None.
        alignpar (:class:`~pypeit.par.pypeitpar.AlignPar`):
            The parameters used for the align traces
        det (:obj:`int`, optional):
            Detector number
        qa_path (:obj:`str`, optional):
            Directory for QA plots
        msbpm (`numpy.ndarray`_, optional):
            Bad pixel mask image

    Attributes:
        steps (:obj:`list`):
            List of the processing steps performed.
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            Relevant spectrograph.
        slits (:class:`~pypeit.slittrace.SlitTraceSet`):
            Slit edge traces.
    """
    def __init__(self, rawalignimg, slits, spectrograph, alignpar, det=1, qa_path=None,
                 msbpm=None):

        # Defaults
        self.spectrograph = spectrograph
        self.PYP_SPEC = spectrograph.name
        # Alignment parameters
        self.alignpar = alignpar

        # Input data
        self.rawalignimg = rawalignimg
        self.slits = slits

        # Optional parameters
        self.bpm = rawalignimg.fullmask.bpm if msbpm is None else msbpm
        self.qa_path = qa_path
        self.det = det

        # Attributes unique to this object
        self._alignprof = None

        # Completed steps
        self.steps = []

    @property
    def nslits(self):
        """
        Return the number of slits.  Pulled directly from :attr:`slits`, if it exists.
        """
        return 0 if self.slits is None else self.slits.nslits

    @property
    def nspec(self):
        """
        Return the number of spectral elements.  Pulled directly from :attr:`slits`, if it exists.
        """
        return self.rawalignimg.shape[0] if self.slits is None else self.slits.nspec

    def build_traces(self, show_peaks=False, debug=False):
        """
        Main routine to generate the align profile traces in all slits

        Args:
             show_peaks (bool, optional):
               Generate QA showing peaks identified by alignment profile tracing
             show_trace (bool, optional):
               Generate QA showing traces identified. Requires an open ginga RC modules window.
               Launch with ``ginga --modules=RC,SlitWavelength &``
             debug (bool, optional):
               Debug the alignment tracing algorithm

        Returns:
            dict:  self.align_dict
        """
        # Generate slits
        slitid_img_init = self.slits.slit_img(initial=True)
        left, right, _ = self.slits.select_edges(initial=True)
        align_prof = dict({})
        # Go through the slits
        for slit_idx, slit_spat in enumerate(self.slits.spat_id):
            specobj_dict = {'SLITID': slit_idx, 'DET': self.rawalignimg.detector.name,
                            'OBJTYPE': "align_profile", 'PYPELINE': self.spectrograph.pypeline}
            msgs.info("Fitting alignment traces in slit {0:d}".format(slit_idx))
            align_traces = findobj_skymask.objs_in_slit(
                self.rawalignimg.image, self.rawalignimg.ivar, slitid_img_init == slit_spat,
                left[:, slit_idx], right[:, slit_idx],
                ncoeff=self.alignpar['trace_npoly'],
                specobj_dict=specobj_dict, snr_thresh=self.alignpar['snr_thresh'],
                show_peaks=show_peaks, show_fits=False,
                trim_edg=self.alignpar['trim_edge'],
                nperslit=len(self.alignpar['locations']))
            if len(align_traces) != len(self.alignpar['locations']):
                # Align tracing has failed for this slit
                msgs.error("Alignment tracing has failed on slit {0:d}".format(slit_idx))
            align_prof['{0:d}'.format(slit_idx)] = align_traces.copy()

        # Steps
        self.steps.append(inspect.stack()[0][3])

        # Return
        return align_prof

    def generate_traces(self, align_prof):
        """
        Generate a dictionary containing all of the information from the profile fitting

        Parameters
        ----------
        align_prof : :obj:`dict`
            Dictionary of SpecObjs classes (one for each slit)

        Returns
        -------
        align_traces : `numpy.ndarray`_
            Spatial traces (3D array of shape [nspec, ntraces, nslits])
        """
        nbars = len(self.alignpar['locations'])
        # Generate an array containing the centroid of all bars
        align_traces = np.zeros((self.nspec, nbars, self.nslits))
        for sl in range(self.nslits):
            sls = '{0:d}'.format(sl)
            for bar in range(nbars):
                if align_prof[sls][bar].SLITID != sl:
                    msgs.error("Alignment profiling failed to generate traces")
                align_traces[:, bar, sl] = align_prof[sls][bar].TRACE_SPAT
        return align_traces

    def run(self, show=False):
        """
        Main driver for alignment profile tracing

        Args:
            show (bool, optional):
                Show the alignment traces?

        Returns:
            :class:`pypeit.alignframe.Alignments`:

        """
        ###############
        # Fill up the calibrations and generate QA
        align_prof = self.build_traces()

        # Prepare the dictionary items for the data container
        align_traces = self.generate_traces(align_prof)

        if show:
            show_alignment(self.rawalignimg.image, align_traces=align_traces, slits=self.slits)

        # Build the alignment calibration frame
        align = Alignments(alignframe=self.rawalignimg.image, nspec=self.nspec,
                           nalign=align_traces.shape[1], nslits=self.nslits, traces=align_traces,
                           PYP_SPEC=self.PYP_SPEC, spat_id=self.slits.spat_id)
        # Copy the internals from the processed alignment image
        align.copy_calib_internals(self.rawalignimg)
        # Return
        return align


def show_alignment(alignframe, align_traces=None, slits=None, clear=False):
    """
    Show one of the class internals

    Parameters
    ----------

    alignframe : `numpy.ndarray`_
        Image to be plotted (i.e. the align frame)
    align_traces : list, optional
        The align traces
    slits : :class:`pypeit.slittrace.SlitTraceSet`, optional
        properties of the slits, including traces.
    clear : bool, optional
        Clear the plotting window in ginga?

    Returns
    -------

    """
    display.connect_to_ginga(raise_err=True, allow_new=True)
    ch_name = 'alignment'
    viewer, channel = display.show_image(alignframe, chname=ch_name, clear=clear, wcs_match=False)

    # Display the slit edges
    if slits is not None and viewer is not None:
        left, right, mask = slits.select_edges()
        display.show_slits(viewer, channel, left, right)

    # Display the alignment traces
    if align_traces is not None and viewer is not None:
        for bar in range(align_traces.shape[1]):
            for slt in range(align_traces.shape[2]):
                # Alternate the colors of the slits
                color = 'orange'
                if slt%2 == 0:
                    color = 'magenta'
                # Display the trace
                display.show_trace(viewer, channel, align_traces[:, bar, slt], trc_name="",
                                   color=color)


class AlignmentSplines:
    def __init__(self, traces, locations, tilts):
        """Convenience class to build and transform between detector pixel coordinates and
        WCS pixel coordinates (i.e. constant wavelength and spatial position).

        Parameters
        ----------
        traces : `numpy.ndarray`
            3D array containing the alignments (traces) of the slits. This allows different
            slices to be aligned correctly. Ideally, this variable will be assigned the
            value of alignments.traces. However, this can also be assigned the
            left and right slit edges:

            .. code-block:: python

                traces = np.append(left.reshape((left.shape[0], 1, left.shape[1])),
                       right.reshape((left.shape[0], 1, left.shape[1])), axis=1)
                # In this case, you should set the locations argument to
                locations=np.array([0,1])

        locations : `numpy.ndarray`_, list
            locations along the slit of the alignment traces. Must
            be a 1D array of the same length as alignments.traces.shape[1]
        tilts : `numpy.ndarray`
            Spectral tilts.
        """
        # Perform some checks
        msgs.work("Spatial flexure is not currently implemented for the astrometric alignment")
        if type(locations) is list:
            locations = np.array(locations)
        if locations.size != traces.shape[1]:
            msgs.error("The size of locations must be the same as traces.shape[1]")
        # Store the relevant input
        self.traces = traces
        self.locations = locations
        self.tilts = tilts
        self.nspec, self.ntrace, self.nslit = traces.shape
        # Transform between detector pixels and the locations:
        self.spl_loc = self.nslit * [self.nspec*[None]]  # Splines - map (x,y) pixels ==> tilts
        self.spl_slen = self.nslit * [None]  # Splines - map y pixel ==> slit length
        self.spl_transform = self.nslit * [None]  # Splines - map x,y pixel ==> offset in pixels from the central trace
        self.spl_fulltilts = RegularGridInterpolator((np.arange(tilts.shape[0]), np.arange(tilts.shape[1])),
                                                     tilts * (self.nspec - 1), method='linear')
        self.build_splines()

    def build_splines(self):
        """
        Build the interpolation transforms for each slit
        """
        spldict = dict(kind='linear', bounds_error=False, fill_value="extrapolate")
        ycoord = np.arange(self.nspec)
        for sl in range(self.nslit):
            msgs.info("Calculating astrometric transform of slit {0:d}/{1:d}".format(sl+1, self.nslit))
            xlr, tlr = np.zeros((self.nspec, 2)), np.zeros((self.nspec, 2))
            eval_trim = 2  # This evaluates the slit length inside the actual slit edges, due to edge effects.
            for sp in range(self.nspec):
                # Calculate x coordinate at the slit edges, and the spectral tilts at those locations
                xlr[sp, :] = interp1d(self.locations, self.traces[sp,:,sl], **spldict)([0.0, 1.0])
                tmptilt = self.spl_fulltilts([[sp, xlr[sp,0] + eval_trim],
                                              [sp, xlr[sp,0] + eval_trim+1],
                                              [sp, xlr[sp,1] - eval_trim],
                                              [sp, xlr[sp,1] - eval_trim-1]])
                tlr[sp, :] = [tmptilt[0]-eval_trim*(tmptilt[1]-tmptilt[0]),
                              tmptilt[2]-eval_trim*(tmptilt[3]-tmptilt[2])]
                # pseudo-2D alignments -> locations
                self.spl_loc[sl][sp] = interp1d(self.traces[sp,:,sl], self.locations, **spldict)
            # For a given tilt value, get the (x,y) coordinate of the right edge
            tilt_ypos = interp1d(tlr[:,1], ycoord, **spldict)
            ypos_xpos = interp1d(ycoord, xlr[:,1], **spldict)
            ypos = tilt_ypos(tlr[:,0])
            xpos = ypos_xpos(ypos)
            # Calculate the slitlength from (xlr, y), (xpos, ypos) coordinates
            slitlen = np.sqrt((ypos - np.arange(ypos.size)) ** 2 + (xpos - xlr[:, 0]) ** 2)
            self.spl_slen[sl] = interp1d(ycoord, slitlen, **spldict)  # The tilt used to calculate the slit length corresponds to the left edge, so use ycoord for the first argument
            # Construct a 2D Regular grid that interpolates over all coordinates
            xcoord = np.arange(np.floor(np.min(xlr)), np.ceil(np.max(xlr))+1, 1.0)
            out_transform = np.zeros((self.nspec, xcoord.size))
            for sp in range(self.nspec):
                out_transform[sp,:] = (self.spl_loc[sl][sp](xcoord) - 0.5) * self.spl_slen[sl](sp)
            self.spl_transform[sl] = RegularGridInterpolator((ycoord, xcoord), out_transform, method='linear',
                                                             bounds_error=False, fill_value=None) # This will extrapolate
            # TODO :: Remove these notes...
            # We now have everything we need to calculate the location and tilt of every pixel in the image.
            # evalpos = (self.spl_loc[sl][ypixels](xpixels) - 0.5) * self.spl_slen[sl](ypixels)
            # wcs.wcs_pix2world(slitID, evalpos, tilts[onslit_init] * (nspec - 1), 0)

    def transform(self, slitnum, spatpix, specpix):
        """
        Convenience function to return the spatial offset in pixels
        from the spatial center of the slit.

        Parameters
        ----------
        slitnum : `int`
            Slit number
        spatpix : `numpy.ndarray`
            Detector pixel coordinate (spatial direction)
        specpix : `numpy.ndarray`
            Detector pixel coordinate (spectral direction)

        Returns
        -------
        spl_transform : `numpy.ndarray`
            The spatial offset (measured in pixels) from the center of the slit.
        """
        return self.spl_transform[slitnum]((specpix, spatpix))
