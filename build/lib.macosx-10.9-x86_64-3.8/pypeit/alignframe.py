"""
Module for generating an Alignment image to map constant spatial locations

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import inspect
import numpy as np
from IPython import embed

from pypeit.display import display
from pypeit.core import extract
from pypeit import datamodel, msgs


class Alignments(datamodel.DataContainer):
    """
    Simple DataContainer for the alignment output

    All of the items in the datamodel are required for instantiation,
      although they can be None (but shouldn't be)

    """
    minimum_version = '1.1.0'
    version = '1.1.0'

    # I/O
    output_to_disk = None  # This writes all items that are not None
    hdu_prefix = None

    # Master fun
    master_type = 'Alignment'
    master_file_format = 'fits'

    datamodel = {'alignframe': dict(otype=np.ndarray, atype=np.floating,
                                    descr='Processed, combined alignment frames'),
                 'nspec': dict(otype=int, descr='The number of spectral elements'),
                 'nalign': dict(otype=int, descr='Number of alignment traces in each slit'),
                 'nslits': dict(otype=int, descr='The number of slits'),
                 'traces': dict(otype=np.ndarray, atype=np.floating,
                                descr='Traces of the alignment frame'),
                 'PYP_SPEC': dict(otype=str, descr='PypeIt spectrograph name'),
                 'spat_id': dict(otype=np.ndarray, atype=np.integer, descr='Slit spat_id ')}

    def __init__(self, alignframe=None, nspec=None, nalign=None, nslits=None,
                 traces=None, PYP_SPEC=None, spat_id=None):
        # Parse
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        d = dict([(k,values[k]) for k in args[1:]])
        # Setup the DataContainer
        datamodel.DataContainer.__init__(self, d=d)

    def _init_internals(self):
        # Master stuff
        self.master_key = None
        self.master_dir = None

    def _validate(self):
        # TBC - need to check that all alignment traces have been correctly traced
        pass

    def is_synced(self, slits):
        """
        Confirm the slits in the alignment are the same as that in SlitTraceSet

        Barfs if not

        Args:
            slits (:class:`pypeit.slittrace.SlitTraceSet`):

        """
        if not np.array_equal(self.spat_id, slits.spat_id):
            msgs.error("Your alignment solutions are out of sync with your slits.  Remove Masters and start from scratch")

    def show(self, slits=None):
        """
        Simple wrapper to show_alignment()

        Parameters
        ----------

        Returns:

        """
        # Show
        show_alignment(self.alignframe, align_traces=self.traces, slits=slits)


class TraceAlignment(object):
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
        binning (:obj:`str`, optional):
            Detector binning in comma separated numbers for the
            spectral and spatial binning.
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
    version = '1.0.0'

    # Frametype is a class attribute
    frametype = 'alignment'
    """
    Frame type designation.
    """

    master_type = 'Alignment'
    master_file_format = 'fits'

    def __init__(self, rawalignimg, slits, spectrograph, alignpar, det=1,
                 binning=None, qa_path=None, msbpm=None):

        # Defaults
        self.spectrograph = spectrograph
        self.PYP_SPEC = spectrograph.name
        self.binning = binning
        # Alignment parameters
        self.alignpar = alignpar

        # Input data
        self.rawalignimg = rawalignimg
        self.slits = slits

        # Optional parameters
        self.bpm = rawalignimg.mask if msbpm is None else msbpm
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
               Generate QA showing traces identified. Requires an open ginga RC modules window
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
            specobj_dict = {'SLITID': slit_idx, 'DET': self.det, 'OBJTYPE': "align_profile",
                            'PYPELINE': self.spectrograph.pypeline}
            msgs.info("Fitting alignment traces in slit {0:d}".format(slit_idx))
            align_traces, _ = extract.objfind(
                self.rawalignimg.image, slitid_img_init == slit_spat,
                left[:, slit_idx], right[:, slit_idx],
                ir_redux=False, ncoeff=self.alignpar['trace_npoly'],
                specobj_dict=specobj_dict, sig_thresh=self.alignpar['sig_thresh'],
                show_peaks=show_peaks, show_fits=False,
                trim_edg=self.alignpar['trim_edge'],
                cont_fit=False, npoly_cont=0,
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

        return Alignments(alignframe=self.rawalignimg.image,
                          nspec=self.nspec,
                          nalign=align_traces.shape[1],
                          nslits=self.nslits,
                          traces=align_traces,
                          PYP_SPEC=self.PYP_SPEC,
                          spat_id=self.slits.spat_id)


def show_alignment(alignframe, align_traces=None, slits=None, clear=False):
    """
    Show one of the class internals

    Parameters
    ----------

    alignframe : `numpy.ndarray`_
        Image to be plotted (i.e. the master align frame)
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
                display.show_trace(viewer, channel, align_traces[:, bar, slt], trc_name="", color=color)
