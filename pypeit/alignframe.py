"""
Module for generating an Alignment image to map constant spatial locations

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import os
import inspect
import numpy as np
from IPython import embed

from pypeit import ginga, msgs
from pypeit.core import extract
from pypeit.par.pypeitpar import AlignPar
from pypeit import bspline
from pypeit import datamodel
from pypeit import masterframe


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

    datamodel = {
        'alignframe':  dict(otype=np.ndarray, atype=np.floating, desc='Processed, combined alignment frames'),
        'traces': dict(otype=np.ndarray, atype=np.floating, desc='Traces of the alignment frame'),
        'PYP_SPEC': dict(otype=str, desc='PypeIt spectrograph name'),
        'nalign': dict(otype=int, desc='Number of alignment traces in each slit'),
        'spat_id': dict(otype=np.ndarray, atype=np.integer, desc='Slit spat_id '),
    }

    def __init__(self, alignframe=None, traces=None, nalign=None, PYP_SPEC=None, spat_id=None):
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


    @classmethod
    def _parse(cls, hdu, ext=None, transpose_table_arrays=False, debug=False,
               hdu_prefix=None):
        # Grab everything but the bspline's
        _d, dm_version_passed, dm_type_passed = super(FlatImages, cls)._parse(hdu)
        # Now the bsplines
        list_of_bsplines = []
        spat_ids = []
        for ihdu in hdu:
            if 'BSPLINE' in ihdu.name:
                ibspl = bspline.bspline.from_hdu(ihdu)
                if ibspl.version != bspline.bspline.version:
                    msgs.warn("Your bspline is out of date!!")
                list_of_bsplines.append(ibspl)
                # Grab SPAT_ID for checking
                i0 = ihdu.name.find('ID-')
                i1 = ihdu.name.find('_BSP')
                spat_ids.append(int(ihdu.name[i0+3:i1]))
        # Check
        if spat_ids != _d['spat_id'].tolist():
            msgs.error("Bad parsing of the MasterFlat")
        # Finish
        _d['spat_bsplines'] = np.asarray(list_of_bsplines)
        return _d, dm_version_passed, dm_type_passed

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
        rawalignimg (:class:`pypeit.images.pypeitimage.PypeItImage` or None):
            Align image, created by the AlignFrame class
        slits (:class:`pypeit.slittrace.SlitTraceSet`, None):
            Slit edges
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph` or None):
            The `Spectrograph` instance that sets the
            instrument used to take the observations.  Used to set
            :attr:`spectrograph`.
        alignpar (:class:`pypeit.par.pypeitpar.AlignPar`):
            The parameters used for the align traces
        det (int, optional): Detector number
        qa_path (str, optional):  For QA
        msbpm (ndarray, optional): Bad pixel mask image

    Attributes:
        frametype : str
            Hard-coded to 'align_prof'
        steps : list
            List of the processing steps performed
    """
    # Frametype is a class attribute
    version = '1.0.0'
    frametype = 'alignment'
    master_type = 'Alignment'
    master_file_format = 'fits'

    def __init__(self, rawalignimg, slits, spectrograph, alignpar, det=1,
                 binning=None, qa_path=None, msbpm=None):

        # Defaults
        self.spectrograph = spectrograph
        self.PYP_SPEC = spectrograph.spectrograph
        self.binning = binning
        # Alignment parameters
        self.alignpar = alignpar

        # Input data
        self.rawalignimg = rawalignimg
        self.slits = slits

        # Check synced
        self.rawalignimg.is_synced(self.slits)

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

    def build_traces(self, show_peaks=False, show_trace=False, debug=False):
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
                self.rawalignimg.image, slitid_img_init == slit_idx,
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

        if show_trace:
            show_alignment(self.rawalignimg.image, align_traces=align_prof, slits=self.slits)

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
        profile_info : list
            A list of the alignment profile information. The elements of the list include:
            [spatial traces, TODO :: INCLUDE ALL ELEMENTS NEEDED HERE FOR THE DATA CONTAINER]
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

    def run(self, show_trace=False):
        """
        Main driver for alignment profile tracing

        Args:
            skip_QA (bool, optional):
                Skip the QA?

        Returns:
            dict, ndarray:  wv_calib dict and maskslits bool array

        """
        ###############
        # Fill up the calibrations and generate QA
        align_prof = self.build_traces(show_trace=show_trace)

        # Prepare the dictionary items for the data container
        align_traces = self.generate_traces(align_prof)

        return Alignments(# TODO :: FILL THIS IN)












def show_alignment(alignframe, align_traces=None, slits=None, clear=False):
    """
    Show one of the class internals

    Parameters
    ----------

    attr : str
        image - plot the master align frame
    image : ndarray
        Image to be plotted (i.e. the master align frame)
    align_traces : list
        The align traces
    chname : str
        The channel name sent to ginga
    slits : :class:`pypeit.slittrace.SlitTraceSet`, optional
        properties of the slits, including traces.
    clear : bool
        Clear the plotting window in ginga?

    Returns
    -------

    """
    ginga.connect_to_ginga(raise_err=True, allow_new=True)
    ch_name = 'alignment'
    viewer, channel = ginga.show_image(alignframe, chname=ch_name, clear=clear, wcs_match=False)

    # Display the slit edges
    if slits is not None and viewer is not None:
        left, right, mask = slits.select_edges()
        ginga.show_slits(viewer, channel, left, right)

    # Display the alignment traces
    if align_traces is not None and viewer is not None:
        for spec in align_traces:
            color = 'magenta' if spec.hand_extract_flag else 'orange'
            ginga.show_trace(viewer, channel, spec.TRACE_SPAT, trc_name="", color=color)
    return
