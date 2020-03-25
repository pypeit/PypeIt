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
from pypeit import masterframe
from pypeit.par import pypeitpar
from pypeit.images import pypeitimage
from pypeit.core import procimg, extract, arc, pixels, load

from pypeit import io
from pypeit.core import extract, arc, save, load


class Alignment():
    """
    Class to guide the determination of the alignment traces

    Args:
        msalign (:class:`pypeit.images.pypeitimage.PypeItImage` or None):
            Align image, created by the AlignFrame class
        slits (:class:`pypeit.slittrace.SlitTraceSet`, None):
            Slit edges
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph` or None):
            The `Spectrograph` instance that sets the
            instrument used to take the observations.  Used to set
            :attr:`spectrograph`.
        par (:class:`pypeit.par.pypeitpar.AlignPar` or None):
            The parameters used for the align traces
        det (int, optional): Detector number
        master_key (str, optional)
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
    file_format = 'fits'

    def __init__(self, msalign, slits, spectrograph, par, det=1,
                 binning=None, master_key=None,
                 qa_path=None, msbpm=None):

        # Required parameters (but can be None)
        self.msalign = msalign
        self.slits = slits
        self.spectrograph = spectrograph
        self.PYP_SPEC = spectrograph.spectrograph
        self.par = par
        self.master_key = master_key
        self.binning = binning

        # Optional parameters
        self.bpm = msalign.mask if msbpm is None else msbpm
        self.qa_path = qa_path
        self.det = det

        # Attributes
        self.viewer, self.channel = None, None
        self.steps = []  # steps executed

        # Get the non-linear count level
        self.nonlinear_counts = 1e10 if self.spectrograph is None \
            else self.spectrograph.nonlinear_counts(self.msalign.detector)

        # --------------------------------------------------------------
        # Set the slitmask and slit boundary related attributes that the
        # code needs for execution. This also deals with alignframes that
        # have a different binning then the images used to defined
        # the slits
        if self.slits is not None and self.msalign is not None:
            # NOTE: This uses the internal definition of `pad`
            self.slitmask_science = self.slits.slit_img()
            gpm = (self.bpm == 0) if self.bpm is not None \
                else np.ones_like(self.slitmask_science, dtype=bool)
            self.shape_science = self.slitmask_science.shape
            self.shape_align = self.msalign.image.shape
            self.nslits = self.slits.nslits
            self.slit_left, self.slit_right = self.slits.select_edges()
            self.slitcen = 0.5 * (self.slit_left + self.slit_right)
            self.slitmask = arc.resize_mask2arc(self.shape_align, self.slitmask_science)
            self.gpm = (arc.resize_mask2arc(self.shape_align, gpm)) & (self.msalign.image < self.nonlinear_counts)
        else:
            self.slitmask_science = None
            self.shape_science = None
            self.shape_align = None
            self.nslits = 0
            self.slit_left = None
            self.slit_right = None
            self.slitcen = None
            self.slitmask = None
            self.gpm = None

    def build_traces(self, show_peaks=False, show_trace=False, debug=False):
        """
        Main routine to generate the align profile traces in all slits

        Args:
             show_peaks (bool, optional):
               Generate QA showing peaks identified by object finding
             show_trace (bool, optional):
               Generate QA  showing traces identified. Requires an open ginga RC modules window
             debug (bool, optional):

        Returns:
            dict:  self.align_dict
        """
        align_prof = dict({})
        # Prepare the plotting canvas
        if show_trace:
            self.show('image', image=self.msalign.image, chname='align_traces', slits=True)
        # Go through the slits
        for sl in range(self.nslits):
            specobj_dict = {'SLITID': sl, 'DET': self.det, 'OBJTYPE': "align_profile",
                            'PYPELINE': self.spectrograph.pypeline}
            msgs.info("Fitting alignment traces in slit {0:d}".format(sl))
            align_traces, _ = extract.objfind(
                self.msalign.image, self.slitmask == sl,
                self.slit_left[:, sl], self.slit_right[:, sl],
                ir_redux=False, ncoeff=self.par['trace_npoly'],
                specobj_dict=specobj_dict, sig_thresh=self.par['sig_thresh'],
                show_peaks=show_peaks, show_fits=False,
                trim_edg=self.par['trim_edge'],
                cont_fit=False, npoly_cont=0,
                nperslit=len(self.par['locations']))
            if len(align_traces) != len(self.par['locations']):
                # Align tracing has failed for this slit
                msgs.warn("Alignment tracing has failed on slit {0:d}".format(sl))
            if show_trace:
                self.show('overplot', chname='align_traces', align_traces=align_traces, slits=False)
            align_prof['{0:d}'.format(sl)] = align_traces.copy()

        align_dict = self.generate_dict(align_prof)

        # Steps
        self.steps.append(inspect.stack()[0][3])

        # Return
        return align_dict

    def generate_dict(self, align_prof):
        """
        Generate a dictionary containing all of the information from the profile fitting

        Args:
            align_prof (:obj:`dict`):
                Dictionary of SpecObjs classes (one for each slit)

        Returns:
            dict:  align_dict
        """
        nbars = len(self.par['locations'])
        nspec, nslits = self.slit_left.shape
        # Generate an array containing the centroid of all bars
        alignprof = np.zeros((nspec, nbars, nslits))
        for sl in range(nslits):
            sls = '{0:d}'.format(sl)
            for bar in range(nbars):
                if align_prof[sls][bar].SLITID != sl:
                    msgs.error("Alignment profiling failed to generate dictionary")
                    alignprof[:, bar, sl] = align_prof[sls][bar].TRACE_SPAT
        # Return the profile information as a single dictionary
        return dict(alignments=alignprof)

    def save(self, master_key, master_dir, outfile=None, overwrite=True):
        """
        Save the alignment traces to a master frame.

        Args:
            outfile (:obj:`str`, optional):
                Name of the output file.  Defaults to :attr:`master_file_path`.
            overwrite (:obj:`bool`, optional):
                Overwrite any existing file.
        """
        _outfile = outfile
        # Check if it exists
        if os.path.exists(_outfile) and not overwrite:
            msgs.warn('Master file exists: {0}'.format(_outfile) + msgs.newline()
                      + 'Set overwrite=True to overwrite it.')
            return

        # Report and save
        # First do the header
        if master_key is not None:
            prihdr = masterframe.build_master_header(self, master_key, master_dir)
        else:
            prihdr = io.initialize_header()
        #   - Add the binning
        prihdr['BINNING'] = (self.binning, 'PypeIt: Binning')
        #   - Add the detector number
        prihdr['DET'] = (self.det, 'PypeIt: Detector')
        #   - Add the tracing parameters
        self.par.to_header(prihdr)

        # Set the data and extension names
        data = [self.align_dict['alignments']]
        extnames = ['ALIGNMENTS']
        # Write the output to a fits file
        save.write_fits(prihdr, data, _outfile, extnames=extnames)
        msgs.info('Master frame written to {0}'.format(_outfile))

    def load(self, ifile=None):
        """
        Load the profiles of the align frame.

        Args:
            ifile (:obj:`str`, optional):
                Name of the master frame file.  Defaults to
                :attr:`master_file_path`.

        Returns:
            dict or None: self.align_dict
        """
        msgs.info('Loading Master frame: {0}'.format(ifile))
        # Load
        extnames = ['ALIGNMENTS']
        *data, head0 = load.load_multiext_fits(ifile, extnames)

        # Fill the dict
        self.align_dict = {}
        keys = []
        for k in keys:
            self.align_dict[k] = head0[k.upper()]
        # Data
        for ii, ext in enumerate(extnames):
            self.align_dict[ext.lower()] = data[ii]
        # Return
        return self.align_dict

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
        self.align_dict = self.build_traces(show_trace=show_trace)

        # Pack up
        self.align_dict['steps'] = self.steps
        sv_par = self.par.data.copy()
        self.align_dict['par'] = sv_par

        return self.align_dict

    def show(self, attr, image=None, align_traces=None,
             chname=None, slits=False, clear=False):
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
        slits : bool
            Overplot the slit edges?
        clear : bool
            Clear the plotting window in ginga?

        Returns:

        """
        if attr == 'image':
            ch_name = chname if chname is not None else 'align_traces'
            self.viewer, self.channel = ginga.show_image(image, chname=ch_name, clear=clear, wcs_match=False)
        elif attr == 'overplot':
            pass
        else:
            msgs.warn("Not an option for show")

        if align_traces is not None and self.viewer is not None:
            for spec in align_traces:
                color = 'magenta' if spec.hand_extract_flag else 'orange'
                ginga.show_trace(self.viewer, self.channel, spec.TRACE_SPAT, trc_name="", color=color)

        if slits:
            if self.slits is not None and self.viewer is not None:
                ginga.show_slits(self.viewer, self.channel, self.slit_left, self.slit_right)
        return

    def __repr__(self):
        # Generate sets string
        txt = '<{:s}: '.format(self.__class__.__name__)
        if len(self.steps) > 0:
            txt += ' steps: ['
            for step in self.steps:
                txt += '{:s}, '.format(step)
            txt = txt[:-2] + ']'  # Trim the trailing comma
        txt += '>'
        return txt
