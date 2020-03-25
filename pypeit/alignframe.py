"""
Module for generating an Alignment image to map constant spatial locations

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import os
import inspect
import numpy as np
from IPython import embed

import astropy.io.fits as fits

from pypeit import masterframe, io
from pypeit import ginga, msgs
from pypeit.core import extract, arc, load
from pypeit.par.pypeitpar import AlignPar


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
        self._alignprof = None
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
        """
        nbars = len(self.par['locations'])
        nspec, nslits = self.slit_left.shape
        # Generate an array containing the centroid of all bars
        self._alignprof = np.zeros((nspec, nbars, nslits))
        for sl in range(nslits):
            sls = '{0:d}'.format(sl)
            for bar in range(nbars):
                if align_prof[sls][bar].SLITID != sl:
                    msgs.error("Alignment profiling failed to generate dictionary")
                    self._alignprof[:, bar, sl] = align_prof[sls][bar].TRACE_SPAT
        # Return the profile information as a single dictionary
        return

    def save(self, outfile, overwrite=True, checksum=True, float_dtype='float32',
             master_dir=None, master_key=None):
        """
        Save the alignment traces to a master frame.

        Args:
            outfile (:obj:`str`):
                Name for the output file.  Defaults to
                :attr:`master_file_path`.
            overwrite (:obj:`bool`, optional):
                Overwrite any existing file.
            checksum (:obj:`bool`, optional):
                Passed to `astropy.io.fits.HDUList.writeto`_ to add
                the DATASUM and CHECKSUM keywords fits header(s).
            float_dtype (:obj:`str`, optional):
                Convert floating-point data to this data type before
                writing.  Default is 32-bit precision.
        """

        # Check if it exists
        if os.path.exists(outfile) and not overwrite:
            msgs.error('Master file exists: {0}'.format(outfile) + msgs.newline()
                       + 'Set overwrite=True to overwrite it.')

        # Report name before starting to write it
        msgs.info('Writing master frame to {0}'.format(outfile))

        # Build the primary header
        #   - Initialize with basic metadata
        # Header
        if master_key is not None:
            prihdr = masterframe.build_master_header(self, master_key, master_dir)
        else:
            prihdr = io.initialize_header()
        #   - Add the binning
        prihdr['BINNING'] = (self.binning, 'PypeIt: Binning')
        #   - Add the detector number
        prihdr['DET'] = (self.det, 'PypeIt: Detector')
        #   - Add the alignment parameters
        self.par.to_header(prihdr)

        # Determine if the file should be compressed
        compress = False
        if outfile.split('.')[-1] == 'gz':
            outfile = outfile[:outfile.rfind('.')]
            compress = True

        # First check if a trace is available
        if self._alignprof is None:
            msgs.error("No alignment information available")
        # Build the list of extensions to write
        hdu = fits.HDUList([fits.PrimaryHDU(header=prihdr),
                            fits.ImageHDU(data=self._alignprof.astype(float_dtype), name='ALIGNMENTS')])
        # Add detector
        hdu += self.msalign.detector.to_hdu()

        # Write the fits file; note not everything is written.
        hdu.writeto(outfile, overwrite=True, checksum=checksum)

        # Compress the file if the output filename has a '.gz'
        # extension; this is slow but still faster than if you have
        # astropy.io.fits do it directly
        if compress:
            msgs.info('Compressing file to: {0}.gz'.format(outfile))
            io.compress_file(outfile, overwrite=overwrite)

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

        # Set the data
        self._alignprof = data[0]

        # Check the parset agrees with the one on file
        par = AlignPar.from_header(head0)
        if self.par.data != par.data:
            # TODO: The above inequality works for non-nested ParSets,
            # but will need to be more careful for nested ones, or just
            # avoid writing nested ParSets to headers...
            msgs.error('This master frame was generated using different parameter values!')
        return self._alignprof

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
        self.build_traces(show_trace=show_trace)

        return self._alignprof

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
