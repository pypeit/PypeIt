"""
Module for generating a Bar image to map constant spatial locations

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import os
import pdb
import copy
import inspect
import numpy as np
from IPython import embed

from pypeit import ginga, msgs, edgetrace
from pypeit import masterframe
from pypeit.par import pypeitpar
from pypeit.images import calibrationimage
from pypeit.images import pypeitimage
from pypeit.core import procimg, extract, arc, pixels, save, load


class BarFrame(calibrationimage.CalibrationImage, masterframe.MasterFrame):
    """
    Class to generate/load the bar image

    This class is primarily designed to generate a Bar frame to map constant
    spatial locations.  It also contains I/O methods for the Master frames
    of PypeIt.  The build_master() method will return a simple command
    (str) if that is the specified parameter (`par['useframe']`).

    Args:
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph used to take the data.

        files (:obj:`list`, optional):
            List of filenames to process.
        det (:obj:`int`, optional):
            The 1-indexed detector number to process.
        par (:class:`pypeit.par.pypeitpar.FrameGroupPar`, optional):
            The parameters used to process the frames.  If None, set
            to::

                pypeitpar.FrameGroupPar('bar')

        master_key (:obj:`str`, optional):
            The string identifier for the instrument configuration.  See
            :class:`pypeit.masterframe.MasterFrame`.
        master_dir (:obj:`str`, optional):
            Path to master frames
        reuse_masters (:obj:`bool`, optional):
            Load from disk if possible
    """

    # Frame type is a class attribute
    frametype = 'bar'
    master_type = 'Bar'

    @classmethod
    def from_master_file(cls, master_file, par=None):
        """
        Instantiate from a master file

        Args:
            master_file (str):
            par (:class:`pypeit.par.pypeitpar.FrameGroupPar`, optional):

        Returns:
            barframe.BarFrame:
                The PypeItImage is loaded into self.pypeitImage

        """
        # Spectrograph
        spectrograph, extras = masterframe.items_from_master_file(master_file)
        head0 = extras[0]
        # Master info
        master_dir = head0['MSTRDIR']
        master_key = head0['MSTRKEY']
        # Instantiate
        slf = cls(spectrograph, par=par, master_dir=master_dir, master_key=master_key,
                  reuse_masters=True)
        slf.pypeitImage = slf.load(ifile=master_file)
        # Return
        return slf

    # Keep order same as processimages (or else!)
    def __init__(self, spectrograph, files=None, det=1, par=None, master_key=None,
                 master_dir=None, reuse_masters=False, msbias=None):

        # Parameters unique to this Object
        self.msbias = msbias

        # Parameters
        self.par = pypeitpar.FrameGroupPar(self.frametype) if par is None else par

        # Start us up
        calibrationimage.CalibrationImage.__init__(self, spectrograph, det, self.par['process'], files=files)

        # MasterFrames: Specifically pass the ProcessImages-constructed
        # spectrograph even though it really only needs the string name
        masterframe.MasterFrame.__init__(self, self.master_type, master_dir=master_dir,
                                         master_key=master_key, reuse_masters=reuse_masters)
        # Process steps
        self.process_steps = procimg.init_process_steps(self.msbias, self.par['process'])
        self.process_steps += ['trim']
        self.process_steps += ['orient']
        self.process_steps += ['apply_gain']

    def save(self, outfile=None, overwrite=True):
        """
        Save the master bar data.

        Args:
            outfile (:obj:`str`, optional):
                Name for the output file.  Defaults to
                :attr:`file_path`.
            overwrite (:obj:`bool`, optional):
                Overwrite any existing file.
        """
        _outfile = self.master_file_path if outfile is None else outfile
        # Check if it exists
        if os.path.exists(_outfile) and not overwrite:
            msgs.warn('Master file exists: {0}'.format(_outfile) + msgs.newline()
                      + 'Set overwrite=True to overwrite it.')
            return
        #
        hdr = self.build_master_header(steps=self.process_steps, raw_files=self.file_list)
        self.pypeitImage.write(_outfile, hdr=hdr, iext='BAR')
        msgs.info('Master frame written to {0}'.format(_outfile))

    def load(self, ifile=None):
        """
        Load the bar frame data from a saved master frame.

        Args:
            ifile (:obj:`str`, optional):
                Name of the master frame file.  Defaults to
                :attr:`file_path`.
            return_header (:obj:`bool`, optional):
                Return the header

        Returns:
            Returns a `numpy.ndarray`_ with the bar master frame image.
            Also returns the primary header, if requested.
        """
        # Check on whether to reuse and whether the file exists
        master_file = self.chk_load_master(ifile)
        if master_file is None:
            return
        # Load it up
        self.pypeitImage = pypeitimage.PypeItImage.from_file(master_file)
        return self.pypeitImage


class BarProfile(masterframe.MasterFrame):
    """
    Class to guide the determination of the bar traces

    Args:
        msbar (:class:`pypeit.images.pypeitimage.PypeItImage` or None):
            Bar image, created by the BarFrame class
        tslits_dict (dict or None):  TraceSlits dict
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph` or None):
            The `Spectrograph` instance that sets the
            instrument used to take the observations.  Used to set
            :attr:`spectrograph`.
        par (:class:`pypeit.par.pypeitpar.BarPar` or None):
            The parameters used for the bar traces
        det (int, optional): Detector number
        master_key (str, optional)
        master_dir (str, optional): Path to master frames
        reuse_masters (bool, optional):  Load from disk if possible
        qa_path (str, optional):  For QA
        msbpm (ndarray, optional): Bad pixel mask image

    Attributes:
        frametype : str
            Hard-coded to 'bar_prof'
        steps : list
            List of the processing steps performed
    """
    # Frametype is a class attribute
    frametype = 'bar_prof'
    master_type = 'BarProfile'

    def __init__(self, msbar, tslits_dict, spectrograph, par, det=1,
                 binning=None, master_key=None, master_dir=None, reuse_masters=False,
                 qa_path=None, msbpm=None):

        # MasterFrame
        masterframe.MasterFrame.__init__(self, self.master_type, master_dir=master_dir,
                                         master_key=master_key, file_format='fits',
                                         reuse_masters=reuse_masters)

        # Required parameters (but can be None)
        self.msbar = msbar
        self.tslits_dict = tslits_dict
        self.spectrograph = spectrograph
        self.par = par
        self.binning = binning

        # Optional parameters
        self.bpm = msbar.mask if msbpm is None else msbpm
        self.qa_path = qa_path
        self.det = det
        self.master_key = master_key

        # Attributes
        self.viewer, self.channel = None, None
        self.steps = []  # steps executed

        # Get the non-linear count level
        self.nonlinear_counts = 1e10 if self.spectrograph is None \
            else self.spectrograph.nonlinear_counts(self.det)

        # --------------------------------------------------------------
        # Set the slitmask and slit boundary related attributes that the
        # code needs for execution. This also deals with barframes that
        # have a different binning then the images used to defined
        # the slits
        if self.tslits_dict is not None and self.msbar is not None:
            self.slitmask_science = pixels.tslits2mask(self.tslits_dict)
            gpm = self.bpm == 0 if self.bpm is not None \
                else np.ones_like(self.slitmask_science, dtype=bool)
            self.shape_science = self.slitmask_science.shape
            self.shape_bar = self.msbar.image.shape
            self.nslits = self.tslits_dict['slit_left'].shape[1]
            self.slit_left = arc.resize_slits2arc(self.shape_bar, self.shape_science,
                                                  self.tslits_dict['slit_left'])
            self.slit_righ = arc.resize_slits2arc(self.shape_bar, self.shape_science,
                                                  self.tslits_dict['slit_righ'])
            self.slitcen = arc.resize_slits2arc(self.shape_bar, self.shape_science,
                                                self.tslits_dict['slitcen'])
            self.slitmask = arc.resize_mask2arc(self.shape_bar, self.slitmask_science)
            self.gpm = arc.resize_mask2arc(self.shape_bar, gpm)
            self.gpm &= self.msbar.image < self.nonlinear_counts
            self.slit_spat_pos = edgetrace.slit_spat_pos(self.tslits_dict['slit_left'],
                                                         self.tslits_dict['slit_righ'],
                                                         self.tslits_dict['nspec'],
                                                         self.tslits_dict['nspat'])
        else:
            self.slitmask_science = None
            self.shape_science = None
            self.shape_bar = None
            self.nslits = 0
            self.slit_left = None
            self.slit_righ = None
            self.slitcen = None
            self.slitmask = None
            self.gpm = None
            self.nonlinear_counts = None

    def build_traces(self, show_peaks=False, show_trace=False, debug=False):
        """
        Main routine to generate the bar profile traces in all slits

        Args:
             show_peaks (bool, optional):
               Generate QA showing peaks identified by object finding
             show_trace (bool, optional):
               Generate QA  showing traces identified. Requires an open ginga RC modules window
             debug (bool, optional):

        Returns:
            dict:  self.bar_dict
        """
        bar_prof = dict({})
        nslits = self.tslits_dict['slit_left'].shape[1]
        # Prepare the plotting canvas
        if show_trace:
            self.show('image', image=self.msbar.image, chname='bar_traces', slits=True)
        # Go through the slits
        for sl in range(nslits):
            specobj_dict = {'setup': "unknown", 'slitid': sl,
                            'det': self.det, 'objtype': "bar_profile", 'pypeline': "MultiSlit"}
            msgs.info("Fitting bar traces in slit {0:d}".format(sl))
            bar_traces, _ = extract.objfind(
                self.msbar.image, self.slitmask == sl,
                self.tslits_dict['slit_left'][:, sl],
                self.tslits_dict['slit_righ'][:, sl],
                ir_redux=False, ncoeff=self.par['trace_npoly'],
                specobj_dict=specobj_dict, sig_thresh=self.par['sig_thresh'],
                show_peaks=show_peaks, show_fits=False,
                trim_edg=self.par['trim_edge'],
                cont_fit=False, npoly_cont=0,
                nperslit=len(self.par['locations']))
            if len(bar_traces) != len(self.par['locations']):
                # Bar tracing has failed for this slit
                msgs.warn("Bar tracing has failed on slit {0:d}".format(sl))
            if show_trace:
                self.show('overplot', chname='bar_traces', bar_traces=bar_traces, slits=False)
            bar_prof['{0:d}'.format(sl)] = bar_traces.copy()

        bar_dict = self.generate_dict(bar_prof)

        # Steps
        self.steps.append(inspect.stack()[0][3])

        # Return
        return bar_dict

    def generate_dict(self, bar_prof):
        """
        Generate a dictionary containing all of the information from the profile fitting

        Args:
            bar_prof (:obj:`dict`):
                Dictionary of SpecObjs classes (one for each slit)

        Returns:
            dict:  bar_dict
        """
        bar_dict = dict({})
        nbars = len(self.par['locations'])
        nspec, nslits = self.tslits_dict['slit_left'].shape
        # Generate an array containing the centroid of all bars
        barprof = np.zeros((nspec, nbars, nslits))
        for sl in range(nslits):
            sls = '{0:d}'.format(sl)
            for bar in range(nbars):
                if bar_prof[sls][bar].SLITID != sl:
                    msgs.error("Bar profiling failed to generate dictionary")
                barprof[:, bar, sl] = bar_prof[sls][bar].TRACE_SPAT
        # Put all of the info into a single dictionary
        bar_dict = dict(bar_profiles=barprof)
        return bar_dict

    def save(self, outfile=None, overwrite=True):
        """
        Save the bar traces to a master frame.

        Args:
            outfile (:obj:`str`, optional):
                Name of the output file.  Defaults to :attr:`master_file_path`.
            overwrite (:obj:`bool`, optional):
                Overwrite any existing file.
        """
        _outfile = self.master_file_path if outfile is None else outfile
        # Check if it exists
        if os.path.exists(_outfile) and not overwrite:
            msgs.warn('Master file exists: {0}'.format(_outfile) + msgs.newline()
                      + 'Set overwrite=True to overwrite it.')
            return

        # Report and save
        prihdr = self.build_master_header(steps=self.steps)
        #   - Add the binning
        prihdr['BINNING'] = (self.binning, 'PypeIt: Binning')
        #   - Add the detector number
        prihdr['DET'] = (self.det, 'PypeIt: Detector')
        #   - Add the tracing parameters
        self.par.to_header(prihdr)

        # Set the data and extension names
        data = [self.bar_dict['bar_profiles']]
        extnames = ['BAR_PROFILES']
        # Write the output to a fits file
        save.write_fits(prihdr, data, _outfile, extnames=extnames)
        msgs.info('Master frame written to {0}'.format(_outfile))

    def load(self, ifile=None):
        """
        Load the profiles of the bar frame.

        Args:
            ifile (:obj:`str`, optional):
                Name of the master frame file.  Defaults to
                :attr:`master_file_path`.

        Returns:
            dict or None: self.bar_dict
        """
        # Check on whether to reuse and whether the file exists
        master_file = self.chk_load_master(ifile)
        if master_file is None:
            return
        msgs.info('Loading Master frame: {0}'.format(master_file))
        # Load
        extnames = ['BAR_PROFILES']
        *data, head0 = load.load_multiext_fits(master_file, extnames)

        # Fill the dict
        self.bar_dict = {}
        keys = []
        for k in keys:
            self.bar_dict[k] = head0[k.upper()]
        # Data
        for ii, ext in enumerate(extnames):
            self.bar_dict[ext.lower()] = data[ii]
        # Return
        return self.bar_dict

    def run(self, show_trace=False, skip_QA=False, debug=False):
        """
        Main driver for bar profile tracing

        Args:
            skip_QA : bool, optional

        Returns:
            dict, ndarray:  wv_calib dict and maskslits bool array

        """
        ###############
        # Fill up the calibrations and generate QA
        self.bar_dict = self.build_traces(show_trace=show_trace)

        # Pack up
        self.bar_dict['steps'] = self.steps
        sv_par = self.par.data.copy()
        self.bar_dict['par'] = sv_par

        return self.bar_dict

    def show(self, attr, image=None, bar_traces=None,
             chname=None, slits=False, clear=False):
        """
        Show one of the class internals

        Args:
            attr : str
                image - plot the master bar frame
            image : ndarray
                Image to be plotted (i.e. the master bar frame)
            bar_traces : list
                The bar traces
            chname : str
                The channel name sent to ginga
            slits : bool
                Overplot the slit edges?
            clear : bool
                Clear the plotting window in ginga?

        Returns:

        """
        if attr == 'image':
            ch_name = chname if chname is not None else 'bar_traces'
            self.viewer, self.channel = ginga.show_image(image, chname=ch_name, clear=clear, wcs_match=False)
        elif attr == 'overplot':
            pass
        else:
            msgs.warn("Not an option for show")

        if bar_traces is not None and self.viewer is not None:
            for spec in bar_traces:
                color = 'magenta' if spec.hand_extract_flag else 'orange'
                ginga.show_trace(self.viewer, self.channel, spec.TRACE_SPAT, trc_name="", color=color)

        if slits:
            if self.tslits_dict is not None and self.viewer is not None:
                slit_ids = [edgetrace.get_slitid(image.shape,
                                                 self.tslits_dict['slit_left'],
                                                 self.tslits_dict['slit_righ'], ii)[0]
                            for ii in range(self.tslits_dict['slit_left'].shape[1])]
                ginga.show_slits(self.viewer, self.channel, self.tslits_dict['slit_left'],
                                 self.tslits_dict['slit_righ'], slit_ids)
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
