from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import time
from abc import ABCMeta
import os
import datetime
import numpy as np
from collections import OrderedDict

from astropy.time import Time
from pypeit import msgs
from pypeit import pypeitsetup
from pypeit import calibrations
from pypeit import scienceimage
from pypeit import specobjs
from pypeit import fluxspec
from pypeit.core import paths
from pypeit.core import qa
from pypeit.core import pypsetup
from pypeit.core import wave
from pypeit.core import save
from pypeit.core import load
from pypeit.spectrographs.util import load_spectrograph
from pypeit.scripts import run_pypeit
from pypeit.par.util import make_pypeit_file, parse_pypeit_file

from pypeit import debugger

class PypeIt(object):
    """
    This class is designed to run PypeIt

    .. todo::
        Improve docstring...

    Args:
        spectrograph (:obj:`str`,
            :class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            The string or `Spectrograph` instance that sets the
            instrument used to take the observations.  Used to set
            :attr:`spectrograph`.
        setups_path (:obj:`str`, optional):
            Path for files related to all setups
        verbosity (:obj:`int`, optional):
            Verbosity level of system output.  Can be::
                - 0: No output
                - 1: Minimal output (default)
                - 2: All output
        overwrite (:obj:`bool`, optional):
            Flag to overwrite any existing files/directories.
        logname (:obj:`str`, optional):
            The name of an ascii log file with the details of the
            reduction.
        show: (:obj:`bool`, optional):
            Show reduction steps via plots (which will block further
            execution until clicked on) and outputs to ginga. Requires
            remote control ginga session via "ginga --modules=RC &"

    Attributes:
        pypeit_file (:obj:`str`):
            Name of the pypeit file to read.  Pypit files have a specific
            set of valid formats. A description can be found `here`_
            (include doc link).
    """
    __metaclass__ = ABCMeta

    def __init__(self, spectrograph, setups_path=None, verbosity=2,
                 overwrite=True, logname=None, show=False):

        # Spectrometer class
        self.spectrograph = load_spectrograph(spectrograph)

        # Init
        self.verbosity = verbosity
        self.overwrite = overwrite
        self.setups_path = os.getcwd() if setups_path is None else setups_path

        # Internals
        self.pypeit_file = None
        self.logname = logname
        self.setup_pypeit_file = None
        self.redux_path = None
        self.show=show
        self.setup = None
        self.det = None
        self.sci_ID = None

        self.tstart = None
        self.fitstbl = None
        self.par = None
        self.caliBrate = None
        self.basename = None
        self.sciI = None
        self.obstime = None
        self.pypeitSetup = None
        self.setup_dict = None


    def build_setup_files(self, files_root, extension='.fits'):
        """
        Generate the setup files for PypeIt from a list of input files

        Args:
            files_root: str
              Root name for the files to be reduced including full path

        Returns:

        """

        # Record the starting time
        self.tstart = time.time()

        pargs, sort_dir, self.setup_pypeit_file \
                = self._make_setup_pypeit_file(files_root, extension=extension)
        self._setup(self.setup_pypeit_file, setup_only=True, calibration_check=False,
                    sort_dir=sort_dir)

        self.print_end_time()

    def build_custom_pypeitfiles(self):
        """
        Build the custom PypeIt files, one per unique instrument configuration
        Each is put in a custom folder parallel to the folder of setups files

        Returns:

        """

        msgs.reset(verbosity=2)

        # Read master file
        cfg_lines, data_files, frametype, usrdata, setups \
                = parse_pypeit_file(self.setup_pypeit_file)
        sorted_file = os.path.splitext(self.setup_pypeit_file)[0]+'.sorted'
        sub_cfg_lines = cfg_lines[0:2]

        # Get paths
        paths = []
        for data_file in data_files:
            islsh = data_file.rfind('/')
            path = data_file[:islsh+1]
            if path not in paths:
                paths.append(path)

        # Generate .pypeit files and sub-folders
        all_setups, all_setuplines, all_setupfiles = pypsetup.load_sorted(sorted_file)
        for setup, setup_lines, sorted_files in zip(all_setups, all_setuplines, all_setupfiles):
            root = self.spectrograph.spectrograph+'_setup_'
            # cfg_lines
            cfg_lines = sub_cfg_lines
            cfg_lines += ['    sortroot = {0}'.format(root + setup)]
            # Make the dir
            newdir = os.path.join(self.setups_path, root+setup)
            if not os.path.exists(newdir):
                os.mkdir(newdir)
            # Now the file
            pypeit_file = os.path.join(newdir, root+setup+'.pypeit')
            # Modify parlines
            for kk in range(len(cfg_lines)):
                if 'sortroot' in cfg_lines[kk]:
                    cfg_lines[kk] = '    sortroot = {0}'.format(root+setup)

            make_pypeit_file(pypeit_file, self.spectrograph.spectrograph, [], cfg_lines=cfg_lines,
                             setup_lines=setup_lines, sorted_files=sorted_files, paths=paths)
            print("Wrote {:s}".format(pypeit_file))

    def build_qa(self):
        """
        Generate QA wrappers

        Returns:

        """
        qa.gen_mf_html(self.pypeit_file)
        qa.gen_exp_html()


    def calibrate_one(self, sci_ID, det):
        """
        Calibrate a science exposure / detector pair

        Args:
            sci_ID: int
              binary flag indicating the science frame
            det: int
              detector number

        Returns:

        """
        # Setup
        self.setup, self.setup_dict = pypsetup.instr_setup(sci_ID, det, self.fitstbl,
                                                           setup_dict=self.setup_dict,
                                                           must_exist=True)
        # Setup
        self.caliBrate.reset(self.setup, det, sci_ID, self.par['calibrations'])
        # Run em
        self.caliBrate.run_the_steps()

        msgs.info("Successful Calibration!")


    def _chk_for_std(self):
        # Can only reduce these frames if the mask is the same
        std_idx = self.fitstbl.find_frames('standard', sci_ID=self.sci_ID, index=True)
        if len(std_idx) > 0:
            std_idx = std_idx[0]
            return std_idx
        else:
            msgs.info("No standard star associated with this science frame")
            return -1

    def reduce_all(self, reuse_masters=False):
        """
        Reduce all of the science exposures
        Generate all needed calibration files

        Args:
            reuse_masters (:obj:`bool`, optional):
                Use the master frames if available (same as setting
                par['calibrations']['masters'] = 'reuse'.

        Returns:

        """

        self.tstart = time.time()
        self.std_dict = {}
        # Science IDs are in a binary system: 1,2,4,8, etc.
        all_sci_ID = self.fitstbl['sci_ID'][self.fitstbl.find_frames('science')]
        numsci = len(all_sci_ID)
        basenames = [None]*numsci  # For fluxing at the very end

        # Check par
        required = ['rdx', 'calibrations', 'scienceframe', 'scienceimage', 'flexure', 'fluxcalib']
        can_be_None = ['flexure', 'fluxcalib']
        self.par.validate_keys(required=required, can_be_None=can_be_None)

        # Save
        for kk,sci_ID in enumerate(all_sci_ID):
            sci_dict = self.reduce_exposure(sci_ID, reuse_masters=reuse_masters)
            scidx = self.fitstbl.find_frames('science', sci_ID=sci_ID, index=True)[0]
            self.save_exposure(scidx, sci_dict, self.basename)
            basenames[kk] = self.basename

        # Standard stars
        for std_idx in self.std_dict.keys():
            # Basename
            ikey = list(self.std_dict[std_idx].keys())[0]  # Any will do, so take the first
            std_spec_objs = self.save_exposure(std_idx, self.std_dict[std_idx],
                                               self.std_dict[std_idx][ikey]['basename'])

        # Flux?
        if self.par['fluxcalib'] is None or len(self.std_dict) == 0:
            msgs.info('Flux calibration is not performed.')
        elif self.par['fluxcalib'] is None and len(self.std_dict) > 0:
            msgs.info('Flux calibration parameters not provided.  Standards not used.')
        else:
            # Standard star (is this a calibration, e.g. goes above?)
            msgs.info("Taking one star per detector mosaic")
            msgs.info("Waited until very end to work on it")
            msgs.warn("You should probably consider using the pypeit_flux_spec script anyhow...")

            # Get the sensitivity function
            if self.par['fluxcalib']['sensfunc'] is None:
                # Take the first standard
                std_idx = list(self.std_dict.keys())[0]
                # Build the list of stdobjs
                #all_std_objs = []
                #for det in self.std_dict[std_idx].keys():
                #    all_std_objs += self.std_dict[std_idx][det]['specobjs']
                # Need the Header for RA/DEC
                std_header = {}
                for key in ['ra', 'dec', 'airmass', 'exptime']:
                    std_header[key.upper()] = self.fitstbl[std_idx][key]
                # Go
                FxSpec = fluxspec.FluxSpec(std_specobjs=std_spec_objs.specobjs, spectrograph=self.spectrograph,
                                           setup=self.setup, master_dir=self.caliBrate.master_dir, std_header=std_header, mode=self.par['calibrations']['masters'])
                sens_dict = FxSpec.get_sens_dict(self.fitstbl[std_idx])
            else:
                # User provided it
                FxSpec = fluxspec.FluxSpec(sens_file=self.par['fluxcalib']['sensfunc'],
                                           spectrograph=self.spectrograph, master_dir=self.caliBrate.master_dir,
                                           mode=self.par['calibrations']['masters'])
                sens_dict = FxSpec.sens_dict

            # Apply the flux calibration
            msgs.info("Fluxing with {:s}".format(sens_dict['std_name']))
            save_format = 'fits'
            for kk, sci_ID in enumerate(all_sci_ID):
                # Load from disk (we zero'd out the object to free memory)
                if save_format == 'fits':
                    sci_spec1d_file = os.path.join(self.par['rdx']['scidir'],
                                                   'spec1d_{:s}.fits'.format(basenames[kk]))

                # Load
                sci_specobjs, sci_header = load.load_specobj(sci_spec1d_file)
                # TODO: (KBW) I'm wary of this kind of approach.  We want
                # FluxSpec to check that its internals make sense and this
                # bypasses any of that checking.
                FxSpec.sci_specobjs = sci_specobjs
                FxSpec.sci_header = sci_header

                # Flux
                FxSpec.flux_science()
                # Over-write
                FxSpec.write_science(sci_spec1d_file)

        # Finish
        self.print_end_time()

    def reduce_exposure(self, sci_ID, reuse_masters=False):
        """
        Reduce a single science exposure

        Args:
            sci_ID: int
              binary flag indicating the science frame
            reuse_masters: bool, optional
              Reuse MasterFrame files (where available)


        Returns:
            sci_dict: dict
              dict containing the primary outputs of extraction

        """
        self.sci_ID = sci_ID

        # Insist on re-using MasterFrames where applicable
        if reuse_masters:
            self.par['calibrations']['masters'] = 'reuse'

        sci_dict = OrderedDict()  # This needs to be ordered
        sci_dict['meta'] = {}
        sci_dict['meta']['vel_corr'] = 0.
        #
        scidx = self.fitstbl.find_frames('science', sci_ID=sci_ID, index=True)[0]
        msgs.info("Reducing file {0:s}, target {1:s}".format(self.fitstbl['filename'][scidx],
                                                             self.fitstbl['target'][scidx]))

        # Loop on Detectors
        for kk in range(self.spectrograph.ndet):
            det = kk + 1  # Detectors indexed from 1
            self.det = det
            if self.par['rdx']['detnum'] is not None:
                detnum = [self.par['rdx']['detnum']] if isinstance(self.par['rdx']['detnum'],int) else self.par['rdx']['detnum']
                if det not in map(int, detnum):
                    msgs.warn("Skipping detector {:d}".format(det))
                    continue
                else:
                    msgs.warn("Restricting the reduction to detector {:d}".format(det))
            # Setup
            msgs.info("Working on detector {0}".format(det))
            sci_dict[det] = {}

            # Calibrate
            self.calibrate_one(sci_ID, det)

            # Init ScienceImage class
            self.init_one_science(sci_ID, det)
            # Extract
            sciimg, sciivar, skymodel, objmodel, ivarmodel, outmask, sobjs, vel_corr = self._extract_one()

            # Save for outputing (after all detectors are done)
            sci_dict[det]['sciimg'] = sciimg
            sci_dict[det]['sciivar'] = sciivar
            sci_dict[det]['skymodel'] = skymodel
            sci_dict[det]['objmodel'] = objmodel
            sci_dict[det]['ivarmodel'] = ivarmodel
            sci_dict[det]['outmask'] = outmask
            sci_dict[det]['specobjs'] = sobjs   #utils.unravel_specobjs([specobjs])
            if vel_corr is not None:
                sci_dict['meta']['vel_corr'] = vel_corr

            # Standard star
            # TODO -- Make this more modular
            self.std_idx = self._chk_for_std()
            if self.std_idx is not -1:
                self._extract_std()

        # Return
        return sci_dict

    def save_exposure(self, sidx, s_dict, basename, only_1d=False):
        """
        Save the outputs from extraction for a given exposure

        Args:
            sidx: int
              Index of the exposure to save
            s_dict: dict
              dict containing the primary outputs of extraction
            only_1d: bool
              Save only the 1D spectra?

        Returns:

        """
        # Build the final list of specobjs and vel_corr
        all_specobjs = specobjs.SpecObjs()

        vel_corr = 0.  # This will not be set for Standard stars, which is fine
        for key in s_dict:
            if key in ['meta']:
                vel_corr = s_dict['meta']['vel_corr']
                continue
            #
            try:
                all_specobjs.add_sobj(s_dict[key]['specobjs'])
            except KeyError:  # No object extracted
                continue

        if len(all_specobjs) == 0:
            msgs.warn('No objects to save!')
            return

        # Write 1D spectra
        save_format = 'fits'
        if save_format == 'fits':
            outfile = os.path.join(self.par['rdx']['redux_path'], self.par['rdx']['scidir'],
                                   'spec1d_{:s}.fits'.format(basename))
            helio_dict = dict(refframe='pixel'
            if self.caliBrate.par['wavelengths']['reference'] == 'pixel'
            else self.caliBrate.par['wavelengths']['frame'],
                              vel_correction=vel_corr)
            # Did the user re-run a single detector?
            save.save_1d_spectra_fits(all_specobjs, self.fitstbl[sidx], outfile,
                                      helio_dict=helio_dict, telescope=self.spectrograph.telescope,
                                      update_det=self.par['rdx']['detnum'])
        #        elif save_format == 'hdf5':
        #            debugger.set_trace()  # NEEDS REFACTORING
        #            arsave.save_1d_spectra_hdf5(None)
        else:
            msgs.error(save_format + ' is not a recognized output format!')
        # 1D only?
        if only_1d:
            return
        # Obj info
        save.save_obj_info(all_specobjs, self.fitstbl, self.spectrograph, basename,
                           os.path.join(self.par['rdx']['redux_path'], self.par['rdx']['scidir']))
        # Write 2D images for the Science Frame
        rawfile = os.path.join(self.fitstbl['directory'][sidx], self.fitstbl['filename'][sidx]) # Need the rawfile header
        save.save_2d_images(s_dict, rawfile, self.spectrograph.primary_hdrext,
                            self.setup, self.caliBrate.master_dir,
                            os.path.join(self.par['rdx']['redux_path'], self.par['rdx']['scidir']),
                            basename, update_det=self.par['rdx']['detnum'])
        return all_specobjs

    def _extract_one(self):
        """
        Dummy method for object extraction

        Returns:

        """
        assert False

    def _extract_std(self):
        """
        Dummy method for std extraction

        Returns:

        """
        assert False

# This is no longer required

#    def _init_calibrations(self):
#        """
#        Instantiate the Calibrations class
#        Returns:
#
#        """
#        # TODO -- Need to make save_masters and write_qa optional
#        # Init calib dict
#        self.caliBrate \
#                = calibrations.MultiSlitCalibrations(self.fitstbl, spectrograph=self.spectrograph,
#                                                     par=self.par['calibrations'],
#                                                     redux_path=self.par['rdx']['redux_path'],
#                                                     save_masters=True, write_qa=True,
#                                                     show=self.show)



    def init_one_science(self, sci_ID, det):
        """
        Instantiate ScienceImage class and run the first step with it

        Args:
            sci_ID: int
              binary flag indicating the science frame
            det: int
              detector index

        Returns:
            self.obstime : Time
            self.basename : str
        """
        self.sci_ID = sci_ID
        self.det = det

        sci_image_files = self.fitstbl.find_frame_files('science', sci_ID=sci_ID)
        scidx = self.fitstbl.find_frames('science', sci_ID=sci_ID, index=True)[0]
        self.sciI = scienceimage.ScienceImage(self.spectrograph, sci_image_files, det=det,
                                              binning = self.fitstbl['binning'][scidx],
                                              objtype='science', scidx=scidx, setup=self.setup,
                                              par=self.par['scienceimage'],
                                              frame_par=self.par['scienceframe'])
        msgs.sciexp = self.sciI  # For QA on crash

        # Names and time
        self.obstime, self.basename = self.init_time_names(self.fitstbl, scidx)
        # Return
        return self.obstime, self.basename  # For fluxing

    # Move out of the scienceimage class to here. It is more appropriate here where fitstable exists
    def init_time_names(self, fitstbl, scidx):
        """
        Setup the basename (for file output mainly)
        and time objects (for heliocentric)

        Parameters
        ----------
        camera : str
          Taken from settings['mosaic']['camera']
        timeunit : str
          mjd

        Returns
        -------
        self.time : Time
        self.basename : str

        """

        timeunit = self.spectrograph.timeunit
        camera = self.spectrograph.camera

        self.fitstbl = fitstbl

        # TODO: Given that we just read the file header to get the
        # datasec_img in the init function above, I don't see why I need
        # to access the fits table for exptime and binning. This
        # information is also in the headers. By simply pulling the
        # stuff from the header, we would remove the fitstbl entirely.
        # Another option would be to put the datasec_img stuff in the
        # fitstbl for each detector
        self.exptime = self.fitstbl['exptime'][scidx]
        try:
            self.binning = self.fitstbl['binning'][scidx]
        except:
            self.binning = (1,1)

        # This should have been set when we construct the fitstbl
        try:
            tval = Time(fitstbl['time'][scidx], format='mjd')#'%Y-%m-%dT%H:%M:%S.%f')
        except:
            debugger.set_trace()

        # Time
        tiso = Time(tval, format='isot')#'%Y-%m-%dT%H:%M:%S.%f')
        dtime = datetime.datetime.strptime(tiso.value, '%Y-%m-%dT%H:%M:%S.%f')
        self.time = tval
        # Basename
        self.inst_name = camera
        self.target_name = self.fitstbl['target'][scidx].replace(" ", "")
        self.basename = self.target_name+'_'+self.inst_name+'_'+ \
                         datetime.datetime.strftime(dtime, '%Y%b%dT') + \
                         tiso.value.split("T")[1].replace(':','')
        # Return
        return self.time, self.basename


    def init_setup(self, pypeit_file, redux_path=None, calibration_check=True):
        """
        Prepare to run redux on a setup

        Args:
            pypeit_file: str
            redux_path: str, optional
            calibration_check (:obj:`bool`, optional):
                Only check that the calibration frames are appropriately
                setup and exist on disk.  Pypit is expected to execute
                in a way that ends after this class is fully
                instantiated such that the user can inspect the results
                before proceeding. 

        Returns:

        """
        self.pypeit_file = pypeit_file

        # This loads the file and sets the following internals:
        #  self.par
        #  self.fitstbl
        #  self.setup_dict
        self._setup(self.pypeit_file, calibration_check=calibration_check)

        # Make the output directories
        if redux_path is not None:
            self.par['rdx']['redux_path'] = redux_path
        else:
            self.par['rdx']['redux_path'] = os.getcwd()
        msgs.info("Setting reduction path to {:s}".format(self.par['rdx']['redux_path']))
        paths.make_dirs(self.spectrograph.spectrograph, self.par['calibrations']['caldir'],
                        self.par['rdx']['scidir'], self.par['rdx']['qadir'],
                        overwrite=self.overwrite, redux_path=self.par['rdx']['redux_path'])
        # Instantiate Calibration class
        self.caliBrate \
                = calibrations.MultiSlitCalibrations(self.fitstbl, spectrograph=self.spectrograph,
                                                     par=self.par['calibrations'],
                                                     redux_path=self.par['rdx']['redux_path'],
                                                     save_masters=True, write_qa=True,
                                                     show=self.show)
        #self._init_calibrations()

    def _make_setup_pypeit_file(self, files_root, extension='.fits', overwrite=False):
        """
        Generate a single PypeIt File for a setup

        Args:
            files_root: str
            extension: str, optional
              Extension of data files
            overwrite: bool, optional

        Returns:
            pargs: ArgParse
            outdir: str
            pypeit_file: str

        """

        # setup_files dir
        outdir = os.path.join(self.setups_path, 'setup_files')
        msgs.info('Setup files will be written to: {0}'.format(outdir))
        if not os.path.isdir(outdir):
            os.mkdir(outdir)

        # Generate a dummy .pypeit file
        date = str(datetime.date.today().strftime('%Y-%b-%d'))
        root = self.spectrograph.spectrograph+'_'+date
        pypeit_file = os.path.join(outdir, root+'.pypeit')

        # Generate
        dfname = '{:s}/*{:s}*'.format(files_root, extension) \
                    if os.path.isdir(files_root) else '{:s}*{:s}*'.format(files_root, extension)
        # configuration lines
        cfg_lines = ['[rdx]']
        cfg_lines += ['    spectrograph = {0}'.format(self.spectrograph.spectrograph)]
        cfg_lines += ['    sortroot = {0}'.format(root)]
        make_pypeit_file(pypeit_file, self.spectrograph, [dfname], cfg_lines=cfg_lines,
                         setup_mode=True)
        msgs.info('Wrote template pypeit file: {0}'.format(pypeit_file))

        # Parser
        pinp = [pypeit_file, '-p', '-r {0}'.format(root) ]
        if overwrite:
            pinp += ['-o']
        pargs = run_pypeit.parser(pinp)
        # Return
        return pargs, outdir, pypeit_file

    def msgs_reset(self):
        """
        Reset the msgs object

        Returns:

        """

        # Reset the global logger
        msgs.reset(log=self.logname, verbosity=self.verbosity)
        msgs.pypeit_file = self.pypeit_file

    def print_end_time(self):
        """
        Print the elapsed time

        Returns:

        """
        # Capture the end time and print it to user
        tend = time.time()
        codetime = tend-self.tstart
        if codetime < 60.0:
            msgs.info('Execution time: {0:.2f}s'.format(codetime))
        elif codetime/60.0 < 60.0:
            mns = int(codetime/60.0)
            scs = codetime - 60.0*mns
            msgs.info('Execution time: {0:d}m {1:.2f}s'.format(mns, scs))
        else:
            hrs = int(codetime/3600.0)
            mns = int(60.0*(codetime/3600.0 - hrs))
            scs = codetime - 60.0*mns - 3600.0*hrs
            msgs.info('Execution time: {0:d}h {1:d}m {2:.2f}s'.format(hrs, mns, scs))

    # TODO: deprecate use_header_id, not propagated anyway
    def _setup(self, pypeit_file, setup_only=False, calibration_check=False, use_header_id=False,
               sort_dir=None):
        """
        Setup PypeIt
           Check files for all calibrations
           Generate the FITS table + write to disk

        Args:
            pypeit_file (str):
            setup_only (:obj:`bool`, optional):
                Only this setup will be performed.  Pypit is expected to
                execute in a way that ends after this class is fully
                instantiated such that the user can inspect the results
                before proceeding.  This has the effect of providing
                more output describing the success of the setup and how
                to proceed, and provides warnings (instead of errors)
                for issues that may cause the reduction itself to fail.

            calibration_check (bool, optional):
                Only check that the calibration frames are appropriately
                setup and exist on disk.  Pypit is expected to execute in a
                way that ends after this class is fully instantiated such
                that the user can inspect the results before proceeding.
            use_header_id (:obj:`bool`, optional):
                Allow setup to use the frame types drawn from the file
                headers using the instrument specific keywords.
            sort_dir (:obj:`str`, optional):
                The directory to put the '.sorted' file.

        Returns:

        """

        # Msgs
        self.msgs_reset()

        # Perform the setup
        self.pypeitSetup = pypeitsetup.PypeItSetup.from_pypeit_file(pypeit_file)
        self.par, _, self.fitstbl, self.setup_dict \
                = self.pypeitSetup.run(setup_only=setup_only, calibration_check=calibration_check,
                                       use_header_id=use_header_id, sort_dir=sort_dir)
        # Write the fits table
        # This is now done in run()
        #self.pypeitSetup.write_metadata(setup_only=setup_only, sort_dir=sort_dir)

    def show_science(self):
        """
        Simple print of science frames

        Returns:

        """
        indx = self.fitstbl.find_frames('science')
        print(self.fitstbl[['target','ra','dec','exptime','dispname','sci_ID']][indx])

    def __repr__(self):
        # Generate sets string
        return '<{:s}: pypeit_file={}>'.format(self.__class__.__name__, self.pypeit_file)


class MultiSlit(PypeIt):
    """
    Child of PypeIt for Multislit and Longslit reductions

    """
    def __init__(self, spectrograph, **kwargs):
        super(MultiSlit, self).__init__(spectrograph, **kwargs)

    def _extract_one(self, std=False):
        """
        Extract a single exposure/detector pair

        sci_ID and det need to have been set internally prior to calling this method

        Returns:
            sciimg
            sciivar
            skymodel
            objmodel
            ivarmodel
            outmask
            sobjs
            vel_corr

        """

        # Standard star specific
        if std:
            # Dict
            msgs.info("Processing standard star")
            if self.std_idx in self.std_dict.keys():
                if self.det in self.std_dict[self.std_idx].keys():
                    return
            else:
                self.std_dict[self.std_idx] = {}
            # Files
            std_image_files = self.fitstbl.find_frame_files('standard', sci_ID=self.sci_ID)
            if self.par['calibrations']['standardframe'] is None:
                msgs.warn('No standard frame parameters provided.  Using default parameters.')

            # Instantiate for the Standard
            # TODO: Uses the same trace and extraction parameter sets used for the science
            # frames.  Should these be different for the standards?
            self.stdI = scienceimage.ScienceImage(self.spectrograph, file_list=std_image_files,
                                          frame_par=self.par['calibrations']['standardframe'],
                                          det=self.det, binning = self.fitstbl['binning'][self.std_idx],
                                          setup=self.setup, scidx=self.std_idx,
                                          objtype='standard', par=self.par['scienceimage'])
            # Names and time
            _, self.std_basename = self.init_time_names(self.fitstbl,self.std_idx)
            #
            sciI = self.stdI
        else:
            sciI = self.sciI

        # Process images (includes inverse variance image, rn2 image,
        # and CR mask)
        sciimg, sciivar, rn2img, crmask \
                = sciI.process(self.caliBrate.msbias, self.caliBrate.mspixflatnrm,
                               self.caliBrate.msbpm, illum_flat=self.caliBrate.msillumflat,
                               apply_gain=True, trim=self.caliBrate.par['trim'], show=self.show)

        # Object finding, first pass on frame without sky subtraction
        maskslits = self.caliBrate.maskslits.copy()
        if not std:
            sobjs_obj0, nobj0 = sciI.find_objects(self.caliBrate.tslits_dict, skysub=False,
                                                   maskslits=maskslits)

        # Global sky subtraction, first pass. Uses skymask from object
        # finding
        global_sky0 = sciI.global_skysub(self.caliBrate.tslits_dict,
                                         self.caliBrate.tilts_dict['tilts'],
                                         use_skymask=True,maskslits=maskslits, show=self.show)

        # Object finding, second pass on frame *with* sky subtraction.
        # Show here if requested
        sobjs_obj, nobj = sciI.find_objects(self.caliBrate.tslits_dict, skysub=True,
                                            maskslits=maskslits, show_peaks=self.show)

        if std:
            if nobj == 0:
                msgs.warn('No objects to extract for standard frame' + msgs.newline()
                      + self.fitstbl['filename'][self.sciI.scidx])
                return
            # Extract
            skymodel, objmodel, ivarmodel, outmask, sobjs = self.stdI.local_skysub_extract(sobjs_obj,
                self.caliBrate.mswave, maskslits=maskslits,show_profile=self.show, show=self.show)

            # Save for fluxing and output later
            self.std_dict[self.std_idx][self.det] = {}
            self.std_dict[self.std_idx][self.det]['basename'] = self.std_basename
            self.std_dict[self.std_idx][self.det]['specobjs'] = sobjs
            # Done
            return

        # If there are objects, do 2nd round of global_skysub,
        # local_skysub_extract, flexure, geo_motion
        vel_corr = None
        if nobj > 0:
            # Global sky subtraction second pass. Uses skymask from object finding
            global_sky = sciI.global_skysub(self.caliBrate.tslits_dict,
                                            self.caliBrate.tilts_dict['tilts'], use_skymask=True,
                                            maskslits=maskslits, show=self.show)

            skymodel, objmodel, ivarmodel, outmask, sobjs \
                    = sciI.local_skysub_extract(sobjs_obj, self.caliBrate.mswave, maskslits=maskslits,
                                                show_profile=self.show, show=self.show)

            # Flexure correction?
            if self.par['flexure']['method'] != 'skip':
                sky_file, sky_spectrum = self.spectrograph.archive_sky_spectrum()
                flex_list = wave.flexure_obj(sobjs, maskslits, self.par['flexure']['method'],
                                             sky_spectrum, sky_file=sky_file,
                                             mxshft=self.par['flexure']['maxshift'])
                # QA
                wave.flexure_qa(sobjs, maskslits, self.basename, self.det, flex_list,
                                out_dir=self.par['rdx']['redux_path'])

            # Helio
            # Correct Earth's motion
            # vel_corr = -999999.9
            if (self.caliBrate.par['wavelengths']['frame'] in ['heliocentric', 'barycentric']) \
                    and (self.caliBrate.par['wavelengths']['reference'] != 'pixel'):
                if sobjs is not None:
                    msgs.info("Performing a {0} correction".format(
                                                    self.caliBrate.par['wavelengths']['frame']))

                    vel, vel_corr \
                            = wave.geomotion_correct(sobjs, maskslits, self.fitstbl, self.sciI.scidx,
                                                     self.obstime,
                                                     self.spectrograph.telescope['longitude'],
                                                     self.spectrograph.telescope['latitude'],
                                                     self.spectrograph.telescope['elevation'],
                                                     self.caliBrate.par['wavelengths']['frame'])
                else:
                    msgs.info('There are no objects on detector {0} to perform a '.format(self.det)
                              + '{1} correction'.format(self.caliBrate.par['wavelengths']['frame']))
            else:
                msgs.info('A wavelength reference-frame correction will not be performed.')

        else:
            msgs.warn('No objects to extract for science frame' + msgs.newline()
                      + self.fitstbl['filename'][self.sciI.scidx])
            # set to first pass global sky
            skymodel = global_sky0
            objmodel = np.zeros_like(sciimg)
            # Set to sciivar. Could create a model but what is the point?
            ivarmodel = np.copy(sciivar)
            # Set to inmask in case on objects were found
            outmask = sciI.mask
            # empty specobjs object from object finding
            sobjs = sobjs_obj

        return sciimg, sciivar, skymodel, objmodel, ivarmodel, outmask, sobjs, vel_corr

    def _extract_std(self):
        self._extract_one(std=True)

    # THESE ARENT USED YET BUT WE SHOULD CONSIDER IT
    @staticmethod
    def default_sci_find_obj_steps():
        return ['process', 'find_objects', 'global_skysub', 'find_objects']

    @staticmethod
    def default_std_find_obj_steps():
        return ['process', 'global_skysub', 'find_objects']



class Echelle(PypeIt):
    """
    Child of PypeIt for Multislit and Longslit reductions

    """
    def __init__(self, spectrograph, **kwargs):
        super(Echelle, self).__init__(spectrograph, **kwargs)

    # JFH This should replace the init_one_science above in Multislit, i.e. it should be in the PypeIt parent
    # class
    def init_one_science(self, sci_ID, frametype, det):
         """
         Instantiate ScienceImage class and run the first step with it

         Args:
             sci_ID: int
               binary flag indicating the science frame
             det: int
               detector index

         Returns:
             self.obstime : Time
             self.basename : str
         """
         self.sci_ID = sci_ID
         self.det = det

         sci_image_files = self.fitstbl.find_frame_files(frametype, sci_ID=sci_ID)
         scidx = self.fitstbl.find_frames(frametype, sci_ID=sci_ID, index=True)[0]
         self.sciI = scienceimage.ScienceImage(self.spectrograph, sci_image_files, det=det,
                                               binning=self.fitstbl['binning'][scidx],
                                               objtype=frametype, scidx=scidx, setup=self.setup,
                                               par=self.par['scienceimage'],
                                               frame_par=self.par['scienceframe'])
         msgs.sciexp = self.sciI  # For QA on crash

         # Names and time
         self.obstime, self.basename = self.init_time_names(self.fitstbl,scidx)
         # Return
         return self.obstime, self.basename  # For fluxing

    # JFH This should be moved to PypeIt
    def reduce_all(self, reuse_masters=False):
        """
        Reduce all of the science exposures
        Generate all needed calibration files

        Args:
            reuse_masters (:obj:`bool`, optional):
                Use the master frames if available (same as setting
                par['calibrations']['masters'] = 'reuse'.

        Returns:

        """

        self.tstart = time.time()

        # Science IDs are in a binary system: 1,2,4,8, etc.
        all_sci_ID = self.fitstbl['sci_ID'][self.fitstbl.find_frames('science')]
        numsci = len(all_sci_ID)

        # Grab the standards
        all_std_ID = self.fitstbl['sci_ID'][self.fitstbl.find_frames('standard')]
        numstd = len(all_std_ID)

        # Check par
        required = ['rdx', 'calibrations', 'scienceframe', 'scienceimage', 'flexure', 'fluxcalib']
        can_be_None = ['flexure', 'fluxcalib']
        self.par.validate_keys(required=required, can_be_None=can_be_None)


        # Reduce the standards first
        for kk, std_ID in enumerate(all_std_ID):
            std_dict = self.reduce_exposure(std_ID, 'standard', reuse_masters=reuse_masters)
            stddx = self.fitstbl.find_frames('standard', sci_ID=std_ID, index=True)[0]
            self.save_exposure(stddx, std_dict, self.basename)

        # Save
        for kk,sci_ID in enumerate(all_sci_ID):
            sci_dict = self.reduce_exposure(sci_ID, 'science', reuse_masters=reuse_masters)
            scidx = self.fitstbl.find_frames('science', sci_ID=sci_ID, index=True)[0]
            self.save_exposure(scidx, sci_dict, self.basename)

        # Finish
        self.print_end_time()

    # JFH This simpler reduce_exposure should be moved to PypeIt
    def reduce_exposure(self, sci_ID, frametype, reuse_masters=False):
        """
        Reduce a single science exposure

        Args:
            sci_ID: int
              binary flag indicating the science frame
            reuse_masters: bool, optional
              Reuse MasterFrame files (where available)


        Returns:
            sci_dict: dict
              dict containing the primary outputs of extraction

        """
        self.sci_ID = sci_ID

        # Insist on re-using MasterFrames where applicable
        if reuse_masters:
            self.par['calibrations']['masters'] = 'reuse'

        sci_dict = OrderedDict()  # This needs to be ordered
        sci_dict['meta'] = {}
        sci_dict['meta']['vel_corr'] = 0.
        #
        scidx = self.fitstbl.find_frames(frametype, sci_ID=sci_ID, index=True)[0]
        msgs.info("Reducing file {0:s}, target {1:s}".format(self.fitstbl['filename'][scidx],
                                                             self.fitstbl['target'][scidx]))
        # Loop on Detectors
        for kk in range(self.spectrograph.ndet):
            det = kk + 1  # Detectors indexed from 1
            self.det = det
            if self.par['rdx']['detnum'] is not None:
                detnum = [self.par['rdx']['detnum']] if isinstance(self.par['rdx']['detnum'],int) else self.par['rdx']['detnum']
                if det not in map(int, detnum):
                    msgs.warn("Skipping detector {:d}".format(det))
                    continue
                else:
                    msgs.warn("Restricting the reduction to detector {:d}".format(det))
            # Setup
            msgs.info("Working on detector {0}".format(det))
            sci_dict[det] = {}

            # Calibrate
            self.calibrate_one(sci_ID, det)

            # Init ScienceImage class
            self.sci_ID = sci_ID
            self.det = det
            self.init_one_science(sci_ID, frametype, det)
            # Extract

            # ToDO make this a method load_std_trace(). Not yet implemented
            # Does a standard exist in the fitstbl? If so grab it since we need the trace. In the future this should
            # use calibgroup matching criteria
            if (frametype == 'science') and np.any(self.fitstbl.find_frames('standard',sci_ID=sci_ID)):
                stddx = self.fitstbl.find_frames('standard',sci_ID=sci_ID, index=True)[0]
                _, std_basename = self.init_time_names(self.fitstbl, stddx)
                std_outfile = os.path.join(self.par['rdx']['redux_path'], self.par['rdx']['scidir'],
                                       'spec1d_{:s}.fits'.format(std_basename))
                pass
                #std_trace = load.load_std_trace(std_outfile)

            sciimg, sciivar, skymodel, objmodel, ivarmodel, outmask, sobjs, vel_corr = self._extract_one()

            # Save for outputing (after all detectors are done)
            sci_dict[det]['sciimg'] = sciimg
            sci_dict[det]['sciivar'] = sciivar
            sci_dict[det]['skymodel'] = skymodel
            sci_dict[det]['objmodel'] = objmodel
            sci_dict[det]['ivarmodel'] = ivarmodel
            sci_dict[det]['outmask'] = outmask
            sci_dict[det]['specobjs'] = sobjs   #utils.unravel_specobjs([specobjs])
            if vel_corr is not None:
                sci_dict['meta']['vel_corr'] = vel_corr

        # Return
        return sci_dict

    # JFH This is the beginning of the echelle class control flow
    def _extract_one(self):
        """
        Extract a single exposure/detector pair

        sci_ID and det need to have been set internally prior to calling this method

        Returns:
            sciimg
            sciivar
            skymodel
            objmodel
            ivarmodel
            outmask
            sobjs
            vel_corr

        """

        sciI = self.sciI

        # Process images (includes inverse variance image, rn2 image,
        # and CR mask)
        sciimg, sciivar, rn2img, crmask \
                = sciI.process(self.caliBrate.msbias, self.caliBrate.mspixflatnrm,
                               self.caliBrate.msbpm, illum_flat=self.caliBrate.msillumflat,
                               apply_gain=True, trim=self.caliBrate.par['trim'], show=self.show)

        # Object finding, first pass on frame without sky subtraction
        maskslits = self.caliBrate.maskslits.copy()

        # Do one iteration of object finding, and sky subtract to get initial sky model
        initial_sky = sciI.get_init_sky(self.caliBrate.tslits_dict, self.caliBrate.tilts_dict['tilts'], show = self.show)

        sobjs_ech, nobj = sciI.get_ech_objects(self.caliBrate.tslits_dict, show=self.show)


        # If there are objects, do 2nd round of global_skysub,
        # local_skysub_extract, flexure, geo_motion
        vel_corr = None
        if nobj > 0:
            skymodel, objmodel, ivarmodel, outmask, sobjs \
                    = sciI.local_skysub_extract(sobjs_ech, self.caliBrate.mswave, maskslits=maskslits,
                                                show_profile=self.show, show=self.show)

            # Flexure correction?
            if self.par['flexure']['method'] != 'skip':
                sky_file, sky_spectrum = self.spectrograph.archive_sky_spectrum()
                flex_list = wave.flexure_obj(sobjs, maskslits, self.par['flexure']['method'],
                                             sky_spectrum, sky_file=sky_file,
                                             mxshft=self.par['flexure']['maxshift'])
                # QA
                wave.flexure_qa(sobjs, maskslits, self.basename, self.det, flex_list,
                                out_dir=self.par['rdx']['redux_path'])

            # Helio
            # Correct Earth's motion
            # vel_corr = -999999.9
            if (self.caliBrate.par['wavelengths']['frame'] in ['heliocentric', 'barycentric']) \
                    and (self.caliBrate.par['wavelengths']['reference'] != 'pixel'):
                if sobjs is not None:
                    msgs.info("Performing a {0} correction".format(
                                                    self.caliBrate.par['wavelengths']['frame']))

                    vel, vel_corr \
                            = wave.geomotion_correct(sobjs, maskslits, self.fitstbl, self.sciI.scidx,
                                                     self.obstime,
                                                     self.spectrograph.telescope['longitude'],
                                                     self.spectrograph.telescope['latitude'],
                                                     self.spectrograph.telescope['elevation'],
                                                     self.caliBrate.par['wavelengths']['frame'])
                else:
                    msgs.info('There are no objects on detector {0} to perform a '.format(self.det)
                              + '{1} correction'.format(self.caliBrate.par['wavelengths']['frame']))
            else:
                msgs.info('A wavelength reference-frame correction will not be performed.')

        else:
            msgs.warn('No objects to extract for science frame' + msgs.newline()
                      + self.fitstbl['filename'][self.sciI.scidx])
            # set to first pass global sky
            skymodel = global_sky0
            objmodel = np.zeros_like(sciimg)
            # Set to sciivar. Could create a model but what is the point?
            ivarmodel = np.copy(sciivar)
            # Set to inmask in case on objects were found
            outmask = sciI.mask
            # empty specobjs object from object finding
            sobjs = sobjs_obj

        return sciimg, sciivar, skymodel, objmodel, ivarmodel, outmask, sobjs, vel_corr

    def _extract_std(self):
        self._extract_one(std=True)

    # THESE ARENT USED YET BUT WE SHOULD CONSIDER IT
    @staticmethod
    def default_sci_find_obj_steps():
        return ['process', 'find_objects', 'global_skysub', 'find_objects']

    @staticmethod
    def default_std_find_obj_steps():
        return ['process', 'global_skysub', 'find_objects']


def instantiate_me(spectrograph, **kwargs):
    """
    Instantiate the PypeIt subclass appropriate for the provided
    spectrograph.

    The class must be subclassed from PypeIt.  See :class:`PypeIt` for
    the description of the valid keyword arguments.

    Args:
        spectrograph
            (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            The instrument used to collect the data to be reduced.

    Returns:
        :class:`PypeIt`: One of the classes with :class:`PypeIt` as its
        base.
    """
    indx = [ c.__name__ == spectrograph.pypeline for c in PypeIt.__subclasses__() ]
    if not np.any(indx):
        msgs.error('Pipeline {0} is not defined!'.format(spectrograph.pypeline))
    return PypeIt.__subclasses__()[np.where(indx)[0][0]](spectrograph, **kwargs)

