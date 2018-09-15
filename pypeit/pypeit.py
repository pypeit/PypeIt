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

from pypeit import msgs
from pypeit.core import fsort
from pypeit.core import qa
from pypeit.par.util import make_pypeit_file, parse_pypeit_file

from pypeit import pypeitsetup
from pypeit.scripts import run_pypeit
from pypeit.core import pypsetup
from pypeit import calibrations
from pypeit import scienceimage
from pypeit.core import wave
from pypeit import specobjs
from pypeit.core import save

from pypeit import debugger


class PypeIt(object):
    """
    This class is designed to run PypeIt

    .. todo::
        Improve docstring...

    Parameters
    ----------
    spectrograph: Spectrograph
    setups_path: str, optional
      Path for files related to all setups
    verbosity: int, optional
      2=verbose
    overwrite: bool, optional
      Overwrite files
    logname: str, optional
      filename for log file
    show: bool, optional
      Show QA on screen along the way

    Attributes
    ----------
        pypeit_file (:obj:`str`):
            Name of the pypeit file to read.  Pypit files have a specific
            set of valid formats. A description can be found `here`_
            (include doc link).

    Inherited Attributes
    --------------------
    """
    __metaclass__ = ABCMeta

    def __init__(self, spectrograph, setups_path=None, verbosity=2,
                 overwrite=True, logname=None, show=False):

        # Init
        self.spectrograph = spectrograph
        self.verbosity = verbosity
        self.overwrite = overwrite
        if setups_path is None:
            self.setups_path = os.getcwd()
        else:
            self.setups_path = setups_path

        # Internals
        self.pypeit_file = None
        self.logname = logname
        self.setup_pypeit_file = None
        self.redux_path = None
        self.show=show
        self.setup = None
        self.det = None
        self.sci_ID = None


    def build_setup_files(self, files_root):
        """
        Generate the setup files for PypeIt from a list of input files

        Args:
            files_root: str
              Root name for the files to be reduced including full path

        Returns:

        """

        # Record the starting time
        self.tstart = time.time()

        pargs, sort_dir, self.setup_pypeit_file = self._make_setup_pypeit_file(files_root)
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
        cfg_lines, data_files, frametype, setups = parse_pypeit_file(self.setup_pypeit_file)
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
        Dummy method.  Set in a child

        Args:
            sci_ID:
            det:

        Returns:

        """
        assert False

    def reduce_all(self, reuse_masters=False):
        """
        Reduce all of the science exposures
        Generate all needed calibration files

        Args:
            reuse_masters: bool, optional
              Reuse MasterFrame files (where available)

        Returns:

        """

        self.tstart = time.time()
        std_dict = {}
        all_sci_ID = self.fitstbl['sci_ID'].data[self.fitstbl['science']]  # Binary system: 1,2,4,8, etc.
        numsci = len(all_sci_ID)
        basenames = [None]*numsci  # For fluxing at the very end

        # Check par
        required = [ 'rdx', 'calibrations', 'scienceframe', 'scienceimage', 'flexure', 'fluxcalib' ]
        can_be_None = [ 'flexure', 'fluxcalib' ]
        self.par.validate_keys(required=required, can_be_None=can_be_None)

        for sci_ID in all_sci_ID:
            sci_dict = self.reduce_exposure(sci_ID, reuse_masters=reuse_masters)
            self.save_exposure(sci_ID, sci_dict)

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
        scidx = np.where((self.fitstbl['sci_ID'] == sci_ID) & self.fitstbl['science'])[0][0]
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

        # Return
        return sci_dict

    def save_exposure(self, sci_ID, sci_dict):
        """
        Save the outputs from extraction for a given exposure

        Args:
            sci_ID: int
              binary flag indicating the science frame
            sci_dict: dict
              dict containing the primary outputs of extraction

        Returns:

        """
        self.sci_ID = sci_ID
        scidx = np.where((self.fitstbl['sci_ID'] == self.sci_ID) & self.fitstbl['science'])[0][0]

        # Build the final list of specobjs and vel_corr
        all_specobjs = specobjs.SpecObjs()

        for key in sci_dict:
            if key in ['meta']:
                continue
            #
            try:
                all_specobjs.add_sobj(sci_dict[key]['specobjs'])
            except KeyError:  # No object extracted
                continue

        if len(all_specobjs) == 0:
            msgs.warn('No objects to save!')
            return

        # Write 1D spectra
        save_format = 'fits'
        if save_format == 'fits':
            outfile = os.path.join(self.par['rdx']['redux_path'], self.par['rdx']['scidir'], 'spec1d_{:s}.fits'.format(self.basename))
            helio_dict = dict(refframe='pixel'
            if self.caliBrate.par['wavelengths']['reference'] == 'pixel'
            else self.caliBrate.par['wavelengths']['frame'],
                              vel_correction=sci_dict['meta']['vel_corr'])
            save.save_1d_spectra_fits(all_specobjs, self.fitstbl[scidx], outfile,
                                      helio_dict=helio_dict, telescope=self.spectrograph.telescope)
        #        elif save_format == 'hdf5':
        #            debugger.set_trace()  # NEEDS REFACTORING
        #            arsave.save_1d_spectra_hdf5(None)
        else:
            msgs.error(save_format + ' is not a recognized output format!')
        # Obj info
        save.save_obj_info(all_specobjs, self.fitstbl, self.spectrograph, self.basename,
                           os.path.join(self.par['rdx']['redux_path'], self.par['rdx']['scidir']))
        # Write 2D images for the Science Frame
        save.save_2d_images(sci_dict, self.fitstbl, scidx, self.spectrograph.primary_hdrext,
                            self.setup, self.caliBrate.master_dir,
                            os.path.join(self.par['rdx']['redux_path'], self.par['rdx']['scidir']),
                            self.basename)
        return

    def _extract_one(self):
        """
        Dummy method for extraction

        Returns:

        """
        assert False

    def _init_calibrations(self):
        """
        Dummy method for instantiating a Calibration class

        Returns:

        """
        pass

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

        sci_image_files = fsort.list_of_files(self.fitstbl, 'science', sci_ID)
        scidx = np.where((self.fitstbl['sci_ID'] == sci_ID) & self.fitstbl['science'])[0][0]
        self.sciI = scienceimage.ScienceImage(self.spectrograph, sci_image_files,
                                         det=det, objtype='science',
                                         scidx=scidx, setup=self.setup,
                                         par=self.par['scienceimage'],
                                         frame_par=self.par['scienceframe'])
        msgs.sciexp = self.sciI  # For QA on crash

        # Names and time
        self.obstime, self.basename = self.sciI.init_time_names(self.fitstbl)
        # Return
        return self.obstime, self.basename  # For fluxing


    def init_setup(self, pypeit_file, redux_path=None, calibration_check=True):
        """
        Prepare to run redux on a setup

        Args:
            pypeit_file: str
            redux_path: str, optional
            calibration_check: bool, optional
              Check calibrations

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
        fsort.make_dirs(self.spectrograph.spectrograph, self.par['calibrations']['caldir'],
                        self.par['rdx']['scidir'], self.par['rdx']['qadir'], overwrite=self.overwrite,
                        redux_path=self.par['rdx']['redux_path'])
        # Instantiate Calibration class
        self._init_calibrations()

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
        pypeit_file = outdir+'/'+root+'.pypeit'

        # Generate
        dfname = "{:s}*{:s}*".format(files_root, extension)
        # configuration lines
        cfg_lines = ['[rdx]']
        cfg_lines += ['    spectrograph = {0}'.format(self.spectrograph.spectrograph)]
        cfg_lines += ['    sortroot = {0}'.format(root)]
        make_pypeit_file(pypeit_file, self.spectrograph, [dfname], cfg_lines=cfg_lines, setup_mode=True)
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
            msgs.info('Data reduction execution time: {0:.2f}s'.format(codetime))
        elif codetime/60.0 < 60.0:
            mns = int(codetime/60.0)
            scs = codetime - 60.0*mns
            msgs.info('Data reduction execution time: {0:d}m {1:.2f}s'.format(mns, scs))
        else:
            hrs = int(codetime/3600.0)
            mns = int(60.0*(codetime/3600.0 - hrs))
            scs = codetime - 60.0*mns - 3600.0*hrs
            msgs.info('Data reduction execution time: {0:d}h {1:d}m {2:.2f}s'.format(hrs, mns, scs))


    def _setup(self, pypeit_file, setup_only=False, calibration_check=False, use_header_frametype=False, sort_dir=None):
        """
        Setup PypeIt
           Check files for all calibrations
           Generate the FITS table + write to disk

        Args:
            pypeit_file (str):
            setup_only (bool, optional):
                Only this setup will be performed.  Pypit is expected to
                execute in a way that ends after this class is fully
                instantiated such that the user can inspect the results
                before proceeding.  This has the effect of providing more
                output describing the success of the setup and how to
                proceed, and provides warnings (instead of errors) for
                issues that may cause the reduction itself to fail.
            calibration_check (bool, optional):
                Only check that the calibration frames are appropriately
                setup and exist on disk.  Pypit is expected to execute in a
                way that ends after this class is fully instantiated such
                that the user can inspect the results before proceeding.
            use_header_frametype (bool, optional):
                Allow setup to use the frame types drawn from the file
                headers using the instrument specific keywords.
            sort_dir (str, optional):
                The directory to put the '.sorted' file.

        Returns:

        """

        # Msgs
        self.msgs_reset()

        # Perform the setup
        self.pypeitSetup = pypeitsetup.PypeItSetup.from_pypeit_file(pypeit_file)
        self.par, _, self.fitstbl, self.setup_dict = self.pypeitSetup.run(setup_only=setup_only,
                                                           calibration_check=calibration_check,
                                                           use_header_frametype=use_header_frametype,
                                                           sort_dir=sort_dir)
        # Write the fits table
        self.pypeitSetup.write_fitstbl()

    def show_science(self):
        """
        Simple print of science frames

        Returns:

        """
        print(self.fitstbl[['target','ra','dec','exptime','dispname','sci_ID']][self.fitstbl['science']])

    def __repr__(self):
        # Generate sets string
        txt = '<{:s}: pypeit_file={}'.format(self.__class__.__name__,
                                                          self.pypeit_file)
        txt += '>'
        return txt


class MultiSlit(PypeIt):
    """
    Child of PypeIt for Multislit and Longslit reductions

    """
    def __init__(self, spectrograph, **kwargs):
        PypeIt.__init__(self, spectrograph, **kwargs)

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
        self.setup = pypsetup.instr_setup(sci_ID, det, self.fitstbl, self.setup_dict,
                                     self.spectrograph.detector[det-1]['numamplifiers'],
                                     must_exist=True)
        # Setup
        self.caliBrate.reset(self.setup, det, sci_ID, self.par['calibrations'])
        # Run em
        self.caliBrate.run_the_steps()

        msgs.info("Successful Calibration!")

    def _init_calibrations(self):
        """
        Instantiate the Calibrations class

        Returns:

        """
        # TODO -- Need to make save_masters and write_qa optional
        # Init calib dict
        self.caliBrate = calibrations.MultiSlitCalibrations(
            self.fitstbl, spectrograph=self.spectrograph,
            par=self.par['calibrations'],
            redux_path=self.par['rdx']['redux_path'],
            save_masters=True, write_qa=True, show=self.show)


    def _extract_one(self):
        """
        Extract a single exposure/detector pair
        sci_ID and det need to have been set internally prior

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
        # TODO -- Turn the following stream into a recipe like in Calibrations
        # TODO -- Should check the Calibs were done already
        scidx = np.where((self.fitstbl['sci_ID'] == self.sci_ID) & self.fitstbl['science'])[0][0]

        # Process images (includes inverse variance image, rn2 image, and CR mask)
        sciimg, sciivar, rn2img, crmask = self.sciI.process(
            self.caliBrate.msbias, self.caliBrate.mspixflatnrm, self.caliBrate.msbpm,
            illum_flat=self.caliBrate.msillumflat, apply_gain=True, trim=self.caliBrate.par['trim'],
            show=self.show)

        # Object finding, first pass on frame without sky subtraction
        maskslits = self.caliBrate.maskslits.copy()
        sobjs_obj0, nobj0 = self.sciI.find_objects(self.caliBrate.tslits_dict, skysub=False,
                                                   maskslits=maskslits)

        # Global sky subtraction, first pass. Uses skymask from object finding
        global_sky0 = self.sciI.global_skysub(self.caliBrate.tslits_dict, self.caliBrate.mstilts,
                                              use_skymask=True, maskslits=maskslits, show=self.show)

        # Object finding, second pass on frame *with* sky subtraction. Show here if requested
        sobjs_obj, nobj = self.sciI.find_objects(self.caliBrate.tslits_dict, skysub=True,
                                                 maskslits=maskslits, show_peaks=self.show)

        # If there are objects, do 2nd round of global_skysub, local_skysub_extract, flexure, geo_motion
        vel_corr = None
        if nobj > 0:
            # Global sky subtraction second pass. Uses skymask from object finding
            global_sky = self.sciI.global_skysub(self.caliBrate.tslits_dict, self.caliBrate.mstilts,
                                                 use_skymask=True, maskslits=maskslits, show=self.show)

            skymodel, objmodel, ivarmodel, outmask, sobjs = self.sciI.local_skysub_extract(self.caliBrate.mswave, maskslits=maskslits,
                                                                                      show_profile=self.show, show=self.show)


            # Flexure correction?
            if self.par['flexure'] is not None and self.par['flexure']['method'] is not None:
                sky_file, sky_spectrum = self.spectrograph.archive_sky_spectrum()
                flex_list = wave.flexure_obj(sobjs, maskslits, self.par['flexure']['method'],
                                             sky_spectrum, sky_file=sky_file,
                                             mxshft=self.par['flexure']['maxshift'])
                # QA
                wave.flexure_qa(sobjs, maskslits, self.basename, self.det, flex_list, out_dir=self.par['rdx']['redux_path'])

            # Helio
            # Correct Earth's motion
            # vel_corr = -999999.9
            if (self.caliBrate.par['wavelengths']['frame'] in ['heliocentric', 'barycentric']) and \
                    (self.caliBrate.par['wavelengths']['reference'] != 'pixel'):
                if sobjs is not None:
                    msgs.info("Performing a {0} correction".format(
                        self.caliBrate.par['wavelengths']['frame']))

                    vel, vel_corr = wave.geomotion_correct(sobjs, maskslits, self.fitstbl, scidx,
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
                      + self.fitstbl['filename'][scidx])
            skymodel = global_sky0  # set to first pass global sky
            objmodel = np.zeros_like(sciimg)
            ivarmodel = np.copy(sciivar)  # Set to sciivar. Could create a model but what is the point?
            outmask = self.sciI.bitmask  # Set to inmask in case on objects were found
            sobjs = sobjs_obj  # empty specobjs object from object finding

        return sciimg, sciivar, skymodel, objmodel, ivarmodel, outmask, sobjs, vel_corr


def instantiate_me(spectrograph, **kwargs):
    """
    Simple wrapper for grabbing the right PypeIt class

    Args:
        name: str
          Allowed options are MultiSlit
        spectrograph: Spectrograph
        **kwargs:

    Returns:

    """
    if spectrograph.pypeit_class() == 'MultiSlit':
        pypeIt = MultiSlit(spectrograph, **kwargs)
    else:
        msgs.error("NOT READY FOR THIS TYPE OF REDUX")
    # Return
    return pypeIt
