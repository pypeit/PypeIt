"""
Script for performing 2d coadds of PypeIt data.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import os
import glob
from collections import OrderedDict

from IPython import embed

import numpy as np

from astropy.io import fits

from pypeit import par, msgs, io
from pypeit import coadd2d
from pypeit import io
from pypeit import specobjs
from pypeit import spec2dobj
from pypeit.pypeit import PypeIt
from pypeit.scripts import scriptbase
from pypeit.spectrographs.util import load_spectrograph


class CoAdd2DSpec(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Coadd 2D spectra produced by PypeIt',
                                    width=width)

        parser.add_argument("--file", type=str, default=None, help="File to guide 2d coadds")
        parser.add_argument('--det', default=None, type=int,
                            help="Only coadd data from this detector (1-indexed)")
        parser.add_argument("--obj", type=str, default=None,
                            help="Object name in lieu of extension, e.g if the spec2d files are "
                                 "named 'spec2d_J1234+5678_GNIRS_2017Mar31T085412.181.fits' then "
                                 "use --obj J1234+5678")
        parser.add_argument("--show", default=False, action="store_true",
                            help="Show the reduction steps. Equivalent to the -s option when "
                                 "running pypeit.")
        parser.add_argument("--debug_offsets", default=False, action="store_true",
                            help="Show QA plots useful for debugging automatic offset "
                                 "determination")
        parser.add_argument("--peaks", default=False, action="store_true",
                            help="Show the peaks found by the object finding algorithm.")
        parser.add_argument("--basename", type=str, default=None,
                            help="Basename of files to save the parameters, spec1d, and spec2d")
        parser.add_argument('--spec_samp_fact', default=1.0, type=float,
                            help="Make the wavelength grid finer (spec_samp_fact < 1.0) or "
                                 "coarser (spec_samp_fact > 1.0) by this sampling factor, i.e. "
                                 "units of spec_samp_fact are pixels.")
        parser.add_argument('--spat_samp_fact', default=1.0, type=float,
                            help="Make the spatial grid finer (spat_samp_fact < 1.0) or coarser "
                                 "(spat_samp_fact > 1.0) by this sampling factor, i.e. units of "
                                 "spat_samp_fact are pixels.")
        parser.add_argument("--debug", default=False, action="store_true", help="show debug plots?")

        # TODO implement an option to only do certian slits
        #parser.add_argument("--only_slits", type=list, default=None,
        #                    help="Only coadd the following slits")

        #parser.add_argument("--wave_method", type=str, default=None,
        #                    help="Wavelength method for wavelength grid. If not set, code will "
        #                         "use linear for Multislit and log10 for Echelle")
        #parser.add_argument("--std", default=False, action="store_true",
        #                    help="This is a standard star reduction.")

        return parser

    @staticmethod
    def main(args):
        """ Executes 2d coadding
        """
        msgs.warn('PATH =' + os.getcwd())
        # Load the file
        if args.file is not None:
            spectrograph_name, config_lines, spec2d_files \
                    = io.read_spec2d_file(args.file, filetype="coadd2d")
            spectrograph = load_spectrograph(spectrograph_name)

            # Parameters
            # TODO: Shouldn't this reinstantiate the same parameters used in
            # the PypeIt run that extracted the objects?  Why are we not
            # just passing the pypeit file?
            # JFH: The reason is that the coadd2dfile may want different reduction parameters
            spectrograph_def_par = spectrograph.default_pypeit_par()
            parset = par.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_def_par.to_config(),
                                                    merge_with=config_lines)
        elif args.obj is not None:
            # TODO: We should probably be reading the pypeit file and using
            # those parameters here rather than using the default parset.
            
            # TODO: This needs to define the science path
            spec2d_files = glob.glob('./Science/spec2d_*' + args.obj + '*')
            head0 = fits.getheader(spec2d_files[0])
            spectrograph_name = head0['PYP_SPEC']
            spectrograph = load_spectrograph(spectrograph_name)
            parset = spectrograph.default_pypeit_par()
        else:
            msgs.error('You must provide either a coadd2d file (--file) or an object name (--obj)')

        # Update with configuration specific parameters (which requires science
        # file) and initialize spectrograph
        spectrograph_cfg_lines = spectrograph.config_specific_par(spec2d_files[0]).to_config()
        parset = par.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_cfg_lines,
                                              merge_with=parset.to_config())

        # If detector was passed as an argument override whatever was in the coadd2d_file
        if args.det is not None:
            msgs.info("Restricting reductions to detector={}".format(args.det))
            parset['rdx']['detnum'] = int(args.det)

        # Get headers (if possible) and base names
        spec1d_files = [files.replace('spec2d', 'spec1d') for files in spec2d_files]
        head1d = None
        for spec1d_file in spec1d_files:
            if os.path.isfile(spec1d_file):
                head1d = fits.getheader(spec1d_file)
                break
        if head1d is None:
            msgs.warn("No 1D spectra so am generating a dummy header for output")
            head1d = io.initialize_header()

        head2d = fits.getheader(spec2d_files[0])
        if args.basename is None:
            #TODO Fix this, currently does not work if target names have - or _
            filename_first = os.path.basename(spec2d_files[0])
            filename_last = os.path.basename(spec2d_files[-1])
            prefix_first = (filename_first.split('_')[1]).split('-')[0]
            prefix_last = (filename_last.split('_')[1]).split('-')[0]
            objname = (filename_first.split('-')[1]).split('_')[0]
            basename = '{:s}-{:s}-{:s}'.format(prefix_first,prefix_last,objname)
        else:
            basename = args.basename

        # Write the par to disk
        par_outfile = basename+'_coadd2d.par'
        print("Writing the parameters to {}".format(par_outfile))
        parset.to_config(par_outfile)

        # Now run the coadds

        skysub_mode = head2d['SKYSUB']
        findobj_mode = head2d['FINDOBJ']
        ir_redux = True if 'DIFF' in skysub_mode else False
        find_negative = True if 'NEG' in findobj_mode else False

        # Print status message
        msgs_string = 'Reducing target {:s}'.format(basename) + msgs.newline()
        msgs_string += 'Coadding frame sky-subtraced with {:s}'.format(skysub_mode)
        msgs_string += 'Searching for objects that are {:s}'.format(findobj_mode)
        msgs_string += msgs.newline() + 'Combining frames in 2d coadd:' + msgs.newline()
        for file in spec2d_files:
            msgs_string += '{0:s}'.format(os.path.basename(file)) + msgs.newline()
        msgs.info(msgs_string)

        # TODO: This needs to be added to the parameter list for rdx
        redux_path = os.getcwd()
        master_dirname = os.path.basename(head2d['PYPMFDIR']) + '_coadd'
        master_dir = os.path.join(redux_path, master_dirname)

        # Make the new Master dir
        if not os.path.isdir(master_dir):
            msgs.info('Creating directory for Master output: {0}'.format(master_dir))
            os.makedirs(master_dir)

        # Instantiate the sci_dict
        sci_dict = OrderedDict()  # This needs to be ordered
        sci_dict['meta'] = {}
        sci_dict['meta']['vel_corr'] = 0.
        sci_dict['meta']['ir_redux'] = ir_redux
        sci_dict['meta']['find_negative'] = find_negative


        # Find the detectors to reduce
        detectors = PypeIt.select_detectors(detnum=parset['rdx']['detnum'], ndet=spectrograph.ndet)
        if len(detectors) != spectrograph.ndet:
            msgs.warn('Not reducing detectors: {0}'.format(' '.join([str(d) for d in
            set(np.arange(spectrograph.ndet) + 1) - set(detectors)])))

        # Loop on detectors
        for det in detectors:
            msgs.info("Working on detector {0}".format(det))
            sci_dict[det] = {}

            # Instantiate Coadd2d
            coadd = coadd2d.CoAdd2D.get_instance(spec2d_files, spectrograph, parset, det=det,
                                                 offsets=parset['coadd2d']['offsets'],
                                                 weights=parset['coadd2d']['weights'],
                                                 spec_samp_fact=args.spec_samp_fact,
                                                 spat_samp_fact=args.spat_samp_fact,
                                                 ir_redux=ir_redux, find_negative=find_negative,
                                                 debug_offsets=args.debug_offsets,
                                                 debug=args.debug)

            # TODO Add this stuff to a run method in coadd2d
            # Coadd the slits
            coadd_dict_list = coadd.coadd(only_slits=None) # TODO implement only_slits later
            # Create the pseudo images
            pseudo_dict = coadd.create_pseudo_image(coadd_dict_list)
            # Reduce
            msgs.info('Running the extraction')

            # TODO -- This should mirror what is in pypeit.extract_one
            
            # TODO -- JFH :: This ought to return a Spec2DObj and SpecObjs which
            # would be slurped into AllSpec2DObj and all_specobsj, as below.
            
            # TODO -- JFH -- Check that the slits we are using are correct

            sci_dict[det]['sciimg'], sci_dict[det]['sciivar'], sci_dict[det]['skymodel'], \
            sci_dict[det]['objmodel'], sci_dict[det]['ivarmodel'], sci_dict[det]['outmask'], \
            sci_dict[det]['specobjs'], sci_dict[det]['detector'], sci_dict[det]['slits'], \
            sci_dict[det]['tilts'], sci_dict[det]['waveimg'] \
                    = coadd.reduce(pseudo_dict, show=args.show, show_peaks=args.peaks)

            # Save pseudo image master files
            #coadd.save_masters()

        # Make the new Science dir
        # TODO: This needs to be defined by the user
        scipath = os.path.join(redux_path, 'Science_coadd')
        if not os.path.isdir(scipath):
            msgs.info('Creating directory for Science output: {0}'.format(scipath))
            os.makedirs(scipath)

        # THE FOLLOWING MIMICS THE CODE IN pypeit.save_exposure()

        # TODO -- These lines should be above once reduce() passes back something sensible
        all_specobjs = specobjs.SpecObjs()
        for det in detectors:
            all_specobjs.add_sobj(sci_dict[det]['specobjs'])

        # Write
        outfile1d = os.path.join(scipath, 'spec1d_{:s}.fits'.format(basename))
        subheader = spectrograph.subheader_for_spec(head2d, head2d)
        all_specobjs.write_to_fits(subheader, outfile1d)

        # 2D spectra
        # TODO -- These lines should be above once reduce() passes back something sensible
        all_spec2d = spec2dobj.AllSpec2DObj()
        all_spec2d['meta']['ir_redux'] = ir_redux
        all_spec2d['meta']['find_negative'] = find_negative
        for det in detectors:
            all_spec2d[det] = spec2dobj.Spec2DObj(det=det,
                                                  sciimg=sci_dict[det]['sciimg'],
                                                  ivarraw=sci_dict[det]['sciivar'],
                                                  skymodel=sci_dict[det]['skymodel'],
                                                  objmodel=sci_dict[det]['objmodel'],
                                                  ivarmodel=sci_dict[det]['ivarmodel'],
                                                  scaleimg=np.array([1.0], dtype=np.float),
                                                  bpmmask=sci_dict[det]['outmask'],
                                                  detector=sci_dict[det]['detector'],
                                                  slits=sci_dict[det]['slits'],
                                                  waveimg=sci_dict[det]['waveimg'],
                                                  tilts=sci_dict[det]['tilts'],
                                                  sci_spat_flexure=None,
                                                  sci_spec_flexure=None,
                                                  vel_corr=None,
                                                  vel_type=None)
        # Build header
        outfile2d = os.path.join(scipath, 'spec2d_{:s}.fits'.format(basename))
        pri_hdr = all_spec2d.build_primary_hdr(head2d, spectrograph,
                                               subheader=subheader,
                                               # TODO -- JFH :: Decide if we need any of these
                                               redux_path=None, master_key_dict=None,
                                               master_dir=None)
        # Write
        all_spec2d.write_to_fits(outfile2d, pri_hdr=pri_hdr)


