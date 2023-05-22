"""
Script for performing 2d coadds of PypeIt data.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
from pypeit.scripts import scriptbase

class CoAdd2DSpec(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Coadd 2D spectra produced by PypeIt',
                                    width=width)

        parser.add_argument('coadd2d_file', type=str, default=None,
                            help='File to guide 2d coadds')

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

        parser.add_argument("--debug", default=False, action="store_true", help="show debug plots?")
        parser.add_argument('-v', '--verbosity', type=int, default=1,
                            help='Verbosity level between 0 [none] and 2 [all]. Default: 1. '
                                 'Level 2 writes a log with filename coadd_2dspec_YYYYMMDD-HHMM.log')

        # TODO: Make spec_samp_fact and spat_samp_fact parameters in CoAdd2DPar,
        # and then move these to setup_coadd2d.py
        parser.add_argument('--spec_samp_fact', default=1.0, type=float,
                            help="Make the wavelength grid finer (spec_samp_fact < 1.0) or "
                                 "coarser (spec_samp_fact > 1.0) by this sampling factor, i.e. "
                                 "units of spec_samp_fact are pixels.")
        parser.add_argument('--spat_samp_fact', default=1.0, type=float,
                            help="Make the spatial grid finer (spat_samp_fact < 1.0) or coarser "
                                 "(spat_samp_fact > 1.0) by this sampling factor, i.e. units of "
                                 "spat_samp_fact are pixels.")

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

        from pathlib import Path
        import os
        import glob
        import copy
        from collections import OrderedDict

        from IPython import embed

        import numpy as np

        from astropy.io import fits

        from pypeit import msgs
        from pypeit import coadd2d
        from pypeit import inputfiles
        from pypeit import specobjs
        from pypeit import spec2dobj
        from pypeit.spectrographs.util import load_spectrograph

        # Set the verbosity, and create a logfile if verbosity == 2
        msgs.set_logfile_and_verbosity('coadd_2dspec', args.verbosity)

        # Load the file
        coadd2dFile = inputfiles.Coadd2DFile.from_file(args.coadd2d_file)
        spectrograph, par, _ = coadd2dFile.get_pypeitpar()

        # Check some of the parameters
        # TODO Heliocentric for coadd2d needs to be thought through. Currently turning it off.
        if par['calibrations']['wavelengths']['refframe'] != 'observed':
            msgs.warn('Wavelength reference frame shift (e.g., heliocentric correction) not yet '
                      'fully developed.  Ignoring input and setting "refframe = observed".')
            par['calibrations']['wavelengths']['refframe'] = 'observed'
        # TODO Flexure correction for coadd2d needs to be thought through. Currently turning it off.
        if par['flexure']['spec_method'] != 'skip':
            msgs.warn('Spectroscopic flexure correction not yet fully developed.  Skipping.')
            par['flexure']['spec_method'] = 'skip'
        # TODO This is currently the default for 2d coadds, but we need a way to toggle it on/off
        if not par['reduce']['findobj']['skip_skysub']:
            msgs.warn('Must skip sky subtraction when finding objects (i.e., sky should have '
                      'been subtracted during primary reduction procedure).  Skipping.')
            par['reduce']['findobj']['skip_skysub'] = True

        # Get the files
        spec2d_files = coadd2dFile.filenames

        # Get the paths
        coadd_scidir, qa_path = map(lambda x : Path(x).resolve(),
                                    coadd2d.CoAdd2D.output_paths(spec2d_files, par))

        # Get the output basename
        head2d = fits.getheader(spec2d_files[0])
        basename = coadd2d.CoAdd2D.default_basename(spec2d_files) \
                        if args.basename is None else args.basename

        # Write the par to disk
        par_outfile = coadd_scidir.parent / f'{basename}_coadd2d.par'
        print(f'Writing full parameter set to {par_outfile}.')
        par.to_config(par_outfile, exclude_defaults=True, include_descr=False)

        # Now run the coadds
        bkg_redux = head2d['SKYSUB'] == 'DIFF'
        find_negative = head2d['FINDOBJ'] == 'POS_NEG'

        # Print status message
        msgs_string = f'Reducing target {basename}' + msgs.newline()
        msgs_string += f"Coadding frame sky-subtracted with {head2d['SKYSUB']}" + msgs.newline()
        msgs_string += f"Searching for objects that are {head2d['FINDOBJ']}" + msgs.newline()
        msgs_string += 'Combining frames in 2d coadd:' + msgs.newline()
        for f, file in enumerate(spec2d_files):
            msgs_string += f'Exp {f}: {Path(file).name}' + msgs.newline()
        msgs.info(msgs_string)

        # Instantiate the sci_dict
        # TODO Why do we need this sci_dict at all?? JFH
        sci_dict = OrderedDict()  # This needs to be ordered
        sci_dict['meta'] = {}
        sci_dict['meta']['vel_corr'] = 0.
        sci_dict['meta']['bkg_redux'] = bkg_redux
        sci_dict['meta']['find_negative'] = find_negative

        # Find the detectors to reduce
        # TODO: Allow slitspatnum to be specified?  E.g.:
#        detectors = spectrograph.select_detectors(
#                subset=par['rdx']['detnum'] if par['rdx']['slitspatnum'] is None else
#                        par['rdx']['slitspatnum'])
        detectors = spectrograph.select_detectors(subset=par['rdx']['detnum'])
        msgs.info(f'Detectors to work on: {detectors}')

        # container for specobjs
        all_specobjs = specobjs.SpecObjs()
        # container for spec2dobj
        all_spec2d = spec2dobj.AllSpec2DObj()
        # set some meta
        all_spec2d['meta']['bkg_redux'] = bkg_redux
        all_spec2d['meta']['find_negative'] = find_negative

        # Loop on detectors
        for det in detectors:
            msgs.info("Working on detector {0}".format(det))

            # Instantiate Coadd2d
            coadd = coadd2d.CoAdd2D.get_instance(spec2d_files, spectrograph, par, det=det,
                                                 offsets=par['coadd2d']['offsets'],
                                                 weights=par['coadd2d']['weights'],
                                                 spec_samp_fact=args.spec_samp_fact,
                                                 spat_samp_fact=args.spat_samp_fact,
                                                 bkg_redux=bkg_redux, find_negative=find_negative,
                                                 debug_offsets=args.debug_offsets,
                                                 debug=args.debug)

            # TODO Add this stuff to a run method in coadd2d
            # Coadd the slits
            coadd_dict_list = coadd.coadd(only_slits=par['coadd2d']['only_slits'])
            # Create the pseudo images
            pseudo_dict = coadd.create_pseudo_image(coadd_dict_list)
            # Reduce
            msgs.info('Running the extraction')

            # TODO -- This should mirror what is in pypeit.extract_one
            
            # TODO -- JFH :: This ought to return a Spec2DObj and SpecObjs which
            # would be slurped into AllSpec2DObj and all_specobsj, as below.
            
            # TODO -- JFH -- Check that the slits we are using are correct

            sci_dict[coadd.detname] = {}
            sci_dict[coadd.detname]['sciimg'], sci_dict[coadd.detname]['sciivar'], \
                sci_dict[coadd.detname]['skymodel'], sci_dict[coadd.detname]['objmodel'], \
                sci_dict[coadd.detname]['ivarmodel'], sci_dict[coadd.detname]['outmask'], \
                sci_dict[coadd.detname]['specobjs'], sci_dict[coadd.detname]['detector'], \
                sci_dict[coadd.detname]['slits'], sci_dict[coadd.detname]['tilts'], \
                sci_dict[coadd.detname]['waveimg'] \
                    = coadd.reduce(pseudo_dict, show=args.show, show_peaks=args.peaks, basename=basename)

            # Tack on detector (similarly to pypeit.extract_one)
            for sobj in sci_dict[coadd.detname]['specobjs']:
                sobj.DETECTOR = sci_dict[coadd.detname]['detector']

            # fill the specobjs container
            all_specobjs.add_sobj(sci_dict[coadd.detname]['specobjs'])

            # fill the spec2dobj container but first ...
            # pull out maskdef_designtab from sci_dict[det]['slits']
            maskdef_designtab = sci_dict[coadd.detname]['slits'].maskdef_designtab
            slits = copy.deepcopy(sci_dict[coadd.detname]['slits'])
            slits.maskdef_designtab = None
            # fill up
            all_spec2d[coadd.detname] = spec2dobj.Spec2DObj(sciimg=sci_dict[coadd.detname]['sciimg'],
                                                          ivarraw=sci_dict[coadd.detname]['sciivar'],
                                                          skymodel=sci_dict[coadd.detname]['skymodel'],
                                                          objmodel=sci_dict[coadd.detname]['objmodel'],
                                                          ivarmodel=sci_dict[coadd.detname]['ivarmodel'],
                                                          scaleimg=np.array([1.0], dtype=float),
                                                          bpmmask=sci_dict[coadd.detname]['outmask'],
                                                          detector=sci_dict[coadd.detname]['detector'],
                                                          slits=slits,
                                                          wavesol=None,
                                                          waveimg=sci_dict[coadd.detname]['waveimg'],
                                                          tilts=sci_dict[coadd.detname]['tilts'],
                                                          sci_spat_flexure=None,
                                                          sci_spec_flexure=None,
                                                          vel_corr=None,
                                                          vel_type=None,
                                                          maskdef_designtab=maskdef_designtab)

        # SAVE TO DISK

        # THE FOLLOWING MIMICS THE CODE IN pypeit.save_exposure()
        subheader = spectrograph.subheader_for_spec(head2d, head2d)
        # Write spec1D
        if all_specobjs.nobj > 0:
            outfile1d = coadd_scidir / f'spec1d_{basename}.fits'
            all_specobjs.write_to_fits(subheader, outfile1d)

            # Info
            outfiletxt = coadd_scidir / f'spec1d_{basename}.txt'
            sobjs = specobjs.SpecObjs.from_fitsfile(outfile1d, chk_version=False)
            sobjs.write_info(outfiletxt, spectrograph.pypeline)

        # Build header for spec2d
        outfile2d = coadd_scidir / f'spec2d_{basename}.fits'
        pri_hdr = all_spec2d.build_primary_hdr(head2d, spectrograph,
                                               subheader=subheader,
                                               # TODO -- JFH :: Decide if we need any of these
                                               redux_path=None)
        # Write spec2d
        all_spec2d.write_to_fits(outfile2d, pri_hdr=pri_hdr)


