"""
Script to run to a single calibration step for an input frame

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase

class ReducebyStep(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):

        parser = super().get_parser(
            description='Run one of the PypeIt reduction steps on a single frame (and detector)',
            width=width, formatter=scriptbase.SmartFormatter
        )
        parser.add_argument(
            'pypeit_file', type=str, help='PypeIt reduction file (must have .pypeit extension)'
        )
        parser.add_argument(
            'frame', type=str,
            help='Raw science/standard frame to reduce as listed in your PypeIt file, e.g. '
                 'b28.fits.gz.')
        parser.add_argument(
            'step', type=str, choices=['process', 'findobj', 'extract'],
            help='Reduction step to perform.  Must be "process" to perform basic image processing '
                 '(bias subtraction, field flattening, etc), "findobj" to perform object '
                 'detection and sky subtraction, or "extract" to extract 1D spectra.'
        )
        parser.add_argument(
            '--det', default=None, type=str,
            help='Single detector number or Mosaic tuple. The Mosaic tuple must include the parentheses '
                 'and be provided as a string, e.g. "(1,2)". Required, but the list of options is provided '
                 'if nothing is provided.'
        )
        parser.add_argument(
            '--show', default=False, action='store_true',
            help='Show reduction steps via plots (which will block further execution until '
                 'clicked on) and outputs to ginga. Requires remote control ginga session via '
                 '"ginga --modules=RC,SlitWavelength &"'
        )
        return parser

    @classmethod
    def main(cls, args):

        import numpy as np
        from pathlib import Path

        from pypeit.core import parse
        from pypeit import pypeit
        from pypeit import pypeit_steps
        from pypeit import log
        from pypeit import PypeItError
        from pypeit import outputfiles
        from pypeit.images import pypeitimage
        from pypeit import specobjs
        from pypeit import spec2dobj
        from pypeit import exposure
        from pypeit import slittrace

        from IPython import embed

        # Set a default log file based on the name of the pypeit file, not the
        # name of the script
        # NOTE: This is copied directly from run_pypeit.py
        if args.log_file == 'default':
            _pypeit_file = Path(args.pypeit_file)
            if _pypeit_file.suffix != '.pypeit':
                raise PypeItError('Input file must have a .pypeit extension!')
            args.log_file = _pypeit_file.with_suffix('.log')

        cls.init_log(args)

        # Instantiate the main pipeline reduction object
        pypeIt = pypeit.PypeIt(args.pypeit_file, show=args.show)
        pypeIt.reuse_calibs = True

        # Detector
        detectors = pypeIt.spectrograph.select_detectors()
        mosaics = pypeIt.spectrograph.allowed_mosaics
        if args.det is None:

            print("---------------------------------------------------------------------")
            print("---------------------------------------------------------------------")
            print("---------------------------------------------------------------------")
            print(f"No detector provided (--det). Choose from one of these: {detectors}.")
            if len(mosaics) > 0:
                print("")
                print(f"This instrument also supports the following mosaics: {mosaics}")
                print('To reduce a mosaic, provide the mosaic as the detectors tuple, e.g. --det "(1,2)"')
                print("")
            print("---------------------------------------------------------------------")
            return
        else:
            det = pypeIt.spectrograph.select_detectors(subset=parse.eval_detectors(args.det))
            if len(det) > 1:
                raise PypeItError("The input --det must be a single detector or mosaic.")
            det = det[0]
        # detector name
        det_name = pypeIt.spectrograph.get_det_name(det)

        # Find the frame
        mt_row = pypeIt.fitstbl['filename'] == args.frame
        if np.sum(mt_row) != 1:
            raise PypeItError(f"Frame {args.frame} not found or not unique")
        frame = int(np.where(mt_row)[0][0])
        calib_IDs = pypeIt.fitstbl.find_frame_calib_groups(frame)
        if len(calib_IDs) != 1:
            raise PypeItError(
                f"Frame {args.frame} is a calibration frame.  This script is for science/standard "
                "frames only"
            )
        calib_ID = calib_IDs[0]

        # Sci metadata
        objtype_out, calib_key, obstime, basename, binning = pypeit_steps.get_sci_metadata(
            pypeIt.spectrograph, pypeIt.fitstbl, frame, det)

        # Find all the frames
        comb_id = pypeIt.fitstbl['comb_id'][frame]
        frames, bg_frames = pypeIt.fitstbl.get_frames_from_combid(comb_id)

        # Intermediate filenames
        sci_filename = outputfiles.intermediate_filename('sciImg', basename, det_name)
        bkg_filename = outputfiles.intermediate_filename('bkgImg', basename, det_name)
        sky_filename = outputfiles.intermediate_filename('Sky', basename,  det_name)
        bkgredux_sky_filename = outputfiles.intermediate_filename('BkgReduxSky', basename, det_name)
        spec1d_filename = outputfiles.intermediate_filename('spec1d', basename, 'all')
        slits_filename = outputfiles.intermediate_filename('slits', basename, det_name)

        # Prep for background subtraction and finding negative traces
        has_bg, bkg_redux, find_negative = pypeit_steps.set_bkg_negative(
            pypeIt.fitstbl, pypeIt.par, bg_frames)

        # Process?
        if args.step == 'process':
            pypeit_steps.process_one_det(
                pypeIt.spectrograph, pypeIt.fitstbl, pypeIt.par,
                frames, det, calib_ID, pypeIt.calibrations_path,
                bg_frames=bg_frames, sci_outfile=sci_filename,
                bkg_outfile=bkg_filename)

            # All done
            return

        # Find Objects
        if args.step == 'findobj':

            # Load intermediate frames needed for finding objects
            log.info(f'Loading images for detector {det}')
            sciImg = pypeitimage.PypeItImage.from_file(sci_filename)
            if bg_frames is not None and len(bg_frames) > 0:
                bkg_redux_sciimg = pypeitimage.PypeItImage.from_file(bkg_filename)
            else:
                bkg_redux_sciimg = None

            # Load up the standard star spec1d file if it exists
            if objtype_out == 'science':
                is_standard = pypeIt.fitstbl.find_frames('standard')
                frame_indx = np.arange(len(pypeIt.fitstbl))
                try:
                    std_outfile = outputfiles.get_std_outfile(pypeIt.fitstbl, pypeIt.par, 
                                                          frame_indx[is_standard])
                except PypeItError:
                    log.warning(
                        'No reduced standard star spec1d file found for this science frame, but '
                        'one was expected because it is in your PypeIt file.\n  Continuing '
                        'without standard star information.'
                    )
                    std_outfile = None
            else:
                std_outfile = None

            # #####################################
            # find objects + initial sky subtraction
            initial_sky, sobjs_obj_find, objFind = pypeit_steps.findobj_on_det(
                sciImg, pypeIt.spectrograph, pypeIt.fitstbl, pypeIt.par,
                frames, calib_ID, det, pypeIt.calibrations_path,
                bkg_redux=bkg_redux,
                find_negative=find_negative,
                std_outfile=std_outfile, show=args.show)

            # #####################################
            # slitmask stuff
            if pypeIt.par['reduce']['slitmask']['assign_obj']:
                detname = sciImg.detector.name
                sciImg_dict = {detname: sciImg}
                slits, sobjs_obj_find = exposure.adjust_for_slitmask(
                    sciImg_dict,  pypeIt.spectrograph, pypeIt.fitstbl, pypeIt.par,
                    frames[0], sobjs_obj_find, [objFind.slits])
                objFind.slits = slits[0]

            # #####################################
            # final sky subtraction
            final_global_sky, bkg_redux_global_sky, objFind = \
                pypeit_steps.finalize_sky_det(pypeIt.spectrograph, pypeIt.fitstbl, pypeIt.par, frames[0],
                                              det, objFind, initial_sky, sobjs_obj_find,
                                              bkg_redux_sciimg=bkg_redux_sciimg, bkg_redux=bkg_redux, show=args.show)

            # update Slits
            _slits = objFind.slits
            flagged_slits = np.where(objFind.reduce_bpm)[0]
            if len(flagged_slits) > 0:
                _slits.mask[flagged_slits] = \
                    _slits.bitmask.turn_on(_slits.mask[flagged_slits], 'BADSKYSUB')

            # Update the sciImg with the scaleImg information
            sciImg.rel_scaleImg = objFind.scaleimg
            # and the global spectral flexure shift
            sciImg.flex_shift = objFind.slitshift

            # Write
            # sobjs object found
            sobjs_obj_find.write_to_fits({}, spec1d_filename)
            log.info(f'Wrote intermediate spec1d file with objects found to {spec1d_filename}')
            # final sky image
            skyimg = pypeitimage.PypeItImage(final_global_sky)
            if not sky_filename.parent.is_dir():
                sky_filename.parent.mkdir()
            skyimg.to_file(sky_filename, overwrite=True)
            log.info(f'Wrote final sky image to {sky_filename}')
            # bkg_redux sky image
            if bkg_redux_global_sky is not None:
                bkgredux_skyimg = pypeitimage.PypeItImage(bkg_redux_global_sky)
                bkgredux_skyimg.to_file(bkgredux_sky_filename, overwrite=True)
                log.info(f'Wrote bkg_redux final sky image to {bkgredux_sky_filename}')
            # slits
            _slits.to_file(slits_filename, overwrite=True)
            log.info(f'Wrote intermediate slits to {slits_filename}')
            # updated sciImg
            sciImg.to_file(sci_filename, overwrite=True)
            log.info(f'Wrote updated science image to {sci_filename}')


        # Extract?
        if args.step == 'extract':
            # Load intermediate frames needed for the extraction
            log.info(f'Loading images for detector {det}')
            sciImg = pypeitimage.PypeItImage.from_file(sci_filename)

            # sky images
            log.info(f'Loading sky image for detector {det}')
            if not sky_filename.is_file():
                raise PypeItError(f'Sky image {sky_filename} not found!')
            skyimg = pypeitimage.PypeItImage.from_file(sky_filename)
            skyimg = skyimg.image
            if bkgredux_sky_filename.is_file():
                log.info(f'Loading bkg_redux sky image for detector {det}')
                bkg_redux_skyimg = pypeitimage.PypeItImage.from_file(bkgredux_sky_filename)
                bkg_redux_skyimg = bkg_redux_skyimg.image
            else:
                bkg_redux_skyimg = None
            # specobjs from findobj
            log.info(f'Loading spec1d file for detector {det}')
            if not spec1d_filename.is_file():
                raise PypeItError(f'spec1d file {spec1d_filename} not found!')
            specobjs_objfind = specobjs.SpecObjs.from_fitsfile(spec1d_filename)
            # slits
            log.info(f'Loading slits for detector {det}')
            if not slits_filename.is_file():
                raise PypeItError(f'Slits file {slits_filename} not found!')
            calib_slits = slittrace.SlitTraceSet.from_file(slits_filename)

            # Container for Spec2DObj
            all_spec2d = spec2dobj.AllSpec2DObj()
            all_spec2d['meta']['bkg_redux'] = bkg_redux
            all_spec2d['meta']['find_negative'] = find_negative

            # Grab the objects
            detname = sciImg.detector.name
            if specobjs_objfind.nobj > 0:
                specobjs_on_det = specobjs_objfind[specobjs_objfind.DET == detname]
            else:
                specobjs_on_det = specobjs_objfind

            # Extract
            all_spec2d[detname], sobjs_extract = pypeit_steps.extract_det(
                pypeIt.spectrograph, pypeIt.fitstbl, pypeIt.par, frames, det,
                calib_ID, pypeIt.calibrations_path,
                sciImg, skyimg, specobjs_on_det,
                calib_slits, bkg_redux_final_sky=bkg_redux_skyimg,
                bkg_redux=bkg_redux, find_negative=find_negative)

            # TODO: Should we add the calibration associations to the SpecObjs object as done in the main run?

            # Save it
            exposure.save_exposure(pypeIt.spectrograph, pypeIt.fitstbl,
                                   pypeIt.par, frame, all_spec2d, sobjs_extract,
                                   pypeIt.calibrations_path,
                                   in_update_det=det)
