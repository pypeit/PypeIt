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
            'step', type=str,
            help='Reduction step to perform.  Must be "process" to perform basic image processing '
                 '(bias subtraction, field flattening, etc), "findobj" to perform object '
                 'detection and initial sky subtraction, or "extract" to extract 1D spectra.'
        )
        parser.add_argument(
            '--det', default=None, type=str,
            help='Detector number or Mosaic tuple. Required, but the list of options is provided '
                 'if nothing is provided.'
        )
        parser.add_argument(
            '--show', default=False, action='store_true',
            help='Show reduction steps via plots (which will block further execution until '
                 'clicked on) and outputs to ginga. Requires remote control ginga session via '
                 '"ginga --modules=RC,SlitWavelength &"'
        )
        return parser

    @staticmethod
    def main(args):

        import numpy as np
        from pathlib import Path

        from pypeit import pypeit
        from pypeit import pypeit_steps
        from pypeit import msgs
        from pypeit import pypmsgs
        from pypeit import outputfiles
        from pypeit.images import pypeitimage
        from pypeit import specobjs
        from pypeit import spec2dobj
        from pypeit import exposure

        from IPython import embed

        # Load options from command line
        pypeit_file = Path(args.pypeit_file).absolute()
        logname = pypeit_file.parent / f'{pypeit_file.stem}.log'

        # Instantiate the main pipeline reduction object
        pypeIt = pypeit.PypeIt(args.pypeit_file, logname=logname, show=args.show) 
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
                print("To reduce a mosaic, provide the mosaic as the detector, e.g. --det '1,2'")
                print("")
            print("---------------------------------------------------------------------")
            return
        else: 
            det = pypeIt.spectrograph.select_detectors(subset=args.det)
            if len(det) > 1:
                msgs.error("The input --det must be a single detector or mosaic.")
        # detector name
        det_name = pypeIt.spectrograph.get_det_name(det)

        # Find the frame
        mt_row = pypeIt.fitstbl['filename'] == args.frame
        if np.sum(mt_row) != 1:
            msgs.error(f"Frame {args.frame} not found or not unique")
        frame = int(np.where(mt_row)[0][0])
        calib_IDs = pypeIt.fitstbl.find_frame_calib_groups(frame)
        if len(calib_IDs) != 1:
            msgs.error(f"Frame {args.frame} is a calibration frame.  This script is for science/standard frames only")
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
        initsky_filename = outputfiles.intermediate_filename('initSky', basename,  det_name)
        spec1d_filename = outputfiles.intermediate_filename('spec1d', basename, 'all')

        # Prep for background subtraction and finding negative traces
        has_bg, bkg_redux, find_negative = pypeit_steps.set_bkg_negative(
            pypeIt.fitstbl, pypeIt.par, bg_frames)

        # Process?
        if args.step == 'process':
            sciImg, bkg_redux_sciimg = pypeit_steps.process_one_det(
                pypeIt.spectrograph, pypeIt.fitstbl, pypeIt.par,
                frames, det, calib_ID, pypeIt.calibrations_path,
                bg_frames=bg_frames, sci_outfile=sci_filename,
                bkg_outfile=bkg_filename)

            # All done
            return

        msgs.info(f'Loading images for detector {det}')
        sciImg = pypeitimage.PypeItImage.from_file(sci_filename)
        if bg_frames is not None and len(bg_frames) > 0:
            bkg_redux_sciimg = pypeitimage.PypeItImage.from_file(bkg_filename)
        else:
            bkg_redux_sciimg = None

        # Find Objects
        if args.step == 'findobj':

            # Load up the standard star spec1d file if it exists
            if objtype_out == 'science':
                is_standard = pypeIt.fitstbl.find_frames('standard')
                frame_indx = np.arange(len(pypeIt.fitstbl))
                try:
                    std_outfile = outputfiles.get_std_outfile(pypeIt.fitstbl, pypeIt.par, 
                                                          frame_indx[is_standard])
                except pypmsgs.PypeItError:
                    msgs.warn('No reduced standard star spec1d file found for this science frame, but one was expected because it is in your PypeIt file.\n'+\
                        'Continuing without standard star information.')
                    std_outfile = None
            else:
                std_outfile = None                                                    

            # Do it
            initial_sky, sobjs_obj, _ = pypeit_steps.findobj_on_det(
                sciImg, pypeIt.spectrograph, pypeIt.fitstbl, pypeIt.par,
                frames, calib_ID, det, pypeIt.calibrations_path,
                bkg_redux=bkg_redux,
                find_negative=find_negative,
                std_outfile=std_outfile, show=args.show)

            # TODO -- write from findobj_on_det?

            # Write
            init_pypeit = pypeitimage.PypeItImage(initial_sky)
            if not initsky_filename.parent.is_dir():
                initsky_filename.parent.mkdir()
            init_pypeit.to_file(initsky_filename, overwrite=True)
            #
            sobjs_obj.write_to_fits({}, spec1d_filename)

            # All done
            return

        # Load
        msgs.info(f'Loading initial sky for detector {det}')
        tmp = pypeitimage.PypeItImage.from_file(initsky_filename)
        initial_sky = tmp.image
        #
        specobjs_objfind = specobjs.SpecObjs.from_fitsfile(spec1d_filename)

        # TODO -- Add objs in here!! -- what did I mean by this?
            
        # Extract?
        if args.step == 'extract':
            this_calib_silts = None

            # Slitmask?
            if pypeIt.par['reduce']['slitmask']['assign_obj']:
                sciImg_dict = {}
                calib_slits = []
                # TODO - ask Debora to check this..
                if isinstance(det, tuple):
                    for mosaic in mosaics:
                        # Science images
                        sci_filename = outputfiles.intermediate_filename('sciImg', basename, 
                                    pypeIt.spectrograph.get_det_name(mosaic))
                        if not sci_filename.is_file():
                            msgs.warn(f"Science image {sci_filename} not found for adjusting for slitmask.  Skipping mosaic {mosaic}.")
                            continue
                        sciImg_dict[mosaic] = pypeitimage.PypeItImage.from_file(sci_filename)
                        # Grab the calibrations
                        caliBrate = pypeit_steps.load_calibrations_for_frame(
                            pypeIt.spectrograph, pypeIt.fitstbl, pypeIt.par, frames[0], mosaic, 
                            calib_ID, pypeIt.calibrations_path)
                        calib_slits.append(caliBrate.slits)
                else:
                    msgs.error("Need to implement this!")

                # Run me
                calib_slits = exposure.adjust_for_slitmask(
                    sciImg_dict, pypeIt.spectrograph, pypeIt.fitstbl, pypeIt.par, 
                    frames[0], pypeIt.fitstbl['binning'][frames[0]],
                    specobjs_objfind, calib_slits)
                # Update this_calib_silts
                idx = mosaics.index(det)
                this_calib_silts = calib_slits[idx]

            # Proceed with extraction

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
                sciImg, bkg_redux_sciimg, initial_sky, specobjs_on_det,
                bkg_redux=bkg_redux,
                find_negative=find_negative,
                calib_slits=this_calib_silts)

            # Save it
            exposure.save_exposure(pypeIt.spectrograph, pypeIt.fitstbl,
                                   pypeIt.par, frame, all_spec2d, sobjs_extract,
                                   pypeIt.calibrations_path,
                                   in_update_det=det)