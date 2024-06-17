"""
Launch the pypeit_identify GUI tool.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import argparse
from IPython import embed

from pypeit.scripts import scriptbase

class Identify(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Launch PypeIt pypeit_identify tool, display extracted '
                                                'Arc, and load linelist.', width=width)
        parser.add_argument('arc_file', type=str, default=None, help='PypeIt Arc file')
        parser.add_argument('slits_file', type=str, default=None, help='PypeIt Slits file')
        parser.add_argument("--lamps", type=str,
                            help="Comma separated list of calibration lamps (no spaces)")
        parser.add_argument('-s', '--solution', default=False, action='store_true',
                            help="Load a wavelength solution from the arc_file (if it exists)")
        parser.add_argument("--wmin", type=float, default=3000.0, help="Minimum wavelength range")
        parser.add_argument("--wmax", type=float, default=50000.0, help="Maximum wavelength range")
        parser.add_argument("--slits", type=str, default='0',
                            help="Which slit to load for wavelength calibration. " 
                            "Format should be [0,1,...] for multiple slits, 0 for only one slit. "
                            "If creating a new WaveCalib with the -n flag, this is not necessary.")
        parser.add_argument('-m', '--multi', default=False, action = 'store_true',
                            help="Set this flag to create wavelength solutions for muliple slits")
        parser.add_argument('-n', '--new_sol', default=False, action = 'store_true',
                            help="Set this flag to construct a new WaveCalib file, rather than using the exising one")
        parser.add_argument("--det", type=int, default=1, help="Detector index")
        parser.add_argument("--rmstol", type=float, default=0.1, help="RMS tolerance")
        parser.add_argument("--fwhm", type=float, default=4., help="FWHM for line finding")
        parser.add_argument("--sigdetect", type=float, help="sigma detection for line finding")
        parser.add_argument("--pixtol", type=float, default=0.1,
                            help="Pixel tolerance for Auto IDs")
        parser.add_argument('--test', default=False, action='store_true',
                            help=argparse.SUPPRESS)
        parser.add_argument("--linear", default=False, action="store_true",
                            help="Show the spectrum in linear (rather than log) scale")
        parser.add_argument('--force_save', default=False, action='store_true',
                            help="Save the solutions, despite the RMS")
        parser.add_argument('--rescale_resid', default=False, action='store_true',
                            help="Rescale the residual plot to include all points?")
        parser.add_argument('-v', '--verbosity', type=int, default=1,
                            help='Verbosity level between 0 [none] and 2 [all]. Default: 1. '
                                 'Level 2 writes a log with filename identify_YYYYMMDD-HHMM.log')
        parser.add_argument('--try_old', default=False, action='store_true',
                            help='Attempt to load old datamodel versions.  A crash may ensue..')
        return parser

    @staticmethod
    def main(args):

        import os
        import json

        import numpy as np
                
        from pypeit import msgs
        from pypeit.spectrographs.util import load_spectrograph
        from pypeit.core.gui.identify import Identify
        from pypeit.wavecalib import BuildWaveCalib, WaveCalib
        from pypeit import slittrace
        from pypeit.images.buildimage import ArcImage
        from pypeit.core.wavecal import autoid
        from linetools.utils import jsonify


        chk_version = not args.try_old

        # Set the verbosity, and create a logfile if verbosity == 2
        msgs.set_logfile_and_verbosity('identify', args.verbosity)

        # Load the Arc file
        msarc = ArcImage.from_file(args.arc_file, chk_version=chk_version)

        # Load the spectrograph
        spec = load_spectrograph(msarc.PYP_SPEC)
        par = spec.default_pypeit_par()['calibrations']['wavelengths']

        # Get the lamp list
        if args.lamps is None:
            lamps = par['lamps']
            if lamps is None or lamps == ['use_header']:
                msgs.error('Cannot determine the lamps; use --lamps argument')
        else:
            lamps = args.lamps.split(",")
        par['lamps'] = lamps

        # Load the slits
        slits = slittrace.SlitTraceSet.from_file(args.slits_file, chk_version=chk_version)
        # Reset the mask
        slits.mask = slits.mask_init

        msgs.info('Loading in Solution if desired and exists')
        # Check if a solution exists
        solnname = WaveCalib.construct_file_name(msarc.calib_key, calib_dir=msarc.calib_dir)
        wv_calib = WaveCalib.from_file(solnname, chk_version=chk_version) \
                        if os.path.exists(solnname) and args.solution else None
        
        # Load the calibration frame (if it exists and is desired).  Bad-pixel mask
        # set to any flagged pixel in Arc.
        wavecal = BuildWaveCalib(msarc, slits, spec, par, lamps, det=args.det,
                                msbpm=msarc.select_flag())

        # If we are dealing with a multi-trace solution
        if args.multi:

            # Obtain a list of good slits
            ok_mask_idx = np.where(np.logical_not(wavecal.wvc_bpm))[0]

            # print to screen the slit widths if maskdef_designtab is available
            if slits.maskdef_designtab is not None:
                msgs.info("Slit widths (arcsec): {}".format(np.round(slits.maskdef_designtab['SLITWID'].data, 2)))

            # Generate a map of the instrumental spectral FWHM
            # TODO nsample should be a parameter
            fwhm_map = autoid.map_fwhm(wavecal.msarc.image, wavecal.gpm, wavecal.slits_left, wavecal.slits_right, wavecal.slitmask,
                                    nsample=10, slit_bpm=wavecal.wvc_bpm, specord=wavecal.par['fwhm_spec_order'],
                                    spatord=wavecal.par['fwhm_spat_order'])
            # Calculate the typical spectral FWHM down the centre of the slit
            measured_fwhms = np.zeros(slits.nslits, dtype=object)
            for islit in range(slits.nslits):
                if islit not in ok_mask_idx:
                    continue
                # Measure the spectral FWHM (in pixels) at the midpoint of the slit
                # (i.e. the midpoint in both the spectral and spatial directions)
                measured_fwhms[islit] = fwhm_map[islit].eval(wavecal.msarc.image.shape[0]//2, 0.5)

            # Save for redo's
            wavecal.measured_fwhms = measured_fwhms
            if args.new_sol:
                sv_par = par.data.copy()
                j_par = jsonify(sv_par)
                strpar = json.dumps(j_par)

                wv_calib = WaveCalib(wv_fits=None,
                                    fwhm_map=fwhm_map,
                                    arc_spectra=np.array(wavecal.extract_arcs()[0]),
                                    nslits=wavecal.slits.nslits,
                                    spat_ids=wavecal.slits.spat_id,
                                    PYP_SPEC=wavecal.spectrograph.name,
                                    lamps=','.join(wavecal.lamps),
                                    strpar = strpar)



            #iterate over each slit to add to the wv_calib
            if args.new_sol:
                slits_inds = np.arange(slits.nslits)
            else:
                if args.slits == 'all':
                    slits_inds = np.arange(slits.nslits)
                else:
                    slits_inds = np.array(list(slits.strip('[]').split(",")), dtype=int)
            fits_dicts = []
            specdata = []
            wv_fits_arr = []
            lines_pix_arr = []
            lines_wav_arr = []
            lines_fit_ord = []
            custom_wav = []
            custom_wav_ind = []
            for slit_val in slits_inds:
                # Load the calibration frame (if it exists and is desired).  Bad-pixel mask
                # set to any flagged pixel in Arc.
                wv_calib_slit = None
                if wv_calib is not None:
                    if not np.any(wv_calib.wv_fits):
                        wv_calib_slit = None
                    else:
                        if not wv_calib.wv_fits[slit_val]['pypeitfit']:
                            wv_calib_slit = None
                        else:
                            wv_calib_slit = wv_calib

                arccen, arc_maskslit = wavecal.extract_arcs(slitIDs=[slit_val])

                # Launch the identify window
                # TODO -- REMOVE THIS HACK
                try:
                    nonlinear_counts = msarc.detector.nonlinear_counts()
                except AttributeError:
                    nonlinear_counts = None
                arcfitter = Identify.initialise(arccen, lamps, slits, slit=int(slit_val), par=par,
                                                wv_calib_all=wv_calib_slit, wavelim=[args.wmin, args.wmax],
                                                nonlinear_counts=nonlinear_counts,
                                                pxtoler=args.pixtol, test=args.test, 
                                                fwhm=args.fwhm,
                                                sigdetect=args.sigdetect,
                                                specname=spec.name,
                                                y_log=not args.linear,
                                                rescale_resid=args.rescale_resid)

                # If testing, return now
                if args.test:
                    return arcfitter, msarc
                final_fit = arcfitter.get_results()
                fits_dicts.append(arcfitter._fitdict)
                specdata.append(arccen[:,slit_val])
                # Build here to avoid circular import
                #  Note:  This needs to be duplicated in test_scripts.py
                # Wavecalib (wanted when dealing with multiple detectors, eg. GMOS)
                if 'WaveFit' in arcfitter._fitdict.keys():

                    lines_pix_arr.append(arcfitter._fitdict['WaveFit']['pixel_fit'])
                    lines_wav_arr.append(arcfitter._fitdict['WaveFit']['wave_fit'])
                    lines_fit_ord.append(arcfitter._fitdict['WaveFit']['pypeitfit']['order'])
                    if np.any(wv_calib.wv_fits):
                        wv_calib.wv_fits[slit_val] = arcfitter._fitdict['WaveFit']
                        wv_calib.wv_fits[slit_val].fwhm = measured_fwhms[slit_val]
                    else:
                        wv_fits_arr.append(arcfitter._fitdict['WaveFit'])
                        wv_fits_arr[-1].fwhm = measured_fwhms[slit_val]
                    waveCalib = wv_calib

            
                else:
                    if np.any(wv_calib.wv_fits):
                        wv_calib.wv_fits[slit_val] = None
                    else: 
                        wv_fits_arr.append('None')
                    waveCalib = None
                    custom_wav_q = ''
                    while custom_wav_q != 'y' and custom_wav_q != 'n':
                        custom_wav_q = input('No solution made! Do you want to input an estimated linear solution? y/[n]: ')
                    if custom_wav_q == 'y':
                        while True:
                            try:
                                min_wav_str = input('Please enter the desired minimum wavelength: ')
                                min_wav = float(min_wav_str)
                            except ValueError:
                                print("Sorry, try that again...")
                                #better try again... Return to the start of the loop
                                continue
                            try:
                                max_wav_str = input('Please enter the desired maximum wavelength: ')
                                max_wav = float(max_wav_str)

                            except ValueError:
                                print("Sorry, try that again...")
                                #better try again... Return to the start of the loop
                                continue
                            else:
                                #wavelengths were successfully parsed!
                                #we're ready to exit the loop.
                                break
                        nspec = np.shape(arccen)[0]
                        wav_x = np.linspace(min_wav, max_wav, nspec)
                        custom_wav.append(wav_x)
                        custom_wav_ind.append(slit_val)
                        # sample the new solution to make a fake IDs list:
                        lines_pix_arr.append(np.linspace(0, nspec-1, 8)[1:-1])
                        lines_wav_arr.append(np.linspace(min_wav, max_wav, 8)[1:-1])
                        lines_fit_ord.append([1])

            if not np.any(wv_calib.wv_fits):
                wv_calib.wv_fits = np.array(wv_fits_arr)
                #wv_calib.to_file()
            if not args.new_sol:
                wv_calib.wv_fit2d=None
            if args.new_sol:
                wv_calib.copy_calib_internals(msarc)

        # If we just want the normal one-trace output
        else:
            arccen, arc_maskslit = wavecal.extract_arcs(slitIDs=[int(args.slits)])
            # Launch the identify window
            # TODO -- REMOVE THIS HACK
            try:
                nonlinear_counts = msarc.detector.nonlinear_counts()
            except AttributeError:
                nonlinear_counts = None
            arcfitter = Identify.initialise(arccen, lamps, slits, slit=int(args.slits), par=par,
                                            wv_calib_all=wv_calib, wavelim=[args.wmin, args.wmax],
                                            nonlinear_counts=nonlinear_counts,
                                            pxtoler=args.pixtol, test=args.test, 
                                            fwhm=args.fwhm,
                                            sigdetect=args.sigdetect,
                                            specname=spec.name,
                                            y_log=not args.linear,
                                            rescale_resid=args.rescale_resid)

            # If testing, return now
            if args.test:
                return arcfitter, msarc
            final_fit = arcfitter.get_results()

            # Build here to avoid circular import
            #  Note:  This needs to be duplicated in test_scripts.py
            # Wavecalib (wanted when dealing with multiple detectors, eg. GMOS)
            if 'WaveFit' in arcfitter._fitdict.keys():
                waveCalib = WaveCalib(nslits=1, wv_fits=np.atleast_1d(arcfitter._fitdict['WaveFit']),
                                    arc_spectra=np.atleast_2d(arcfitter.specdata).T,
                                    spat_ids=np.atleast_1d(int(arcfitter._spatid)),
                                    PYP_SPEC=msarc.PYP_SPEC, lamps=','.join(lamps))
                waveCalib.copy_calib_internals(msarc)
            else:
                waveCalib = None

            fits_dicts = None
            specdata = None
            slits = None 
            lines_pix_arr = None
            lines_wav_arr = None
            lines_fit_ord = None 
            custom_wav = None 
            custom_wav_ind = None 
        # TODO: Make the following more elegant:
        # fill lines with dummy values to make this work
        # Ask the user if they wish to store the result in PypeIt calibrations
        arcfitter.store_solution(final_fit, slits.binspec,
                                wvcalib=wv_calib,
                                rmstol=args.rmstol,
                                force_save=args.force_save, 
                                multi = args.multi, fits_dicts = fits_dicts,
                                specdata = np.array(specdata),
                                slits = slits,
                                lines_pix_arr = lines_pix_arr,
                                lines_wav_arr = lines_wav_arr,
                                lines_fit_ord = np.array(lines_fit_ord),
                                custom_wav = np.array(custom_wav),
                                custom_wav_ind = np.array(custom_wav_ind) )
            


