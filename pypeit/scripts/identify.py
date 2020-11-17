#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
Launch the identify GUI tool.
"""

from IPython import embed


def parse_args(options=None, return_parser=False):
    import argparse

    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Launch PypeIt identify tool, display extracted '
                                                 'MasterArc, and load linelist.'
                                                 'Run above the Masters/ folder',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('arc_file', type=str, default=None, help='PypeIt MasterArc file')
    parser.add_argument('slits_file', type=str, default=None, help='PypeIt MasterSlits file')
    parser.add_argument("--lamps", type=str, help="Comma separated list of calibration lamps (no spaces)")
    parser.add_argument('-s', '--solution', default=False, action='store_true',
                        help="Load a wavelength solution from the arc_file (if it exists)")
    parser.add_argument("--wmin", type=float, default=3000.0, help="Minimum wavelength range")
    parser.add_argument("--wmax", type=float, default=10000.0, help="Maximum wavelength range")
    parser.add_argument("--slit", type=int, default=0, help="Which slit to load for wavelength calibration")
    parser.add_argument("--det", type=int, default=1, help="Detector index")
    parser.add_argument("--rmstol", type=float, default=0.1, help="RMS tolerance")
    parser.add_argument("--fwhm", type=float, default=4., help="FWHM for line finding")
    parser.add_argument("--pixtol", type=float, default=0.1, help="Pixel tolerance for Auto IDs")
    parser.add_argument('--test', default=False, action='store_true',
                        help="Unit tests?")
    parser.add_argument('--force_save', default=False, action='store_true',
                        help="Save the solutions, despite the RMS")

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):

    import os
    import sys
    import numpy as np
    from pypeit import masterframe
    from pypeit.spectrographs.util import load_spectrograph
    from pypeit.core.gui.identify import Identify
    from pypeit.core.wavecal import waveio
    from pypeit.wavecalib import BuildWaveCalib, WaveCalib
    from pypeit import slittrace
    from pypeit.images.buildimage import ArcImage

    # Load the MasterArc file
    if os.path.exists(args.arc_file):
        arcfil = args.arc_file
    else:
        try:
            arcfil = "Masters/{0:s}".format(args.arc_file)
        except FileNotFoundError:
            print("Could not find MasterArc file.")
            sys.exit()
    msarc = ArcImage.from_file(arcfil)

    mdir = msarc.head0['MSTRDIR']
    mkey = msarc.head0['MSTRKEY']

    # Load the spectrograph
    specname = msarc.head0['PYP_SPEC']
    spec = load_spectrograph(specname)
    par = spec.default_pypeit_par()['calibrations']['wavelengths']

    # Get the lamp list
    if args.lamps is None:
        lamplist = par['lamps']
        if lamplist is None:
            print("ERROR :: Cannot determine the lamps")
            sys.exit()
    else:
        lamplist = args.lamps.split(",")
    par['lamps'] = lamplist

    # Load the slits
    slits = slittrace.SlitTraceSet.from_file(args.slits_file)
    # Reset the mask
    slits.mask = slits.mask_init

    # Check if a solution exists
    solnname = os.path.join(mdir, masterframe.construct_file_name(WaveCalib, mkey))
    wv_calib = waveio.load_wavelength_calibration(solnname) if os.path.exists(solnname) and args.solution else None

    # Load the MasterFrame (if it exists and is desired)?
    wavecal = BuildWaveCalib(msarc, slits, spec, par, binspectral=slits.binspec, det=args.det,
                        master_key=mkey, msbpm=msarc.fullmask)
    arccen, arc_maskslit = wavecal.extract_arcs(slitIDs=[args.slit])

    # Launch the identify window
    arcfitter = Identify.initialise(arccen, slits, slit=int(args.slit), par=par, wv_calib_all=wv_calib,
                                    wavelim=[args.wmin, args.wmax],
                                    nonlinear_counts=spec.nonlinear_counts(msarc.detector),
                                    pxtoler=args.pixtol, test=args.test, fwhm=args.fwhm)
    # Testing?
    if args.test:
        return arcfitter
    final_fit = arcfitter.get_results()

    # Build here to avoid circular import
    #  Note:  This needs to be duplicated in test_scripts.py
    # Wavecalib (wanted when dealing with multiple detectors, eg. GMOS)
    if 'WaveFit' in arcfitter._fitdict.keys():
        waveCalib = WaveCalib(nslits=1, wv_fits=np.atleast_1d(arcfitter._fitdict['WaveFit']),
                                    arc_spectra=np.atleast_2d(arcfitter.specdata).T,
                                    spat_ids=np.atleast_1d(arcfitter._slit),
                                    PYP_SPEC=specname,
                                    )
    else:
        waveCalib = None

    # Ask the user if they wish to store the result in PypeIt calibrations
    arcfitter.store_solution(final_fit, mdir, slits.binspec,
                             wvcalib=waveCalib,
                             rmstol=args.rmstol,
                             force_save=args.force_save)
