#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
Launch the identify GUI tool.
"""

from IPython import embed


def parser(options=None):
    import argparse

    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Launch PypeIt identify tool, display extracted'
                                                 'MasterArc, and load linelist.'
                                                 'Run above the Masters/ folder',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('arc_file', type=str, default=None, help='PypeIt MasterArc file')
    parser.add_argument('slits_file', type=str, default=None, help='PypeIt MasterSlits file')
    parser.add_argument("--lamps", type=str, help="Comma separated list of calibration lamps (no spaces)")
    parser.add_argument("--wmin", default=3000.0, help="Minimum wavelength range",
                        action="store_true")
    parser.add_argument("--wmax", default=10000.0, help="Maximum wavelength range",
                        action="store_true")
    parser.add_argument("--slit", default=0, help="Slit number to wavelength calibrate",
                        action="store_true")
    parser.add_argument("--det", default=1, help="Detector index", action="store_true")
    parser.add_argument("--rmstol", default=0.1, help="RMS tolerance", action="store_true")

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):

    import os
    import sys
    import astropy.io.fits as fits
    from pypeit import masterframe
    from pypeit.spectrographs.util import load_spectrograph
    from pypeit.core import parse
    from pypeit.core import gui
    from pypeit.core.wavecal import waveio, templates
    from pypeit.wavecalib import WaveCalib
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
    if args.lamps == '':
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
    wv_calib = waveio.load_wavelength_calibration(solnname) if os.path.exists(solnname) else None

    # Load the MasterFrame (if it exists and is desired)?
    wavecal = WaveCalib(msarc, slits, spec, par, binspectral=slits.binspec, det=args.det,
                        master_key=mkey, msbpm=msarc.fullmask)
    arccen, arc_maskslit = wavecal.extract_arcs()

    # Launch the identify window
    arcfitter = gui.identify.initialise(arccen, slit=int(args.slit), par=par, wv_calib_all=wv_calib,
                                        wavelim=[args.wmin, args.wmax],
                                        nonlinear_counts=spec.nonlinear_counts(msarc.detector))
    final_fit = arcfitter.get_results()

    # Ask the user if they wish to store the result in PypeIt calibrations
    if final_fit['rms'] < args.rmstol:
        ans = ''
        while ans != 'y' and ans != 'n':
            ans = input("Would you like to store this wavelength solution in the archive? (y/n): ")
        if ans == 'y':
            gratname = fits.getheader(msarc.head0['F1'])[spec.meta['dispname']['card']].replace("/", "_")
            dispangl = "UNKNOWN"
            templates.pypeit_identify_record(final_fit, slits.binspec, specname, gratname, dispangl)
            print("Your wavelength solution has been stored")
            print("Please consider sending your solution to the PypeIt team!")
    else:
        print("Final fit RMS: {0:0.3f} is larger than the allowed tolerance: {1:0.3f}".format(final_fit['rms'], args.rmstol))
        print("Set the variable --rmstol on the command line to allow a more flexible RMS tolerance")
        ans = ''
        while ans != 'y' and ans != 'n':
            ans = input("Would you like to store the line IDs? (y/n): ")
        if ans == 'y':
            arcfitter.save_IDs()
