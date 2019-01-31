#!/usr/bin/env python

"""
Script for fluxing PYPEIT 1d spectra
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse

# pypeit_flux_spec sensfunc --std_file=spec1d_Feige66_KASTb_2015May20T041246.96.fits  --instr=shane_kast_blue --sensfunc_file=tmp.yaml
# pypeit_flux_spec flux --sci_file=spec1d_J1217p3905_KASTb_2015May20T045733.56.fits --sensfunc_file=tmp.yaml --flux_file=tmp.fits
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# pypeit_flux_spec sensfunc --std_file=spec1d_G191b2b_KASTr_2015Jan23T024320.75.fits  --instr=shane_kast_red_ret --sensfunc_file=tmp.yaml --telluric
# pypeit_flux_spec flux     --sci_file=spec1d_J0025-0312_KASTr_2015gen23T025323.85.fits --sensfunc_file=tmp.yaml --flux_file=tmp.fits
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Echelle examples:
## Generate sensfunc
# pypeit_flux_spec sensfunc keck_nires --std_file=spec1d_HIP13917_V8p6_NIRES_2018Oct01T094225.598.fits
#         --sensfunc_file=spec1d_HIP13917_V8p6_NIRES.yaml --telluric --echelle --star_type A0 --star_mag 8.6 --debug
## flux calibrate your science.
# pypeit_flux_spec flux keck_nires --sci_file=spec1d_J0252-0503_NIRES_2018Oct01T100254.698.fits
#         --sensfunc_file=spec1d_HIP13917_V8p6_NIRES.yaml
#         --flux_file=spec1d_J0252-0503_NIRES_2018Oct01T100254.698_flux.fits --echelle


def parser(options=None):
    parser = argparse.ArgumentParser(description='Parse')
    parser.add_argument("steps", type=str, help="Steps to perform [sensfunc,flux]")
    parser.add_argument("spectrograph", type=str, help="Spectrograph")
    parser.add_argument("--std_file", type=str, help="File containing the standard 1d spectrum")
    parser.add_argument("--std_obj", type=str, help="Standard star identifier, e.g. O479-S5009-D01-I0023")
    parser.add_argument("--sci_file", type=str, help="File containing the science 1d spectra")
    parser.add_argument("--sensfunc_file", type=str, help="YAML file containing the sensitivity function (input or output)")
    parser.add_argument("--flux_file", type=str, help="Output filename for fluxed science spectra")
    parser.add_argument("--plot", default=False, action="store_true", help="Show the sensitivity function?")
    parser.add_argument("--multi_det", type=str, help="Multiple detectors (e.g. 3,7 for DEIMOS)")
    parser.add_argument("--telluric", default=False, action="store_true", help="Correct for telluric absorptions?")
    parser.add_argument("--echelle", default=False, action="store_true", help="Echelle spectra?")
    parser.add_argument("--star_type", type=str, help="Which type of yur telluric/standard star?")
    parser.add_argument("--star_mag", type=float, help="What's the V-band magnitude of your telluric/standard star?")
    parser.add_argument("--debug", default=False, action="store_true", help="show debug plots?")
    parser.add_argument("--param_file", type=str, help="Parameter file for fluxing.. (config)")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args, unit_test=False):
    """ Runs fluxing steps
    """
    import pdb

    from pypeit import fluxspec
    from pypeit.spectrographs.util import load_spectrograph
    from pypeit.par import pypeitpar
    from pypeit.par import util

    # Parse the steps
    steps = args.steps.split(',')

    spectrograph = load_spectrograph(args.spectrograph)
    spectrograph_def_par = spectrograph.default_pypeit_par()

    # Extra parameters?
    if args.param_file is not None:
        lines = util._read_pypeit_file_lines(args.param_file)
        par = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_def_par.to_config(), merge_with=list(lines))
    else:
        par = spectrograph_def_par


    # Multi detector?
    if args.multi_det is not None:
        multi_det = [int(det) for det in args.multi_det.split(',')]
    else:
        multi_det = None

    # Checks
    if 'sensfunc' in steps:
        #if args.instr is None:
        #    raise IOError("You must set the instrument to generate the sensfunc")
        if args.std_file is None:
            raise IOError("You must input a spec1d file of the standard to generate the sensfunc")
        if args.sensfunc_file is None:
            raise IOError("You must give the output filename in --sensfunc_file to generate the sensfunc")
    if 'flux' in steps:
        if args.sci_file is None:
            raise IOError("You must input a spec1d file of the science spectra to flux them")
        if args.flux_file is None:
            raise IOError("You must give the output filename in --flux_file for the fluxed spectra")
        if 'sensfunc' not in steps:
            if args.sensfunc_file is None:
                raise IOError("You must input a sensfunc in -sensfunc_file to flux.")

    # Instantiate
    if 'sensfunc' in steps:
        sfile = None  # Need to create it, not load it
    else:
        sfile = args.sensfunc_file

    if args.echelle:
        FxSpec = fluxspec.EchFluxSpec(spectrograph,
                                      std_spec1d_file=args.std_file,
                                   sci_spec1d_file=args.sci_file,
                                   telluric=args.telluric,
                                   sens_file=sfile,
                                   star_type=args.star_type,
                                   star_mag=args.star_mag,
                                   debug = args.debug)
    else:
        FxSpec = fluxspec.FluxSpec(spectrograph, par['fluxcalib'],
                                   std_spec1d_file=args.std_file,
                                   sci_spec1d_file=args.sci_file,
                                   telluric=args.telluric,
                                   sens_file=sfile,
                                   multi_det=multi_det,
                                   debug = args.debug)
    # Step through
    if 'sensfunc' in steps:
        # For echelle, the code will deal with the standard star in the ech_fluxspec.py
        if not args.echelle:
            # Find the star automatically?
            if args.std_obj is None:
                _ = FxSpec.find_standard()
            else:
                _ = FxSpec._set_std_obj(args.std_obj)
        # Sensitivity
        _ = FxSpec.generate_sensfunc()
        # Output
        _ = FxSpec.save_master(FxSpec.sens_dict, outfile=args.sensfunc_file)
        # Show
        if args.plot:
            FxSpec.show_sensfunc()

    if 'flux' in steps:
        FxSpec.flux_science()
        FxSpec.write_science(args.flux_file)



