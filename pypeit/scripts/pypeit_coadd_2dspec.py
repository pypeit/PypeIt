#!/usr/bin/env python

"""
Script for fluxing PYPEIT 1d spectra
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from configobj import ConfigObj
import numpy as np
from pypeit import par, msgs
from pypeit.spectrographs.util import load_spectrograph
import argparse

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Echelle examples:
## Generate sensfunc
# pypeit_flux_spec sensfunc keck_nires --std_file=spec1d_HIP13917_V8p6_NIRES_2018Oct01T094225.598.fits
#         --sensfunc_file=spec1d_HIP13917_V8p6_NIRES.yaml --telluric --echelle --star_type A0 --star_mag 8.6 --debug
## flux calibrate your science.
# pypeit_flux_spec flux keck_nires --sci_file=spec1d_J0252-0503_NIRES_2018Oct01T100254.698.fits
#         --sensfunc_file=spec1d_HIP13917_V8p6_NIRES.yaml
#         --flux_file=spec1d_J0252-0503_NIRES_2018Oct01T100254.698_flux.fits --echelle


def read_fluxfile(ifile):
    """
    Read a PypeIt flux file, akin to a standard PypeIt file

    The top is a config block that sets ParSet parameters
      The spectrograph is required

    Args:
        ifile: str
          Name of the flux file

    Returns:
        spectrograph: Spectrograph
        cfg_lines: list
          Config lines to modify ParSet values
        flux_dict: dict
          Contains spec1d_files and flux_files
          Empty if no flux block is specified

    """
    # Read in the pypeit reduction file
    msgs.info('Loading the fluxcalib file')
    lines = par.util._read_pypeit_file_lines(ifile)
    is_config = np.ones(len(lines), dtype=bool)


    # Parse the fluxing block
    flux_dict = {}
    s, e = par.util._find_pypeit_block(lines, 'flux')
    if s >= 0 and e < 0:
        msgs.error("Missing 'flux end' in {0}".format(ifile))
    elif (s < 0) or (s==e):
        msgs.warn("No flux block, you must be making the sensfunc only..")
    else:
        flux_dict['spec1d_files'] = []
        flux_dict['flux_files'] = []
        for line in lines[s:e]:
            prs = line.split(' ')
            flux_dict['spec1d_files'].append(prs[0])
            flux_dict['flux_files'].append(prs[1])
        is_config[s-1:e+1] = False

    # Construct config to get spectrograph
    cfg_lines = list(lines[is_config])
    cfg = ConfigObj(cfg_lines)
    spectrograph_name = cfg['rdx']['spectrograph']
    spectrograph = load_spectrograph(spectrograph_name)

    # Return
    return spectrograph, cfg_lines, flux_dict

def parser(options=None):
    parser = argparse.ArgumentParser(description='Parse')
    parser.add_argument("flux_file", type=str, help="File to guide fluxing process")
    parser.add_argument("--debug", default=False, action="store_true", help="show debug plots?")
    parser.add_argument("--plot", default=False, action="store_true", help="Show the sensitivity function?")
    parser.add_argument("--par_outfile", default='fluxing.par', action="store_true", help="Output to save the parameters")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args, unit_test=False):
    """ Runs fluxing steps
    """
    import pdb
    import os
    import numpy as np

    from pypeit import fluxspec
    from pypeit.core import flux
    from pypeit.par import pypeitpar

    # Load the file
    spectrograph, config_lines, flux_dict = read_fluxfile(args.flux_file)

    # Parameters
    spectrograph_def_par = spectrograph.default_pypeit_par()
    par = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_def_par.to_config(),
                                             merge_with=config_lines)

    if unit_test:
        path = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Science')
        par['fluxcalib']['std_file'] = os.path.join(path, par['fluxcalib']['std_file'])
        for kk, spec1d_file, flux_file in zip(np.arange(len(flux_dict['spec1d_files'])), flux_dict['spec1d_files'], flux_dict['flux_files']):
            flux_dict['spec1d_files'][kk] = os.path.join(path, spec1d_file)
            flux_dict['flux_files'][kk] = os.path.join(path, flux_file)

    # Write the par to disk
    print("Writing the parameters to {}".format(args.par_outfile))
    par.to_config(args.par_outfile)

    # Instantiate
    if spectrograph.pypeline == 'Echelle':
        # THIS MAY BE BROKEN
        FxSpec = fluxspec.EchFluxSpec(spectrograph, par['fluxcalib'], debug=args.debug)
    else:
        FxSpec = fluxspec.FluxSpec(spectrograph, par['fluxcalib'], debug=args.debug)

    # Generate sensfunc??
    if par['fluxcalib']['std_file'] is not None:
        # Load standard
        FxSpec.load_objs(par['fluxcalib']['std_file'], std=True)
        # For echelle, the code will deal with the standard star in the ech_fluxspec.py
        if not spectrograph.pypeline == 'Echelle':
            # Find the star
            _ = FxSpec.find_standard()
        # Sensitivity
        _ = FxSpec.generate_sensfunc()
        # Output
        _ = FxSpec.save_sens_dict(FxSpec.sens_dict, outfile=par['fluxcalib']['sensfunc'])
        # Show
        if args.plot:
            FxSpec.show_sensfunc()

    # Flux?
    if len(flux_dict) > 0:
        for spec1d_file, flux_file in zip(flux_dict['spec1d_files'], flux_dict['flux_files']):
            FxSpec.flux_science(spec1d_file)
            FxSpec.write_science(flux_file)





