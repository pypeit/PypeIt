#!/usr/bin/env python

"""
Script for performing 2d coadds of PypeIt data.
"""
from configobj import ConfigObj
from astropy.io import fits
from pypeit import par, msgs
from pypeit.core import coadd2d
from pypeit.core import save
from pypeit.spectrographs.util import load_spectrograph
from collections import OrderedDict
import argparse
import glob
import os
import numpy as np
from pypeit.par import pypeitpar


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Echelle examples:
## Generate sensfunc
# pypeit_flux_spec sensfunc keck_nires --std_file=spec1d_HIP13917_V8p6_NIRES_2018Oct01T094225.598.fits
#         --sensfunc_file=spec1d_HIP13917_V8p6_NIRES.yaml --telluric --echelle --star_type A0 --star_mag 8.6 --debug
## flux calibrate your science.
# pypeit_flux_spec flux keck_nires --sci_file=spec1d_J0252-0503_NIRES_2018Oct01T100254.698.fits
#         --sensfunc_file=spec1d_HIP13917_V8p6_NIRES.yaml
#         --flux_file=spec1d_J0252-0503_NIRES_2018Oct01T100254.698_flux.fits --echelle


def read_coadd2d_file(ifile):
    """
    Read a PypeIt coadd2d file, akin to a standard PypeIt file

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
    msgs.info('Loading the coadd2d file')
    lines = par.util._read_pypeit_file_lines(ifile)
    is_config = np.ones(len(lines), dtype=bool)

    # Parse the fluxing block
    spec2d_files = []
    s, e = par.util._find_pypeit_block(lines, 'coadd2d')
    if s >= 0 and e < 0:
        msgs.error("Missing 'coadd2d end' in {0}".format(ifile))
    else:
        for line in lines[s:e]:
            prs = line.split(' ')
            spec2d_files.append(os.path.join('Science/', os.path.basename(prs[0])))
        is_config[s-1:e+1] = False
    #elif (s < 0) or (s==e):
    #    msgs.warn("No spec2d file block, you must have passed --obj as an argument..")

    # Construct config to get spectrograph
    cfg_lines = list(lines[is_config])
    cfg = ConfigObj(cfg_lines)
    spectrograph_name = cfg['rdx']['spectrograph']
    spectrograph = load_spectrograph(spectrograph_name)

    # Return
    return spectrograph, cfg_lines, spec2d_files


# TODO this is an exact copy of the same function in pypeit class. We should probably make it a utility in par??
def select_detectors(par, spectrograph):
    """
    Return the 1-indexed list of detectors to reduce.

    Returns:
    list:  List of detectors to be reduced

    """
    if par['rdx']['detnum'] is None:
        return np.arange(spectrograph.ndet)+1
    return [par['rdx']['detnum']] if isinstance(par['rdx']['detnum'], int) else par['rdx']['detnum']


def parser(options=None):
    parser = argparse.ArgumentParser(description='Parse')
    parser.add_argument("--file", type=str, default=None, help="File to guide 2d coadds")
    parser.add_argument('--det', default=1, type=int, help="Only coadd this detector number")
    parser.add_argument("--obj", type=str, default=None,
                        help="Object name in lieu of extension, e.g if the spec2d files are named "
                             "'spec2d_J1234+5678_GNIRS_2017Mar31T085412.181.fits. then obj=J1234+5678")
    parser.add_argument("--show", default=False, action="store_true",
                        help="Show the reduction steps. Equivalent to the -s option when running pypeit.")
    parser.add_argument("--peaks", default=False, action="store_true",
                        help="Show the peaks found by the object finding algorithm.")
    parser.add_argument("--par_outfile", default='coadd2d.par', action="store_true",
                        help="Output file to save the parameters")
    parser.add_argument("--debug", default=False, action="store_true", help="show debug plots?")


    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args):
    """ Executes 2d coadding
    """

    # Load the file
    if args.file is not None:
        spectrograph, config_lines, spec2d_files = read_coadd2d_file(args.file)
        # Parameters
        spectrograph_def_par = spectrograph.default_pypeit_par()
        par = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_def_par.to_config(),
                                                 merge_with=config_lines)
    elif args.obj is not None:
        spec2d_files = glob.glob('./Science/spec2d_' + args.obj + '*')
        head0 = fits.getheader(spec2d_files[0])
        spectrograph_name = head0['SPECTROG']
        spectrograph = load_spectrograph(spectrograph_name)
        par = spectrograph.default_pypeit_par()
    else:
        msgs.error('You must either input a coadd2d file with --file or an object name with --obj')

    # If detector was passed as an argument override whatever was in the coadd2d_file
    if args.det is not None:
        msgs.info("Restricting reductions to detector={}".format(args.det))
        par['rdx']['detnum'] = int(args.det)

    # Write the par to disk
    print("Writing the parameters to {}".format(args.par_outfile))
    par.to_config(args.par_outfile)

    # Now run the coadds
    spec1d_files = [files.replace('spec2d', 'spec1d') for files in spec2d_files]
    head1d = fits.getheader(spec1d_files[0])
    head2d = fits.getheader(spec2d_files[0])
    filename = os.path.basename(spec2d_files[0])
    basename = filename.split('_')[1]

    skysub_mode = head2d['SKYSUB']
    ir_redux = True if 'DIFF' in skysub_mode else False

    # Print status message
    msgs_string = 'Reducing target {:s}'.format(basename) + msgs.newline()
    msgs_string += 'Performing coadd of frames reduce with {:s} imaging'.format(skysub_mode) + msgs.newline()
    msgs_string += 'Combining frames in 2d coadd:' + msgs.newline()
    for file in spec2d_files:
        msgs_string += '{0:s}'.format(os.path.basename(file)) + msgs.newline()
    msgs.info(msgs_string)

    redux_path = './'
    master_dirname = os.path.basename(head2d['PYPMFDIR'])+'/'
    master_dir = os.path.join(redux_path,os.path.normpath(master_dirname) + '_coadd/')

    # Make the new master dir and Science dir
    if os.path.exists(master_dir):
        msgs.info("The following directory already exists:"+msgs.newline()+master_dir)
    else:
        os.mkdir(master_dir)

    # Instantiate the sci_dict
    sci_dict = OrderedDict()  # This needs to be ordered
    sci_dict['meta'] = {}
    sci_dict['meta']['vel_corr'] = 0.
    sci_dict['meta']['ir_redux'] = ir_redux

    # Find the detectors to reduce
    detectors = select_detectors(par, spectrograph)
    if len(detectors) != spectrograph.ndet:
        msgs.warn('Not reducing detectors: {0}'.format(' '.join([str(d) for d in
        set(np.arange(spectrograph.ndet)) - set(detectors)])))

    # Loop on detectors
    for det in detectors:
        msgs.info("Working on detector {0}".format(det))
        sci_dict[det] = {}

        # Read in the images stacks and other clibration/meta data for this detector
        stack_dict = coadd2d.load_coadd2d_stacks(spec2d_files, det)

        sci_dict[det]['sciimg'], sci_dict[det]['sciivar'], sci_dict[det]['skymodel'], \
        sci_dict[det]['objmodel'], sci_dict[det]['ivarmodel'], sci_dict[det]['outmask'], \
        sci_dict[det]['specobjs'] = coadd2d.extract_coadd2d(stack_dict, master_dir, ir_redux=ir_redux, par=par,
                                                            show=args.show, show_peaks=args.peaks)

    # Make the science directory, and write outputs to disk
    master_key_dict = stack_dict['master_key_dict']
    scipath = redux_path + 'Science_coadd'
    if os.path.exists(scipath):
        msgs.info("The following directory already exists:" + msgs.newline() + scipath)
    else:
        os.mkdir(scipath)
    save.save_all(sci_dict, master_key_dict, master_dir, spectrograph, head1d, head2d, scipath, basename)



    # I don't think I need this below
    #try:
    #    mjd = head1d['mjd']  # recorded as 'mjd' in fitstbl
    #except KeyError:
    #    mjd = head1d['MJD-OBS']
    #obstime = time.Time(mjd, format='mjd')


