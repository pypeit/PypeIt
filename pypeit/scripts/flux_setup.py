#!/usr/bin/env python
import argparse
import os,time
import numpy as np
from astropy.io import fits
from astropy.table import Table
from pypeit import msgs
from pypeit import io
from pypeit.par.util import make_pypeit_file


class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)

def parse_args(options=None, return_parser=False):
    parser = argparse.ArgumentParser(description='Setup to perform flux calibration',
                                     formatter_class=SmartFormatter)
    parser.add_argument("sci_path", type=str, help="Path for Science folder")
    parser.add_argument("--objmodel", type=str, default='qso', choices=['qso', 'star', 'poly'],
                        help="R|Science object model used in the telluric fitting.\n"
                        "The options are:\n"
                        "\n"
                        "    qso  = For quasars. You might need to set redshift, bal_wv_min_mx in the tell file.\n"
                        "\n"
                        "    star  = For stars. You need to set star_type, star_ra, star_dec, and star_mag in the tell_file.\n"
                        "\n"
                        "    poly = For other type object, You might need to set fit_wv_min_mx, \n"
                        "           and norder in the tell_file."
                        )

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):
    """
      This setups PypeIt files for fluxing, coadding and telluric corrections.
      It will produce three files named as your_spectragraph.flux, your_spectragraph.coadd1d,
      and your_spectragraph.tell
    """
    allfiles = os.listdir(args.sci_path)
    allfiles = np.sort(allfiles)
    spec1dfiles = []
    spec2dfiles = []
    spec1dinfos = []
    for ifile in allfiles:
        if ('spec1d' in ifile) and ('.fits' in ifile):
            spec1dfiles.append(ifile)
        elif ('spec2d' in ifile) and ('.fits' in ifile):
            spec2dfiles.append(ifile)
        elif ('spec1d' in ifile) and ('.txt' in ifile):
            spec1dinfos.append(ifile)
        else:
            msgs.warn('{:} is not a standard PypeIt output.'.format(ifile))
    if len(spec2dfiles) > len(spec1dfiles):
        msgs.warn('The following exposures do not have 1D extractions:')
        for ii in range(len(spec2dfiles)):
            if not os.path.exists(os.path.join(args.sci_path, spec2dfiles[ii].replace('spec2d','spec1d'))):
                msgs.info('\t {:}'.format(spec2dfiles[ii]))

    if len(spec1dfiles) > 0:
        par = io.fits_open(os.path.join(args.sci_path, spec1dfiles[0]))

        ## fluxing pypeit file
        spectrograph = par[0].header['PYP_SPEC']
        pypeline = par[0].header['PYPELINE']
        flux_file = '{:}.flux'.format(spectrograph)
        cfg_lines = ['[fluxcalib]']
        cfg_lines += ['  extinct_correct = False # Set to True if your SENSFUNC derived with the UVIS algorithm\n']
        cfg_lines += ['# Please add your SENSFUNC file name below before running pypeit_flux_calib']
        make_pypeit_file(flux_file, spectrograph, spec1dfiles, cfg_lines=cfg_lines, setup_mode=True)
        fin = open(flux_file, "rt")
        data = fin.read()
        data = data.replace('spec1d_', os.path.join(args.sci_path,'spec1d_'))
        data = data.replace('data', 'flux')
        fin.close()
        fin = open(flux_file, "wt")
        fin.write(data)
        fin.close()

        ## coadd1d pypeit file
        coadd1d_file = '{:}.coadd1d'.format(spectrograph)
        cfg_lines = ['[coadd1d]']
        cfg_lines += ['  coaddfile = YOUR_OUTPUT_FILE_NAME # Please set your output file name']
        cfg_lines += ['  sensfuncfile = YOUR_SENSFUNC_FILE # Please set your SENSFUNC file name']
        if pypeline == 'Echelle':
            cfg_lines += ['  wave_method = velocity # creates a uniformly space grid in log10(lambda)\n']
        else:
            cfg_lines += ['  wave_method = linear # creates a uniformly space grid in lambda\n']

        cfg_lines += ['# This file includes all extracted objects. You need to figure out which object you want to \n'+\
                      '# coadd before running pypeit_coadd_1dspec!!!']
        spec1d_info = []
        for ii in range(len(spec1dfiles)):
            meta_tbl = Table.read(os.path.join(args.sci_path, spec1dfiles[ii]).replace('.fits', '.txt'),
                                  format='ascii.fixed_width')
            _, indx = np.unique(meta_tbl['name'],return_index=True)
            objects = meta_tbl[indx]
            for jj in range(len(objects)):
                spec1d_info.append(spec1dfiles[ii] + ' '+ objects['name'][jj])
        make_pypeit_file(coadd1d_file, spectrograph, spec1d_info, cfg_lines=cfg_lines, setup_mode=True)
        fin = open(coadd1d_file, "rt")
        data = fin.read()
        data = data.replace('spec1d_', os.path.join(args.sci_path,'spec1d_'))
        data = data.replace('data', 'coadd1d')
        fin.close()
        fin = open(coadd1d_file, "wt")
        fin.write(data)
        fin.close()

        ## tellfit pypeit file
        tellfit_file = '{:}.tell'.format(spectrograph)
        cfg_lines = ['[tellfit]']
        if args.objmodel == 'qso':
            cfg_lines += ['  objmodel = qso']
            cfg_lines += ['  redshift = 0.0']
            cfg_lines += ['  bal_wv_min_max = 10000.,11000.']
        elif args.objmodel == 'star':
            cfg_lines += ['  objmodel = star']
            cfg_lines += ['  star_type = A0']
            cfg_lines += ['  star_mag = 0.0']
        elif args.objmodel == 'poly':
            cfg_lines += ['  objmodel = poly']
            cfg_lines += ['  polyorder = 5']
            cfg_lines += ['  fit_wv_min_max = 17000.0,22000.0']

        with open(tellfit_file, 'w') as f:
            f.write('# Auto-generated PypeIt file\n')
            f.write('# {0}\n'.format(time.strftime("%a %d %b %Y %H:%M:%S", time.localtime())))
            f.write("\n")
            f.write("# User-defined execution parameters\n")
            f.write("# This is only an example. Make sure to change the following parameters accordingly.\n")
            f.write('\n'.join(cfg_lines))
            f.write('\n')
            f.write('\n')
        msgs.info('PypeIt file written to: {0}'.format(tellfit_file))


