"""
Setup files for flux calibration.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import os
import time

import numpy as np

from astropy.table import Table

from pypeit import msgs
from pypeit import io
from pypeit.scripts import scriptbase
from pypeit import inputfiles

# TODO -- We need a test of this script

class FluxSetup(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Setup to perform flux calibration',
                                    width=width, formatter=scriptbase.SmartFormatter)
        parser.add_argument("sci_path", type=str, help="Path for Science folder")

        parser.add_argument("--objmodel", type=str, default='qso', choices=['qso', 'star', 'poly'],
                            help='R|science object model used in the telluric fitting. The '
                                 'options are:\n\n'
                                 'qso = For quasars. You might need to set redshift, '
                                 'bal_wv_min_max in the tell file.\n'
                                 '\n'
                                 'star = For stars. You need to set star_type, star_ra, star_dec, '
                                 'and star_mag in the tell_file.\n'
                                 '\n'
                                 'poly = For other type object, You might need to set '
                                 'fit_wv_min_max, and norder in the tell_file.\n'
                                 '\n')
        return parser

    @staticmethod
    def main(args):
        """
        This setups PypeIt input files for fluxing, coadding, and telluric
        corrections.  It will produce three files named as
        your_spectragraph.flux, your_spectragraph.coadd1d, and
        your_spectragraph.tell.
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
                if not os.path.exists(os.path.join(args.sci_path,
                                                   spec2dfiles[ii].replace('spec2d','spec1d'))):
                    msgs.info('\t {:}'.format(spec2dfiles[ii]))

        if len(spec1dfiles) > 0:
            par = io.fits_open(os.path.join(
                args.sci_path, spec1dfiles[0]))

            ## fluxing pypeit file
            spectrograph = par[0].header['PYP_SPEC']
            pypeline = par[0].header['PYPELINE']

            # Build the bits and pieces
            cfg_lines = ['[fluxcalib]']
            cfg_lines += ['  extinct_correct = False # Set to True if your SENSFUNC derived with the UVIS algorithm\n']
            cfg_lines += ['# Please add your SENSFUNC file name below before running pypeit_flux_calib']
            data = Table()
            data['filename'] = spec1dfiles
            data['sensfile'] = ''

            # Instantiate
            fluxFile = inputfiles.FluxFile(
                config=cfg_lines,
                file_paths = [args.sci_path], 
                data_table=data)
            # Write
            flux_file = f'{spectrograph}.flux'
            fluxFile.write(flux_file)

            ## coadd1d pypeit file
            cfg_lines = ['[coadd1d]']
            cfg_lines += ['  coaddfile = YOUR_OUTPUT_FILE_NAME # Please set your output file name']
            cfg_lines += ['  sensfuncfile = YOUR_SENSFUNC_FILE # Please set your SENSFUNC file name. Only required for Echelle']
            if pypeline == 'Echelle':
                cfg_lines += ['  wave_method = velocity # creates a uniformly space grid in log10(lambda)\n']
            else:
                cfg_lines += ['  wave_method = linear # creates a uniformly space grid in lambda\n']

            cfg_lines += ['# This file includes all extracted objects. You need to figure out which object you want to \n'+\
                          '# coadd before running pypeit_coadd_1dspec!!!']


            all_specfiles, all_obj = [], []
            for ii in range(len(spec1dfiles)):
                meta_tbl = Table.read(os.path.join(args.sci_path, spec1dfiles[ii]).replace('.fits', '.txt'),
                                      format='ascii.fixed_width')
                _, indx = np.unique(meta_tbl['name'],return_index=True)
                objects = meta_tbl[indx]
                for jj in range(len(objects)):
                    all_specfiles.append(spec1dfiles[ii])
                    all_obj.append(objects['name'][jj])
            data = Table()
            data['filename'] = all_specfiles
            data['obj_id'] = all_obj
            # Instantiate
            coadd1dFile = inputfiles.Coadd1DFile(
                config=cfg_lines,
                file_paths = [args.sci_path], 
                data_table=data)
            # Write
            coadd1d_file = '{:}.coadd1d'.format(spectrograph)
            coadd1dFile.write(coadd1d_file)

            ## tellfit pypeit file
            cfg_lines = ['[telluric]']
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
            # Instantiate
            tellFile = inputfiles.TelluricFile(
                config=cfg_lines)
            # Write
            tellfit_file = f'{spectrograph}.tell'
            tellFile.write(tellfit_file)



def write_input_file(outfile, files, 
                     cfg_lines=None, file_block='data', 
                     paths=None):
    """
    Generate a default PypeIt file

    Args:
        outfile (str): Name of output file to be generated
        spectrograph (str):  Name of spectrograph
        files (list):  List of data files -- essentially Deprecated
        file_block (str): Type of file block to generate ['data', 'flux']
        cfg_lines (list, optional):  List of configuration lines for parameters
        paths (list, optional): List of paths for slurping data files
    """
    # Error checking
    if not isinstance(files, list):
        msgs.error("files needs to be a list")

    # Defaults
    if cfg_lines is None:
        _cfg_lines = ['[rdx]']
        _cfg_lines += ['    spectrograph = {0}'.format(spectrograph)]
    else:
        _cfg_lines = list(cfg_lines)

    # Here we go
    with open(outfile, 'w') as f:
        f.write('# Auto-generated PypeIt input file using PypeIt version: \n')
        #f.write('# {0}\n'.format(time.strftime("%a %d %b %Y %H:%M:%S",time.localtime())))
        f.write('# {0}\n'.format(time.strftime("%Y-%m-%d",time.localtime())))
        f.write("\n")
        f.write("# User-defined execution parameters\n")
        f.write('\n'.join(_cfg_lines))
        f.write('\n')
        f.write('\n')

        # Data
        if file_block == 'data':
            f.write("# Read in the data\n")
        elif file_block == 'flux':
            f.write("# Files to flux \n")
        elif file_block == 'coadd1d':
            f.write("# Files to coadd \n")
        else:
            msgs.error("Bad file block name")
        f.write(f"{file_block} read\n")
        # Old school
        for datafile in files:
            f.write(' '+datafile+'\n')
        # paths and Setupfiles
        if paths is not None:
            for path in paths:
                f.write(' path '+path+'\n')
        #if sorted_files is not None:
        #    f.write('\n'.join(sorted_files))
        #    f.write('\n')
        f.write(f"{file_block} end\n")
        f.write("\n")

    msgs.info('PypeIt input file written to: {0}'.format(outfile))


