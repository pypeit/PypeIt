"""
Setup files for flux calibration.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import os
from pathlib import Path

import numpy as np

from astropy.table import Table

from pypeit import msgs
from pypeit import io
from pypeit.scripts import scriptbase
from pypeit import inputfiles
from pypeit.spectrographs.util import load_spectrograph

def match_spec1ds_to_sensfuncs(spectrograph_name, spec1dfiles, sensfiles):
    """
    This needs a docstring
    """
    result_map = {}
    spectrograph = load_spectrograph(spectrograph_name)

    # Read configurations of each sensfile
    sens_configs = []
    for sensfile in sensfiles:
        with io.fits_open(sensfile) as hdul:
            header = hdul[0].header

        sens_configs.append({config_key: header[config_key] for config_key in spectrograph.configuration_keys() if config_key in header})

    # Read the configurations of each spec1d file, and try to find a matching sensfunc file
    for spec1dfile in spec1dfiles:        
        with io.fits_open(spec1dfile) as hdul:
            header = hdul[0].header

        spec1d_config = {config_key: header[config_key] for config_key in spectrograph.configuration_keys() if config_key in header}
        for i, sens_config in enumerate(sens_configs):
            if spectrograph.same_configuration([spec1d_config, sens_config],check_keys=False):
                result_map[spec1dfile.name] = sensfiles[i]
                break

    return result_map

class FluxSetup(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Setup configuration files to perform flux calibration, 1D coadding, and telluric correction.',
                                    width=width, formatter=scriptbase.SmartFormatter)
        parser.add_argument("paths", type=str, nargs='+', help="One or more paths for Science folders or sensitivity functions. Sensitivity functions must start with 'sens_' to be detected.")
        parser.add_argument("--name", type=str, default=None, help="The base name to use for the output files. Defaults to the instrument name is used.")
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
        name.flux, name.coadd1d, and name.tell. "name" defaults to the 
        spectrograph name but can be overriden on the command line.

        """
        allfiles = []
        for path in args.paths:
            allfiles += Path(path).iterdir()
        spec1dfiles = []
        spec2dfiles = []
        spec1dinfos = []
        unique_paths = set()
        sensfiles = []
        for ifile in allfiles:
            if ('spec1d' in ifile.name) and ('.fits' in ifile.name):
                spec1dfiles.append(ifile)
                unique_paths.add(str(ifile.parent))
            elif ('spec2d' in ifile.name) and ('.fits' in ifile.name):
                spec2dfiles.append(ifile)
            elif ('spec1d' in ifile.name) and ('.txt' in ifile.name):
                spec1dinfos.append(ifile)
            elif ifile.name.startswith('sens_') and ('.fits' in ifile.name):
                sensfiles.append(ifile)
                unique_paths.add(str(ifile.parent))
            else:
                msgs.info('{:} is not a standard PypeIt output, skipping.'.format(ifile))
        if len(spec2dfiles) > len(spec1dfiles):
            msgs.warn('The following exposures do not have 1D extractions:')
            for ii in range(len(spec2dfiles)):
                if (spec2dfiles[ii].parent / spec2dfiles[ii].name.replace("spec2d", "spec1d")).exists():
                    msgs.info('\t {:}'.format(spec2dfiles[ii]))

        if len(spec1dfiles) > 0:
            with io.fits_open(spec1dfiles[0]) as hdul:

                # Get basic configuration info from first spec1d
                spectrograph = hdul[0].header['PYP_SPEC']
                pypeline = hdul[0].header['PYPELINE']

            if args.name is None:
                output_basename = spectrograph
            else:
                output_basename = args.name            

            # Determine how to map sensitivity functions to spec1d files
            sensfile_mapping=match_spec1ds_to_sensfuncs(spectrograph, spec1dfiles, sensfiles)

            ## fluxing pypeit file
            # Build the bits and pieces
            cfg_lines = ['[fluxcalib]']
            cfg_lines += ['  extinct_correct = False # Set to True if your SENSFUNC derived with the UVIS algorithm\n']
            cfg_lines += ['# Please add your SENSFUNC file name below before running pypeit_flux_calib']
            data = Table()
            data['filename'] = [x.name for x in spec1dfiles]
            data['sensfile'] = ['' if x.name not in sensfile_mapping else sensfile_mapping[x.name].name for x in spec1dfiles]

            # Instantiate
            fluxFile = inputfiles.FluxFile(
                config=cfg_lines,
                file_paths = unique_paths, 
                data_table=data)
            # Write
            flux_file = f'{output_basename}.flux'
            fluxFile.write(flux_file)

            ## coadd1d pypeit file
            cfg_lines = ['[coadd1d]']
            cfg_lines += ['  coaddfile = YOUR_OUTPUT_FILE_NAME # Please set your output file name']
            if pypeline == 'Echelle':
                cfg_lines += ['  wave_method = velocity # creates a uniformly space grid in log10(lambda)\n']
            else:
                cfg_lines += ['  wave_method = linear # creates a uniformly space grid in lambda\n']

            cfg_lines += ['# This file includes all extracted objects. You need to figure out which object you want to \n'+\
                          '# coadd before running pypeit_coadd_1dspec!!!']
            if pypeline == 'Echelle':
                cfg_lines += ['# For Echelle spectrographs, please double check the sensfunc file and setup id\n']


            all_specfiles, all_obj = [], []
            for ii in range(len(spec1dfiles)):
                txtinfofile = spec1dfiles[ii].parent / (spec1dfiles[ii].stem + ".txt")
                meta_tbl = Table.read(txtinfofile,
                                      format='ascii.fixed_width')
                _, indx = np.unique(meta_tbl['name'],return_index=True)
                objects = meta_tbl[indx]
                for jj in range(len(objects)):
                    all_specfiles.append(spec1dfiles[ii])
                    all_obj.append(objects['name'][jj])
            data = Table()
            data['filename'] = [str(x.name) for x in all_specfiles]
            data['obj_id'] = all_obj
            if pypeline == 'Echelle':

                if len(sensfiles) > 1:
                    # If there are multiple sensfunc files, try to set sensibile values
                    # for the 'sensfile' and 'setup_id' columns
                    all_sensfiles = []
                    all_setup_ids = []
                    for spec1d in data['filename']:
                        if spec1d in sensfile_mapping:
                            sensfile = sensfile_mapping[spec1d]
                            # Use the index of the sensfile in sensfiles as the setup id,
                            # converted to a letter. Hopefully we don't get more than 26
                            setup_id = chr(ord('A') + sensfiles.index(sensfile))
                        else:
                            sensfile = 'SENSFUNC FILE'
                            setup_id = 'default'

                        all_sensfiles.append(sensfile.name)
                        all_setup_ids.append(setup_id)
                    data['sensfile'] = all_sensfiles
                    data['setup_id'] = all_setup_ids

                else:
                    # Just use one default sensfunc file and one setup.
                    if len(sensfiles) == 1:
                        default_sensfile = sensfiles[0].name # Use the first sensfunc and only sensfunc
                    else:
                        default_sensfile = 'SENSFUNC FILE' # Use a dummy sensfunc filename

                    data['sensfile'] = [default_sensfile] + ([''] * (len(all_obj)-1))
                    data['setup_id'] = ['A'] + ([''] * (len(all_obj)-1))


                
            # Instantiate
            coadd1dFile = inputfiles.Coadd1DFile(
                config=cfg_lines,
                file_paths = unique_paths, 
                data_table=data)
            # Write
            coadd1d_file = '{:}.coadd1d'.format(output_basename)
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
            tellfit_file = f'{output_basename}.tell'
            tellFile.write(tellfit_file)


