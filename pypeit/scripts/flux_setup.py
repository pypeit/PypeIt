"""
Setup files for flux calibration.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import os
from pathlib import Path
from datetime import datetime, timedelta

import numpy as np

from astropy.table import Table
from astropy.time import Time

from pypeit import log
from pypeit import io
from pypeit.scripts import scriptbase
from pypeit import inputfiles
from pypeit.spectrographs.util import load_spectrograph
from IPython import embed


def parse_date_from_filename(path):
    name_wo_ext = path.stem
    datetime_str = name_wo_ext.split("_")[-1]
    try:
        result = datetime.fromisoformat(datetime_str)
    except Exception as e:
        log.warning(f"Failed to parse date from {path.name}")
        raise
    return result

def read_config_from_pypeit_fits(spectrograph, pypeit_fits_file):

    with io.fits_open(pypeit_fits_file) as hdul:
        header = hdul[0].header
        return {config_key: header[config_key] for config_key in spectrograph.configuration_keys() if config_key in header}



def read_metadata(spectrograph, pypeit_file_names, spec1dfiles):

    colnames = ['dirname', 'filename', 'frametype', 'mjd', 'setup_id'] + spectrograph.configuration_keys()
    data = []

    pypeit_files = [inputfiles.PypeItFile.from_file(file) for file in pypeit_file_names]

    # Find the PypeIt files for these spec1ds
    for spec1dfile in spec1dfiles:
        found = False
        for i, pf in enumerate(pypeit_files):
            for row in pf.data:
                if Path(row['filename']).name.split('.fits')[0] in spec1dfile.name:
                    spec1d_config = read_config_from_pypeit_fits(spectrograph, spec1dfile)
                    data.append([str(spec1dfile.parent),spec1dfile.name,str(row['frametype']),str(row['mjd']), chr(ord('A') + i)] + [spec1d_config[key] for key in spectrograph.configuration_keys()])
                    found = True
                    break
            if found:
                break
        if not found:
            raise ValueError(f"Failed to find .pypeit file metadata for {spec1dfile.name}")

    metadata_table = Table(rows=data, names=colnames)
    return metadata_table

def match_spec1ds_to_sensfuncs(args, spectrograph, metadata, sensfiles):
    """
    Match sensitivity functions to spec1d files for fluxing.

    Args:
        spectrograph_name (str): Name of the spectrograph being used
        spec1dfiles (list of Path): list of names of spec1d files to match to sensitivity function files.
        specobjs_list (list of :obj:`pypeit.specobjs.SpecObjs`): List of SpecObjs objects corresponding to each spec1d file.
        sensfiles (list of Path): Lis tof names of sensitivity funciton files to match to spec1d files.

    Return (list of tuple):
        A 1-to-1 map of spec1d filenames to sensitivity function files. Each file is returend as a tuple of its
        filename and setup_id. If a sensitivity function file could not be found for a spec1d file, that spec1d file
        is not returned in the results and a warning is logged.
    """
    spec1d_to_sensfile = {}
    # Read observing configurations of each sensfunc file
    sens_configs = []
    for sensfile in sensfiles:
        sens_configs.append(read_config_from_pypeit_fits(spectrograph, sensfile))

    # Read the configurations of each spec1d file, and try to find a matching sensfunc file.
    for metadata_row in metadata:

        if args.skip_standards:
            if metadata_row['frametype'] == 'standard':
                continue

        matching_sensfiles = []
        spec1d_config = {config_key: metadata_row[config_key] for config_key in spectrograph.configuration_keys() if config_key in metadata.colnames}

        # Go through all of the sensfunc configs and find a match for the spec1d
        for j, sens_config in enumerate(sens_configs):

            if spectrograph.same_configuration([spec1d_config, sens_config],check_keys=False):
                matching_sensfiles.append(sensfiles[j])
                
        if len(matching_sensfiles) == 0:
            log.warning(f"{metadata_row['filename']} does not have any matching sensitivity functions.")
        elif len(matching_sensfiles) == 1:
            spec1d_to_sensfile[metadata_row['filename']] = matching_sensfiles[0]
        else:
            # Pick between multiple sensitivity functions by using the closest date
            
            try:
                # Convert mjd from metadata to a python datetime
                spec1d_date = Time(metadata_row['mjd'],format='mjd').to_datetime()
            
                # Get time differences between the dates of each each sens file and the spec1d file
                time_differences = [(sensfile, abs(parse_date_from_filename(sensfile)-spec1d_date)) for sensfile in matching_sensfiles]
                
                # Select the file with the smallest time difference
                spec1d_to_sensfile[metadata_row['filename']] = min(time_differences, key=lambda x: x[1])[0]

            except Exception as e:
                log.warning(f"Could not compare sensfunc files by date. Using first found sensfile for {metadata_row['filename']}")
                spec1d_to_sensfile[metadata_row['filename']] = matching_sensfiles[0]
            
    return spec1d_to_sensfile

class FluxSetup(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Setup configuration files to perform flux calibration, 1D coadding, and telluric correction.',
                                    width=width, formatter=scriptbase.SmartFormatter)
        parser.add_argument("paths", type=str, nargs='+', help="One or more paths for Science folders or sensitivity functions. Sensitivity functions must start with 'sens_' to be detected.")
        parser.add_argument("--recursive","-r", default=False, action="store_true", help="Whether to recursively search subdirectories of the given paths.")
        parser.add_argument("--name", type=str, default=None, help="The base name to use for the output files. Defaults to the instrument name is used.")
        parser.add_argument("--coadd_output", "-o", type=str, default="YOUR_OUTPUT_FILE_NAME", help="Output file for the 1D coadfding, placed into the coadd1d configuration file.")
        parser.add_argument("--skip_standards", default=False, action="store_true", help="If set, the flux calibration file will not contain any standard star observations.")
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

    @classmethod
    def main(cls, args):
        """
        This setups PypeIt input files for fluxing, coadding, and telluric
        corrections.  It will produce three files named as
        name.flux, name.coadd1d, and name.tell. "name" defaults to the 
        spectrograph name but can be overriden on the command line.

        """
        # Initialize the log
        cls.init_log(args)
        allpaths = []
        for path in args.paths:
            allpaths.append(Path(path))
        spec1dfiles = []
        spec2dfiles = []
        spec1dinfos = []
        pypeit_file_names = []
        unique_paths = set()
        sensfiles = []
        for path in allpaths:
            for ifile in path.iterdir():
                if ifile.is_dir():
                    if args.recursive:
                        allpaths.append(ifile)
                elif ifile.name.startswith('spec1d') and ifile.name.endswith('.fits'):
                    spec1dfiles.append(ifile)
                    unique_paths.add(str(ifile.parent))
                elif ifile.name.startswith('spec2d') and ifile.name.endswith('.fits'):
                    spec2dfiles.append(ifile)
                elif ifile.name.startswith('spec1d') and ifile.name.endswith('.txt'):
                    spec1dinfos.append(ifile)
                elif ifile.name.startswith('sens_') and ('.fits' in ifile.name):
                    sensfiles.append(ifile)
                    unique_paths.add(str(ifile.parent))
                elif ifile.suffix == ".pypeit":
                    pypeit_file_names.append(ifile)
                else:
                    # log.info('{:} is not a standard PypeIt output, skipping.'.format(ifile))
                    continue

        if len(spec2dfiles) > len(spec1dfiles):
            log.warning('The following exposures do not have 1D extractions:')
            for ii in range(len(spec2dfiles)):
                if not (spec2dfiles[ii].parent / spec2dfiles[ii].name.replace("spec2d", "spec1d")).exists():
                    log.info('\t {:}'.format(spec2dfiles[ii]))

        if len(spec1dfiles) > 0:
            with io.fits_open(str(spec1dfiles[0])) as hdul:
                spectrograph_name = hdul[0].header['PYP_SPEC']
            spectrograph = load_spectrograph(spectrograph_name)
            pypeline = spectrograph.pypeline

            metadata_table = read_metadata(spectrograph, pypeit_file_names, spec1dfiles)
           
            # Get basic configuration info from first spec1d

            if args.name is None:
                output_basename = spectrograph_name
            else:
                output_basename = args.name            

            # Determine how to map sensitivity functions to spec1d files
            spec1d_to_sensfile_map=match_spec1ds_to_sensfuncs(args, spectrograph, metadata_table, sensfiles)

            ## fluxing pypeit file
            # Build the bits and pieces
            cfg_lines = ['[fluxcalib]']
            cfg_lines += ['  extinct_correct = False # Set to True if your SENSFUNC derived with the UVIS algorithm\n']
            cfg_lines += ['# Please add your SENSFUNC file name below before running pypeit_flux_calib']
            data = Table()
            # Data list, excluding standard stars if requested and using a blank sensfile name if no sensfile was found
            data['filename'] = [metadata_row['filename'] for metadata_row in metadata_table if (args.skip_standards is False or metadata_row['frametype']=='science')]
            data['sensfile'] = ['' if spec1dname not in spec1d_to_sensfile_map else spec1d_to_sensfile_map[spec1dname].name for spec1dname in data['filename']]

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
            cfg_lines += [f'  coaddfile = {args.coadd_output} # Please set your output file name']
            if pypeline == 'Echelle':
                cfg_lines += ['  wave_method = velocity # creates a uniformly space grid in log10(lambda)\n']
            else:
                cfg_lines += ['  wave_method = linear # creates a uniformly space grid in lambda\n']

            cfg_lines += ['# This file includes all extracted objects. You need to figure out which object you want to \n'+\
                          '# coadd before running pypeit_coadd_1dspec!!!']
            if pypeline == 'Echelle':
                cfg_lines += ['# For Echelle spectrographs, please double check the sensfunc file and setup id\n']


            all_specfiles, all_obj, all_setup_ids = [], [], []
            metadata_table.sort('filename')
            for metadata_row in metadata_table:
                if 'standard' in metadata_row['frametype']:
                    # Skip coadding standard frames
                    continue
                txtinfofile = Path(metadata_row['dirname'], Path(metadata_row['filename']).stem + ".txt")
                meta_tbl = Table.read(txtinfofile,
                                      format='ascii.fixed_width')
                _, indx = np.unique(meta_tbl['name'],return_index=True)
                objects = meta_tbl[indx]
                for jj in range(len(objects)):
                    all_specfiles.append(metadata_row['filename'])
                    all_setup_ids.append(metadata_row['setup_id'])
                    all_obj.append(objects['name'][jj])
            data = Table()
            data['filename'] = all_specfiles
            data['obj_id'] = all_obj
            if pypeline == 'Echelle':

                if len(sensfiles) > 1:
                    # If there are multiple sensfunc files, try to set sensibile values
                    # for the 'sensfile' column
                    all_sensfiles = []
                    for spec1d in data['filename']:
                        if spec1d in spec1d_to_sensfile_map:
                            sensfile = spec1d_to_sensfile_map[spec1d]
                        else:
                            sensfile = Path('SENSFUNC FILE')

                        all_sensfiles.append(sensfile.name)
                    data['sensfile'] = all_sensfiles
                    data['setup_id'] = all_setup_ids

                else:
                    # Just use one default sensfunc file
                    if len(sensfiles) == 1:
                        default_sensfile = sensfiles[0].name # Use the first sensfunc and only sensfunc
                    else:
                        default_sensfile = 'SENSFUNC FILE' # Use a dummy sensfunc filename

                    data['sensfile'] = [default_sensfile] + ([''] * (len(all_obj)-1))


                
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


