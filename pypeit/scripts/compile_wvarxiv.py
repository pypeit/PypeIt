"""
This script enables the user to convert a MasterWaveCalib wavelength solution fits file
into a PypeIt arxiv solution that can be used with the full_template method.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import time
from pypeit import msgs
from pypeit import par
from pypeit import inputfiles
from pypeit import utils
from pypeit.scripts import scriptbase


class CompileWVarxiv(scriptbase.ScriptBase):
    """
    A class for compiling a set of wxarxiv solutions from Identify into a single fits file.

    Args:
        wvarxiv_folder (str): Location of the WVarxiv files.
        instrument (str): Name of the instrument (e.g., keck_lris_blue, keck_deimos, gemini_gmos_s_ham).
        grating (str): Instrument grating name (e.g., B600, R400, 600_10000).
        append (bool, optional): Append to an existing file for this instrument. Defaults to False.

    Example:
        parser = WvarxivCompile.get_parser()
        args = parser.parse_args()
        WvarxivCompile.main(args)
    """    

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Read in a set of wxarxiv solutions from Identify and compile them into a '
                                                'single fits file to be used with the reidentify method.', width=width)
        parser.add_argument('wvarxiv_folder', type = str, default=None, help='Location of the WVarxiv files')
        parser.add_argument('instrument', type = str, default=None, help='Name of instrument. e.g. keck_lris_blue, keck_deimos, gemini_gmos_south_ham')
        parser.add_argument('grating', type=str, help="Instrument grating name. E.g. b600, r400, 600_10000.")
        parser.add_argument('--append', default=False, action='store_true', help='Append to an existing file for this instrument.')

        return parser

    @staticmethod
    def main(args):
        from astropy.table import Table, join
        from importlib_resources import files as imres_files
        import glob, os

        # Read in the wvarxiv files
        assert os.path.isdir(args.wvarxiv_folder), 'The wvarxiv_folder does not exist'
        wavarxiv_files = glob.glob(args.wvarxiv_folder + '/*.fits')

        assert len(wavarxiv_files) > 0, 'No wvarxiv fits files found in the folder.'

        # Generate the fits file
        array_len = len(Table.read(wavarxiv_files[0])['wave'].data)

        reid_table = Table(names=('wave', 'flux', 'order'), dtype=(f'({array_len},)>f8', f'({array_len},)>f8', '>i8'))

        # Loop over all wvarxiv files and merge them into a single table
        for fitsfile in wavarxiv_files:
            tab = Table.read(fitsfile)
            # Make array length errors clear to the user.
            assert len(tab['wave'].data) == array_len, 'The wvarxiv arrays are not the same length. Check the fits files and ensure these correspond to the same grating and instrument.'

            # Add as a row to the table
            reid_table.add_row([tab['wave'].data, tab['flux'].data, 0])
        
        # Write to file
        out_path = imres_files('pypeit').joinpath('data', 'arc_lines', 'reid_arxiv').joinpath(f'{args.instrument}_{args.grating}_compiled.fits')

        # Does a file already exist?
        if out_path.exists() and not args.append:
            msgs.error(f'File {out_path} already exists. Use --append to overwrite the file and add your new solutions to the existing ones.')
        # What if user asks to append solutions?
        elif out_path.exists() and args.append:
            old_table = Table.read(out_path)
            old_array_len = len(old_table['wave'][0].data)
            if old_array_len != array_len:
                msgs.error(f'The old file has an array length of {old_array_len} while the new files have an array length of {array_len}. Cannot merge these files.')
            else:
                reid_table = join(old_table, reid_table)
                reid_table.write(out_path, format='fits', overwrite=args.append)
                msgs.info(f'Wrote the compiled wvarxiv file to {out_path}.')
        
        # If the file does not exist, just write it out
        else:
            reid_table.write(out_path, format='fits')
            msgs.info(f'Wrote the compiled wvarxiv file to {out_path}.')
