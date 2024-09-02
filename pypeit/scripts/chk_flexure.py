"""
This script displays the flexure (spatial or spectral) applied to the science data.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


class ChkFlexure(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Print QA on flexure to the screen',
                                    width=width)

        parser.add_argument('input_file', type=str, nargs='+', help='One or more PypeIt spec2d or spec1d file')
        inp = parser.add_mutually_exclusive_group(required=True)
        inp.add_argument('--spec', default=False, action='store_true', help='Check the spectral flexure')
        inp.add_argument('--spat', default=False, action='store_true', help='Check the spatial flexure')
        parser.add_argument('--try_old', default=False, action='store_true',
                            help='Attempt to load old datamodel versions.  A crash may ensue..')
        return parser

    @staticmethod
    def main(args):

        from astropy.io import fits
        from pypeit import msgs
        from pypeit.core import flexure

        chk_version = not args.try_old

        # Loop over the input files
        for in_file in args.input_file:

            msgs.info(f'Checking flexure for file: {in_file}')

            # What kind of file are we??
            hdul = fits.open(in_file)
            head0 = hdul[0].header
            file_type = None
            if 'PYP_CLS' in head0.keys() and head0['PYP_CLS'].strip() == 'AllSpec2DObj':
                file_type = 'spec2d'
            elif 'DMODCLS' in head0.keys() and head0['DMODCLS'].strip() == 'SpecObjs':
                file_type = 'spec1d'
            else:
                msgs.error("Bad file type input!")

            # Check the flexure
            flexure.flexure_diagnostic(in_file, file_type=file_type, flexure_type='spat' if args.spat else 'spec',
                                       chk_version=chk_version)

            #  space between files for clarity
            print('')






