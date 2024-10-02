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

        parser.add_argument('input_file', type=str, nargs='+',
                            help='One or more PypeIt spec2d or spec1d file')
        inp = parser.add_mutually_exclusive_group(required=True)
        inp.add_argument('--spec', default=False, action='store_true',
                         help='Check the spectral flexure')
        inp.add_argument('--spat', default=False, action='store_true',
                         help='Check the spatial flexure')
        parser.add_argument('--try_old', default=False, action='store_true',
                            help='Attempt to load old datamodel versions.  A crash may ensue..')
        return parser

    @staticmethod
    def main(args):

        from IPython import embed
        from astropy.io import fits
        from pypeit import msgs
        from pypeit import specobjs
        from pypeit import spec2dobj

        chk_version = not args.try_old
        flexure_type = 'spat' if args.spat else 'spec'

        # Loop over the input files
        for in_file in args.input_file:

            msgs.info(f'Checking fluxure for file: {in_file}')

            # What kind of file are we??
            hdul = fits.open(in_file)
            head0 = hdul[0].header

            if 'PYP_CLS' in head0.keys() and head0['PYP_CLS'].strip() == 'AllSpec2DObj':
                # load the spec2d file
                allspec2D = spec2dobj.AllSpec2DObj.from_fits(in_file, chk_version=chk_version)
                allspec2D.flexure_diagnostics(flexure_type=flexure_type)
            elif 'DMODCLS' in head0.keys() and head0['DMODCLS'].strip() == 'SpecObjs':
                if flexure_type == 'spat':
                    msgs.error("Spat flexure not available in the spec1d file, try with a "
                               "spec2d file")
                # load the spec1d file
                sobjs = specobjs.SpecObjs.from_fitsfile(in_file, chk_version=chk_version)
                sobjs.flexure_diagnostics()
            else:
                msgs.error("Bad file type input!")

            #  space between files for clarity
            print('')

