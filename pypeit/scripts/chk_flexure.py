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
        # parser.add_argument('type', type=str,
        #                     help='Type of flexure to check. Options are: spat, spec')
        parser.add_argument('--try_old', default=False, action='store_true',
                            help='Attempt to load old datamodel versions.  A crash may ensue..')
        return parser

    @staticmethod
    def main(args):

        from IPython import embed
        import numpy as np
        from astropy.io import fits
        from astropy.table import Table
        from pypeit import spec2dobj, specobjs, msgs

        chk_version = not args.try_old

        # tables to return
        return_tables = []

        # Loop over the input files
        for in_file in args.input_file:

            print(f'\nCheck fluxure for file: {in_file}')

            # What kind of file are we??
            hdul = fits.open(in_file)
            head0 = hdul[0].header
            file_type = None
            if 'PYP_CLS' in head0.keys() and head0['PYP_CLS'].strip() == 'AllSpec2DObj':
                file_type = 'AllSpec2D'
            elif 'DMODCLS' in head0.keys() and head0['DMODCLS'].strip() == 'SpecObjs':
                file_type = 'SpecObjs'
            else:
                msgs.error("Bad file type input!")

            if file_type == 'AllSpec2D':
                # load the spec2d file
                allspec2D = spec2dobj.AllSpec2DObj.from_fits(in_file, chk_version=chk_version)
                # Loop on Detectors
                for det in allspec2D.detectors:
                    print('')
                    print('='*50 + f'{det:^7}' + '='*51)
                    # get and print the spectral flexure
                    if args.spec:
                        spec_flex = allspec2D[det].sci_spec_flexure
                        spec_flex.rename_column('sci_spec_flexure', 'global_spec_shift')
                        if np.all(spec_flex['global_spec_shift'] != None):
                            spec_flex['global_spec_shift'].format = '0.3f'
                        # print the table
                        spec_flex.pprint_all()
                        # return the table
                        return_tables.append(spec_flex)
                    # get and print the spatial flexure
                    if args.spat:
                        spat_flex = allspec2D[det].sci_spat_flexure
                        # print the value
                        print(f'Spat shift: {spat_flex}')
                        # return the value
                        return_tables.append(spat_flex)
            elif file_type == 'SpecObjs':
                # no spat flexure in spec1d file
                if args.spat:
                    msgs.error("Spat flexure not available in the spec1d file, try with a spec2d file")
                # load the spec1d file
                sobjs = specobjs.SpecObjs.from_fitsfile(in_file, chk_version=chk_version)
                spec_flex = Table()
                spec_flex['NAME'] = sobjs.NAME
                spec_flex['global_spec_shift'] = sobjs.FLEX_SHIFT_GLOBAL
                if np.all(spec_flex['global_spec_shift'] != None):
                    spec_flex['global_spec_shift'].format = '0.3f'
                spec_flex['local_spec_shift'] = sobjs.FLEX_SHIFT_LOCAL
                if np.all(spec_flex['local_spec_shift'] != None):
                    spec_flex['local_spec_shift'].format = '0.3f'
                spec_flex['total_spec_shift'] = sobjs.FLEX_SHIFT_TOTAL
                if np.all(spec_flex['total_spec_shift'] != None):
                    spec_flex['total_spec_shift'].format = '0.3f'
                # print the table
                spec_flex.pprint_all()
                # return the table
                return_tables.append(spec_flex)

        return return_tables






