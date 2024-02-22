"""
This script parses info from one or more SlitTraceSet objects

.. include:: ../include/links.rst
"""

from IPython import embed

from pypeit.scripts import scriptbase

from pypeit import slittrace
from pypeit import spec2dobj
from pypeit import msgs

from astropy.table import Table
from astropy.io import fits


def print_slits(slits):
    """
    Print slits info
    Args:
        slits (:class:`~pypeit.slittrace.SlitTraceSet`):
            Slits object

    """
    diag_table = Table()
    diag_table['SpatID'] = slits.spat_id
    if slits.pypeline == 'Echelle':
        diag_table['Order'] = slits.ech_order
    if slits.maskdef_id is not None:
        diag_table['MaskID'] = slits.maskdef_id
        diag_table['MaskOFF (pix)'] = [0 if slits.maskdef_offset is None else slits.maskdef_offset]
        diag_table['MaskOFF (pix)'].format = '{:.2f}'

    # get flags
    # bitmask
    bitmask = slittrace.SlitTraceBitMask()
    allflags = ['None']*slits.nslits
    for slit_idx in range(slits.nslits):
        this_flags = []
        if slits.mask[slit_idx] != 0:
            for key in bitmask.keys():
                if bitmask.flagged(slits.mask[slit_idx], flag=key):
                    this_flags += [key]
            allflags[slit_idx] = ', '.join(this_flags)

    diag_table['Flags'] = allflags

    # print table
    diag_table.pprint_all(align='>')


class ParseSlits(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Print info on slits from a input file',
                                    width=width)
        parser.add_argument('input_file', type=str, help='Either a spec2D or Slits filename')
        parser.add_argument('--try_old', default=False, action='store_true',
                            help='Attempt to load old datamodel versions.  A crash may ensue..')
        return parser


    @staticmethod
    def main(args):

        chk_version = not args.try_old

        # What kind of file are we??
        hdul = fits.open(args.input_file)
        head0 = hdul[0].header
        head1 = hdul[1].header
        file_type = None
        if 'MSTRTYP' in head0.keys() and head0['MSTRTYP'].strip() == 'Slits':
            file_type = 'Slits'
        elif 'CALIBTYP' in head1.keys() and head1['CALIBTYP'].strip() == 'Slits':
            file_type = 'Slits'
        elif 'PYP_CLS' in head0.keys() and head0['PYP_CLS'].strip() == 'AllSpec2DObj':
            file_type = 'AllSpec2D'
        else:
            msgs.error("Bad file type input!")

        if file_type == 'Slits':
            slits = slittrace.SlitTraceSet.from_file(args.input_file, chk_version=chk_version)
            print('')
            print_slits(slits)

        elif file_type == 'AllSpec2D':
            allspec2D = spec2dobj.AllSpec2DObj.from_fits(args.input_file, chk_version=chk_version)
            # Loop on Detectors
            for det in allspec2D.detectors:
                print('')
                print('='*30 + f'{det:^7}' + '='*30)
                spec2Dobj = allspec2D[det]
                print_slits(spec2Dobj.slits)
        else:
            msgs.error("Bad file type input!  Must be a Slits calibration frame or a spec2d file.")

