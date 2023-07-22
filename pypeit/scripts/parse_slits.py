"""
This script parses info from one or more SlitTraceSet objects

.. include:: ../include/links.rst
"""

from IPython import embed

from pypeit.scripts import scriptbase

from pypeit import slittrace
from pypeit import msgs

from IPython import embed


def print_slits(slits):
    # bitmask
    bitmask = slittrace.SlitTraceBitMask()
    if slits.pypeline  in ['MultiSlit', 'IFU']:
        slitord_id = slits.spat_id
        slit_label = 'SpatID'
    elif slits.pypeline == 'Echelle':
        slitord_id = slits.slitord_id
        slit_label = 'Order'
    else:
        msgs.error('Not ready for this pypeline: {0}'.format(slits.pypeline))
    print(f'{slit_label:<8} {"MaskID":<8} {"MaskOFF (pix)":<14} {"Flags":<20}')
    # TODO JFH No need to print out the MaskID and MaskOFF for echelle
    for slit_idx, slit_spat in enumerate(slitord_id):
        maskdefID = 0 if slits.maskdef_id is None else slits.maskdef_id[slit_idx]
        maskoff = 0 if slits.maskdef_offset is None else slits.maskdef_offset
        # Flags
        flags = []
        if slits.mask[slit_idx] == 0:
            flags += ['None']
        else:
            for key in bitmask.keys():
                if bitmask.flagged(slits.mask[slit_idx], key):
                    flags += [key]
        print('{0:<8} {1:<8} {2:<14} {3:<20}'.format(f'{slit_spat:04d}', f'{maskdefID:04d}',
                                                     f'{maskoff:.2f}', ', '.join(flags)))


class ParseSlits(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Print info on slits from a input file',
                                    width=width)
        parser.add_argument('input_file', type=str, help='Either a spec2D or Slits filename')
        return parser


    @staticmethod
    def main(pargs):

        from astropy.io import fits
        from IPython import embed

        from pypeit import spec2dobj


        # What kind of file are we??
        hdul = fits.open(pargs.input_file)
        head0 = hdul[0].header
        head1 = hdul[1].header
        if 'MSTRTYP' in head0.keys() and head0['MSTRTYP'].strip() == 'Slits':
            file_type = 'Slits'
        elif 'CALIBTYP' in head1.keys() and head1['CALIBTYP'].strip() == 'Slits':
            file_type = 'Slits'
        elif 'PYP_CLS' in head0.keys() and head0['PYP_CLS'].strip() == 'AllSpec2DObj':
            file_type = 'AllSpec2D'
        else:
            raise IOError("Bad file type input!")

        if file_type == 'Slits':
            try:
                slits = slittrace.SlitTraceSet.from_file(pargs.input_file, chk_version=False)
            # TODO: Should this specify the type of exception to pass?
            except:
                pass
            else:
                print_slits(slits)
                return

        try:
            allspec2D = spec2dobj.AllSpec2DObj.from_fits(pargs.input_file, chk_version=False)
            # TODO: Should this specify the type of exception to pass?
        except:
            pass
        else:
            # Loop on Detectors
            for det in allspec2D.detectors:
                print('='*30 + f'{det:^7}' + '='*30)
                spec2Dobj = allspec2D[det]
                print_slits(spec2Dobj.slits)
            return
        
        raise IOError("Bad file type input!  Must be a Slits calibration frame or a spec2d file.")

