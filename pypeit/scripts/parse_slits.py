"""
This script parses info from one or more SlitTraceSet objects

.. include:: ../include/links.rst
"""

from IPython import embed

from pypeit.scripts import scriptbase

from pypeit import slittrace


def print_slits(slits):
    # bitmask
    bitmask = slittrace.SlitTraceBitMask()
    print(f'{"SpatID":<8} {"MaskID":<8} {"MaskOFF (pix)":<14} {"Flags":<20}')
    for slit_idx, slit_spat in enumerate(slits.spat_id):
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

        try:
            slits = slittrace.SlitTraceSet.from_file(pargs.input_file)#, chk_version=False)
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

