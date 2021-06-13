"""
This script displays the flat images in an RC Ginga window.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


class Chk2DSlits(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Print info on slits from a spec2D file',
                                    width=width)
        parser.add_argument('spec2d_file', type=str, help='spec2D filename')
        return parser

    @staticmethod
    def main(args):

        from IPython import embed

        from pypeit import slittrace
        from pypeit import spec2dobj

        # bitmask
        bitmask = slittrace.SlitTraceBitMask()
        # Load
        allspec2D = spec2dobj.AllSpec2DObj.from_fits(args.spec2d_file, chk_version=False)
        # Loop on Detectors
        for det in allspec2D.detectors:
            print("================ DET {:02d} ======================".format(det))
            spec2Dobj = allspec2D[det]
            print("SpatID  MaskID  MaskOFF(pix)  Flags")
            for slit_idx, slit_spat in enumerate(spec2Dobj.slits.spat_id):
                maskdefID = 0 if spec2Dobj.slits.maskdef_id is None \
                                else spec2Dobj.slits.maskdef_id[slit_idx]
                maskoff = 0 if spec2Dobj.slits.mask_median_off is None \
                                else spec2Dobj.slits.mask_median_off
                line = '{:04d}    {:04d}    {:.2f}   '.format(slit_spat, maskdefID, maskoff)
                # Flags
                flags = []
                if spec2Dobj.slits.mask[slit_idx] == 0:
                    flags += ['None']
                else:
                    for key in bitmask.keys():
                        if bitmask.flagged(spec2Dobj.slits.mask[slit_idx], key):
                            flags += [key]
                # Finish
                sflag = ', '
                line += '    '+sflag.join(flags)
                print(line)


