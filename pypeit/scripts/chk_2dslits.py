"""
This script displays the flat images in an RC Ginga window.
"""

def parse_args(options=None, return_parser=False):
    import argparse
    parser = argparse.ArgumentParser(description='Print info on slits from a spec2D file',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('spec2d_file', type=str, help='spec2D filename')
    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


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
            maskdefID = 0 if spec2Dobj.slits.maskdef_id is None else spec2Dobj.slits.maskdef_id[slit_idx]
            maskoff = 0 if spec2Dobj.slits.mask_median_off is None else spec2Dobj.slits.mask_median_off
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


def entry_point():
    main(parse_args())


if __name__ == '__main__':
    entry_point()
