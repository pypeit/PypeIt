"""
This script displays the flat images in an RC Ginga window.
"""
import argparse

from pypeit import slittrace
from pypeit import msgs
from IPython import embed
from pypeit import spec2dobj


def parser(options=None):
    parser = argparse.ArgumentParser(description='Print info on slits from a spec2D file',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('spec2d_file', type=str, help='spec2D filename')
    return parser.parse_args() if options is None else parser.parse_args(options)


def main(pargs):
    # bitmask
    bitmask = slittrace.SlitTraceBitMask()
    # Load
    allspec2D = spec2dobj.AllSpec2DObj.from_fits(pargs.spec2d_file)
    # Loop on Detectors
    for det in allspec2D.detectors:
        print("================ DET {:02d} ======================".format(det))
        spec2Dobj = allspec2D[det]
        print("SpatID  MaskID  Flags")
        for slit_idx, slit_spat in enumerate(spec2Dobj.slit_info[:,0]):
            line = '{:04d}    {:04d}'.format(slit_spat, spec2Dobj.slit_info[slit_idx,2])
            # Flags
            flags = []
            if spec2Dobj.slit_info[slit_idx, 1] == 0:
                flags += ['None']
            else:
                for key in bitmask.keys():
                    if bitmask.flagged(spec2Dobj.slit_info[slit_idx, 1], key):
                        flags += [key]
            # Finish
            sflag = ', '
            line += '    '+sflag.join(flags)
            print(line)
