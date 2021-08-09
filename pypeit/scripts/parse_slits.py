"""
This script parses info from one or more SlitTraceSet objects
"""

from pypeit import slittrace


def parse_args(options=None, return_parser=False):
    import argparse
    parser = argparse.ArgumentParser(description='Print info on slits from a spec2D file',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input_file', type=str, help='Either a spec2D or MasterSlits filename')
    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def print_slits(slits):
    # bitmask
    bitmask = slittrace.SlitTraceBitMask()
    print("SpatID  MaskID  MaskOFF(pix)  Flags")
    for slit_idx, slit_spat in enumerate(slits.spat_id):
        maskdefID = 0 if slits.maskdef_id is None else slits.maskdef_id[slit_idx]
        maskoff = 0 if slits.mask_median_off is None else slits.mask_median_off
        line = '{:04d}    {:04d}    {:.2f}   '.format(slit_spat, maskdefID, maskoff)
        # Flags
        flags = []
        if slits.mask[slit_idx] == 0:
            flags += ['None']
        else:
            for key in bitmask.keys():
                if bitmask.flagged(slits.mask[slit_idx], key):
                    flags += [key]
        # Finish
        sflag = ', '
        line += '    ' + sflag.join(flags)
        print(line)


def main(pargs):

    from astropy.io import fits
    from IPython import embed

    from pypeit import spec2dobj


    # What kind of file are we??
    hdul = fits.open(pargs.input_file)
    head0 = hdul[0].header
    if 'MSTRTYP' in head0.keys() and head0['MSTRTYP'].strip() == 'Slits':
        file_type = 'MasterSlits'
    # TODO -- Replace the next elif with checking PYP_CLS in head0 someday
    elif 'DMODCLS' in hdul[1].header and hdul[1].header['DMODCLS'] == 'Spec2DObj':
        file_type = 'AllSpec2D'
    else:
        raise IOError("Bad file type input!")

    if file_type == 'MasterSlits':
        slits = slittrace.SlitTraceSet.from_file(pargs.input_file)#, chk_version=False)
        print_slits(slits)
    elif file_type == 'AllSpec2D':
        # Load
        allspec2D = spec2dobj.AllSpec2DObj.from_fits(pargs.input_file, chk_version=False)
        # Loop on Detectors
        for det in allspec2D.detectors:
            print("================ DET {:02d} ======================".format(det))
            spec2Dobj = allspec2D[det]
            print_slits(spec2Dobj.slits)

def entry_point():
    main(parse_args())


if __name__ == '__main__':
    entry_point()