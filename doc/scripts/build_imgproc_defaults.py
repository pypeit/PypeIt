"""
Construct an rst table with the dependencies
"""

from IPython import embed
import os
from pkg_resources import resource_filename

import numpy

from pypeit.utils import string_table
from pypeit.par.pypeitpar import ProcessImagesPar, PypeItPar
from pypeit.spectrographs import spectrograph_classes


def write_imgproc_def_table(ofile, spec=None):

    # List of parameters to include in the table
    # NOTE: These are ordered according to their use in
    # pypeit.images.rawimage.RawImage.process.
    par_list = ['apply_gain',
                'use_pattern',
                'empirical_rn',
                'use_overscan',
                'trim',
                'orient',
                'use_biasimage',
                'use_darkimage',
                'spat_flexure_correct',
                'use_pixelflat',
                'use_illumflat',
                'use_specillum',
                'shot_noise',
                'noise_floor',
                'mask_cr']

    # NOTE: These are ordered according to their use in
    # pypeit.calibrations.Calibrations (although there are slight
    # differences between MultiSlit/Echelle and IFU).  And I've skipped
    # the pinholeframe.
    frame_list = ['biasframe',
                  'darkframe',
                  'traceframe',
                  'arcframe',
                  'tiltframe',
                  'alignframe',
                  'pixelflatframe',
                  'illumflatframe',
                  'skyframe',
                  'standardframe',
                  'scienceframe']

    procpar = ProcessImagesPar()
    par = PypeItPar() if spec is None else spec.default_pypeit_par()

    data_table = numpy.empty((len(par_list)+1, len(frame_list)+2), dtype=object)
    data_table[0,:] = ['Parameter', 'Default'] \
                        + [f'``{t}``'.replace('frame','') for t in frame_list]
    # Parameter names and defaults
    for i,p in enumerate(par_list):
        data_table[i+1,0] = f'``{p}``'
        data_table[i+1,1] = f'``{procpar[p]}``'
    # Frame-dependent parameter defaults
    for j,t in enumerate(frame_list):
        _par = par[t]['process'] if t == 'scienceframe' else par['calibrations'][t]['process']
        for i,p in enumerate(par_list):
            data_table[i+1,j+2] = '' if _par[p] == procpar[p] else f'``{_par[p]}``'

    lines = string_table(data_table, delimeter='rst')
    with open(ofile, 'w') as f:
        f.write(lines)
    print('Wrote: {}'.format(ofile))


def main():
    output_root = os.path.join(os.path.split(os.path.abspath(resource_filename('pypeit', '')))[0],
                               'doc', 'include')
    if not os.path.isdir(output_root):
        raise NotADirectoryError(f'{output_root} does not exist!')

    ofile = os.path.join(output_root, 'imgproc_defaults_table.rst')
    write_imgproc_def_table(ofile)

#    allspec = spectrograph_classes()
#    for key, spec_c in allspec.items():
#        ofile = os.path.join(output_root, f'imgproc_{key}_table.rst')
#        write_imgproc_def_table(ofile, spec=spec_c())

if __name__ == '__main__':
    main()

