#!/usr/bin/env python3

"""
Dynamically build the rst documentation for the Calibration Images
"""

import os
import time
import numpy

from pkg_resources import resource_filename
from pypeit.utils import to_string, string_table
from pypeit.images import buildimage
from pypeit.flatfield import FlatImages
from pypeit.wavetilts import WaveTilts
from pypeit.wavecalib import WaveCalib

from IPython import embed

def link_string(p):
    return '`{0} Keywords`_'.format(type(p).__name__)

#-----------------------------------------------------------------------------

def build_tbl(obj):

    data_model = obj.full_datamodel(include_children=False)
    keys = list(data_model.keys())
    keys.sort()

    data_table = numpy.empty((len(keys)+1, 4), dtype=object)
    data_table[0,:] = ['HDU Name', 'Obj Type', 'Array Type', 'Description']
    alternate_keys = []
    for i,k in enumerate(keys):
        # Key
        # Rename?
        _k = k.upper()
        if obj.hdu_prefix is not None:
            _k = obj.hdu_prefix+_k
        alternate_keys.append(_k)
        data_table[i+1,0] = to_string(_k, use_repr=False, verbatim=True)
        # Object Type
        if isinstance(data_model[k]['otype'], (list,tuple)):
            data_table[i+1,1] = ', '.join([t.__name__ for t in data_model[k]['otype']])
        else:
            data_table[i+1,1] = data_model[k]['otype'].__name__
        # Array type
        if 'atype' in data_model[k].keys():
            data_table[i+1,2] = data_model[k]['atype'].__name__
        else:
            data_table[i+1,2] = ' '
        # Description
        data_table[i+1,3] = to_string(data_model[k]['descr'])

    # Restrict by output_to_disk?
    if obj.output_to_disk is not None:
        keep_rows = [0]
        for _k in obj.output_to_disk:
            keep_rows.append(alternate_keys.index(_k)+1)
    else:
        keep_rows = numpy.arange(len(data_table)).astype(int)

    # Parse
    data_table = data_table[numpy.asarray(keep_rows)]

    tbl_lines = [string_table(data_table, delimeter='rst')]
    return tbl_lines


if __name__ == '__main__':
    t = time.perf_counter()

    # Set the output directory
    output_root = os.path.join(os.path.split(os.path.abspath(resource_filename('pypeit', '')))[0],
                               'doc', 'include')
    for obj in [buildimage.ArcImage, buildimage.BiasImage, buildimage.TiltImage, WaveCalib, WaveTilts,
                FlatImages]:

        ofile = os.path.join(output_root, 'datamodel_{0}.rst'.format(obj.__name__.lower()))

        # Build the Table
        lines = build_tbl(obj)
        lines = [''] + ['Version {:s}'.format(obj.version)] + [''] + lines

        # Finish
        with open(ofile, 'w') as f:
            f.write('\n'.join(lines))

        print('Wrote: {}'.format(ofile))
    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



