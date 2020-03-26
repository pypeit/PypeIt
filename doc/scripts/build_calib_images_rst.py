#!/usr/bin/env python3

"""
Dynamically build the rst documentation for the Calibration Images
"""

import os
import time
import numpy

from pkg_resources import resource_filename
from pypeit.par.parset import ParSet
from pypeit.images import buildimage
from pypeit.flatfield import FlatImages
from pypeit.wavetilts import WaveTilts

from IPython import embed

def link_string(p):
    return '`{0} Keywords`_'.format(type(p).__name__)

#-----------------------------------------------------------------------------

def build_tbl(imgtyp):

    data_model = imgtyp.full_datamodel()
    keys = list(data_model.keys())
    keys.sort()

    data_table = numpy.empty((len(keys)+1, 4), dtype=object)
    data_table[0,:] = ['HDU Name', 'Obj Type', 'Array Type', 'Description']
    alternate_keys = []
    for i,k in enumerate(keys):
        # Key
        # Rename?
        _k = k.upper()
        if imgtyp.hdu_prefix is not None:
            _k = imgtyp.hdu_prefix+_k
        alternate_keys.append(_k)
        data_table[i+1,0] = ParSet._data_string(_k, use_repr=False, verbatum=True)
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
        data_table[i+1,3] = ParSet._data_string(data_model[k]['desc'])

    # Restrict by output_to_disk?
    if imgtyp.output_to_disk is not None:
        keep_rows = [0]
        for _k in imgtyp.output_to_disk:
            keep_rows.append(alternate_keys.index(_k)+1)
    else:
        keep_rows = numpy.arange(len(data_table)).astype(int)

    # Parse
    data_table = data_table[numpy.asarray(keep_rows)]

    tbl_lines = [ParSet._data_table_string(data_table, delimeter='rst')]
    return tbl_lines

if __name__ == '__main__':
    t = time.perf_counter()

    # Read the baseline file that is not changed and must be edited by
    # the person building the documentation as necessary.
    pypeit_root = os.path.dirname(resource_filename('pypeit', ''))
    output_path = os.path.join(pypeit_root, 'doc', 'include')

    for imgtyp,ofile in zip([buildimage.ArcImage,
                             buildimage.BiasImage,
                             buildimage.TiltImage,
                             WaveTilts,
                             FlatImages],
                            [os.path.join(output_path, 'datamodel_arcimage.rst'),
                             os.path.join(output_path, 'datamodel_biasimage.rst'),
                             os.path.join(output_path, 'datamodel_tiltimage.rst'),
                             os.path.join(output_path, 'datamodel_wavetilts.rst'),
                             os.path.join(output_path, 'datamodel_flatimages.rst'),
                              ]):
        # Build the Table
        tbl_lines = build_tbl(imgtyp)
        tbl_lines = [''] + ['Version {:s}'.format(imgtyp.version)] + [''] + tbl_lines

        # Insert lines
        #pos = lines.index(insert_line)
        #if pos < 0:
        #    raise ValueError("Missing insert line!")

        # Finish
        #output_rst = os.path.join(pypeit_root, 'doc', 'calib_images.rst')
        with open(ofile, 'w') as f:
            f.write('\n'.join(tbl_lines))

        print('Wrote: {}'.format(ofile))
    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



