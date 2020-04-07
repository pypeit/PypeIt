#!/usr/bin/env python3

"""
Dynamically build the rst documentation for the specobj and spec2dobj objects
"""

import os
import time
import numpy

from pkg_resources import resource_filename
from pypeit.par.parset import ParSet
from pypeit import specobj
from pypeit import spec2dobj
from pypeit import coadd1d

from IPython import embed

def link_string(p):
    return '`{0} Keywords`_'.format(type(p).__name__)

#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = time.perf_counter()

    # Read the baseline file that is not changed and must be edited by
    # the person building the documentation as necessary.
    pypeit_root = os.path.dirname(resource_filename('pypeit', ''))

    for obj, basefile in zip([specobj.SpecObj, spec2dobj.Spec2DObj, coadd1d.OneSpec],
                             ['datamodel_specobj.rst', 'datamodel_spec2dobj.rst',
                              'datamodel_onespec.rst']):
        lines = []
        lines += ['']

        # Start to append the automatically generated documentation
        lines += ['']
        lines += ['Version: {:s}'.format(obj.version)]
        lines += ['']

        data_model = obj.datamodel

        keys = list(data_model.keys())
        keys.sort()

        data_table = numpy.empty((len(data_model)+1, 4), dtype=object)
        data_table[0,:] = ['Obj', 'Obj Type', 'Array Type', 'Description']
        for i,k in enumerate(keys):
            # Key
            data_table[i+1,0] = ParSet._data_string(k, use_repr=False, verbatum=True)
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
            if data_table[i+1,2][-1] == '_':
                data_table[i+1,2] = data_table[i+1,2][:-1]

            # Description
            data_table[i+1,3] = ParSet._data_string(data_model[k]['desc'])

        lines += [ParSet._data_table_string(data_table, delimeter='rst')]

        # Finish
        output_rst = os.path.join(pypeit_root, 'doc', 'include', basefile)
        with open(output_rst, 'w') as f:
            f.write('\n'.join(lines))

        print('Wrote: {}'.format(output_rst))
    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



