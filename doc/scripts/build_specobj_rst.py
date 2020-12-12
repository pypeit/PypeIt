#!/usr/bin/env python3

"""
Dynamically build the rst documentation for the specobj and spec2dobj objects
"""

import os
import time
import numpy

from pkg_resources import resource_filename
from pypeit.utils import to_string, string_table
from pypeit import specobj
from pypeit import spec2dobj
from pypeit import coadd1d

from IPython import embed

def link_string(p):
    return '`{0} Keywords`_'.format(type(p).__name__)

#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = time.perf_counter()

    # Set the output directory
    output_root = os.path.join(os.path.split(os.path.abspath(resource_filename('pypeit', '')))[0],
                               'doc', 'include')

    # Iterate through all the specobj classes
    for obj in [specobj.SpecObj, spec2dobj.Spec2DObj, coadd1d.OneSpec]:

        ofile = os.path.join(output_root, 'datamodel_{0}.rst'.format(obj.__name__.lower()))

        lines = []
        lines += ['']

        # Start to append the automatically generated documentation
        lines += ['']
        lines += ['Version: {:s}'.format(obj.version)]
        lines += ['']

        keys = list(obj.datamodel.keys())
        keys.sort()

        data_table = numpy.empty((len(obj.datamodel)+1, 4), dtype=object)
        data_table[0,:] = ['Obj Key', 'Obj Type', 'Array Type', 'Description']
        for i,k in enumerate(keys):
            # Key
            data_table[i+1,0] = to_string(k, use_repr=False, verbatim=True)
            # Object Type
            if isinstance(obj.datamodel[k]['otype'], (list,tuple)):
                data_table[i+1,1] = ', '.join([t.__name__ for t in obj.datamodel[k]['otype']])
            else:
                data_table[i+1,1] = obj.datamodel[k]['otype'].__name__
            # Array type
            if 'atype' in obj.datamodel[k].keys():
                data_table[i+1,2] = obj.datamodel[k]['atype'].__name__
            else:
                data_table[i+1,2] = ' '
            if data_table[i+1,2][-1] == '_':
                data_table[i+1,2] = data_table[i+1,2][:-1]

            # Description
            data_table[i+1,3] = to_string(obj.datamodel[k]['descr'])

        lines += [string_table(data_table, delimeter='rst')]

        # Finish
        with open(ofile, 'w') as f:
            f.write('\n'.join(lines))

        print('Wrote: {}'.format(ofile))
    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



