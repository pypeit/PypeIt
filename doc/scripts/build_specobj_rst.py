"""
Dynamically build the rst documentation for the specobj and spec2dobj objects
"""

from importlib import resources

import numpy

from pypeit.utils import to_string
from rst_table import format_rst_table
from pypeit import specobj
from pypeit import spec2dobj
from pypeit import coadd1d

from IPython import embed

def write_specobj_table(obj, output_root):
    """Write the datamodel table for one extracted-spectrum container."""
    ofile = output_root / f'datamodel_{obj.__name__.lower()}.rst'

    lines = ['', '', f'Version: {obj.version}', '']

    keys = sorted(obj.datamodel)
    data_table = numpy.empty((len(keys)+1, 4), dtype=object)
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
            if isinstance(obj.datamodel[k]['atype'], tuple):
                data_table[i+1,2] = ','.join(['.'.join([a.__module__, a.__name__])
                                                if a.__module__ == 'numpy' else a.__name__
                                                for a in obj.datamodel[k]['atype']])
            else:
                data_table[i+1,2] = obj.datamodel[k]['atype'].__name__
        else:
            data_table[i+1,2] = ' '
        if data_table[i+1,2][-1] == '_':
            data_table[i+1,2] = data_table[i+1,2][:-1]

        # Description
        data_table[i+1,3] = to_string(obj.datamodel[k]['descr'])

    lines += [format_rst_table(data_table)]

    with open(ofile, 'w') as f:
        f.write('\n'.join(lines))
    print(f'Wrote: {ofile}')


if __name__ == '__main__':
    output_root = resources.files('pypeit').parent / 'doc' / 'include'
    for obj in [specobj.SpecObj, spec2dobj.Spec2DObj, coadd1d.OneSpec]:
        write_specobj_table(obj, output_root)



