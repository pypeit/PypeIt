"""
Dynamically build table listing available spectrographs.
"""

from importlib import resources

import astropy.table

from pypeit.utils import to_string
from pypeit.spectrographs import spectrograph_classes

from IPython import embed

#-----------------------------------------------------------------------------

def write_spec_table(path):
    ofile = path / 'spectrographs_table.rst'

    data_table = []
    for cls in spectrograph_classes().values():

        data_table.append(
            {
                '``PypeIt`` Name': cls.name,
                '``PypeIt`` Class': ':class:`~' + cls.__module__ + '.' + cls.__name__ + '`',
                'Telescope': cls.telescope['name'],
                'Camera': cls.camera,
                'URL': '' if cls.url is None else f'`Link <{cls.url}>`__',
                'Pipeline': cls.pypeline,
                'Supported': to_string(cls.supported),
                'QL Tested': to_string(cls.ql_supported),
                'Comments': '' if cls.comment is None else cls.comment
            }
        )

    tbl = astropy.table.Table(data_table)

    tbl.write(ofile, format="ascii.rst", overwrite=True)
    print(f'Wrote: {ofile}')


if __name__ == '__main__':
    output_root = resources.files('pypeit').parent / 'doc' / 'include'
    if not output_root.is_dir():
        output_root.mkdir(parents=True)

    write_spec_table(output_root)
