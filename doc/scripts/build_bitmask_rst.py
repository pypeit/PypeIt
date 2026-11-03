"""
Dynamically build tables for the bitmasks.
"""

from importlib import resources

import astropy.table

from pypeit.utils import to_string

from IPython import embed

#-----------------------------------------------------------------------------

def write_bitmask_table(obj, path):
    bm = obj()
    ofile = path / f'{obj.__name__.lower()}_table.rst'

    data_table = []
    for k in bm.bits.keys():
        data_table.append(
            {
                'Bit Name': k,
                'Bit Number': to_string(bm.bits[k]),
                'Decimal Value': to_string(int(2**bm.bits[k])),
                'Description': to_string(bm.descr[bm.bits[k]]),
            }
        )

    tbl = astropy.table.Table(data_table)

    tbl.write(ofile, format="ascii.rst", overwrite=True)
    print(f'Wrote: {ofile}')


if __name__ == '__main__':
    path = resources.files('pypeit').parent / 'doc' / 'include'
    if not path.is_dir():
        path.mkdir(parents=True)

    from pypeit.images.imagebitmask import ImageBitMask
    write_bitmask_table(ImageBitMask, path)
