"""
Dynamically build tables listing available wavelength calibration line list and
reid_arxiv spectra.
"""

from importlib import resources
import time
import datetime

# TODO: datetime.UTC is not defined in python 3.10.  Remove this when we decide
# to no longer support it.
try:
    __UTC__ = datetime.UTC
except AttributeError as e:
    __UTC__ = datetime.timezone.utc

import numpy
import astropy.table

from pypeit.utils import string_table
from pypeit import dataPaths

from IPython import embed


def write_linelist_table(output_root):

    files = sorted(dataPaths.linelist.glob('*_lines.dat'))
    nfiles = len(files)

    lamp = numpy.empty(nfiles, dtype=object)
    wv_rng = numpy.empty(nfiles, dtype=object)
    mod = numpy.empty(nfiles, dtype=object)

    for i,f in enumerate(files):
        lamp[i] = f.name.replace('_lines.dat', '')

        tbl = astropy.table.Table.read(f, format='ascii.fixed_width', comment='#')
        wv_rng[i] = f'{numpy.floor(numpy.amin(tbl['wave'])):.0f} - ' \
                    f'{numpy.ceil(numpy.amax(tbl['wave'])):.0f}'
        mod[i] = datetime.datetime.fromtimestamp(f.stat().st_mtime, tz=__UTC__).strftime('%Y-%m-%d')

    nrows = int(numpy.ceil(nfiles / 2))
    data_table = numpy.full((nrows+1, 7), '', dtype=object)
    data_table[0,:] = ['Lamp', 'Range (Å)', 'Last Mod', ' ',
                       'Lamp', 'Range (Å)', 'Last Mod']

    data_table[1:,0] = lamp[:nrows]
    data_table[1:,1] = wv_rng[:nrows]
    data_table[1:,2] = mod[:nrows]

    data_table[1:nfiles-nrows+1,4] = lamp[nrows:]
    data_table[1:nfiles-nrows+1,5] = wv_rng[nrows:]
    data_table[1:nfiles-nrows+1,6] = mod[nrows:]

    ofile = output_root / 'linelist_table.rst'
    lines = string_table(data_table, delimeter='rst')
    with open(ofile, 'w') as f:
        f.write(lines)
    print(f'Wrote: {ofile}')


def write_reid_arxiv_table(output_root):

    # TODO: Pull wavelength range (and resolution?) from files

    reid_dir = dataPaths.reid_arxiv

    # Find all the files
    files = sorted(list(reid_dir.glob('*.fits')) + list(reid_dir.glob('*.json')))

    # Read the summary file
    tbl = astropy.table.Table.read(reid_dir / 'summary.txt', format='ascii', delimiter='|') 

    # Check that all the files are in the table and vice versa
    names = numpy.array([f.name for f in files])
    indx = numpy.isin(names, tbl['File'])
    if not numpy.all(indx):
        raise ValueError('Files found in reid_arxiv directory that are *not* in the summary.txt '
                         f'file.  Include summary entries for the following (or delete them): '
                         f'{names[numpy.logical_not(indx)]}')
    indx = numpy.isin(tbl['File'], names)
    if not numpy.all(indx):
        raise ValueError('Files found in the summary.txt file that are *not* in the reid_arxiv '
                         'directory.  Remove summary entries for the following: '
                         f'{tbl["File"][numpy.logical_not(indx)]}')
    
    nrow = len(tbl)
    ncol = len(tbl.keys())
    data_table = numpy.full((nrow+1, ncol), '', dtype=object)
    data_table[0,:] = tbl.keys()
    for i,key in enumerate(tbl.keys()):
        try:
            data_table[1:,i] = tbl[key].data.astype(str).filled(' ')
        except:
            data_table[1:,i] = tbl[key].data.astype(str)

    ofile = output_root / 'reid_arxiv_table.rst'
    lines = string_table(data_table, delimeter='rst')
    with open(ofile, 'w') as f:
        f.write(lines)
    print(f'Wrote: {ofile}')


def main():
    t = time.perf_counter()

    output_root = resources.files('pypeit').parent / 'doc' / 'include'
    if not output_root.is_dir():
        output_root.mkdir(parents=True)

    write_linelist_table(output_root)
    write_reid_arxiv_table(output_root)

    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))


if __name__ == '__main__':
    main()


