
from importlib import resources
import os
import subprocess

import astropy.table

from IPython import embed


def tellgrid_table(root, ofile):

    os.environ['TZ'] = 'UTC'
    try:
        result = subprocess.run(
            [
                'aws', '--endpoint', 'https://s3-west.nrp-nautilus.io', 's3', 'ls',
                f's3://pypeit/telluric/atm_grids/{root}', '--human-readable'
            ], capture_output=True, text=True)
    except Exception:
        print(
            f'Exception raised by subprocess call to `aws`, used to list the available {root} '
            'telluric files.  The associated telluric file table will not be updated!'
        )
        raise

    # Check for success
    if result.returncode != 0:
        raise ValueError(
            f'Return code from `aws` system call was {result.returncode}, not 0.  The stderr '
            f'output was {"empty" if len(result.stderr) == 0 else result.stderr}.'
        )

    # Construct and write the table
    data = [l.split() for l in result.stdout.split('\n')[:-1]]
    if not all(len(d) == 5 for d in data):
        raise ValueError(
            'Unexpected output from `aws` call.  Splitting each entry should yield 5 components.'
        )
    
    data_table = []
    for d in data:
        data_table.append(
            {
                'File': d[4],
                'Size': ' '.join(d[2:4]),
                'Last Modified (UTC)': ' '.join(d[:2]),
            }
        )

    tbl = astropy.table.Table(data_table)

    tbl.write(ofile, format="ascii.rst", overwrite=True)
    print(f'Wrote: {ofile}')

def main():
    output_root = resources.files('pypeit').parent / 'doc' / 'include'
    root = 'TelFit'
    ofile = output_root / f'{root}_files.rst'
    tellgrid_table(root, ofile)
    root = 'TellPCA'
    ofile = output_root / f'{root}_files.rst'
    tellgrid_table(root, ofile)


if __name__ == '__main__':
    main()
