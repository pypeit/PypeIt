"""
Convert the contributors.csv list into an rst table to be included in the team
doc.  This script *can* be run every time the docs are updated.
"""
import ast
from importlib import resources

from astropy.table import Table
from IPython import embed

def main():
    doc_root = resources.files('pypeit').parent / 'doc'
    contributors = doc_root / 'scripts' / 'contributors.csv'
    out_file = doc_root / 'include' / 'contributors.rst'

    tbl = Table.read(
        contributors, format='ascii.no_header', names=['GitHub handle', 'Name', 'valid'],
        delimiter=',',
    )
    tbl['valid'] = [ast.literal_eval(v) for v in tbl['valid']]
    tbl = tbl[tbl['valid']]['Name', 'GitHub handle']
    tbl['GitHub handle'] = [f'@{login}' for login in tbl['GitHub handle']]

    tbl.write(out_file, format='ascii.rst', overwrite=True)


if __name__ == '__main__':
    main()
