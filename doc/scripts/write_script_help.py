"""
Dynamically build the rst documentation with the script help text.
"""

import os
import time
import importlib
from pkg_resources import resource_filename

from pypeit.scripts import script_classes


#-----------------------------------------------------------------------------

def write_help(script_cls, opath, width=80):
    exe = script_cls.name()
    ofile = os.path.join(opath, '{0}.rst'.format(exe))
    lines = ['.. code-block:: console', '']
    lines += ['    $ {0} -h'.format(exe)]
    parser = script_cls.get_parser(width=80)
    parser.prog = exe
    lines += ['    ' + l for l in parser.format_help().split('\n')]
    print('Writing: {0}'.format(ofile))
    with open(ofile, 'w') as f:
        f.write('\n'.join(lines))

if __name__ == '__main__':
    t = time.perf_counter()

    pypeit_root = os.path.dirname(resource_filename('pypeit', ''))
    path = os.path.join(pypeit_root, 'doc', 'help')
    if not os.path.isdir(path):
        os.makedirs(path)

    # Get the list of script names and script classes
    scr_clss = script_classes()

    for name, script_cls in scr_clss.items():
        write_help(script_cls, path)

    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))


