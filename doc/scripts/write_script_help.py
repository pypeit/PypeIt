"""
Dynamically build the rst documentation with the script help text.
"""

import os
import time
from pkg_resources import resource_filename


#-----------------------------------------------------------------------------

def write_help(script_mod, opath, prepend_pypeit=False):
    exe = script_mod.__name__.split('.')[-1]
    if prepend_pypeit:
        exe = 'pypeit_'+exe
    ofile = os.path.join(opath, '{0}.rst'.format(exe))
    lines = ['.. code-block:: console', '']
    lines += ['    $ {0} -h'.format(exe)]
    parser = script_mod.parse_args(return_parser=True)
    parser.prog = exe
    lines += ['    ' + l for l in parser.format_help().split('\n')]
    with open(ofile, 'w') as f:
        f.write('\n'.join(lines))

if __name__ == '__main__':
    t = time.perf_counter()

    pypeit_root = os.path.dirname(resource_filename('pypeit', ''))
    path = os.path.join(pypeit_root, 'doc', 'help')
    if not os.path.isdir(path):
        os.makedirs(path)

    from pypeit.scripts import show_2dspec
    write_help(show_2dspec, path, prepend_pypeit=True)

    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))


