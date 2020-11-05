"""
Dynamically build the rst documentation with the script help text.
"""

import os
import time
import importlib
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
    print('Writing: {0}'.format(ofile))
    with open(ofile, 'w') as f:
        f.write('\n'.join(lines))

if __name__ == '__main__':
    t = time.perf_counter()

    pypeit_root = os.path.dirname(resource_filename('pypeit', ''))
    path = os.path.join(pypeit_root, 'doc', 'help')
    if not os.path.isdir(path):
        os.makedirs(path)

    # Make a dictionary with all of the script modules and whether or
    # not `pypeit` should be prepended to the name of the executable

    # TODO: There might be a smarter way of getting all the attributes
    # of a module...
    scr_mod = {s:True for s in ['arcid_plot', 'chk_2dslits', 'chk_alignments', 'chk_edges',
                                'chk_flats', 'chk_for_calibs', 'chk_tilts', 'coadd_1dspec',
                                'coadd_2dspec', 'coadd_datacube', 'compare_sky', 'find_objects',
                                'flux_calib', 'flux_setup', 'identify', 'lowrdx_pixflat', 
                                'lowrdx_skyspec', 'qa_html', 'ql_keck_mosfire', 'ql_keck_nires',
                                'ql_mos', 'sensfunc', 'setup', 'show_1dspec', 'show_2dspec',
                                'show_arxiv', 'show_wvcalib', 'skysub_regions', 'tellfit',
                                'trace_edges', 'view_fits', 'run_pypeit']}
    scr_mod['run_pypeit'] = False

    for mod,prepend in scr_mod.items():
        write_help(importlib.import_module('pypeit.scripts.{0}'.format(mod)), path,
                   prepend_pypeit=prepend)

    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))


