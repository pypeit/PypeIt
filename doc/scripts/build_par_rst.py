#!/usr/bin/env python3

"""
Dynamically build the rst documentation of the pypit parameters.
"""

import os
import time
import numpy

from pkg_resources import resource_filename
from pypit.par import pypitpar
from pypit.par.parset import ParSet

#-----------------------------------------------------------------------------
#def class_name(p):
#    return '.'.join([type(p).__module__, type(p).__name__])


def link_string(p):
    return '`{0} Keywords`_'.format(type(p).__name__)


def par_hierarchy(p, indent_level=0, key=''):
    indent_step = ' '*indent_level*4
    line_head = '['*indent_level + key + ']'*indent_level
    if len(line_head) > 0:
        line_head = '``' + line_head + '``: '
    lines = [ indent_step + line_head + link_string(p) ]
    lines += [ '' ]

    for k in p.keys():
        if not isinstance(p[k], ParSet):
            continue
        lines += par_hierarchy(p[k], indent_level=indent_level+1, key=k)
    
    return lines

#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = time.clock()

    # Read the baseline file that is not changed and must be edited by
    # the person building the documentation as necessary.
    pypit_root = os.path.dirname(resource_filename('pypit', ''))
    input_base = os.path.join(pypit_root, 'doc', 'scripts', 'base_par_rst.txt')
    with open(input_base, 'r') as f:
        lines = [ l.replace('\n','') for l in f.readlines() ]
    lines += ['']

    # Start to append the automatically generated documentation
    lines += ['Current PypitPar Parameter Hierarchy']
    lines += ['++++++++++++++++++++++++++++++++++++']
    lines += ['']

    p = pypitpar.PypitPar(skysubtract=pypitpar.SkySubtractionPar(),
                          flexure=pypitpar.FlexurePar(),
                          fluxcalib=pypitpar.FluxCalibrationPar())

    lines += par_hierarchy(p)
    lines += ['']
    lines += ['----']
    lines += ['']

    lines += p.to_rst_table()
    lines += ['']

    output_rst = os.path.join(pypit_root, 'doc', 'pypit_par.rst')
    with open(output_rst, 'w') as f:
        f.write('\n'.join(lines))
    
    print('Elapsed time: {0} seconds'.format(time.clock() - t))



