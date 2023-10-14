"""
Dynamically build the rst documentation of the pypeit parameters.
"""

from importlib import resources
import time
import textwrap

from pypeit.par import pypeitpar
from pypeit.par.parset import ParSet
from pypeit.spectrographs.util import load_spectrograph
from pypeit.spectrographs import available_spectrographs

from IPython import embed

#-----------------------------------------------------------------------------
#def class_name(p):
#    return '.'.join([type(p).__module__, type(p).__name__])


def link_string(p):
    return f':ref:`{type(p).__name__.lower()}`'
#    return '`{0} Keywords`_'.format(type(p).__name__)


def par_hierarchy(p, indent_level=0, key=''):
    indent_step = ' '*indent_level*4
    line_head = '['*indent_level + key + ']'*indent_level
    if len(line_head) > 0:
        line_head = '``' + line_head + '``: '
    lines = [ indent_step + line_head + link_string(p) ]
#    lines += [ '' ]

    for k in p.keys():
        if not isinstance(p[k], ParSet):
            continue
        lines += par_hierarchy(p[k], indent_level=indent_level+1, key=k)
    
    return lines

#-----------------------------------------------------------------------------

if __name__ == '__main__':
    t = time.perf_counter()

    # Read the baseline file that is not changed and must be edited by
    # the person building the documentation as necessary.
    pypeit_root = resources.files('pypeit').parent 
    input_base = pypeit_root / 'doc' / 'scripts' / 'base_par.rst'
    with open(input_base, 'r') as f:
        lines = [ l.replace('\n','') for l in f.readlines() ]
    lines += ['']

    # Start to append the automatically generated documentation
    lines += ['Current PypeItPar Parameter Hierarchy']
    lines += ['=====================================']
    lines += ['']

    p = pypeitpar.PypeItPar(flexure=pypeitpar.FlexurePar(),
                            fluxcalib=pypeitpar.FluxCalibratePar())

    lines += ['| '+ l for l in par_hierarchy(p)]
    lines += ['']
    lines += ['----']
    lines += ['']

    lines += p.to_rst_table()
    lines += ['']

    lines += ['.. _instr_par:']
    lines += ['']

    lines += ['Instrument-Specific Default Configuration']
    lines += ['=========================================']
    lines += ['']

    lines += textwrap.wrap('The following provides the changes to the global default parameters '
                           'provided above for each instrument.  That is, if one were to include '
                           'these in the PypeIt file, you would be reproducing the effect of the '
                           '`default_pypeit_par` method specific to each derived '
                           ':class:`~pypeit.spectrographs.spectrograph.Spectrograph` class.', 72)
    lines += ['']

    for spec in available_spectrographs:
        s = load_spectrograph(spec)
        lines += [ f'.. _instr_par-{s.name}:']
        lines += ['']
        lines += [ ' '.join([s.telescope['name'], s.camera, '(``{0}``)'.format(s.name)]) ]
        lines += [ '-'*len(lines[-1]) ]
        lines += [ 'Alterations to the default parameters are:' ]
        lines += ['']
        lines += ['.. code-block:: ini']
        lines += ['']
        sl = s.default_pypeit_par().to_config(include_descr=False, exclude_defaults=True)
        lines += [ '  ' + l for l in sl ]
        lines += ['']
    lines += ['']

    output_rst = pypeit_root / 'doc' / 'pypeit_par.rst'
    with open(output_rst, 'w') as f:
        f.write('\n'.join(lines))
    
    print('Elapsed time: {0} seconds'.format(time.perf_counter() - t))



