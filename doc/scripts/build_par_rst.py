"""
Dynamically build the rst documentation of the pypeit parameters.
"""

from importlib import resources
import textwrap
import warnings

import numpy as np

from pypeit import utils
from pypeit.par import pypeitpar
from pypeit.par.parset import ParSet
from pypeit.spectrographs.util import load_spectrograph, available_spectrographs

from IPython import embed


def class_index(names):
    """
    Build a ``sphinx_design`` grid of cards linking to each ``ParSet``
    class's dropdown section, keyed by its :func:`~pypeit.par.parset.ParSet.to_rst_table`
    ``:name:`` target (the lower-cased class name).

    Parameters
    ----------
    names : list
        The (unsorted, possibly duplicated) list of ``ParSet`` subclass
        names to include in the index.

    Returns
    -------
    list
        The list of rst lines with the grid of cards.
    """
    lines = ['.. grid:: 2 3 4 4']
    lines += ['    :gutter: 2']
    lines += ['']
    for name in sorted(set(names)):
        lines += [f'    .. grid-item-card:: {name}']
        lines += [f'        :link: {name.lower()}']
        lines += ['        :link-type: ref']
        lines += ['']
    return lines


def par_hierarchy(p, indent_level=0, key=''):
    indent_step = ' '*indent_level*4
    line_head = '['*indent_level + key + ']'*indent_level
    lines = [ indent_step + line_head ]

    for k in p.keys():
        if not isinstance(p[k], ParSet):
            continue
        lines += par_hierarchy(p[k], indent_level=indent_level+1, key=k)
    
    return lines

def individual_par_tables():

    pypeit_root = resources.files('pypeit').parent 
    output_root = pypeit_root / 'doc' / 'include'

    parset_subclasses = utils.all_subclasses(ParSet)
    parset_subclasses = np.asarray([
        p for p in parset_subclasses
        if not issubclass(p, pypeitpar.TelescopePar) or p is pypeitpar.TelescopePar
    ])
    srt = np.argsort([cls.__name__ for cls in parset_subclasses])
    for p in parset_subclasses[srt]:
        try:
            def_par = p()
        except:
            warnings.warn(f'Skipping table for ParSet subclass {p.__name__}')
        lines = def_par.to_rst_table(include_keyword_link=False, top_level_only=True)
        ofile = output_root / f'parset_{p.__name__}.rst'
        with open(ofile, 'w') as f:
            f.write(lines)
        print(f'Wrote: {ofile}')

#-----------------------------------------------------------------------------

if __name__ == '__main__':

    # Construct the individual tables for the docstrings of each ParSet
    # subclass.
    individual_par_tables()

    # Now create the full page for the "User-level Parameters"

    # Read the baseline file that is not changed and must be edited by
    # the person building the documentation as necessary.
    pypeit_root = resources.files('pypeit').parent 
    input_base = pypeit_root / 'doc' / 'scripts' / 'base_par.rst'
    with open(input_base, 'r') as f:
        lines = [ l.replace('\n','') for l in f.readlines() ]
    lines += ['']

    # Start to append the automatically generated documentation
    lines += ['Parameter Hierarchy and Definition Tables']
    lines += ['='*len(lines[-1])]
    lines += ['']

    p = pypeitpar.PypeItPar()

    # Build the full set of per-class tables first so that the class names
    # collected along the way (parsets_listed) can be used to build the
    # jump-to-class index below.
    parsets_listed = []
    table_lines = p.to_rst_table(dropdown=True, parsets_listed=parsets_listed)

    lines += class_index(parsets_listed + [type(p).__name__])
    lines += ['']

    lines += ['.. dropdown:: Current PypeItPar Parameter Hierarchy']
    lines += ['    :name: par-hierarchy']
    lines += ['']
    lines += ['    .. code-block:: ini']
    lines += ['']
    lines += ['        ' + l for l in par_hierarchy(p)]
    lines += ['']
    lines += ['----']
    lines += ['']

    lines += table_lines
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
        title = ' '.join([s.telescope['name'], s.camera, f'(``{s.name}``)'])
        lines += [f'.. dropdown:: {title}']
        lines += [f'    :name: instr_par-{s.name}']
        lines += ['']
        lines += ['    .. code-block:: ini']
        lines += ['']
        sl = s.default_pypeit_par().to_config(include_descr=False, exclude_defaults=True)
        lines += ['        ' + l for l in sl]
        lines += ['']
    lines += ['']

    # Automatically open a dropdown when a link (e.g., from the class index
    # above) jumps to it.  This is included here, rather than registered
    # site-wide in conf.py, so the behavior is confined to this page.
    lines += ['.. raw:: html']
    lines += ['']
    lines += ['    <script src="_static/js/dropdown_autoopen.js"></script>']
    lines += ['']

    output_rst = pypeit_root / 'doc' / 'pypeit_par.rst'
    with open(output_rst, 'w') as f:
        f.write('\n'.join(lines))
