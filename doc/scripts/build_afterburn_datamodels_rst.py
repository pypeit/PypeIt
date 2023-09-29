"""
Dynamically build the rst documentation for the Calibration Images
"""

from importlib import resources

import numpy

from pypeit.utils import to_string, string_table
from pypeit.datamodel import DataContainer
from pypeit.images import buildimage
from pypeit.flatfield import FlatImages
from pypeit.edgetrace import EdgeTraceSet
from pypeit.slittrace import SlitTraceSet

from IPython import embed

def link_string(p):
    return '`{0} Keywords`_'.format(type(p).__name__)

#-----------------------------------------------------------------------------

def type_name(t):
    if t is numpy.bool_:
        return 'np.bool'
    return t.__name__


def type_names(types):
    if isinstance(types, (list,tuple)):
        return ', '.join([type_name(t) for t in types])
    return type_name(types)

def column_type(t):
    if t.dtype.type.__name__ == 'str_':
        return 'str'
    if t.dtype.type.__name__ == 'bool_':
        return 'bool'
    return t.dtype.type.__name__


def single_table_datamodel(obj, output_root, ext, descr):

    data_table = numpy.empty((3, 4), dtype=object)
    data_table[0,:] = ['HDU Name', 'HDU Type', 'Data Type', 'Description']
    data_table[1,:] = ['``PRIMARY``', '`astropy.io.fits.PrimaryHDU`_', '...',
                       'Empty data HDU.  Contains basic header information.']
    data_table[2,:] = [f'``{ext}``', '`astropy.io.fits.BinTableHDU`_', '...', descr]

    lines = [''] + ['Version {:s}'.format(obj.version)] + [''] \
                + [string_table(data_table, delimeter='rst')]

    ofile = output_root / f'datamodel_{obj.__name__.lower()}.rst'
    with open(ofile, 'w') as f:
        f.write('\n'.join(lines))
    print('Wrote: {}'.format(ofile))


def sens_datamodel(output_root):

    from pypeit.sensfunc import SensFunc
    from pypeit.core.telluric import Telluric
    from astropy import table

    hdu_table = numpy.empty((7, 4), dtype=object)
    hdu_table[0,:] = ['HDU Name', 'HDU Type', 'Data Type', 'Description']
    hdu_table[1,:] = ['``PRIMARY``', '`astropy.io.fits.PrimaryHDU`_', '...',
                      'Empty data HDU.  Contains basic header information.']
    hdu_table[2,:] = [f'``TELLURIC``', '`astropy.io.fits.BinTableHDU`_', '...',
                      SensFunc.datamodel['telluric']['descr'] + '.  This is identical to the HDU '
                      'extension produced in the file produced by :ref:`pypeit_tellfit`.  Only '
                      'provided when using the IR algorithm.']
    hdu_table[3,:] = [f'``SENS``', '`astropy.io.fits.BinTableHDU`_', '...',
                      SensFunc.datamodel['sens']['descr']]
    hdu_table[4,:] = [f'``WAVE``', '`astropy.io.fits.ImageHDU`_', '...',
                      SensFunc.datamodel['wave']['descr'] + '.  May be combined from many '
                      'detectors; see the ``multi_spec_det`` parameter in :ref:`sensfuncpar`. ']
    hdu_table[5,:] = [f'``ZEROPOINT``', '`astropy.io.fits.ImageHDU`_', '...',
                      SensFunc.datamodel['zeropoint']['descr'] + '.  May be combined from many '
                      'detectors; see the ``multi_spec_det`` parameter in :ref:`sensfuncpar`. ']
    hdu_table[6,:] = [f'``THROUGHPUT``', '`astropy.io.fits.ImageHDU`_', '...',
                      SensFunc.datamodel['throughput']['descr'] + '.  May be combined from many '
                      'detectors; see the ``multi_spec_det`` parameter in :ref:`sensfuncpar`. ']

    telluric = Telluric.empty_model_table(1, 1)
    ncol = len(telluric.keys())
    tell_table = numpy.empty((ncol+1, 3), dtype=object)
    tell_table[0,:] = ['Column', 'Data Type', 'Description']
    for i,key in enumerate(telluric.keys()):
        tell_table[i+1,:] = [f'``{key}``', column_type(telluric[key]), telluric[key].description]

    sens = SensFunc.empty_sensfunc_table(1, 1, 1)
    ncol = len(sens.keys())
    sens_table = numpy.empty((ncol+1, 3), dtype=object)
    sens_table[0,:] = ['Column', 'Data Type', 'Description']
    for i,key in enumerate(sens.keys()):
        sens_table[i+1,:] = [f'``{key}``', column_type(sens[key]), sens[key].description]

    lines = [''] + ['Version {:s}'.format(SensFunc.version)] + [''] \
                + [string_table(hdu_table, delimeter='rst')] \
                + [''] + ['TELLURIC table (if present)'] + [''] \
                + [string_table(tell_table, delimeter='rst')] \
                + [''] + ['SENS table'] + [''] \
                + [string_table(sens_table, delimeter='rst')]

    ofile = output_root / f'datamodel_{SensFunc.__name__.lower()}.rst'
    with open(ofile, 'w') as f:
        f.write('\n'.join(lines))
    print('Wrote: {}'.format(ofile))


def telluric_datamodel(output_root):

    from pypeit.core.telluric import Telluric
    from astropy import table

    hdu_table = numpy.empty((3, 4), dtype=object)
    hdu_table[0,:] = ['HDU Name', 'HDU Type', 'Data Type', 'Description']
    hdu_table[1,:] = ['``PRIMARY``', '`astropy.io.fits.PrimaryHDU`_', '...',
                      'Empty data HDU.  Contains basic header information.']
    hdu_table[2,:] = [f'``TELLURIC``', '`astropy.io.fits.BinTableHDU`_', '...',
                      'Results of the telluric modeling']

    telluric = Telluric.empty_model_table(1, 1)
    ncol = len(telluric.keys())
    tell_table = numpy.empty((ncol+1, 3), dtype=object)
    tell_table[0,:] = ['Column', 'Data Type', 'Description']
    for i,key in enumerate(telluric.keys()):
        tell_table[i+1,:] = [f'``{key}``', column_type(telluric[key]), telluric[key].description]

    lines = [''] + ['Version {:s}'.format(Telluric.version)] + [''] \
                + [string_table(hdu_table, delimeter='rst')] \
                + [''] + ['TELLURIC table'] + [''] \
                + [string_table(tell_table, delimeter='rst')]

    ofile = output_root / f'datamodel_{Telluric.__name__.lower()}.rst'
    with open(ofile, 'w') as f:
        f.write('\n'.join(lines))
    print('Wrote: {}'.format(ofile))


if __name__ == '__main__':
    # Set the output directory
    output_root = resources.files('pypeit').parent / 'doc' / 'include'

    # Sensitivity file
    sens_datamodel(output_root)

    # Telluric file
    telluric_datamodel(output_root)

    from pypeit.onespec import OneSpec
    single_table_datamodel(OneSpec, output_root, 'Spectrum',
                           'Single spectrum data; see :class:`~pypeit.onespec.OneSpec`.')



