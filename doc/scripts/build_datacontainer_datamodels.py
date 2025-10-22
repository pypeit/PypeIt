"""
Dynamically build the rst documentation for the Calibration Images
"""

from importlib import resources

from IPython import embed

import numpy

from pypeit.utils import to_string, string_table
from pypeit import datamodel

def link_string(p):
    return '`{0} Keywords`_'.format(type(p).__name__)

#-----------------------------------------------------------------------------

def type_name(t):
    if issubclass(t, datamodel.DataContainer):
        return f':class:`~{t.__module__}.{t.__name__}`'
    if any([m in t.__module__ for m in ['numpy', 'astropy']]):
        name = 'bool' if t.__name__ == 'bool_' else t.__name__
        return f'`{t.__module__}.{name}`_'
    return t.__name__


def type_names(types):
    if isinstance(types, (list,tuple)):
        return ', '.join([type_name(t) for t in types])
    return type_name(types)


def build_datamodel_tbl(obj):

    data_model = obj.datamodel
    keys = list(data_model.keys())
    keys.sort()

    data_table = numpy.empty((len(keys)+1, 4), dtype=object)
    data_table[0,:] = ['Attribute', 'Type', 'Array Type', 'Description']
    for i,k in enumerate(keys):
        # Key
        data_table[i+1,0] = to_string(k, use_repr=False, verbatim=True)
        # Object Type
        data_table[i+1,1] = type_names(data_model[k]['otype'])
        # Array type
        if 'atype' in data_model[k].keys():
            data_table[i+1,2] = type_names(data_model[k]['atype'])
        else:
            data_table[i+1,2] = ' '
        # Description
        data_table[i+1,3] = to_string(data_model[k]['descr'])

    return [string_table(data_table, delimeter='rst')]


if __name__ == '__main__':

    # Set the output directory
    output_root = resources.files('pypeit').parent / 'doc' / 'include'

    # All DataContainer objects
    # TODO: automate this?
    from pypeit.alignframe import Alignments
    from pypeit.edgetrace import EdgeTraceSet
    from pypeit.flatfield import FlatImages
    from pypeit.manual_extract import ManualExtractionObj
    from pypeit.onespec import OneSpec
    from pypeit.orderstack import OrderStack
    from pypeit.scattlight import ScatteredLight
    from pypeit.sensfunc import SensFunc
    from pypeit.slittrace import SlitTraceSet
    from pypeit.spec2dobj import Spec2DObj
    from pypeit.specobj import SpecObj
    from pypeit.tracepca import TracePCA
    from pypeit.wavecalib import WaveCalib
    from pypeit.wavetilts import WaveTilts
    from pypeit.bspline import bspline
    from pypeit.coadd3d import DataCube
    from pypeit.core.fitting import PypeItFit
    from pypeit.core.flexure import MultiSlitFlexure
    from pypeit.core.telluric import Telluric
    from pypeit.core.wavecal.wv_fitting import WaveFit
    from pypeit.images.detector_container import DetectorContainer
    from pypeit.images.mosaic import Mosaic
    from pypeit.images.pypeitimage import PypeItImage

    from pypeit.images import buildimage

    datacontainers = [Alignments, EdgeTraceSet, FlatImages, ManualExtractionObj, OneSpec, OrderStack, ScatteredLight,
                      SensFunc, SlitTraceSet, Spec2DObj, SpecObj, TracePCA, WaveCalib, WaveTilts,
                      bspline, DataCube, PypeItFit, MultiSlitFlexure, Telluric, DetectorContainer,
                      Mosaic, PypeItImage, WaveFit]

    for obj in datacontainers:

        ofile = output_root / f'class_datamodel_{obj.__name__.lower()}.rst'

        # Build the Table
        lines = [''] + [f'**Version**: {obj.version}'] + [''] + build_datamodel_tbl(obj)

        # Finish
        with open(ofile, 'w') as f:
            f.write('\n'.join(lines))

        print('Wrote: {}'.format(ofile))



