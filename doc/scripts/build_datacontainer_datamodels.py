"""
Dynamically build the rst documentation for the Calibration Images
"""

from importlib import resources
import io

from IPython import embed

import astropy.table

from pypeit.utils import to_string
from pypeit import datamodel

def link_string(p):
    return f'`{type(p).__name__} Keywords`_'

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


def build_datamodel_tbl(obj:datamodel.DataContainer) -> str:

    data_model = obj.datamodel
    keys = list(data_model.keys())
    keys.sort()


    data_table = []
    for k in keys:
        data_table.append(
            {
                'Attribute': to_string(k, use_repr=False, verbatim=True),
                'Type': type_names(data_model[k]['otype']),
                'Array Type': type_names(data_model[k]['atype']) if 'atype' in data_model[k] else ' ',
                'Description': to_string(data_model[k]['descr'])
            }
        )
    tbl  = astropy.table.Table(data_table)
    
    tbl_stream = io.StringIO()
    tbl.write(tbl_stream, format="ascii.rst")
    return tbl_stream.getvalue()

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

    datacontainers = [Alignments, EdgeTraceSet, FlatImages, ManualExtractionObj, OneSpec, OrderStack, ScatteredLight,
                      SensFunc, SlitTraceSet, Spec2DObj, SpecObj, TracePCA, WaveCalib, WaveTilts,
                      bspline, DataCube, PypeItFit, MultiSlitFlexure, Telluric, DetectorContainer,
                      Mosaic, PypeItImage, WaveFit]

    for obj in datacontainers:

        ofile = output_root / f'class_datamodel_{obj.__name__.lower()}.rst'

        # Build the Table
        lines = ['', f'**Version**: {obj.version}', '', build_datamodel_tbl(obj)]

        # Finish
        with open(ofile, 'w', encoding='utf-8') as f_obj:
            f_obj.write('\n'.join(lines))

        print(f'Wrote: {ofile}')



