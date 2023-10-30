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


def basic_pypeitimage_datamodel(obj):
    """
    Build the datamodel tables for the base calibrations images.

    This is not general code for all DataContainer-based outputs.
    """
    data_model = obj.datamodel
    keys = list(data_model.keys())
    keys.sort()

    data_table = numpy.empty((len(keys)+2, 4), dtype=object)
    data_table[0,:] = ['HDU Name', 'HDU Type', 'Data Type', 'Description']
    data_table[1,:] = ['``PRIMARY``', '`astropy.io.fits.PrimaryHDU`_', '...',
                       'Empty data HDU.  Contains basic header information.']
    alternate_keys = []
    for i,k in enumerate(keys):
        j = i + 2
        # Key
        # Rename?
        _k = k.upper()
        if obj.hdu_prefix is not None:
            _k = obj.hdu_prefix+_k
        alternate_keys.append(_k)
        data_table[j,0] = to_string(_k, use_repr=False, verbatim=True)
        if isinstance(data_model[k]['otype'], (list, tuple)) \
                and any([issubclass(otype, DataContainer) for otype in data_model[k]['otype']]):
            data_table[j,1] = '`astropy.io.fits.BinTableHDU`_'
            data_table[j,2] = ''
        elif not isinstance(data_model[k]['otype'], (list, tuple)) \
                and issubclass(data_model[k]['otype'], DataContainer):
            data_table[j,1] = '`astropy.io.fits.BinTableHDU`_'
            data_table[j,2] = ''
        elif data_model[k]['otype'] is numpy.ndarray:
            data_table[j,1] = '`astropy.io.fits.ImageHDU`_'
            data_table[j,2] = type_names(data_model[k]['atype'])
        data_table[j,3] = to_string(data_model[k]['descr'])

    # Restrict by output_to_disk?
    if obj.output_to_disk is not None:
        keep_rows = [0, 1]
        for _k in obj.output_to_disk:
            keep_rows.append(alternate_keys.index(_k)+2)
    else:
        keep_rows = numpy.arange(len(data_table)).astype(int)

    # Parse
    data_table = data_table[numpy.asarray(keep_rows)]
    return [string_table(data_table, delimeter='rst')]

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

def edges_datamodel(output_root):

    from pypeit.edgetrace import EdgeTraceSet

    data_table = numpy.empty((21, 4), dtype=object)
    data_table[0,:] = ['HDU Name', 'HDU Type', 'Data Type', 'Description']
    data_table[1,:] = ['``PRIMARY``', '`astropy.io.fits.PrimaryHDU`_', '...',
                       'Empty data HDU.  Contains basic header information.']
    data_table[2,:] = ['``TRACE_*``', '...', '...', EdgeTraceSet.datamodel['traceimg']['descr']]
    data_table[3,:] = [ '...', '...', '...', '...']
    data_table[4,:] = ['``TRACEBPM``', '`astropy.io.fits.ImageHDU`_', 'uint8',
                       EdgeTraceSet.datamodel['tracebpm']['descr']]
    data_table[5,:] = ['``SOBELSIG``', '`astropy.io.fits.ImageHDU`_', 'float32',
                       EdgeTraceSet.datamodel['sobelsig']['descr']]
    data_table[6,:] = ['``TRACEID``', '`astropy.io.fits.ImageHDU`_', 'int64',
                       EdgeTraceSet.datamodel['traceid']['descr']]
    data_table[7,:] = ['``MASKDEF_ID``', '`astropy.io.fits.ImageHDU`_', 'int64',
                       EdgeTraceSet.datamodel['maskdef_id']['descr']]
    data_table[8,:] = ['``EDGE_CEN``', '`astropy.io.fits.ImageHDU`_', 'float32',
                       EdgeTraceSet.datamodel['edge_cen']['descr']]
    data_table[9,:] = ['``EDGE_ERR``', '`astropy.io.fits.ImageHDU`_', 'float32',
                       EdgeTraceSet.datamodel['edge_err']['descr']]
    data_table[10,:] = ['``EDGE_MSK``', '`astropy.io.fits.ImageHDU`_', 'int32',
                        EdgeTraceSet.datamodel['edge_msk']['descr']]
    data_table[11,:] = ['``EDGE_FIT``', '`astropy.io.fits.ImageHDU`_', 'float32',
                        EdgeTraceSet.datamodel['edge_fit']['descr']]
    data_table[12,:] = ['``PCA``', '`astropy.io.fits.BinTableHDU`_', '...',
                        EdgeTraceSet.datamodel['pca']['descr']]
    data_table[13,:] = ['``PCA_MODEL_?``', '`astropy.io.fits.BinTableHDU`_', '...',
                        'Model parameters for the ?th component of the PCA; see '
                        '``pca_coeffs_model`` in :class:`~pypeit.tracepca.TracePCA` and '
                        ':class:`~pypeit.core.fitting.PypeItFit`.']
    data_table[14,:] = [ '...', '...', '...', '...']
    data_table[15,:] = ['``LEFT_PCA``', '`astropy.io.fits.BinTableHDU`_', '...',
                        EdgeTraceSet.datamodel['left_pca']['descr']]
    data_table[16,:] = ['``LEFT_PCA_MODEL_?``', '`astropy.io.fits.BinTableHDU`_', '...',
                        'Model parameters for the ?th component of the left-edge PCA; see '
                        '``pca_coeffs_model`` in :class:`~pypeit.tracepca.TracePCA` and '
                        ':class:`~pypeit.core.fitting.PypeItFit`.']
    data_table[17,:] = [ '...', '...', '...', '...']
    data_table[18,:] = ['``RIGHT_PCA``', '`astropy.io.fits.BinTableHDU`_', '...',
                        EdgeTraceSet.datamodel['right_pca']['descr']]
    data_table[19,:] = ['``RIGHT_PCA_MODEL_?``', '`astropy.io.fits.BinTableHDU`_', '...',
                        'Model parameters for the ?th component of the right-edge PCA; see '
                        '``RIGHT_pca_coeffs_model`` in :class:`~pypeit.tracepca.TracePCA` and '
                        ':class:`~pypeit.core.fitting.PypeItFit`.']
    data_table[20,:] = [ '...', '...', '...', '...']

    lines = [''] + ['Version {:s}'.format(EdgeTraceSet.version)] + [''] \
                + [string_table(data_table, delimeter='rst')]

    ofile = output_root / f'datamodel_{EdgeTraceSet.__name__.lower()}.rst'
    with open(ofile, 'w') as f:
        f.write('\n'.join(lines))
    print('Wrote: {}'.format(ofile))


def slits_datamodel(output_root):

    from pypeit.edgetrace import EdgeTraceSet
    from pypeit.slittrace import SlitTraceSet
    from astropy import table

    hdu_table = numpy.empty((4, 4), dtype=object)
    hdu_table[0,:] = ['HDU Name', 'HDU Type', 'Data Type', 'Description']
    hdu_table[1,:] = ['``PRIMARY``', '`astropy.io.fits.PrimaryHDU`_', '...',
                      'Empty data HDU.  Contains basic header information.']
    hdu_table[2,:] = [f'``SLITS``', '`astropy.io.fits.BinTableHDU`_', '...',
                      'All data from the :class:`~pypeit.slittrace.SlitTraceSet` datamodel, '
                      'except ``maskdef_designtab`` if present.']
    hdu_table[3,:] = [f'``MASKDEF_DESIGNTAB``', '`astropy.io.fits.BinTableHDU`_', '...',
                      SlitTraceSet.datamodel['maskdef_designtab']['descr']]


    design = EdgeTraceSet.empty_design_table(rows=1)
    objects = EdgeTraceSet.empty_objects_table(rows=1)

    combined = table.join(design, objects)
    combined['TRACEID'].description = design['TRACEID'].description
    ncol = len(combined.keys())

    data_table = numpy.empty((ncol+1, 3), dtype=object)
    data_table[0,:] = ['Column', 'Data Type', 'Description']
    for i,key in enumerate(combined.keys()):
        t = 'str' if combined[key].dtype.type.__name__ == 'str_' \
                else combined[key].dtype.type.__name__
        data_table[i+1,:] = [f'``{key}``', t, combined[key].description]

    lines = [''] + ['Version {:s}'.format(SlitTraceSet.version)] + [''] \
                + [string_table(hdu_table, delimeter='rst')] + [''] \
                + ['MASKDEF_DESIGNTAB content'] + [''] \
                + [string_table(data_table, delimeter='rst')]

    ofile = output_root / f'datamodel_{SlitTraceSet.__name__.lower()}.rst'
    with open(ofile, 'w') as f:
        f.write('\n'.join(lines))
    print('Wrote: {}'.format(ofile))


def wavecalib_datamodel(output_root):

    from pypeit.wavecalib import WaveCalib
    from pypeit.core.wavecal.wv_fitting import WaveFit

    data_table = numpy.empty((7, 4), dtype=object)
    data_table[0,:] = ['HDU Name', 'HDU Type', 'Data Type', 'Description']
    data_table[1,:] = ['``PRIMARY``', '`astropy.io.fits.PrimaryHDU`_', '...',
                       'Empty data HDU.  Contains basic header information.']
    data_table[2,:] = ['``SPAT_IDS``', '`astropy.io.fits.ImageHDU`_',
                       type_name(WaveCalib.datamodel['spat_ids']['atype']),
                       WaveCalib.datamodel['spat_ids']['descr']]
    prefix = WaveFit.hduext_prefix_from_spatid('?')
    data_table[3,:] = [f'``{prefix}_WAVEFIT``', '`astropy.io.fits.BinTableHDU`_', '...',
                       ':class:`~pypeit.core.wavecal.wv_fitting.WaveFit` result for ``slit_id=?``']
    data_table[4,:] = [f'``{prefix}_PYPEITFIT``', '`astropy.io.fits.BinTableHDU`_', '...',
                       ':class:`~pypeit.core.fitting.PypeItFit` element of the '
                       ':class:`~pypeit.core.wavecal.wv_fitting.WaveFit` datamodel for '
                       '``slit_id=?``']
    data_table[5,:] = [ '...', '...', '...', '...']
    data_table[6,:] = ['``ARC_SPECTRA``', '`astropy.io.fits.ImageHDU`_',
                       type_name(WaveCalib.datamodel['arc_spectra']['atype']),
                       WaveCalib.datamodel['arc_spectra']['descr']]


    lines = [''] + ['Version {:s}'.format(WaveCalib.version)] + [''] \
                + [string_table(data_table, delimeter='rst')]

    ofile = output_root / f'datamodel_{WaveCalib.__name__.lower()}.rst'
    with open(ofile, 'w') as f:
        f.write('\n'.join(lines))
    print('Wrote: {}'.format(ofile))


def flatfield_datamodel(output_root):

    from pypeit.flatfield import FlatImages

    data_table = numpy.empty((19, 4), dtype=object)
    data_table[0,:] = ['HDU Name', 'HDU Type', 'Data Type', 'Description']
    data_table[1,:] = ['``PRIMARY``', '`astropy.io.fits.PrimaryHDU`_', '...',
                       'Empty data HDU.  Contains basic header information.']
    data_table[2,:] = ['``PIXELFLAT_RAW``', '`astropy.io.fits.ImageHDU`_', 'float64',
                       FlatImages.datamodel['pixelflat_raw']['descr']]
    data_table[3,:] = ['``PIXELFLAT_NORM``', '`astropy.io.fits.ImageHDU`_', 'float64',
                       FlatImages.datamodel['pixelflat_norm']['descr']]
    data_table[4,:] = ['``PIXELFLAT_MODEL``', '`astropy.io.fits.ImageHDU`_', 'float64',
                       FlatImages.datamodel['pixelflat_model']['descr']]
    data_table[5,:] = ['``PIXELFLAT_SPAT_ID-?-BSPLINE``', '`astropy.io.fits.BinTableHDU`_', '...',
                       'bspline model of the pixelflat for ``spat_id=?``; see '
                       ':class:`~pypeit.bspline.bspline.bspline`']
    data_table[6,:] = ['...', '...', '...', '...']
    data_table[7,:] = ['``PIXELFLAT_SPAT_ID-?-FINECORR``', '`astropy.io.fits.BinTableHDU`_', '...',
                       '2D polynomial fits to the fine correction of the spatial illumination '
                       'profile of the pixelflat for ``slit_id=?``; see '
                       ':class:`~pypeit.core.fitting.PypeItFit`']
    data_table[8,:] = ['...', '...', '...', '...']
    data_table[9,:] = ['``PIXELFLAT_BPM``', '`astropy.io.fits.ImageHDU`_', 'int16',
                       FlatImages.datamodel['pixelflat_bpm']['descr']]
    data_table[10,:] = ['``PIXELFLAT_SPEC_ILLUM``', '`astropy.io.fits.ImageHDU`_', 'float64',
                        FlatImages.datamodel['pixelflat_spec_illum']['descr']]
    data_table[11,:] = ['``PIXELFLAT_WAVEIMG``', '`astropy.io.fits.ImageHDU`_', 'float64',
                        FlatImages.datamodel['pixelflat_waveimg']['descr']]
    data_table[12,:] = ['``ILLUMFLAT_RAW``', '`astropy.io.fits.ImageHDU`_', 'float64',
                        FlatImages.datamodel['illumflat_raw']['descr']]
    data_table[13,:] = ['``ILLUMFLAT_SPAT_ID-?-BSPLINE``', '`astropy.io.fits.BinTableHDU`_', '...',
                        'bspline model of the illumflat for ``spat_id=?``; see '
                        ':class:`~pypeit.bspline.bspline.bspline`']
    data_table[14,:] = ['...', '...', '...', '...']
    data_table[15,:] = ['``ILLUMFLAT_SPAT_ID-?-FINECORR``', '`astropy.io.fits.BinTableHDU`_', '...',
                        '2D polynomial fits to the fine correction of the spatial illumination '
                        'profile of the illumflat for ``slit_id=?``; see '
                        ':class:`~pypeit.core.fitting.PypeItFit`']
    data_table[16,:] = ['...', '...', '...', '...']
    data_table[17,:] = ['``ILLUMFLAT_BPM``', '`astropy.io.fits.ImageHDU`_', 'int16',
                        FlatImages.datamodel['illumflat_bpm']['descr']]
    data_table[18,:] = ['``SPAT_ID``', '`astropy.io.fits.ImageHDU`_', 'int64',
                        FlatImages.datamodel['spat_id']['descr']]

    lines = [''] + ['Version {:s}'.format(FlatImages.version)] + [''] \
                + [string_table(data_table, delimeter='rst')]

    ofile = output_root / f'datamodel_{FlatImages.__name__.lower()}.rst'
    with open(ofile, 'w') as f:
        f.write('\n'.join(lines))
    print('Wrote: {}'.format(ofile))


if __name__ == '__main__':
    # Set the output directory
    output_root = resources.files('pypeit').parent / 'doc' / 'include'

    # Simple datamodels for Arc, Bias, Dark, Tiltimg
    for obj in [buildimage.ArcImage, buildimage.BiasImage, buildimage.DarkImage,
                buildimage.TiltImage]:

        # Build the Table
        lines = [''] + ['Version {:s}'.format(obj.version)] + [''] \
                    + basic_pypeitimage_datamodel(obj)

        ofile = output_root / f'datamodel_{obj.__name__.lower()}.rst'
        with open(ofile, 'w') as f:
            f.write('\n'.join(lines))

        print('Wrote: {}'.format(ofile))


    # All data written to a single extension for Alignments and Tilts
    from pypeit.alignframe import Alignments
    single_table_datamodel(Alignments, output_root, 'ALIGN',
                           'Spatial alignment data; see :class:`~pypeit.alignframe.Alignments`.')
    from pypeit.wavetilts import WaveTilts
    single_table_datamodel(WaveTilts, output_root, 'TILTS',
                           'Tilts data; see :class:`~pypeit.wavetilts.WaveTilts`.')

    # Edges
    edges_datamodel(output_root)

    # Slits
    slits_datamodel(output_root)

    # WaveCalib
    wavecalib_datamodel(output_root)

    # Flat
    flatfield_datamodel(output_root)




