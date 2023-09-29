"""
Construct an rst table with the detector properties
"""

from importlib import resources

import numpy

from pypeit.utils import string_table
from pypeit.spectrographs import spectrograph_classes

from IPython import embed

#def write_detector_datamodel(ofile):
#    from pypeit.images.detector_container import DetectorContainer
#    det_dm = DetectorContainer.datamodel
#    keys = numpy.sort(list(det_dm.keys()))
#    nkeys = len(keys)
#    data_table = numpy.empty((nkeys+1, 3), dtype=object)
#    data_table[0,:] = ['Key', 'Type', 'Description']
#    for i,key in enumerate(keys):
#        data_table[i+1,0] = f'``{key}``'
#        if isinstance(det_dm[key]['otype'], tuple):
#            data_table[i+1,1] = ','.join([t.__name__ for t in det_dm[key]['otype']])
#        else:
#            data_table[i+1,1] = det_dm[key]['otype'].__name__
#        data_table[i+1,2] = det_dm[key]['descr']
#    lines = string_table(data_table, delimeter='rst')
#    with open(ofile, 'w') as f:
#        f.write(lines)
#    print('Wrote: {}'.format(ofile))


def write_detector_table(ofile):

    allspec = spectrograph_classes()
    data_table = [['Instrument', 'Det', 'specaxis', 'specflip', 'spatflip', 'namp', 'gain',
                   'RN', 'darkcurr', 'min', 'sat', 'nonlinear', 'platescale']]
    for key in allspec.keys():
        spec = allspec[key]()
        for i in range(spec.ndet):
            det = spec.get_detector_par(i+1)
            if det is None: # For HIRES
                continue
            dt_row = [f'``{key}``'] if i == 0 else ['...']
            dt_row += [str(i+1), str(det.specaxis), str(det.specflip), str(det.spatflip),
                       '``None``' if det.numamplifiers is None else str(det.numamplifiers),
                       '``None``' if det.gain is None else ', '.join([str(g) for g in det.gain]),
                       '``None``' if det.ronoise is None 
                            else ', '.join([str(r) for r in det.ronoise]),
                       '``None``' if det.darkcurr is None else str(det.darkcurr),
                       f'{det.mincounts:.1e}', str(det.saturation), f'{det.nonlinear:.4f}',
                       f'{det.platescale:.4f}']
            data_table += [dt_row]
        if key == 'vlt_fors2':
            # Get the second detector for VLT-FORS2
            det = spec.get_detector_par(2)
            dt_row = ['...', '2', str(det.specaxis), str(det.specflip), str(det.spatflip),
                       '``None``' if det.numamplifiers is None else str(det.numamplifiers),
                       '``None``' if det.gain is None else ', '.join([str(g) for g in det.gain]),
                       '``None``' if det.ronoise is None 
                            else ', '.join([str(r) for r in det.ronoise]),
                       '``None``' if det.darkcurr is None else str(det.darkcurr),
                       f'{det.mincounts:.1e}', str(det.saturation), f'{det.nonlinear:.4f}',
                       f'{det.platescale:.4f}']
            data_table += [dt_row]

    lines = string_table(numpy.atleast_1d(data_table), delimeter='rst')
    with open(ofile, 'w') as f:
        f.write(lines)
    print('Wrote: {}'.format(ofile))


def main():
    output_root = resources.files('pypeit').parent / 'doc' / 'include'
    if not output_root.is_dir():
        raise NotADirectoryError(f'{output_root} does not exist!')

#    ofile = output_root / 'datamodel_detector.rst'
#    write_detector_datamodel(ofile)

    ofile = output_root / 'inst_detector_table.rst'
    write_detector_table(ofile)


if __name__ == '__main__':
    main()

