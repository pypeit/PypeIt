"""
Construct an rst table with the detector properties
"""

from importlib import resources

import astropy.table

from pypeit.spectrographs import spectrograph_classes

from IPython import embed


def write_detector_table(ofile):

    allspec = spectrograph_classes()

    data_table = []
    for key in allspec:
        # Get the spectrograph
        spec = allspec[key]()
        # Loop through detectors
        for i in range(spec.ndet):
            det = spec.get_detector_par(i+1)
            if det is None: # For HIRES
                continue
            data_table.append(
                {
                    'Instrument': f'``{key}``' if i == 0 else '...',
                    'Det': str(i+1),
                    'specaxis': str(det.specaxis),
                    'specflip': str(det.specflip),
                    'spatflip': str(det.spatflip),
                    'namp': '``None``' if det.numamplifiers is None else str(det.numamplifiers),
                    'gain': '``None``' if det.gain is None else ', '.join([str(g) for g in det.gain]),
                    'RN': '``None``' if det.ronoise is None else ', '.join([str(r) for r in det.ronoise]),
                    'darkcurr': '``None``' if det.darkcurr is None else str(det.darkcurr),
                    'min': f'{det.mincounts:.1e}',
                    'sat': str(det.saturation),
                    'nonlinear': f'{det.nonlinear:.4f}',
                    'platescale': '``None``' if det.platescale is None else f'{det.platescale:.4f}',
                }
            )
        if key == 'vlt_fors2':
            # Get the second detector for VLT-FORS2
            det = spec.get_detector_par(2)
            data_table.append(
                {
                    'Instrument': '...',
                    'Det': '2',
                    'specaxis': str(det.specaxis),
                    'specflip': str(det.specflip),
                    'spatflip': str(det.spatflip),
                    'namp': '``None``' if det.numamplifiers is None else str(det.numamplifiers),
                    'gain': '``None``' if det.gain is None else ', '.join([str(g) for g in det.gain]),
                    'RN': '``None``' if det.ronoise is None else ', '.join([str(r) for r in det.ronoise]),
                    'darkcurr': '``None``' if det.darkcurr is None else str(det.darkcurr),
                    'min': f'{det.mincounts:.1e}',
                    'sat': str(det.saturation),
                    'nonlinear': f'{det.nonlinear:.4f}',
                    'platescale': '``None``' if det.platescale is None else f'{det.platescale:.4f}',
                }
            )

    tbl = astropy.table.Table(data_table)

    tbl.write(ofile, format="ascii.rst", overwrite=True)
    print(f'Wrote: {ofile}')


def main():
    output_root = resources.files('pypeit').parent / 'doc' / 'include'
    if not output_root.is_dir():
        raise NotADirectoryError(f'{output_root} does not exist!')

    ofile = output_root / 'inst_detector_table.rst'
    write_detector_table(ofile)


if __name__ == '__main__':
    main()

