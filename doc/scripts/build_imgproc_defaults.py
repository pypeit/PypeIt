"""
Construct an rst table with the dependencies
"""

from importlib import resources

import astropy.table

from pypeit.par.pypeitpar import ProcessImagesPar, PypeItPar

from IPython import embed


def write_imgproc_def_table(ofile, spec=None):

    # List of parameters to include in the table
    # NOTE: These are ordered according to their use in
    # pypeit.images.rawimage.RawImage.process.
    par_list = ['apply_gain',
                'use_pattern',
                'empirical_rn',
                'use_overscan',
                'trim',
                'orient',
                'use_biasimage',
                'use_darkimage',
                'spat_flexure_correct',
                'use_pixelflat',
                'use_illumflat',
                'use_specillum',
                'shot_noise',
                'noise_floor',
                'mask_cr']

    # NOTE: These are ordered according to their use in
    # pypeit.calibrations.Calibrations (although there are slight
    # differences between MultiSlit/Echelle and IFU).  And I've skipped
    # the pinholeframe.
    frame_list = ['biasframe',
                  'darkframe',
                  'traceframe',
                  'arcframe',
                  'tiltframe',
                  'alignframe',
                  'pixelflatframe',
                  'illumflatframe',
                  'skyframe',
                  'standardframe',
                  'scienceframe']

    procpar = ProcessImagesPar()
    par = PypeItPar() if spec is None else spec.default_pypeit_par()

    data_table = []
    for p in par_list:
        pdict = {'Parameter':f'``{p}``','Default':f'``{procpar[p]}``'}
        for t in frame_list:
            _par = par[t]['process'] if t == 'scienceframe' else par['calibrations'][t]['process']
            pdict.update(
                {f'``{t}``'.replace('frame',''):'' if _par[p] == procpar[p] else f'``{_par[p]}``'
                }
            )
        data_table.append(pdict)

    tbl = astropy.table.Table(data_table)

    tbl.write(ofile, format="ascii.rst", overwrite=True)
    print(f'Wrote: {ofile}')


def main():
    output_root = resources.files('pypeit').parent / 'doc' / 'include'
    if not output_root.is_dir():
        raise NotADirectoryError(f'{output_root} does not exist!')

    ofile = output_root / 'imgproc_defaults_table.rst'
    write_imgproc_def_table(ofile)

if __name__ == '__main__':
    main()

