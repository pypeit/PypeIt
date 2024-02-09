"""
This script runs PypeIt on a pair of NIRES images (A-B)

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


class QLKeckNIRES(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Run PypeIt on an A-B pair of NIRES files',
                                    width=width)
        parser.add_argument('full_rawpath', type=str, help='Full path to the raw files')
        parser.add_argument('fileA', type=str, help='A frame')
        parser.add_argument('fileB', type=str, help='B frame')
        parser.add_argument('-b', '--box_radius', type=float,
                            help='Set the radius for the boxcar extraction')
        return parser

    @staticmethod
    def main(args):

        import os
        import sys
        import numpy as np

        from IPython import embed

        from pypeit import pypeit
        from pypeit import pypeitsetup
        from pypeit.core import framematch
        from pypeit import msgs
        from pypeit import data

        # Setup
        data_files = [os.path.join(args.full_rawpath, args.fileA),
                      os.path.join(args.full_rawpath, args.fileB)]
        ps = pypeitsetup.PypeItSetup(data_files, path='./', spectrograph_name='keck_nires')
        ps.build_fitstbl()
        # TODO -- Get the type_bits from  'science'
        bm = framematch.FrameTypeBitMask()
        file_bits = np.zeros(2, dtype=bm.minimum_dtype())
        file_bits[0] = bm.turn_on(file_bits[0], ['arc', 'science', 'tilt'])
        file_bits[1] = bm.turn_on(file_bits[0], ['arc', 'science', 'tilt'])

        # Setup (this needs to be before set_combination_groups() since it's used there)
        ps.fitstbl['setup'] = 'A'
        ps.fitstbl.set_frame_types(file_bits)
        ps.fitstbl.set_combination_groups()
        # Extras
        # A-B
        ps.fitstbl['bkg_id'] = [2,1]

        # Calibrations
        default_master_dir = os.getenv('QL_MASTERS')
        master_dir = os.path.join(data.Paths.data, 'QL_MASTERS') \
                        if default_master_dir is None else default_master_dir
        master_dir = os.path.join(master_dir, 'NIRES_MASTERS')
        if not os.path.isdir(master_dir):
            msgs.error(f'{master_dir} does not exist!  You must install the QL_MASTERS '
                       'directory; download the data from the PypeIt dev-suite Google Drive and '
                       'either define a QL_MASTERS environmental variable or use the '
                       'pypeit_install_ql_masters script.')

        # Config the run
        cfg_lines = ['[rdx]']
        cfg_lines += ['    spectrograph = {0}'.format('keck_nires')]
        cfg_lines += ['    redux_path = {0}'.format(os.path.join(os.getcwd(),'keck_nires_A'))]
        # Calibrations
        cfg_lines += ['[baseprocess]']
        cfg_lines += ['    use_biasimage = False']
        cfg_lines += ['    use_overscan = False']
        cfg_lines += ['    use_pixelflat = False']
        cfg_lines += ['[calibrations]']
        cfg_lines += ['    master_dir = {0}'.format(master_dir)]
        cfg_lines += ['    raise_chk_error = False']
        cfg_lines += ['[scienceframe]']
        cfg_lines += ['    [[process]]']
        cfg_lines += ['        mask_cr = False']
        cfg_lines += ['[reduce]']
        cfg_lines += ['    [[extraction]]']
        cfg_lines += ['        skip_optimal = True']
        if args.box_radius is not None: # Boxcar radius
            cfg_lines += ['        boxcar_radius = {0}'.format(args.box_radius)]
        cfg_lines += ['    [[findobj]]']
        cfg_lines += ['        skip_second_find = True']

        # Write
        ofiles = ps.fitstbl.write_pypeit(configs='A', write_bkg_pairs=True, cfg_lines=cfg_lines)
        if len(ofiles) > 1:
            msgs.error("Bad things happened..")

        # Instantiate the main pipeline reduction object
        pypeIt = pypeit.PypeIt(ofiles[0], verbosity=2,
                               reuse_masters=True, overwrite=True,
                               logname='nires_proc_AB.log', show=False)
        # Run
        pypeIt.reduce_all()
        msgs.info('Data reduction complete')
        # QA HTML
        msgs.info('Generating QA HTML')
        pypeIt.build_qa()

        return 0


