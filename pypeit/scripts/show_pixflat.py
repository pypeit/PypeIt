"""
Show on a ginga window the archived pixel flat field image

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase
from IPython import embed


class ShowPixFlat(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Show an archived Pixel Flat image in a ginga window.',
                                    width=width)
        parser.add_argument("file", type=str, help="Pixel Flat filename, e.g. pixelflat_keck_lris_blue.fits.gz")
        parser.add_argument('--det', default=None, type=int, nargs='+',
                            help='Detector(s) to show.  If more than one, list the detectors as, e.g. --det 1 2 '
                                 'to show detectors 1 and 2. If not provided, all detectors will be shown.')
        return parser

    @staticmethod
    def main(args):
        import numpy as np
        from pypeit import msgs
        from pypeit import io
        from pypeit.display import display
        from pypeit import dataPaths

        # check if the file exists
        file_path = dataPaths.pixelflat.get_file_path(args.file, return_none=True)
        if file_path is None:
            msgs.error(f'Provided pixelflat file, {args.file} not found. It is not a direct path, '
                       f'a cached file, or a file that can be downloaded from a PypeIt repository.')

        # Load the image
        with io.fits_open(file_path) as hdu:
            # get all the available detectors in the file
            file_dets = [int(h.name.split('-')[0].split('DET')[1]) for h in hdu[1:]]
            # if detectors are provided, check if they are in the file
            if args.det is not None:
                in_file = np.isin(args.det, file_dets)
                # if none of the provided detectors are in the file, raise an error
                if not np.any(in_file):
                    msgs.error(f"Provided detector(s) not found in the file. Available detectors are {file_dets}")
                # if some of the provided detectors are not in the file, warn the user
                elif np.any(np.logical_not(in_file)):
                    det_not_in_file = np.array(args.det)[np.logical_not(in_file)]
                    msgs.warn(f"Detector(s) {det_not_in_file} not found in the file. Available detectors are {file_dets}")

            # show the image
            display.connect_to_ginga(raise_err=True, allow_new=True)
            for h in hdu[1:]:
                det = int(h.name.split('-')[0].split('DET')[1])
                if args.det is not None and det not in args.det:
                    continue
                display.show_image(h.data, chname=h.name, cuts=(0.9, 1.1), clear=False, wcs_match=True)



