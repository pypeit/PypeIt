"""
This script prints a user-friendly description of the bad pixel mask
based on a spec2d file.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import os

import numpy as np

from IPython import embed

from astropy.io import fits
from astropy.stats import sigma_clipped_stats

from pypeit import msgs
from pypeit import specobjs
from pypeit import io
from pypeit import utils
from pypeit import __version__
from pypeit.pypmsgs import PypeItError, PypeItDataModelError

from pypeit.display import display
from pypeit.images.imagebitmask import ImageBitMask
from pypeit.images.detector_container import DetectorContainer
from pypeit import spec2dobj
from pypeit.scripts import scriptbase


class Show2DSpec(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Print out an informative description of a '
                                                'bad pixel masked value. Usually, you should '
                                                'run pypeit_show_2dspec --showmask first to '
                                                'see the bad pixel mask values. Then, call this '
                                                'script with the BPM value that you want to find'
                                                'more information about.',
                                    width=width)

        parser.add_argument('file', type=str, default=None, help='Path to a PypeIt spec2d file')
        parser.add_argument('bpm', type=int, default=None, help='Bad pixel mask value to describe in plain text')
        parser.add_argument('--det', default='1', type=str,
                            help='Detector name or number.  If a number, the name is constructed '
                                 'assuming the reduction is for a single detector.  If a string, '
                                 'it must match the name of the detector object (e.g., DET01 for '
                                 'a detector, MSC01 for a mosaic). Not essential, as the BPM is '
                                 'stored in the header, independently of the detector.')

        return parser

    @staticmethod
    def main(args):

        # List only?
        if args.list:
            io.fits_open(args.file).info()
            return

        # Set the verbosity, and create a logfile if verbosity == 2
        msgs.set_logfile_and_verbosity('print_bpm', args.verbosity)

        # Parse the detector name
        try:
            det = int(args.det)
        except:
            detname = args.det
        else:
            detname = DetectorContainer.get_name(det)

        # Try to read the Spec2DObj using the current datamodel, but allowing
        # for the datamodel version to be different
        try:
            spec2DObj = spec2dobj.Spec2DObj.from_file(args.file, detname, chk_version=False)
        except PypeItDataModelError:
            try:
                # Try to get the pypeit version used to write this file
                file_pypeit_version = fits.getval(args.file, 'VERSPYP', 0)
            except KeyError:
                file_pypeit_version = '*unknown*'
            msgs.warn(f'Your installed version of PypeIt ({__version__}) cannot be used to parse '
                      f'{args.file}, which was reduced using version {file_pypeit_version}.  You '
                      'are strongly encouraged to re-reduce your data using this (or, better yet, '
                      'the most recent) version of PypeIt.  Script will try to parse only the '
                      'relevant bits from the spec2d file and continue (possibly with more '
                      'limited functionality).')
            spec2DObj = None

        if spec2DObj is None:
            # Try to get the relevant elements directly from the fits file
            with io.fits_open(args.file) as hdu:
                names = [h.name for h in hdu]
                _ext = f'{detname}-BPMMASK'
                bpmmask = hdu[_ext].data if _ext in names else None
        else:
            # Use the parsed SpecObjs object
            bpmmask = spec2DObj.bpmmask

        # Finally, print out a message to point users to the online documentation
        msgs.info("Please see the following website for more information:" + msgs.newline() +
                  "https://pypeit.readthedocs.io/en/release/out_masks.html")
