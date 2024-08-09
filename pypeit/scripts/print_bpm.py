"""
This script prints a user-friendly description of the bad pixel mask
based on a spec2d file.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
from IPython import embed

from astropy.io import fits

from pypeit import __version__
from pypeit import msgs, spec2dobj
from pypeit.images.detector_container import DetectorContainer
from pypeit.images.imagebitmask import ImageBitMask
from pypeit.pypmsgs import PypeItDataModelError
from pypeit.scripts import scriptbase


class PrintBPM(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Print out an informative description of a '
                                                'bad pixel masked value. Usually, you should '
                                                'run pypeit_show_2dspec --showmask first to '
                                                'see the bad pixel mask values. Then, call this '
                                                'script with the BPM value that you want to find'
                                                'more information about.',
                                    width=width)

        parser.add_argument('bit', type=int, default=None, help='Bad pixel mask value to describe in plain text')
        parser.add_argument('--file', type=str, default=None, help='PypeIt spec2d file to use for the description'
                                                                 '(optional). If provided, the bitmask contained '
                                                                 'in the spec2d file will be used to describe the '
                                                                 'bad pixel mask value. If not provided, the default '
                                                                 'pypeit bad pixel mask will be used.')
        parser.add_argument('--det', default='1', type=str,
                            help='Detector name or number.  If a number, the name is constructed '
                                 'assuming the reduction is for a single detector.  If a string, '
                                 'it must match the name of the detector object (e.g., DET01 for '
                                 'a detector, MSC01 for a mosaic). This is not required, and the '
                                 'value is acceptable.  Default is 1.')
        return parser

    @staticmethod
    def main(args):

        # Convert the integer bitmask value to a list of binary numbers
        binvals = [int(x) for x in bin(args.bit)[2:]][::-1]

        if args.file is None:
            msgs.info("Using the default PypeIt bad pixel mask.")
            # Generate an Image BitMask object
            bpm = ImageBitMask()
            descr = bpm.descr
        else:
            # Read the spec2d file
            msgs.info("Using the bad pixel mask from the following spec2d file:" + msgs.newline() + f"{args.file}.")
            spec2d_file = args.file

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
                # Generate an Image BitMask object
                msgs.info("Using the default PypeIt bad pixel mask.")
                bpm = ImageBitMask()
                descr = bpm.descr
            else:
                bpm = spec2DObj.bpmmask
                descr = bpm.bitmask.descr

        # Print the description of the bad pixel mask value
        outstr = f"The bad pixel mask value ({args.bit}) corresponds to the following:" \
                 + msgs.newline() + msgs.newline()
        bitkeys = list(bpm.bits.keys())
        # Pad the bit keys with spaces so that they all have the same length
        bitlen = len(max(bitkeys, key=len))
        for i in range(len(binvals)):
            if binvals[i] == 1:
                outstr += f"* {bitkeys[i].ljust(bitlen)} : {descr[i]}" + msgs.newline()

        # Print the message to the user
        msgs.info(outstr)

        # Finally, print out a message to point users to the online documentation
        msgs.info("Please see the following website for more information:" + msgs.newline() +
                  "https://pypeit.readthedocs.io/en/release/out_masks.html")
