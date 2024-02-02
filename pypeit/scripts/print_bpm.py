"""
This script prints a user-friendly description of the bad pixel mask
based on a spec2d file.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
from IPython import embed

from pypeit import msgs
from pypeit.images.imagebitmask import ImageBitMask
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
        return parser

    @staticmethod
    def main(args):

        # Convert the integer bitmask value to a list of binary numbers
        binvals = [int(x) for x in bin(args.bit)[2:]][::-1]

        # Generate an Image BitMask object
        bpm = ImageBitMask()

        # Print the description of the bad pixel mask value
        outstr = "This bad pixel mask value corresponds to the following:" + msgs.newline() + msgs.newline()
        for i in range(len(binvals)):
            if binvals[i] == 1:
                outstr += "* " + bpm.descr[i] + msgs.newline()

        # Print the list of bad pixel mask value descriptions
        msgs.info(outstr)

        # Finally, print out a message to point users to the online documentation
        msgs.info("Please see the following website for more information:" + msgs.newline() +
                  "https://pypeit.readthedocs.io/en/release/out_masks.html")
