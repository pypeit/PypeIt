"""
Script to measure and correct for flexure in multi-slit data.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import os

from IPython import embed

import numpy as np

from pypeit import msgs
from pypeit import inputfiles
from pypeit.spectrographs.util import load_spectrograph
from pypeit.par import pypeitpar
from pypeit.core import flexure
from pypeit.scripts import scriptbase


#def read_flexfile(ifile):
#    """
#    Read a ``PypeIt`` flexure file, akin to a standard ``PypeIt`` file.
#
#    The top is a config block that sets ParSet parameters.
#
#    Args:
#        ifile (:obj:`str`):
#            Name of the flexure file
#
#    Returns:
#        :obj:`tuple`:  Two objects are returned: a :obj:`list` with the
#        configuration entries used to modify the relevant
#        :class:`~pypeit.par.parset.ParSet` parameters and a :obj:`list` with the
#        names of spec1d files to be flexure corrected.
#    """
#    # Read in the pypeit reduction file
#    msgs.info('Loading the flexure file')
#    lines = inputfiles.read_pypeit_file_lines(ifile)
#    is_config = np.ones(len(lines), dtype=bool)
#
#    # Parse the fluxing block
#    spec1dfiles = []
#    objids_in = []
#    s, e = inputfiles.InputFile.find_block(lines, 'flexure')
#    if s >= 0 and e < 0:
#        msgs.error("Missing 'flexure end' in {0}".format(ifile))
#    elif (s < 0) or (s == e):
#        msgs.error(
#            "Missing flexure read block in {0}. Check the input format for the .flex file".format(ifile))
#    else:
#        for ctr, line in enumerate(lines[s:e]):
#            prs = line.split(' ')
#            spec1dfiles.append(prs[0])
#            if len(prs) > 1:
#                msgs.error('Invalid format for .flex file.' + msgs.newline() +
#                           'You must specify only spec1dfiles in the block ')
#        is_config[s-1:e+1] = False
#
#    # Chck the sizes of the inputs
#    nspec = len(spec1dfiles)
#
#    # Construct config to get spectrograph
#    cfg_lines = list(lines[is_config])
#
#    # Return
#    return cfg_lines, spec1dfiles


# TODO: Maybe not a good idea to name this script the same as the
# flexure.MultiSlitFlexure class, but it is technically okay...
class MultiSlitFlexure(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Calculate and apply flexure corrections for 1D '
                                                'spectra produced by PypeIt.',
                                    width=width, formatter=scriptbase.SmartFormatter)
        parser.add_argument("flex_file", type=str,
                            help="R|File to guide flexure corrections for this multi-slit mode."
                                 "  This file must have the following format: \n\n"
                                 "F|flexure read\n"
                                 "F|  filename\n"
                                 "F|  spec1dfile1\n"
                                 "F|  spec1dfile2\n"
                                 "F|     ...    \n"
                                 "F|flexure end\n"
                                 "\n\n")
        parser.add_argument("outroot", type=str,
                            help='Output fileroot for the flexure fits saved as FITS.')
        parser.add_argument("--clobber", default=True,
                            action="store_true", help="Clobber output files")
        parser.add_argument("--debug", default=False,
                            action="store_true", help="show debug plots?")
        return parser

    @staticmethod
    def main(pargs):

        from astropy.io import fits

        # Load the file
        flexFile = inputfiles.FlexureFile.from_file(pargs.flex_file)

        # Read in spectrograph from spec1dfile header
        header = fits.getheader(flexFile.filenames[0])
        spectrograph = load_spectrograph(header['PYP_SPEC'])

        # Parameters
        spectrograph_def_par = spectrograph.default_pypeit_par()
        par = pypeitpar.PypeItPar.from_cfg_lines(
            cfg_lines=spectrograph_def_par.to_config(), 
            merge_with=(flexFile.cfg_lines,))

        # Loop to my loop
        for filename in flexFile.filenames:
            # Instantiate
            mdFlex = flexure.MultiSlitFlexure(s1dfile=filename)
            # Initalize 
            msgs.info("Setup")
            mdFlex.init(spectrograph, par['flexure'])

            # INITIAL SKY LINE STUFF
            msgs.info("Measuring sky lines")
            mdFlex.measure_sky_lines()

            # FIT SURFACES
            msgs.info("Fitting the surface")
            mdFlex.fit_mask_surfaces()

            # Apply
            msgs.info("Applying flexure correction")
            mdFlex.update_fit()

            # REFIT FOR QA PLOTS
            msgs.info("Generate QA")
            mask = header['TARGET'].strip()
            fnames = header['FILENAME'].split('.')
            root = mask+'_'+fnames[2]
            mdFlex.qa_plots('./', root)

            # Write
            msgs.info("Write to disk")
            mdFlex.to_file(pargs.outroot+root+'.fits',
                           overwrite=pargs.clobber)

            # Apply??

        print("All done!!")


