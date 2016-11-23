""" Module for basic utils outside the PYPIT domain, i.e. msgs is not used
"""
from __future__ import absolute_import, division, print_function


def make_pypit_file(pyp_file, spectrograph, files_root,
                    datfil_extension='.fits'):
    """ Generate a default PYPIT settings file

    Parameters
    ----------
    pyp_file : str
      Name of PYPIT file to be generated
    spectrograph : str
    files_root : list
      Path + file root of datafiles
    datfil_extension : str, optional
      Extension of data file, e.g.  '.fits', '.fits.gz'

    Returns
    -------
    Creates a PYPIT File

    """
    # Error checking
    if not isinstance(files_root, list):
        raise IOError("files_root needs to be a list")

    # Here we go
    with open(pyp_file, 'w') as f:
        f.write("# This is a comment line\n")
        f.write("\n")
        f.write("# Change the default settings\n")
        f.write("run ncpus 1\n")
        f.write("run calcheck True\n")  # This is the key line here
        f.write("run spectrograph {:s}\n".format(spectrograph))
        f.write("output overwrite True\n")
        #f.write("output sorted {:s}\n".format(root))
        f.write("\n")
        f.write("# Reduce\n")
        f.write("\n")
        f.write("# Read in the data\n")
        f.write("data read\n")
        for file_root in files_root:
            f.write(" {:s}*{:s}*\n".format(file_root, datfil_extension))
        f.write("data end\n")
        f.write("\n")
        f.write("spect read\n")
        f.write(" pixelflat number 0\n")
        f.write(" arc number 1\n")
        f.write(" slitflat number 0\n")
        f.write(" bias number 0\n")
        f.write(" standard number 0\n")
        f.write("spect end\n")
