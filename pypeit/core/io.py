"""
Provides a set of file I/O routines.
"""
import os
import gzip
import shutil

def compress_file(ifile, overwrite=True, remove_uncompressed=True):
    """
    Compress a file using gzip.
    
    Args:
        ifile (:obj:`str`):
            Name of the file to compress.  Should *not* have a '.gz'
            extension.  The output file has the same name as this input
            file, just with '.gz' appended.
        overwrite (:obj:`bool`, optional):
            Overwrite any existing file with the same name.
        remove_uncompressed(:obj:`bool`, optional):
            Remove the uncompressed file once the compressed version has
            been created.

    Raises:
        ValueError:
            Raised if the input file has a 'gz' extension.
        FileExistsError:
            Raised if the file already exists and `overwrite` is False.
    """
    if ifile.split('.')[-1] == 'gz':
        raise ValueError('File appears to already have been compressed! {0}'.format(ifile))

    ofile = '{0}.gz'.format(ifile)
    if os.path.isfile(ofile) and not clobber:
        raise FileExistsError('File already exists: {0}.\nTo overwrite, set clobber=True.'.format(
                                ofile))

    with open(ifile, 'rb') as f_in:
        with gzip.open(ofile, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    if remove_uncompressed:
        os.remove(ifile)

