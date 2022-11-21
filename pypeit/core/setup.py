""" Simple routines for setting up PypeIt from a set of files """

from genericpath import isfile
import os
import glob

from pypeit import inputfiles

from IPython import embed

def grab_rawfiles(raw_paths:list=None, 
               file_of_files:str=None,
               list_of_files:list=None,
               extension:str='.fits'):
    """ Grab the list of raw files

    Args:
        raw_paths (list, optional): 
            List of paths to raw files
        file_of_files (str, optional): 
            File with list of raw files. PypeIt input file formatted
        list_of_files (list, optional): 
            List of raw files.  These are combined with raw_paths
        extension (str, optional): 
            File extension to search on.  Defaults to '.fits'.

    Returns:
        list: List of raw data filenames with full path
    """

    # Which option>
    if file_of_files is not None: # PypeIt formatted list of files
        iRaw = inputfiles.RawFiles.from_file(file_of_files)
        data_files = iRaw.filenames
    elif list_of_files is not None: # An actual list
        data_files = []
        for raw_path in raw_paths:
            for ifile in list_of_files:
                if os.path.isfile(os.path.join(
                    raw_path, ifile)):
                    data_files.append()
    else: # Search in raw_paths for files with the given extension
        data_files = []
        for raw_path in raw_paths:
            data_files += files_from_extension(raw_path, 
                                               extension)
    # Finish
    return data_files
    

def files_from_extension(raw_path:str,
                         extension:str='fits'):
    """ Grab the list of files with a given extension 

        Args:
            raw_path (str):
                Path to raw files
            extension (str, optional):
                File extension to search on.  Defaults to '.fits'.

        Returns:
            list: List of raw data filenames (sorted) with full path
    """
    # Grab the list of files
    dfname = os.path.join(raw_path, '*{0}*'.format(extension)) \
                if os.path.isdir(raw_path) else '{0}*{1}*'.format(raw_path, extension)
    data_files = glob.glob(dfname)
    data_files.sort()

    # Return
    return data_files

