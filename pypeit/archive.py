"""
Module for abstracting common code used in archiving files and metadata.

This module splits archiving functionality into ArchiveMetadata and ArchiveDir
classes. ArchiveMetadata is responsible for archiving metadata from objects
and ArchiveDir is responsible for archiving files.

What types of files and what metadata is archived is delegated to a
``get_metdadata_func`` function provided by the caller. It's signature is:

    Args:
        object: 
            The object to be archived.

    Returns:
        tuple: metadata, files_to_copy

        metadata (list of list of str):
            The metadata rows read from the object. One object can result
            in multiple rows.  Each row is a list of strings.

        files_to_copy (iterable):
            An iterable of tuples. Each tuple has a src file to copy to the archive
            and a destination pathname for that file in the archive. The destination
            pathname is relative to the archive's root directory.

Below is an example ``get_metadata_func`` function that gets metadata for a single fits file and
archives it in a directory based on observation date::

    import os.path
    from astropy.io import fits
    from astropy.time import Time

    def get_target_metadata(file_info):

        header = fits.getheader(file_info)

        # Determine the path within the archive to store the file. The subdir
        # has the format YYYYMM based on the observation date.
        mjd = header['MJD']
        subdir_name = Time(mjd, format='mjd').strftime("%Y%m")
        dest_pathname = os.path.join(subdir_name, os.path.basename(file_info))

        # Return a single data row with a column for file name, date and target
        data_rows = [ [dest_pathname, mjd, header['TARGET']] ]

        return (data_rows, [(file_info, dest_pathname)])


.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""

import os
import shutil

from astropy.io import ascii
from astropy.table import Table

from pypeit import msgs


class ArchiveMetadata():
    """
    Reads metadata from objects and writes it to a file for archiving. This
    class can be used on its own for saving metadata or passed to an 
    ArchiveDir object for saving files and metadata. Currently the files
    are saved in `ipac <https://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html>`_ 
    format. 


    Args:
        metadata_file (str):
            File the metadata should be written to.

        col_names (list of str): 
            The column names of the metadata

        get_metadata_func (func): 
            Function that reads metadata and file information from the objects
            being archived.

        append (bool):
            If true append new metadata to an existing file. If this is false
            any existing file will be overwritten.

        formats (dict, optional):
            Dictionary mapping column names to a format. The format 
            can be any format accepted by a Astropy Table `Column <https://docs.astropy.org/en/stable/api/astropy.table.Column.html#astropy.table.Column>`_.                       
    """

    def __init__(self, metadata_file, col_names, get_metadata_func, append, formats={}):
        self.metadata_file = metadata_file
        self.col_names = col_names
        self.get_metadata_func = get_metadata_func
        self.formats = formats        
        
        # Load metadata from any pre-existing metadata file.
        # Because astropy Tables are slow at adding rows, we convert 
        # the metadata to a list of lists for performance as adding rows is
        # the primary feature of this class.
        if append and os.path.exists(self.metadata_file):
            table = ascii.read(self.metadata_file)
            self._metadata = [list(row) for row in table]
        else:
            self._metadata = []

    def add(self, item):
        """Adds an item to the ArchiveMetadata.
        
        Args:
            item: The object to add to the archive.

        Returns: 
            iterable: An iterable of tuples, as returned by ``get_metadata_func``. 
            Each tuple has a src file to copy to the archive
            and a destination pathname for that file in the archive. The
            destination pathname is relative to the archive's root directory.
        """        
        (rows, files_to_copy) = self.get_metadata_func(item)

        if rows is not None:
            self._metadata += rows

        return files_to_copy

    def save(self):
        """
        Saves the metadata in this class to IPAC files in the archive directory
        """
        if len(self._metadata) > 0:
            with open(self.metadata_file, 'w') as f:
                t = Table(rows=self._metadata, names=self.col_names)
                
                # Set column formats before writing
                for key in self.formats.keys():
                    t[key].format = self.formats[key]

                ascii.write(t, f, format='ipac')


class ArchiveDir():
    """
    Copies files to a directory for archival purposes.

    One or more ArchiveMetadata objects are associated with this
    ArchiveDir, and those metadata files are also stored in the archive.
    These objects also translate the passed in objects to filenames to
    copy into the archive.
        
    Args:
        archive_root (str): 
            Root directory where the files and metadata will be placed. 
            This will be created if needed.

        metadata_list (:obj:`ArchiveMetadata`): 
            One or more ArchiveMetadata objects responsible for translating
            objects passed to the archive to paths for copying.

    """

    def __init__(self, archive_root, metadata_list, copy_to_archive=True):
        self.archive_root = archive_root
        self.metadata_list = metadata_list        
        self._copy_files = copy_to_archive

    def add(self, items):
        """
        Adds items to the archive.

        Args:
            items (object, or iterable): 
                The item or iterable object of items to add to the archive. The items
                will be passed to a ArchiveMetadata object so that 
                metadata and file names can be read from it.
        """
        # Allow any iterable object, or a single object, by testing the "iter"
        # built in function. If it fails, then it's a single object
        try:
            iter(items)

            # It's iterable, it could also be a str which we want to treat it as a single object
            if isinstance(items, str):
                items = [items]
                
        except TypeError:
            # A single object, wrap it in a list
            items = [items]

        for item in items:
            for metadata in self.metadata_list:
                files_to_copy = metadata.add(item)
                if files_to_copy  is not None:
                    for source_file, dest_file in files_to_copy:
                        self._archive_file(source_file, dest_file)





    def save(self):
        """
        Saves the metadata in this class to IPAC files in the archive directory
        """
        for metadata in self.metadata_list:
            metadata.save()

    def _archive_file(self, orig_file, dest_file):
        """
        Copies a file to the archive directory, if copying
        is enable.

        Args:
            orig_file (str): Path to the file to copy.

            dest_file (str): Relative path to the file to copy.

        Returns:
            str: The full path to the new copy in the archive.
        """

        if self._copy_files is False or orig_file is None:
            return orig_file

        if not os.path.exists(orig_file):
            msgs.error(f'File {orig_file} does not exist')

        full_dest_path = os.path.join(self.archive_root, dest_file)
        os.makedirs(os.path.dirname(full_dest_path), exist_ok=True)

        msgs.info(f'Copying {orig_file} to archive root {self.archive_root}')
        try:
            shutil.copy2(orig_file, full_dest_path)
        except:
            msgs.error(f'Failed to copy {orig_file} to {full_dest_path}')

        return full_dest_path
