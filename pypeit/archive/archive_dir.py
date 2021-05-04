# init:
# takes names of columns
# Takes filename
# Takes whether to overwrite or append
# Takes get_metadata function that takes file location + extra args that are
# Passed to add
# boolean whether to save files

# add
# Calls delegated get_metadata function with file + extra args/kwargs
# Puts row into metadata

# save
# saves metadata

#ArchiveDir
# init takes root dir, list of metadata objects

# add copies file to archive, passes arguments to each metadata object,
# Takes subdirectory to place file

# save calls to metadata

# Not sure I like two metadata/archive stuff....
# Have metadata be seaparate?
# Have onlyh 1 but collate makes a second?
# Combine to one clas...
# 

import os
import shutil

from pypeit import msgs

class ArchiveDir():
    """
    Copies files and metadata to a directory for archival purposes.

    Files are all copied to the top level directory in the archive. 

    If a file originates from KOA the KOAID will be extracted either from the
    ``KOAID`` header keyword or from the ``FILENAME`` header keyword.

    A KOAID has the format: ``II.YYYYMMDD.xxxxx``
        See the `KOA FAQ <https://www2.keck.hawaii.edu/koa/public/faq/koa_faq.php>`_ 
        for more information.

    Metadata is written to two files in the 
    `ipac <https://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html>`_ format. 

    ``by_id_meta.dat`` contains metadata for the spec1d and spec2d files in
    the archive. It is organzied by the id (either KOAID, or file name) of the 
    original science image.

    ``by_object_meta.dat`` contains metadata for the coadded output files.
    This may have multiple rows for each file depending on how many science
    images were coadded. The primary key is a combined key of the source 
    object name, filename, and koaid columns.

    Args:
        archive_root (str): Root directory where the files and metadata will be
        placed. This will be created if needed.

    Attributes:
    by_id_metadata (:obj:`list of :obj:`list of str``): List of metadata rows
        for the spec1d and spec2d files in the archive diredctory. This will be
        loaded from any existing by_id_meta.dat file in the archive_root on 
        object initialization.

    by_object_metadata (:obj:`list of :obj:`list of str``): List of metadata rows
        for the coadd output files in the archive directory. This will be
        loaded from any existing by_id_meta.dat file in the archive_root on 
        object initialization.

    """

    def __init__(self, archive_root, metadata_list, copy_to_archive=True):
        self.archive_root = archive_root
        self.metadata_list = metadata_list        
        self._copy_files = copy_to_archive

    def add(self, items):
        if not isinstance(items, list):
            items = [items]
        for item in items:
            for metadata in self.metadata_list:
                (source_file, dest_file) = metadata.add(item)
                if source_file is not None and dest_file is not None:
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

        if self._copy_files is False:
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
