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

from astropy.io import ascii
from astropy.table import Table

class ArchiveMetadata():
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

    # Header and SpecObj keys for metadata needed for the IPAC files

    def __init__(self, metadata_file, col_names, get_metadata_func, append):
        self.metadata_file = metadata_file
        self.col_names = col_names
        self.get_metadata_func = get_metadata_func        
        
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
        (rows, orig_file, dest_file) = self.get_metadata_func(item)

        if rows is not None:
            self._metadata += rows

        return orig_file, dest_file

    def save(self):
        """
        Saves the metadata in this class to IPAC files in the archive directory
        """
        if len(self._metadata) > 0:
            with open(self.metadata_file, 'w') as f:
                ascii.write(Table(rows=self._metadata, names=self.col_names), f, format='ipac')

