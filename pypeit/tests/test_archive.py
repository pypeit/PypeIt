"""
Module to run tests on the ArchiveDir/ArchiveMetadata script.
"""

import pytest
import os
from functools import partial
import filecmp

from pypeit.archive import ArchiveMetadata, ArchiveDir
from pypeit.tests.tstutils import data_path


class mock_file_info:
    def __init__(self, path, name):
        self.path = path
        self.name = name

_DATA_FILE1 = "spec1d_b27.fits"
_DATA_FILE2 = "spec1d_b28.fits"
def get_simple_metadata(file_info):
    if isinstance(file_info, mock_file_info):
        return (None, None)

    dest_file = os.path.basename(file_info)
    metadata = {_DATA_FILE1: [1, dest_file, 1.234, "2021-01-01"],
                _DATA_FILE2: [2, dest_file, 1.234, "2021-01-02"] }

    return ([metadata[dest_file]], [(file_info, dest_file)])

def get_multirow_metadata(file_info):
    if isinstance(file_info, str):
        return (None, None)

    dest_file = file_info.name
    metadata = {_DATA_FILE1: [[1, dest_file, "part1"],
                                [1, dest_file, "part2"]],
                _DATA_FILE2: [[2, dest_file, "part1"],
                                [2, dest_file, "part2"],
                                [2, dest_file, "part3"]]}

    return (metadata[dest_file], [(os.path.join(file_info.path, file_info.name), dest_file)])

def cmp_files(file1, file2):
    """
    Compre two text files line by line. This function exists because
    out test suite has to work under Windows and Linux, and the filecmp
    package doesn't handle mixed newline styles.

    Returns:
    True if the files are the same, False otherwise.
    """
    file1_lines =[]
    file2_lines =[]

    with open(file1, "r") as f1:
        file1_lines = [x.rstrip() for x in f1.readlines() if x.rstrip() != ""]

    with open(file2, "r") as f2:
        file2_lines = [x.rstrip() for x in f2.readlines() if x.rstrip() != ""]

    if len(file1_lines) == len(file2_lines):
        for i in range(len(file1_lines)):
            if file1_lines[i] != file2_lines[i]:
                return False
    else:
        return False

    return True

def test_archive_meta(tmp_path):
    test_meta_path = str(tmp_path / "test_meta.dat")

    col_names = ["id", "file", "num", "date"]

    orig_file1 = os.path.join("orig_path", _DATA_FILE1)
    dest_file1 = _DATA_FILE1
    orig_file2 = os.path.join("orig_path", _DATA_FILE2)
    dest_file2 = _DATA_FILE2

    expected_rows = [ [1, dest_file1, 1.23, '2021-01-01'],
                      [2, dest_file2, 1.234, '2021-01-02']]

    # Test creating a new file, and adding a new row to it, with formatting
    test_meta1 = ArchiveMetadata(test_meta_path, col_names, get_simple_metadata, True, formats={"num": "%.2f"})

    files_to_add = test_meta1.add(orig_file1)
    assert files_to_add[0][0] == orig_file1
    assert files_to_add[0][1] == dest_file1

    test_meta1.save()

    # Test loading an existing file, and adding a second row to it, without formatting
    test_meta2 = ArchiveMetadata(test_meta_path, col_names, get_simple_metadata, True)
    (orig_file, dest_file) = test_meta2.add(orig_file2)[0]
    assert orig_file == orig_file2
    assert dest_file == dest_file2

    test_meta2.save()

    # Test loading as file, and making sure all previous rows are still there
    test_meta3 = ArchiveMetadata(test_meta_path, col_names, get_simple_metadata, True)

    assert test_meta3._metadata == expected_rows

def test_archive_dir(tmp_path):

    archive_root = str(tmp_path)
    metadata_file1 = os.path.join(archive_root, "by_id_meta.dat")
    metadata_file2 = os.path.join(archive_root, "by_object_meta.dat")
    source_path = data_path("")

    # Test an Archiver with two metadta files, one with multiple rows per item
    col_names1 = ["id", "file", "num", "date"]
    test_meta1 = ArchiveMetadata(metadata_file1, col_names1, get_simple_metadata, True)

    col_names2 = ["id", "file", "name"]
    test_meta2 = ArchiveMetadata(metadata_file2, col_names2, get_multirow_metadata, True)

    test_meta3 = ArchiveMetadata(metadata_file1, col_names1, get_simple_metadata, False)

    # Test an archvier that doesn't copy files.
    archive_dir_nocopy = ArchiveDir(archive_root, [test_meta1, test_meta2], False)
    archive_dir_nocopy.add(os.path.join(source_path, _DATA_FILE1))
    archive_dir_nocopy.save()

    assert not os.path.exists(os.path.join(archive_root, _DATA_FILE1))

    # Test an archiver that does copy. This uses a metadta object that won't
    # append to the metadata from the previous archive_dir_nocopy save
    archive_dir = ArchiveDir(archive_root, [test_meta3, test_meta2], True)

    archive_dir.add(os.path.join(source_path, _DATA_FILE1))
    archive_dir.add(mock_file_info(source_path, _DATA_FILE2))
    archive_dir.save()

    # Verify files were copied
    dest_file1 = tmp_path / _DATA_FILE1
    dest_file2 = tmp_path / _DATA_FILE2
    assert dest_file1.exists()
    assert dest_file2.exists()

    # Verify ipac metadata files are correct
    good_path = data_path("ipac")
    good_file1 = os.path.join(good_path, 'by_id_meta.dat')
    good_file2 = os.path.join(good_path, 'by_object_meta.dat')

    assert cmp_files(good_file1, metadata_file1) is True
    assert cmp_files(good_file2, metadata_file2) is True
