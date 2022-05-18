""" Tests on PypeItFile """

import pytest

def test_read_pypeit_file():
    # Read the PypeIt file
    cfg, data, frametype, usrdata, setups, _ \
            = util.parse_pypeit_file(data_path('example_pypeit_file.pypeit'), file_check=False)