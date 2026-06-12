"""
Module to run tests on the PypeIt Dashboard helpers.

These are lightweight, deterministic unit tests that do not require a display
or a running reduction.  They cover the setup-file parser, the STEP-message
parser, and the custom STEP logging level used to drive the dashboard status
bar.
"""
import io
import logging

import pytest

# qtpy/pyqt6 are dependencies, but skip cleanly if the GUI stack is absent.
pytest.importorskip("qtpy")

from pypeit import log
from pypeit.pkg import logger as pypeit_logger
from pypeit.dashboard.dashboard import (
    parse_step_message,
    parse_pypeit_setup_file,
)


# Minimal modern-format .pypeit file used by the parser test.
PYPEIT_FILE_TEXT = """\
# Auto-generated PypeIt file
[rdx]
spectrograph = shane_kast_blue
[baseprocess]
  use_biasimage = False

# Read in the data
data read
 path /path/to/RAW_DATA/shane_kast_blue/600_4310_d55
|     filename |       frametype | exptime |
| b1.fits.gz   |       arc,tilt  |    15.0 |
| b10.fits.gz  | illumflat,trace |    30.0 |
| b27.fits.gz  |         science |  1200.0 |
| b24.fits.gz  |        standard |   180.0 |
data end
"""


def test_parse_step_message():
    # The logging prefix and unknown tokens must be ignored, and keys/values
    # stripped of surrounding whitespace.
    msg = " STEP | calib_id=1|det=2|calib_step=arc|meta_step=Calibrations"
    fields = parse_step_message(msg)
    assert fields == {
        "calib_id": "1",
        "det": "2",
        "calib_step": "arc",
        "meta_step": "Calibrations",
    }
    # A message with no key=value tokens yields an empty dict.
    assert parse_step_message(" STEP | starting up") == {}


def test_parse_pypeit_setup_file(tmp_path):
    # Regression test for the infinite-loop bug: this must return promptly.
    pfile = tmp_path / "shane_kast_blue_600_4310_d55.pypeit"
    pfile.write_text(PYPEIT_FILE_TEXT)

    spectrograph, raw_path, files, science_file = \
        parse_pypeit_setup_file(str(pfile))

    assert spectrograph == "shane_kast_blue"
    assert raw_path.endswith("shane_kast_blue/600_4310_d55")
    # Four data rows, parsed as (filename, frametype) tuples.
    assert len(files) == 4
    assert ("b1.fits.gz", "arc,tilt") in files
    # The first science frame is returned.
    assert science_file == "b27.fits.gz"


def test_parse_pypeit_setup_file_missing_keywords(tmp_path):
    # A file with no spectrograph/path/data rows must not hang and must
    # return None/empty rather than raising.
    pfile = tmp_path / "empty.pypeit"
    pfile.write_text("# nothing useful here\n# just comments\n")
    spectrograph, raw_path, files, science_file = \
        parse_pypeit_setup_file(str(pfile))
    assert spectrograph is None
    assert raw_path is None
    assert files == []
    assert science_file is None


def test_step_level_registered():
    # The custom level exists and round-trips to the name "STEP".
    assert pypeit_logger.STEP == 15
    assert logging.getLevelName(pypeit_logger.STEP) == "STEP"


def test_log_step_emits_at_step_level():
    # A STEP record is produced with levelname "STEP" and the given message.
    captured = []

    class _Capture(logging.Handler):
        def emit(self, record):
            captured.append((record.levelname, record.getMessage()))

    handler = _Capture()
    log.addHandler(handler)
    try:
        log.step("calib_id=1|det=1|calib_step=bias")
    finally:
        log.removeHandler(handler)

    assert ("STEP", "calib_id=1|det=1|calib_step=bias") in captured


def test_log_step_silent_on_info_console():
    # STEP sits below INFO, so it must not appear on an INFO-level console.
    logst = io.StringIO()
    log.init(level=logging.INFO, stream=logst)
    log.step("calib_id=1|det=1|calib_step=bias")
    assert "calib_id" not in logst.getvalue()
