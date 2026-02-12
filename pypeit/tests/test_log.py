
import io
import logging
from pathlib import Path

from IPython import embed

from pypeit import log
from pypeit.logger import clear_text_color

# TODO:  It's difficult to test the PypeItLogger capturing of the warnings and
# exceptions because pytest overrides them as well.  To test these, we would
# need to create a subprocess and actually run a script.  See
# 
# https://stackoverflow.com/questions/46310034/how-to-test-that-a-custom-excepthook-is-installed-correctly
# 
# Punting for now.

def test_debug():

    # Test the different levels
    logst = io.StringIO()
    log.init(level=logging.DEBUG, stream=logst)

    log.debug('test')
    log.info('test')
    log.warning('test')
    log.error('test')
    log.critical('test')

    msg = clear_text_color(logst.getvalue()).split('\n')[:-1]
    assert len(msg) == 5, f'Incorrect number of stream logs {len(msg)}'
    assert all(['test_log.py:test_debug' in m for m in msg]), \
        'Calling function should be in all messages when using DEBUG level'


def test_info():

    # Test the different levels
    logst = io.StringIO()
    log.init(level=logging.INFO, stream=logst)

    log.debug('test')
    log.info('test')
    log.warning('test')
    log.error('test')
    log.critical('test')

    msg = clear_text_color(logst.getvalue()).split('\n')[:-1]
    assert len(msg) == 4, 'Incorrect number of stream logs (should not print DEBUG message)'
    assert not any(['test_log.py:test_info' in m for m in msg]), \
        'Calling function should not be in messages when using INFO level'


def test_warning():

    # Test the different levels
    logst = io.StringIO()
    log.init(level=logging.WARNING, stream=logst)

    log.debug('test')
    log.info('test')
    log.warning('test')
    log.error('test')
    log.critical('test')

    msg = clear_text_color(logst.getvalue()).split('\n')[:-1]
    assert len(msg) == 3, 'Incorrect number of stream logs (should not print DEBUG/INFO messages)'
    assert not any(['test_log.py:test_warning' in m for m in msg]), \
        'Calling function should not be in messages when using WARNING level'


def test_log_file():

    lf = Path().absolute() / 'test_log.txt'

    # Test the different levels
    logst = io.StringIO()
    log.init(level=logging.DEBUG, stream=logst, log_file=lf)

    log.debug('test')
    log.info('test')
    log.warning('test')
    log.error('test')
    log.critical('test')

    assert lf.is_file(), 'Log file not produced'
    with open(lf, 'r') as f:
        log_lines = f.readlines()

    msg = clear_text_color(logst.getvalue()).split('\n')[:-1]

    assert len(msg) == len(log_lines), 'Number of lines in file should match stream'
    assert all(['test_log.py:test_log_file' in l for l in log_lines]), \
        'Calling function should be included in all file logs'
    
    # Close the log so that we can delete the file
    log.close_file()
    lf.unlink()


def test_log_file_info():

    lf = Path().absolute() / 'test_log.txt'

    # Start the log
    logst = io.StringIO()
    log.init(level=logging.INFO, stream=logst, log_file=lf)
    log.debug('test')
    log.info('test')
    log.warning('test')
    assert lf.is_file(), 'Log file not produced'
    with open(lf, 'r') as f:
        log_lines = f.readlines()

    msg = clear_text_color(logst.getvalue()).split('\n')[:-1]

    assert len(msg) == len(log_lines), 'Number of lines in file should match stream'
    assert all(['test_log.py:test_log_file_info' in l for l in log_lines]), \
        'Calling function should be included in all file logs'
    assert not any(['test_log.py:test_log_file_levels' in m for m in msg]), \
        'Calling function should not be in stream messages when using INFO level'

    # Close the log so that we can delete the file
    log.close_file()
    lf.unlink()


def test_log_file_level_diff():

    lf = Path().absolute() / 'test_log.txt'

    # Start the log
    logst = io.StringIO()
    log.init(level=logging.INFO, stream=logst, log_file=lf, log_file_level=logging.DEBUG)
    log.debug('test')
    log.info('test')
    log.warning('test')
    assert lf.is_file(), 'Log file not produced'
    with open(lf, 'r') as f:
        log_lines = f.readlines()

    msg = clear_text_color(logst.getvalue()).split('\n')[:-1]

    assert len(msg) == 2 and len(log_lines) == 3, \
        'Log file should include all entries, but stream should skip DEBUG message'

    # Close the log so that we can delete the file
    log.close_file()
    lf.unlink()


def test_log_overwrite():

    lf = Path().absolute() / 'test_log.txt'

    # Start the log
    log.init(level=logging.DEBUG, log_file=lf)
    log.debug('test')
    assert lf.is_file(), 'Log file not produced'

    # Reinit, which should restart the file
    log.init(level=logging.DEBUG, log_file=lf)
    log.debug('test')
    assert lf.is_file(), 'Log file not produced'
    with open(lf, 'r') as f:
        log_lines = f.readlines()
    assert len(log_lines) == 1, 'reinitializing the log should overwrite the log file'

    # Close the log so that we can delete the file
    log.close_file()
    lf.unlink()
