# -*- coding: utf-8 -*-
"""
Utility functions for PypIt parameter sets
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import glob
import warnings
import textwrap
import sys
if sys.version > '3':
    long = int
from pkg_resources import resource_filename

try:
    basestring
except NameError:
    basestring = str

import numpy as np

from configobj import ConfigObj

from pypeit import msgs

#-----------------------------------------------------------------------
# Parameter utility functions
#-----------------------------------------------------------------------
# TODO: This should go in a different module, or in __init__
def pypeit_root_directory():
    """
    Get the root directory for the PYPIT source distribution.

    .. todo::
        - Set this in __init__.py

    Returns:
        str: Root directory to PYPIT

    Raises:
        OSError: Raised if `pkg_resources.resource_filename` fails.
    """
    try:
        # Get the directory with the pypeit source code
        code_dir = resource_filename('pypeit', '')
    except:
        # TODO: pypeit should always be installed as a package, so is
        # this try/except block necessary?
        raise OSError('Could not find PYPIT package!')
    # Root directory is one level up from source code
    return os.path.split(code_dir)[0]


def _eval_ignore():
    """Provides a list of strings that should not be evaluated."""
    return [ 'open', 'file', 'dict' ]


def recursive_dict_evaluate(d):
    """
    Recursively run :func:`eval` on each element of the provided
    dictionary.

    A raw read of a configuration file with `ConfigObj` results in a
    dictionary that contains strings or lists of strings.  However, when
    assigning the values for the various ParSets, the `from_dict`
    methods expect the dictionary values to have the appropriate type.
    E.g., the ConfigObj will have something like d['foo'] = '1', when
    the `from_dict` method expects the value to be an integer (d['foo']
    = 1).

    This function tries to evaluate *all* dictionary values, except for
    those listed above in the :func:`_eval_ignore` function.  Any value
    in this list or where::

        eval(d[k]) for k in d.keys()

    raises an exception is returned as the original string.

    This is currently only used in :func:`PypitPar.from_cfg_file`; see
    further comments there.

    Args:
        d (dict):
            Dictionary of values to evaluate

    Returns:
        dict: Identical to input dictionary, but with all string values
        replaced with the result of `eval(d[k])` for all `k` in
        `d.keys()`.
    """
    ignore = _eval_ignore()
    for k in d.keys():
        if isinstance(d[k], dict):
           d[k] = recursive_dict_evaluate(d[k])
        elif isinstance(d[k], list):
            replacement = []
            for v in d[k]:
                if v in ignore:
                    replacement += [ v ]
                else:
                    try:
                        replacement += [ eval(v) ]
                    except:
                        replacement += [ v ]
            d[k] = replacement
        else:
            try:
                d[k] = eval(d[k]) if d[k] not in ignore else d[k]
            except:
                pass

    return d


# NOT USED and causes python2 to barf in Travis
#def merge_two_dicts(a, b):
#    """
#    Merge two dictionaries.
#
#    Values in the second dictionary take precedence over the values in
#    the first.
#
#    Args:
#        a (:obj:`dict`):
#            Dictionary with values to be updated.
#        b (:obj:`dict`):
#            Dictionary to use to replace values in ``a``.
#
#    Returns:
#        :obj:`dict`: Dictionary with the merged values.
#    """
#    try:
#        # This is a python > 3.5 thing...
#        return {**a, **b}
#    except:
#        # Assume we have to take the longer route
#        c = a.copy()
#        c.update(b)
#        return c


def get_parset_list(cfg, pk, parsetclass):
    """
    Create a list of ParSets based on a root keyword for a set of
    defined groups in the configuration file.
    
    For example, the :class:`InstrumentPar` group allows for a list of
    detectors (:class:`DetectorPar`) with keywords like `detector1`,
    `detector2`, etc.  This function parses the provided configuration
    object (`cfg`) to find any sections with `detector` (`pk`) as its
    root.  The remainder of the section name must be able to be
    converted to an integer and the section itself must be able to setup
    an instance of `parsetclass`.  The sections must be number
    sequentially from 1..N.  E.g., the :class:`InstrumentPar`
    configuration file cannot have `dectector1` and `detector3`, but no
    `detector2`.  The call to setup the detectors in the
    :class:`InstrumentPar` is::

        kwargs['detector'] = get_parset_list(cfg, 'detector', DetectorPar)

    Args:
        cfg (:class:`ConfigObj`, :obj:`dict`):
            The top-level configuration that defines a list of
            sub-ParSets.
        pk (str):
            The root of the keywords used to set a list of sub-ParSets.
        parsetclass (:class:`pypeit.par.parset.ParSet`):
            The class used to construct each element in the list of
            parameter subsets.  The class **must** have a `from_dict`
            method that instantiates the
            :class:`pypeit.par.parset.ParSet` based on the provide
            subsection/subdict from cfg.

    Returns:
        list: A list of instances of `parsetclass` parsed from the
        provided configuration data.

    Raises:
        ValueError:
            Raised if the indices of the subsections are not sequential
            and 1-indexed.
    """
    # Get the full list of keys
    k = cfg.keys()

    # Iterate through the list of keys to find the appropriate sub
    # parameter sets and their order.
    par = []
    order = []
    for _k in k:
        if _k == pk and cfg[_k] is None:
            continue
        if pk in _k:
            try:
                # Get the order for this subgroup (e.g., 2 for
                # 'detector2'
                order += [ int(_k.replace(pk,'')) ]
                # And instantiate the parameter set
                par += [ parsetclass.from_dict(cfg[_k]) ]
            except:
                continue

    if len(par) > 0:
        # Make sure the instances are correctly sorted and sequential
        srt = np.argsort(order)
        if np.any(np.array(order)[srt]-1 != np.arange(order[srt[-1]])):
            raise ValueError('Parameter set series must be sequential and 1-indexed.')
        # Return the sorted instances
        return [par[i] for i in srt]

    # No such subsets were defined, so return a null result
    return None


def parset_to_dict(par):
    """
    Convert the provided parset into a dictionary.
    """
    try:
        d = dict(ConfigObj(par.to_config(section_name='tmp'))['tmp'])
    except:
        d = dict(ConfigObj(par.to_config()))
    return recursive_dict_evaluate(d)


#-----------------------------------------------------------------------
# Functions for parsing the input pypeit file
# TODO: Should these go into a different module? PypitSetup?
#-----------------------------------------------------------------------
def _read_pypeit_file_lines(ifile):
    """
    General parser for a pypeit file.

    - Checks that the file exists.
    - Reads all the lines in the file
    - Removes comments, empty lines, and replaces special characters.
    
    Applies to settings, setup, and user-level reduction files.

    Args:
        ifile (str): Name of the file to parse.

    Returns:
        :obj:`numpy.ndarray`: Returns a list of the valid lines in the
        files.
    """
    # Check the files
    if not os.path.isfile(ifile):
        msgs.error('The filename does not exist -' + msgs.newline() + ifile)

    # Read the input lines and replace special characters
    with open(ifile, 'r') as f:
        lines = np.array([l.replace('\t', ' ').replace('\n', ' ').strip() \
                                for l in f.readlines()])
    # Remove empty or fully commented lines
    lines = lines[np.array([ len(l) > 0 and l[0] != '#' for l in lines ])]
    # Remove appended comments and return
    return np.array([ l.split('#')[0] for l in lines ])


def _find_pypeit_block(lines, group):
    """Find the start and end of a pypeit file block."""
    start = -1
    end = -1
    for i, l in enumerate(lines):
        entries = l.split()
        if start < 0 and entries[0] == group and entries[1] == 'read':
            start = i+1
            continue
        if entries[0] == group and entries[1] == 'end':
            end = i
            continue
        if start >= 0 and end >= 0:
            break
    return start, end


def _parse_data_file_name(inp, current_path):
    """Expand the data file name as necessary."""
    out = os.path.expanduser(inp) if inp[0] == '~' else inp
    if current_path is not None:
        out = os.path.join(current_path, out)
    return glob.glob(out)
    

def _read_data_file_names(lines, file_check=True):
    """Read the raw data file format."""
    # Pass through all the lines and:
    #   - Determine if a path is set
    #   - Gather the files to skip, can include wildcards
    #   - Gather the files to read, can include wildcards
    current_path = None
    skip_inp = []
    read_inp = []
    for l in lines:
        _l = l.split(' ')

        if _l[0] == 'skip':
            space_ind = l.index(" ")
            path = l[space_ind + 1:]
            skip_inp += _parse_data_file_name(path, current_path)
            continue

        if _l[0] == 'path':
            space_ind = l.index(" ")
            current_path = l[space_ind + 1:]
            continue

        read_inp += _parse_data_file_name(l, current_path)

    # Remove any repeated lines
    if len(skip_inp) > 0 and len(skip_inp) != len(set(skip_inp)):
        msgs.warn('There are duplicated files to skip.')
        skip_inp = list(set(skip_inp))
    if len(read_inp) > 0 and len(read_inp) != len(set(read_inp)):
        msgs.warn('There are duplicated files to read.')
        read_inp = list(set(read_inp))

    # Remove any files to skip
    for _skip in skip_inp:
        if _skip in read_inp:
            read_inp.remove(_skip)

    # Check that the files exist
    if file_check:
        for f in read_inp:
            if not os.path.isfile(f):
                raise FileNotFoundError('{0} does not exist!'.format(f))

    return read_inp


def _read_pypeit_file_lines(ifile):
    """
    General parser for a pypeit file.

    - Checks that the file exists.
    - Reads all the lines in the file
    - Removes comments, empty lines, and replaces special characters.
    
    Applies to settings, setup, and user-level reduction files.

    Args:
        ifile (str): Name of the file to parse.

    Returns:
        :obj:`numpy.ndarray`: Returns a list of the valid lines in the
        files.
    """
    # Check the files
    if not os.path.isfile(ifile):
        msgs.error('The filename does not exist -' + msgs.newline() + ifile)

    # Read the input lines and replace special characters
    with open(ifile, 'r') as f:
        lines = np.array([l.replace('\t', ' ').replace('\n', ' ').strip() \
                                for l in f.readlines()])
    # Remove empty or fully commented lines
    lines = lines[np.array([ len(l) > 0 and l[0] != '#' for l in lines ])]
    # Remove appended comments and return
    return np.array([ l.split('#')[0] for l in lines ])


def _determine_data_format(lines):
    """Determine how the data is formatted in the pypeit file."""
    return 'table' if len(lines) > 1 and lines[1][0] == '|' else 'raw'


def _read_data_file_table(lines, file_check=True):
    """Read the file table format."""
    space_ind = lines[0].index(" ")
    path = lines[0][space_ind+1:]
    header = np.array([ l.strip() for l in lines[1].split('|') ])

    file_col = np.where(header == 'filename')[0]
    if len(file_col) == 0:
        msgs.error('Table format failure: No \'filename\' column.')
    file_col = file_col[0]

    frame_col = np.where(header == 'frametype')[0]
    if len(frame_col) == 0:
        msgs.error('Table format failure: No \'frametype\' column.')
    frame_col = frame_col[0]

    nfiles = len(lines) - 2
    frametype = {}
    data_files = []
    for i in range(nfiles):
        row_values = lines[i+2].split('|') 
        filename = row_values[file_col].strip()
        frametype[filename] = row_values[frame_col].strip()
        filename = os.path.join(path, filename)
        if os.path.isfile(filename):
            data_files.append(filename)
        elif file_check:
            msgs.error('File does not exist: {0}'.format(filename))
    return data_files, frametype


def _parse_setup_lines(lines):
    """Return a list of the setup names"""
    return [ l.split()[1].strip() for l in lines if 'Setup' in l ]


def parse_pypeit_file(ifile, file_check=True):
    """
    Parse the user-provided .pypeit reduction file.

    Args:
        ifile (:obj:`str`):
            Name of pypeit file
        file_check (:obj:`bool`, optional):
            Check that the files in the pypeit configuration data file
            exist, and fault if they do not.

    Returns:
        :obj:`lists`: Four lists are provided: (1) the list of
        configuration lines, (2) the list of datafiles to read, (3) the
        list of frametypes for each file, and (4) the list of setup
        lines.
    """
    # Read in the pypeit reduction file
    msgs.info('Loading the reduction file')
    lines = _read_pypeit_file_lines(ifile)

    # Used to select the configuration lines: Anything that isn't part
    # of the data or setup blocks is assumed to be part of the
    # configuration
    is_config = np.ones(len(lines), dtype=bool)

    # Parse data block
    s, e = _find_pypeit_block(lines, 'data')
    if s >= 0 and e < 0:
        msgs.error("Missing 'data end' in {0}".format(ifile))
    if s < 0:
        msgs.error("You haven't specified any data!")
    data_format = _determine_data_format(lines[s:e])
    if data_format == 'raw':
        frametype = None
        data_files = _read_data_file_names(lines[s:e], file_check=file_check)
    elif data_format == 'table':
        data_files, frametype = _read_data_file_table(lines[s:e], file_check=file_check)
    is_config[s-1:e+1] = False
    if len(data_files) == 0 and file_check:
        msgs.error('There are no raw data frames' + msgs.newline() +
                   'Perhaps the path to the data is incorrect?')
    else:
        msgs.info('Found {0:d} raw data frames'.format(len(data_files)))

    # Parse the setup block
    s, e = _find_pypeit_block(lines, 'setup')
    if s >= 0 and e < 0:
        msgs.error("Missing 'setup end' in {0}".format(ifile))
    if s < 0:
        setups = []
    else:
        setups = _parse_setup_lines(lines[s:e])
        is_config[s-1:e+1] = False

    msgs.info('Input file loaded successfully')
    return list(lines[is_config]), data_files, frametype, setups


def pypeit_config_lines(ifile):
    """
    Return the config lines from a pypeit file.
    """
    lines = _read_pypeit_file_lines(ifile)

    # Find the config lines, assumed to be everything *except* the lines
    # in the data and setup blocks
    is_config = np.ones(len(lines), dtype=bool)

    s, e = _find_pypeit_block(lines, 'data')
    if s >= 0 and e < 0:
        msgs.error("Missing 'data end' in {0}".format(ifile))
    if not s < 0:
        is_config[s-1:e+1] = False
    
    s, e = _find_pypeit_block(lines, 'setup')
    if s >= 0 and e < 0:
        msgs.error("Missing 'setup end' in {0}".format(ifile))
    if not s < 0:
        is_config[s-1:e+1] = False

    return list(lines[is_config])
    

#def pypeit_file_line_groups(ifile):
#    """
#    Return the config, setup, and data lines from a pypeit file.
#    """
#    lines = _read_pypeit_file_lines(ifile)
#
#    # Find the config lines, assumed to be everything *except* the lines
#    # in the data and setup blocks
#    is_config = numpy.ones(len(lines), dtype=bool)
#
#    s, e = find_setting_block(lines, 'data')
#    if s >= 0 and e < 0:
#        msgs.error("Missing 'data end' in {0}".format(ifile))
#    if not s < 0:
#        is_config[s-1:e+1] = False
#    data_lines = lines[s:e]
#    
#    s, e = find_setting_block(lines, 'setup')
#    if s >= 0 and e < 0:
#        msgs.error("Missing 'setup end' in {0}".format(ifile))
#    if not s < 0:
#        is_config[s-1:e+1] = False
#    setup_lines = lines[s:e]
#
#    return list(lines[is_config]), list(setup_lines), list(data_lines)
    

def make_pypeit_file(pypeit_file, spectrograph, data_files, cfg_lines=None, setup_mode=False,
                    setup_lines=None, sorted_files=None, paths=None):
    """ Generate a default PYPIT file

    Parameters
    ----------
    pyp_file : str
      Name of PYPIT file to be generated
    spectrograph : str
    dfnames : list
      Path + file root of datafiles
      Includes skip files
    parlines : list, optional
      Standard parameter calls
    spclines : list, optional
      Lines related to filetype and calibrations
    setup_script : bool, optional
      Running setup script?
    calcheck : bool, optional
      Run calcheck?

    Returns
    -------
    Creates a PYPIT File

    """
    # Error checking
    if not isinstance(data_files, list):
        raise IOError("files_root needs to be a list")

    # Defaults
    if cfg_lines is None:
        _cfg_lines = ['[rdx]']
        _cfg_lines += ['    spectrograph = {0}'.format(spectrograph)]
    else:
        _cfg_lines = list(cfg_lines)
    if setup_mode:
        _cfg_lines += ['[calibrations]']
        _cfg_lines += ['    [[biasframe]]']
        _cfg_lines += ['        number = 0']
        _cfg_lines += ['    [[pixelflatframe]]']
        _cfg_lines += ['        number = 0']
        _cfg_lines += ['    [[arcframe]]']
        _cfg_lines += ['        number = 1']
        _cfg_lines += ['    [[pinholeframe]]']
        _cfg_lines += ['        number = 0']
        _cfg_lines += ['    [[traceframe]]']
        _cfg_lines += ['        number = 0']
        _cfg_lines += ['    [[standardframe]]']
        _cfg_lines += ['        number = 0']
    else:
        _cfg_lines += ['[calibrations]']
        _cfg_lines += ['    [[arcframe]]']
        _cfg_lines += ['        number = 1']

    # TODO: Clean up and check validity of _cfg_lines by reading it into
    # a ConfigObj?

    # Here we go
    with open(pypeit_file, 'w') as f:
        f.write("# This is a comment line\n")
        f.write("\n")
        f.write("# User-defined execution parameters\n")
        f.write('\n'.join(_cfg_lines))
        f.write('\n')
        f.write('\n')
        if setup_lines is not None:
            f.write("# Setup\n")
            f.write("setup read\n")
            for sline in setup_lines:
                f.write(' '+sline)
            f.write("setup end\n")
            f.write("\n")
        # Data
        f.write("# Read in the data\n")
        f.write("data read\n")
        # Old school
        for datafile in data_files:
            f.write(' '+datafile+'\n')
        # paths and Setupfiles
        if paths is not None:
            for path in paths:
                f.write(' path '+path+'\n')
        if sorted_files is not None:
            for sfile in sorted_files:
                f.write(sfile)
        f.write("data end\n")
        f.write("\n")


