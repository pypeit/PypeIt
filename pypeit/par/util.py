# -*- coding: utf-8 -*-
"""
Utility functions for PypeIt parameter sets
"""
import os
import time
import glob
import warnings
import textwrap
from IPython import embed

import numpy as np

from astropy.table import Table

from configobj import ConfigObj

from pypeit import msgs


#-----------------------------------------------------------------------
# Parameter utility functions
#-----------------------------------------------------------------------
def _eval_ignore():
    """Provides a list of strings that should not be evaluated."""
    return [ 'open', 'file', 'dict', 'list', 'tuple' ]


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

    Args:
        par (ParSet):

    Returns:
        dict: Converted ParSet

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
    """
    Find the PypeIt group block

    Args:
        lines (:obj:`list`):
            List of file lines
        group (:obj:`str`):
            Name of group to parse

    Returns:
        int, int: Starting,ending line of the block;  -1 if not present

    """
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
    """
    Expand the data file name as necessary and
    then search for all data files

    Args:
        inp (str): Path
        current_path (str or None):

    Returns:
        list: Glob list of files in the generated a path

    """
    out = os.path.expanduser(inp) if inp[0] == '~' else inp
    if current_path is not None:
        out = os.path.join(current_path, out)
    return glob.glob(out)
    

def _read_data_file_names(lines, file_check=True):
    """
    Read the raw data file format

    Args:
        lines (list):
        file_check (bool, optional):

    Returns:
        list: List of data file names

    """
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



def _determine_data_format(lines):
    """
    Determine the format of the data block in the .pypeit file.

    The test used in this function is pretty basic.  A table format is
    assumed if the first character in *any* line is `|`.

    Args:
        lines (:obj:`list`):
            The list of lines read from the data block of the pypeit
            file.
    
    Returns:
        str: The syntax of the data files to read::

            'raw': A (list of) file roots to be read or found using
            `glob`.

            'table': ASCII output of an astropy.table.Table

    """
    for l in lines:
        if l[0] == '|':
            return 'table'
    return 'raw'


def _read_data_file_table(lines, file_check=True):
    """
    Read the file table format.
    
    Args:
        lines (:obj:`list`):
            List of lines *within the data* block read from the pypeit
            file.
        file_check (:obj:`bool`, optional):
            Check if the specified data files exist.
    
    Returns:
        list, dict, Table:  Returns the list of data file names, a
        dictionary with the frame types of each file where the key of
        the dictionary is the file name, and a Table with the data
        provided in the pypeit file.  Note that the files listed in the
        first object contain the full path, whereas the file names in
        the frame type dictionary and the data table do not include the
        full path to the file.

    Raise:
        PypeItError:
            Raised if `file_check=True` and any of the specified files
            to not exist, or if the table does not have a 'filename' or
            'frametype' column.
    """

    # Allow for multiple paths
    paths = []
    for l in lines:
        space_ind = l.index(" ")
        if l[:space_ind].strip() != 'path':
            break
        paths += [ l[space_ind+1:] ]

    npaths = len(paths)
    header = [ l.strip() for l in lines[npaths].split('|') ][1:-1]

    # Minimum columns required
    if 'filename' not in header:
        msgs.error('Table format failure: No \'filename\' column.')
    if 'frametype' not in header:
        msgs.error('Table format failure: No \'frametype\' column.')

    # Build the table
    nfiles = len(lines) - npaths - 1
    tbl = np.empty((nfiles, len(header)), dtype=object)

    for i in range(nfiles):
        row = np.array([ l.strip() for l in lines[i+npaths+1].split('|') ])[1:-1]
        if len(row) != tbl.shape[1]:
            raise ValueError('Data and header lines have mismatched columns!')
        tbl[i,:] = row
    data = {}
    for i,key in enumerate(header):
        data[key] = tbl[:,i]
    tbl = Table(data)

    # Build full paths to file and set frame types
    frametype = {}
    data_files = []
    for i in range(nfiles):
        frametype[tbl['filename'][i]] = tbl['frametype'][i]
        for p in paths:
            filename = os.path.join(p, tbl['filename'][i])
            if os.path.isfile(filename):
                break
        data_files.append(filename)
        if not os.path.isfile(filename) and file_check:
            msgs.error('File does not exist: {0}'.format(filename))

    return data_files, frametype, tbl


def _parse_setup_lines(lines):
    """Return a list of the setup names"""
    setups = []
    for l in lines:
        if 'Setup' in l:
            tsetup = l.split()[1].strip()
            # Remove any lingering colon
            if tsetup[-1] == ':':
                setup = tsetup[:-1]
            else:
                setup = tsetup
            setups.append(setup)
    #
    return setups


def parse_pypeit_file(ifile, file_check=True, runtime=False):
    """
    Parse the user-provided .pypeit reduction file.

    Args:
        ifile (:obj:`str`):
            Name of pypeit file
        file_check (:obj:`bool`, optional):
            Check that the files in the pypeit configuration data file
            exist, and fault if they do not.
        runtime (:obj:`bool`, optional):
            Perform additional checks if called to run PypeIt

    Returns:
        5-element tuple containing

        - list:  List of configuration lines,
        - list:  List of datafiles to read,
        - list:  List of frametypes for each file
        - :obj:`astropy.table.Table`:  Table of user supplied info on data files
        - list:  List of setup lines.
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
        usrtbl = None
        data_files = _read_data_file_names(lines[s:e], file_check=file_check)
    elif data_format == 'table':
        data_files, frametype, usrtbl = _read_data_file_table(lines[s:e], file_check=file_check)
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

    # TODO: This should be moved to the PypeIt class
    # Running PypeIt?
    if runtime:
        for key in ['filename', 'frametype']:
            if key not in usrtbl.keys():
                msgs.error("Add {:s} to your PypeIt file before using run_pypeit".format(key))
        # Setup
        if len(setups) != 1:
            msgs.error("Add setup info to your PypeIt file in the setup block!")

    msgs.info('Input file loaded successfully')
    return list(lines[is_config]), data_files, frametype, usrtbl, setups


def pypeit_config_lines(ifile):
    """
    Return the config lines from a PypeIt file.

    Args:
        ifile (str): Name of PypeIt file

    Returns:
        list: List of configuration lines; will be used for ConfigObj

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
    

def make_pypeit_file(pypeit_file, spectrograph, data_files, cfg_lines=None, setup_mode=False,
                     setup_lines=None, sorted_files=None, paths=None):
    """
    Generate a default PypeIt file

    Args:
        pypeit_file (str): Name of PYPIT file to be generated
        spectrograph (str):  Name of spectrograph
        data_files (list):  List of data files -- essentially Deprecated
        cfg_lines (list, optional):  List of configuration lines for parameters
        setup_mode (bool, optional):  If True, 0 out required files for everything except Arc
        setup_lines (list, optional):
        sorted_files (list, optional):
        paths (list, optional): List of paths for slurping data files
    """
    # Error checking
    if not isinstance(data_files, list):
        raise IOError("data_files needs to be a list")

    # Defaults
    if cfg_lines is None:
        _cfg_lines = ['[rdx]']
        _cfg_lines += ['    spectrograph = {0}'.format(spectrograph)]
    else:
        _cfg_lines = list(cfg_lines)

    # TODO: Bring back checks for the appropriate number of calibration
    # frames?

    # TODO: Clean up and check validity of _cfg_lines by reading it into
    # a ConfigObj?

    # Here we go
    with open(pypeit_file, 'w') as f:
        f.write('# Auto-generated PypeIt file\n')
        f.write('# {0}\n'.format(time.strftime("%a %d %b %Y %H:%M:%S",time.localtime())))
        f.write("\n")
        f.write("# User-defined execution parameters\n")
        f.write('\n'.join(_cfg_lines))
        f.write('\n')
        f.write('\n')
        if setup_lines is not None:
            f.write("# Setup\n")
            f.write("setup read\n")
            f.write('\n'.join(setup_lines)+'\n')
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
            f.write('\n'.join(sorted_files))
            f.write('\n')
        f.write("data end\n")
        f.write("\n")

    msgs.info('PypeIt file written to: {0}'.format(pypeit_file))


