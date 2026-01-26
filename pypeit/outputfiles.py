"""
Module for constructing output file names and paths for PypeIt reductions.

.. include:: ../include/links.rst

"""
import numpy as np
from pathlib import Path

from pypeit import msgs

def get_std_outfile(fitstbl, par, standard_frames:list):
    """
    Return the spec1d file name for a reduced standard to use as a tracing
    crutch.

    The file is either constructed using the provided standard frame indices
    or it is directly pulled from the
    :class:`~pypeit.par.pypeitpar.FindObjPar` parameters in :attr:`par`.
    The latter takes precedence.  If more than one row is provided by
    ``standard_frames``, the first index is used.

    Args:
        fitstbl (:class:`pypeit.metadata.PypeItMetaData`):
            The metadata table for the current reduction.
        par (:class:`pypeit.par.pypeitpar.PypeItPar`):
            The parameter set for the current reduction.
        standard_frames (:obj:`list`):
            List of frame indices in ``fitstbl`` corresponding to
            standard star frames.

    Returns:
        :obj:`str`: Full path to the standard spec1d output file to use.
    """
    std_outfile = par['reduce']['findobj']['std_spec1d']
    if std_outfile is not None:
        if not par['reduce']['findobj']['use_std_trace']:
            msgs.error('If you provide a standard star spectrum for tracing, you must set use_std_trace=True')
        elif not Path(std_outfile).absolute().exists():
            msgs.error(f'Provided standard spec1d file does not exist: {std_outfile}')
        return std_outfile

    # TODO: Need to decide how to associate standards with
    # science frames in the case where there is more than one
    # standard associated with a given science frame.  Below, I
    # just use the first standard

    std_frame = None if (len(standard_frames) == 0 or not par['reduce']['findobj']['use_std_trace']) \
        else standard_frames[0]
    # Prepare to load up standard?
    if std_frame is not None:
        std_outfile = spec_output_file(fitstbl, par, std_frame) \
                        if isinstance(std_frame, (int,np.integer)) else None
    if std_outfile is not None and not std_outfile.is_file():
        msgs.error(f'Could not find standard file: {std_outfile}')
    return std_outfile

def intermediate_filename(itype:str, basename:str, det_name:str, 
                          inter_path:str='Intermediate'):
    """
    Construct the intermediate file name for a given type and detector

    Args:
        itype (:obj:`str`):
            Type of intermediate file
            Standard options are: sciImg, bkgImg, initSky, spec1d
        basename (:obj:`str`):
            Basename of the file
        det_name (:obj:`str`):
            Name of the detector
        inter_path (:obj:`str`, optional):
            Path to the intermediate files

    Returns:
        :obj:`str`: The full path to the intermediate file
    """
    return Path(inter_path) / f'{itype}_{basename}_{det_name}.fits'

def science_path(par) -> Path:
    """
    Constructs the path to the science directory based on the provided parameters.

    Args:
        par (:class:`~pypeit.par.pypeitpar.PypeItPar`):
            The parameter set for the reduction, including slitmask and
            object finding parameters.

    Returns:
        Path: The full path to the science directory as a `Path` object.
    """
    return Path(par['rdx']['redux_path']) / par['rdx']['scidir']

def spec_output_file(fitstbl, par, frame:int, twod:bool=False,
                     ext:str='.fits', sci_path:Path=None) -> Path:
    """
    Return the path to the spectral output data file.
    
    Args:
        fitstbl (:class:`pypeit.metadata.PypeItMetaData`):
            The metadata table for the current reduction.
        par (:class:`pypeit.par.pypeitpar.PypeItPar`):
            The parameter set for the current reduction.
        frame (:obj:`int`):
            Frame index from :attr:`fitstbl`.
        twod (:obj:`bool`), optional:
            Name for the 2D output file; 1D file otherwise.
        ext (:obj:`str`, optional):
            Extension for the output file.  Default is '.fits'.
    
    Returns:
        `Path`_: The path for the output file
    """
    if sci_path is None:
        sci_path = science_path(par) 
    # Bits and pieces
    basename = fitstbl.construct_basename(frame)
    # Finish``
    return sci_path / f'spec{"2" if twod else "1"}d_{basename}{ext}'