""" Module for I/O in arclines
"""
import pathlib

import astropy.table
import linetools.utils
import numpy as np

from pypeit import msgs
from pypeit.core.wavecal import defs
from pypeit import data

from IPython import embed


# TODO -- Move this to the WaveCalib object
def load_wavelength_calibration(filename: pathlib.Path) -> dict:
    """
    Load the wavelength calibration data from a file.

    Args:
        filename (:obj:`pathlib.Path`):
            Name of the json file.

    Returns:
        :obj:`dict`: Returns the wavelength calibration dictionary.
        Lists read from the json file are returnes as numpy arrays.
    """
    if not filename.is_file():
        msgs.error(f"File does not exist: {filename}")

    # Force any possible pathlib.Path object to string before `loadjson`
    wv_calib = linetools.utils.loadjson(str(filename))

    # Recast a few items as arrays
    for key in wv_calib.keys():
        if key in ['steps', 'par']:  # This isn't really necessary
            continue
        # Masked slit?
        if wv_calib[key] is None:
            continue
        # Arrays
        for tkey in wv_calib[key].keys():
            if isinstance(wv_calib[key][tkey], list):
                wv_calib[key][tkey] = np.array(wv_calib[key][tkey])

    return wv_calib


def load_template(arxiv_file, det, wvrng=None):
    """
    Load a full template file from disk

    Args:
        arxiv_file: str
        det: int
        wvrng (list, optional):
            min, max wavelength range for the arxiv


    Returns:
        wave: ndarray
        flux: ndarray
        binning: int, Of the template arc spectrum

    """
    # Path already included?
    if pathlib.Path(arxiv_file).name == arxiv_file:
        calibfile, _ = data.get_reid_arxiv_filepath(arxiv_file)
    else:
        calibfile = pathlib.Path(arxiv_file)
    # Read me
    tbl = astropy.table.Table.read(calibfile, format='fits')
    # Parse on detector?
    if 'det' in tbl.keys():
        idx = np.where(tbl['det'].data & 2**det)[0]
    else:
        idx = np.arange(len(tbl)).astype(int)
    tbl_wv = tbl['wave'].data[idx]
    tbl_fx = tbl['flux'].data[idx]

    # Cut down?
    if wvrng is not None:
        gd_wv = (tbl_wv >= wvrng[0]) & (tbl_wv <= wvrng[1])
        tbl_wv = tbl_wv[gd_wv]
        tbl_fx = tbl_fx[gd_wv]

    # Return
    return tbl_wv, tbl_fx, tbl.meta['BINSPEC']


def load_reid_arxiv(arxiv_file):
    """
    Load a REID arxiv file
    Now there are 2 possible formats.  We need to consolidate

    Args:
        arxiv_file (str):

    Returns:
        dict, dict-like:

    """
    # This function allows users to specify their own `reid_arxiv`, in
    #   particular, the output from `pypeit_identify`.
    calibfile, arxiv_fmt = data.get_reid_arxiv_filepath(arxiv_file)

    # This is a hack as it will fail if we change the data model yet again for wavelength solutions
    if arxiv_fmt == 'json':
        wv_calib_arxiv = load_wavelength_calibration(calibfile)
        par = wv_calib_arxiv['par'].copy()
        # Pop out par and steps if they were inserted in this calibration dictionary
        try:
            wv_calib_arxiv.pop('steps')
        except KeyError:
            pass
        try:
            wv_calib_arxiv.pop('par')
        except KeyError:
            pass
    elif arxiv_fmt == 'fits':
        # The following is a bit of a hack too
        par = None
        wv_tbl = astropy.table.Table.read(calibfile, format='fits')
        wv_calib_arxiv = {}
        nrow = wv_tbl['wave'].shape[0]
        for irow in np.arange(nrow):
            wv_calib_arxiv[str(irow)] = {}
            wv_calib_arxiv[str(irow)]['spec'] = wv_tbl['flux'][irow,:]
            wv_calib_arxiv[str(irow)]['wave_soln'] = wv_tbl['wave'][irow,:]
            wv_calib_arxiv[str(irow)]['order'] = wv_tbl['order'][irow]

    else:
        msgs.error(f"Not ready for this `reid_arxiv` extension: {arxiv_fmt}")

    return wv_calib_arxiv, par


def load_line_list(line_file, use_ion=False):
    """

    Parameters
    ----------
    line_file : str
        Full path to line_list or name of ion
    use_ion : bool, optional
        Interpret line_file as an ion, e.g. CuI

    Returns
    -------
    line_list : Table

    """
    line_file = data.get_linelist_filepath(f'{line_file}_lines.dat') if use_ion else \
        data.get_linelist_filepath(line_file)
    return astropy.table.Table.read(line_file, format='ascii.fixed_width', comment='#')


def load_line_lists(lamps, unknown=False, all=False, restrict_on_instr=None):
    """
    Loads a series of line list files

    Parameters
    ----------
    lamps : list
        List of arc lamps to be used for wavelength calibration.
        E.g., ['ArI','NeI','KrI','XeI']
    unknown : bool, optional
        Load the unknown list
    restrict_on_instr : str, optional
        Restrict according to the input spectrograph

    Returns
    -------
    line_list : Table

    """
    # All?
    if all:
        # Search both in the package directory and the PypeIt cache
        line_files = list(data.Paths.linelist.glob('*_lines.dat'))
        line_files.append(data.search_cache('_lines.dat'))
        lamps = []
        for line_file in line_files:
            i0 = line_file.rfind('/')
            i1 = line_file.rfind('_')
            lamps.append(line_file[i0+1:i1])

    msgs.info(f"Arc lamps used: {', '.join(lamps)}")
    # Read standard files
    # NOTE: If one of the `lamps` does not exist, data.get_linelist_filepath()
    #       will exit with msgs.error().
    lists = [load_line_list(data.get_linelist_filepath(f'{lamp}_lines.dat')) for lamp in lamps]
    # Stack
    if len(lists) == 0:
        return None
    line_lists = astropy.table.vstack(lists, join_type='exact')

    # Restrict on the spectrograph?
    if restrict_on_instr is not None:
        instr_dict = defs.instruments()
        gdI = (line_lists['Instr'] & instr_dict[restrict_on_instr]) > 0
        line_lists = line_lists[gdI]

    # Unknown
    if unknown:
        unkn_lines = load_unknown_list(lamps)
        unkn_lines.remove_column('line_flag')  # may wish to have this info
        # Stack
        line_lists = astropy.table.vstack([line_lists, unkn_lines])

    # Return
    return line_lists


def load_tree(polygon=4, numsearch=20):
    """
    Load a KDTree of ThAr patterns that is stored on disk

    Parameters
    ----------
    polygon : int
        Number of sides to the polygon used in pattern matching:

            - polygon=3  -->  trigon (two anchor lines and one floating line)
            - polygon=4  -->  tetragon (two anchor lines and two floating lines)
            - polygon=5  -->  pentagon (two anchor lines and three floating lines)

    numsearch : int
        Number of consecutive detected lines used to generate a pattern.
        For example, if numsearch is 4, then for a trigon, the following
        patterns will be generated (assuming line #1 is the left
        anchor):

            - 1 2 3  (in this case line #3 is the right anchor)
            - 1 2 4  (in this case line #4 is the right anchor)
            - 1 3 4  (in this case line #4 is the right anchor)

    Returns
    -------
    file_load : KDTree instance
        The KDTree containing the patterns
    index : ndarray
        For each pattern in the KDTree, this array stores the
        corresponding index in the linelist
    """

    # TODO: Can we save these as fits files instead?
    # TODO: Please don't use imports within functions
    import pickle
    filename = data.get_linelist_filepath(f'ThAr_patterns_poly{polygon}_search{numsearch}.kdtree')
    fileindx = data.get_linelist_filepath(
        f'ThAr_patterns_poly{polygon}_search{numsearch}.index.npy'
    )
    try:
        with open(filename, "rb", encoding="utf-8") as f_obj:
            file_load = pickle.load(f_obj)
        index = np.load(fileindx)
    except FileNotFoundError:
        msgs.info('The requested KDTree was not found on disk' + msgs.newline() +
                  'please be patient while the ThAr KDTree is built and saved to disk.')
        from pypeit.core.wavecal import kdtree_generator
        file_load, index = kdtree_generator.main(polygon, numsearch=numsearch, verbose=True,
                                                 ret_treeindx=True, outname=filename)

    return file_load, index


def load_unknown_list(lines, unknwn_file=None, all=False):
    """

    Parameters
    ----------
    lines : list
        Restricted lines;  use all=True for all
    unknwn_file : str, optional
    all : bool, optional


    Returns
    -------
    unknwn_lines : Table

    """
    line_dict = defs.lines()
    # Load
    if unknwn_file is None:
        unknwn_file = data.get_linelist_filepath('UNKNWNs.dat')
    line_list = load_line_list(unknwn_file)
    # Cut on input lamps?
    if all:
        return line_list

    # Otherwise
    msk = np.zeros(len(line_list), dtype=bool)
    for line in lines:
        line_flag = line_dict[line]
        match = line_list['line_flag'] % (2*line_flag) >= line_flag
        msk[match] = True
    # Finish
    return line_list[msk]
