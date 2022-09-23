"""
Module containing unused routines from ``pypeit.core.wavecal.waveio``


Created: 05-19-2022, TEB
"""

import glob
import os
import datetime

import numpy as np

from astropy.table import Table, Column, vstack

from pypeit import msgs
from pypeit.core.wavecal import defs, waveio
from pypeit import data

from IPython import embed



def load_by_hand():
    """
    By-hand line list

    Parameters
    ----------
    line_file
    add_path

    Returns
    -------
    byhand : Table

    """
    str_len_dict = defs.str_len()

    src_file = os.path.join(data.Paths.arclines, 'sources', 'by_hand_list.ascii')
    # Read
    line_list = Table.read(src_file, format='ascii.fixed_width', comment='#')
    # Add
    line_list['NIST'] = 1
    # Deal with Instr and Source
    ilist, slist = [], []
    for row in line_list:
        ilist.append(defs.instruments()[row['sInstr']])  # May need to split
        slist.append(row['sSource'])
    line_list['Instr'] = ilist
    line_list['Source'] = np.array(slist, dtype='S{:d}'.format(str_len_dict['Source']))
    # Trim
    return line_list[['ion', 'wave', 'NIST', 'Instr', 'amplitude', 'Source']]


def load_nist(ion):
    """
    Parse a NIST ASCII table.  Note that the long ---- should have been
    commented out and also the few lines at the start.

    Parameters
    ----------
    ion : str
        Name of ion

    Returns
    -------
    tbl : Table
        Table of lines

    """
    import glob
    # Find file
    srch_file = os.path.join(data.Paths.nist, f'{ion}_vacuum.ascii')
    nist_file = glob.glob(srch_file)
    if len(nist_file) == 0:
        raise IOError(f"Cannot find NIST file {srch_file}")
    # Read
    nist_tbl = Table.read(nist_file[0], format='ascii.fixed_width')
    gdrow = nist_tbl['Observed'] > 0.  # Eliminate dummy lines
    nist_tbl = nist_tbl[gdrow]
    # Now unique values only (no duplicates)
    uniq, indices = np.unique(nist_tbl['Observed'],return_index=True)
    nist_tbl = nist_tbl[indices]
    # Deal with Rel
    agdrel = []
    for row in nist_tbl:
        try:
            gdrel = int(row['Rel.'])
        except:
            try:
                gdrel = int(row['Rel.'][:-1])
            except:
                gdrel = 0
        agdrel.append(gdrel)
    agdrel = np.array(agdrel)
    # Remove and add
    nist_tbl.remove_column('Rel.')
    nist_tbl.remove_column('Ritz')
    nist_tbl['RelInt'] = agdrel
    #nist_tbl.add_column(Column([ion]*len(nist_tbl), name='Ion', dtype='S5'))
    nist_tbl.add_column(Column([ion]*len(nist_tbl), name='Ion', dtype='U5'))
    nist_tbl.rename_column('Observed','wave')
    # Return
    return nist_tbl


def load_source_table():
    """
    Load table of arcline sources

    Returns
    -------
    sources : Table

    """
    src_file = os.path.join(data.Paths.arclines, 'sources', 'arcline_sources.ascii')
    # Load
    sources = Table.read(src_file, format='ascii.fixed_width', comment='#')
    # Return
    return sources



#=============================================================================#
# These two routines include reference to the NIST files; the versions in
#  `pypeit.core.wavecal.waveio` no longer contain the NIST sections.
#  TEB: 05-19-2022

def load_line_list(line_file, add_path=False, use_ion=False, NIST=False):
    """

    Parameters
    ----------
    line_file : str
        Full path to line_list or name of ion
    add_path : bool, optional
        Not yet implemented
    NIST : bool, optional
        NIST formatted table?
    use_ion : bool, optional
        Interpret line_file as an ion, e.g. CuI

    Returns
    -------
    line_list : Table

    """
    if NIST:
        if use_ion:
            line_file = os.path.join(data.Paths.nist, f'{line_file}_vacuum.ascii')
        line_list = Table.read(line_file, format='ascii.fixed_width', comment='#')
        # Remove unwanted columns
        tkeys = line_list.keys()
        for badkey in ['Ritz','Acc.','Type','Ei','Lower','Upper','TP','Line']:
            for tkey in tkeys:
                if badkey in tkey:
                    line_list.remove_column(tkey)
        # Relative intensity -- Strip junk off the end
        reli = []
        for imsk, idat in zip(line_list['Rel.'].mask, line_list['Rel.'].data):
            if imsk:
                reli.append(0.)
            else:
                try:
                    reli.append(float(idat))
                except ValueError:
                    try:
                        reli.append(float(idat[:-1]))
                    except ValueError:
                        reli.append(0.)
        line_list.remove_column('Rel.')
        line_list['RelInt'] = reli
        #
        gdrows = line_list['Observed'] > 0.  # Eliminate dummy lines
        line_list = line_list[gdrows]
        line_list.rename_column('Observed','wave')
        # Others
        # Grab ion name
        i0 = line_file.rfind('/')
        i1 = line_file.rfind('_')
        ion = line_file[i0+1:i1]
        line_list.add_column(Column([ion]*len(line_list), name='Ion', dtype='U5'))
        line_list.add_column(Column([1]*len(line_list), name='NIST'))

    else:
        line_file = data.get_linelist_filepath(f'{line_file}_lines.dat') if use_ion else \
            data.get_linelist_filepath(line_file)
        line_list = Table.read(line_file, format='ascii.fixed_width', comment='#')

    # Return
    return line_list


def load_line_lists(lamps, unknown=False, skip=False, all=False, NIST=False,
                    restrict_on_instr=None):
    """
    Loads a series of line list files

    Parameters
    ----------
    lamps : list
        List of arc lamps to be used for wavelength calibration.
        E.g., ['ArI','NeI','KrI','XeI']
    unknown : bool, optional
    skip : bool, optional
        Skip missing line lists (mainly for building)
    NIST : bool, optional
        Load the full NIST linelists
    restrict_on_instr : str, optional
        Restrict according to the input spectrograph

    Returns
    -------
    line_list : Table

    """
    # All?
    if all:
        ### Also search the cache for linelists
        line_files = glob.glob(os.path.join(data.Paths.linelist, '*_lines.dat'))
        lamps = []
        for line_file in line_files:
            i0 = line_file.rfind('/')
            i1 = line_file.rfind('_')
            lamps.append(line_file[i0+1:i1])

    msgs.info(f"Arc lamps used: {', '.join(lamps)}")
    # Read standard files
    lists = []
    for lamp in lamps:
        if NIST:
            line_file = os.path.join(data.Paths.nist, f'{lamp}_vacuum.ascii')
        else:
            line_file = os.path.join(data.Paths.linelist, f'{lamp}_lines.dat')
        if not os.path.isfile(line_file):
            if not skip:
                line_files = glob.glob(os.path.join(data.Paths.linelist, '*_lines.dat'))
                all_list = [os.path.split(ll)[1].replace("_lines.dat", "") for ll in line_files]
                msgs.warn("Input line {:s} is not included in arclines".format(lamp))
                msgs.info("Please choose from the following list:" + msgs.newline() +
                          ",".join(all_list))
                import pdb; pdb.set_trace()
                raise IOError("Cannot continue without list")
        else:
            lists.append(load_line_list(line_file, NIST=NIST))
    # Stack
    if len(lists) == 0:
        return None
    line_lists = vstack(lists, join_type='exact')

    # Restrict on the spectrograph?
    if restrict_on_instr is not None:
        instr_dict = defs.instruments()
        gdI = (line_lists['Instr'] & instr_dict[restrict_on_instr]) > 0
        line_lists = line_lists[gdI]
        
    # Unknown
    if unknown:
        unkn_lines = waveio.load_unknown_list(lamps)
        unkn_lines.remove_column('line_flag')  # may wish to have this info
        # Stack
        line_lists = vstack([line_lists, unkn_lines])

    # Return
    return line_lists


#def load_spectrum(spec_file, index=0):
#    """
#    Load a simple spectrum from input file
#
#    Parameters
#    ----------
#    spec_file : str
#        Possible formats are:
#
#            - .fits --  Assumes simple ndarray in 0 extension
#            - .ascii -- Assumes Table.read(format='ascii') will work with single column
#
#    Returns
#    -------
#
#    """
#    import h5py
#    iext = spec_file.rfind('.')
#    if 'ascii' in spec_file[iext:]:
#        tbl = Table.read(spec_file, format='ascii')
#        key = tbl.keys()[0]
#        spec = tbl[key].data
#    elif 'fits' in spec_file[iext:]:
#        spec = fits.open(spec_file)[0].data
#    elif 'hdf5' in spec_file[iext:]:
#        hdf = h5py.File(spec_file, 'r')
#        if 'arcs' in hdf.keys():
#            print("Taking arc={:d} in this file".format(index))
#            spec = hdf['arcs/'+str(index)+'/spec'].value
#        else:
#            raise IOError("Not ready for this hdf5 file")
#    elif 'json' in spec_file[iext:]:
#        jdict = linetools.utils.loadjson(spec_file)
#        try:
#            spec = np.array(jdict['spec'])
#        except KeyError:
#            raise IOError("spec not in your JSON dict")
#    # Return
#    return spec


#=============================================================================#
# This function exists in `pypeit.data.arc_lines.convert_NIST_to_lists`
def write_line_list(tbl, outfile, overwrite=True):
    """
    Parameters
    ----------
    tbl
    outfile
    overwrite (optional), default=True
    """
    # Format
    tbl['wave'].format = '10.4f'
    # Write
    with open(outfile,'w') as f:
        f.write('# Creation Date: {:s}\n'.format(str(datetime.date.today().strftime('%Y-%m-%d'))))
        tbl.write(f, format='ascii.fixed_width', overwrite=overwrite)
