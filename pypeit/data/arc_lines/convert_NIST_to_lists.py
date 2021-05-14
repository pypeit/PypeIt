"""
Module for generating Arc Line lists
  Should be run where it is located (for now)
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import os
import pdb
import datetime

from pkg_resources import resource_filename
from collections import OrderedDict
from astropy.table import Table, Column

line_path = resource_filename('pypeit', '/data/arc_lines/lists/')
nist_path = resource_filename('pypeit', '/data/arc_lines/NIST/')


def parser(options=None):
    import argparse
    # Parse
    parsefunc = argparse.ArgumentParser(
        description='Build the PypeIt line lists from NIST tables')
    parsefunc.add_argument("-w", "--write", default=False, action='store_true', help="Actually write files?")
    parsefunc.add_argument("--skip_stop", default=False, action='store_true', help="Skip warning stop?")
    parsefunc.add_argument("-r", "--relint", type=float, default=1000.0, help="Set the relative intensity threshold")
    parsefunc.add_argument("line", default='', help="Name of ion")

    if options is None:
        args = parsefunc.parse_args()
    else:
        args = parsefunc.parse_args(options)
    return args


def init_line_list():
    """ Initialize a Table for a linelist
    Rigidly enforces table column formats
    Strings are the most annoying

    Returns
    -------
    init_tbl : Table
      One dummy row
    """
    dummy_src = str('#')*50
    # Arc Line name
    dummy_line = str('#')*8
    #

    # Dict for Table
    idict = OrderedDict()
    idict['ion'] = dummy_line
    idict['wave'] = 0.
    idict['NIST'] = 0
    idict['Instr'] = 0  # Flag for instrument
    idict['amplitude'] = 0
    idict['Source'] = dummy_src

    # Table
    tkeys = idict.keys()
    lst = [[idict[tkey]] for tkey in tkeys]
    init_tbl = Table(lst, names=tkeys)

    # Return
    return init_tbl


def load_line_list(line):
    """
    Parameters
    ----------
    line : str
      Name of ion

    Returns
    -------
    line_list : Table

    """
    line_file = nist_path + '{:s}_vacuum.ascii'.format(line)

    # Check the NIST lines file exists
    if not os.path.isfile(line_file):
        raise IOError("Input line {:s} is not available".format(line))

    line_list = Table.read(line_file, format='ascii.fixed_width', comment='#')

    # Remove unwanted columns
    tkeys = line_list.keys()
    for badkey in ['Ritz', 'Acc.', 'Type', 'Ei', 'Lower', 'Upper', 'TP', 'Line']:
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
    line_list['Rel.'] = reli
    #
    gdrows = line_list['Observed'] > 0.  # Eliminate dummy lines
    line_list = line_list[gdrows]
    line_list.rename_column('Observed', 'wave')
    # Others
    # Grab ion name
    i0 = line_file.rfind('/')
    i1 = line_file.rfind('_')
    ion = line_file[i0+1:i1]
    line_list.add_column(Column([ion]*len(line_list), name='Ion', dtype='U5'))
    line_list.add_column(Column([1]*len(line_list), name='NIST'))
    return line_list


def main(args=None):
    """ This script convert an input NIST table into a line list that can be used by PypeIt

    Parameters

    ----------
    args

    Returns
    -------

    """
    # Grab arguments
    pargs = parser(options=args)
    line = pargs.line
    relIntThreshold = pargs.relint

    print("=============================================================")
    print("This script is for EXPERTS ONLY")
    print("Continue only if you know what you are doing")
    print("Otherwise exit")
    print("p.s.  You need to remove the files you wish to re-build")
    print("=============================================================")
    if not pargs.skip_stop:
        pdb.set_trace()

    # Load the NIST ThAr list
    llist = load_line_list(line)
    # ['wave', 'Aki', 'Rel.', 'Ion', 'NIST']

    # Generate a table
    linelist = init_line_list()

    # now add all NIST lines
    nlines = llist['Ion'].size
    for ll in range(nlines):
        if llist['Rel.'][ll] > relIntThreshold:
            linelist.add_row([llist['Ion'][ll], llist['wave'][ll], 1, 0, llist['Rel.'][ll], 'NIST'])
        if ll+1 % 100 == 0:
            print(ll+1, '/', nlines)
    # Remove the first dummy row
    linelist.remove_row(0)

    # Finally, sort the list by increasing wavelength
    linelist.sort('wave')

    # Write?
    if not pargs.write:
        print("=============================================================")
        print("Rerun with --write if you are happy with what you see.")
        print("=============================================================")
        return

    # Write the table to disk
    outfile = line_path + '{:s}_lines.dat'.format(line)
    write_line_list(linelist, outfile)
    return


def write_line_list(tbl, outfile):
    """
    Parameters
    ----------
    tbl
    outfile
    """
    # Format
    tbl['wave'].format = '10.4f'
    # Write
    with open(outfile, 'w') as f:
        f.write('# Creation Date: {:s}\n'.format(str(datetime.date.today().strftime('%Y-%m-%d'))))
        tbl.write(f, format='ascii.fixed_width')


if __name__ == '__main__':
    main()
