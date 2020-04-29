#!/usr/bin/env python
import argparse
import os
import numpy as np
from astropy.io import fits
from astropy.table import Table
from pypeit import msgs


def parser(options=None):
    parser = argparse.ArgumentParser(description='Parse')
    parser.add_argument("sci_path", type=str, help="Path for Science folder")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args):
    """
      Prepare fluxing, coadd1d files.
    """
    allfiles = os.listdir(args.sci_path)
    allfiles = np.sort(allfiles)
    spec1dfiles = []
    spec2dfiles = []
    spec1dinfos = []
    for ifile in allfiles:
        if ('spec1d' in ifile) and ('.fits' in ifile):
            spec1dfiles.append(ifile)
        elif ('spec2d' in ifile) and ('.fits' in ifile):
            spec2dfiles.append(ifile)
        elif ('spec1d' in ifile) and ('.txt' in ifile):
            spec1dinfos.append(ifile)
        else:
            msgs.warn('{:} is not a standard PypeIt output.'.format(ifile))
    if len(spec2dfiles) > len(spec1dfiles):
        msgs.warn('The following exposures do not have 1D extractions:')
        for ii in range(len(spec2dfiles)):
            if not os.path.exists(os.path.join(args.sci_path, spec2dfiles[ii].replace('spec2d','spec1d'))):
                msgs.info('\t {:}'.format(spec2dfiles[ii]))

    if len(spec1dfiles) > 0:
        par = fits.open(os.path.join(args.sci_path, spec1dfiles[0]))

        ## fluxing pypeit file
        f1 = open('{:}.flux'.format(par[0].header['PYP_SPEC']), 'w')
        print('# User-defined fluxing parameters', file=f1)
        print('[fluxcalib]', file=f1)
        print('  extinct_correct = False # Set to True if your SENSFUNC derived with the UVIS algorithm', file=f1)
        print('\nflux read', file=f1)

        ## coadd1d pypeit file
        f2 = open('{:}.coadd1d'.format(par[0].header['PYP_SPEC']), 'w')
        print('# User-defined coadd1d parameters', file=f2)
        print('[coadd1d]', file=f2)
        print('  coaddfile = YOUR_OUTPUT_FILE_NAME', file=f2)
        print('  sensfuncfile = YOUR_SENSFUNC_FILE', file=f2)
        print('\n# This file includes all extracted objects. You need to figure out which object you want to \n' + \
              '# coadd before running pypeit_coadd_1dspec!!!\n', file=f2)
        print('coadd1d read', file=f2)

        ## tellfit pypeit file
        f3 = open('{:}.tell'.format(par[0].header['PYP_SPEC']), 'w')
        print('# User-defined tellfit parameters. Please only use one of the following object models.', file=f3)
        print('\n# Algorithm for Quasars', file=f3)
        print('[tellfit]', file=f3)
        print('  objmodel = qso', file=f3)
        print('  redshift = 7.0', file=f3)
        print('\n# Algorithm for Stars', file=f3)
        print('[tellfit]', file=f3)
        print('  objmodel = star', file=f3)
        print('  star_type = A0', file=f3)
        print('  star_mag = 8.0', file=f3)
        print('\n# Algorithm for Other objects', file=f3)
        print('[tellfit]', file=f3)
        print('  objmodel = poly', file=f3)
        print('  polyorder = 3', file=f3)
        print('  fit_region_mask = 9000.0,9500.', file=f3)
        f3.close()

        for ii in range(len(spec1dfiles)):
            if ii == 0:
                print('  ' + os.path.join(args.sci_path, spec1dfiles[ii]) + ' YOUR_SENSFUNC_FILE', file=f1)
            else:
                print('  ' + os.path.join(args.sci_path, spec1dfiles[ii]), file=f1)

            meta_tbl = Table.read(os.path.join(args.sci_path, spec1dfiles[ii]).replace('.fits', '.txt'),
                                  format='ascii.fixed_width')
            objects = np.unique(meta_tbl['name'])
            for jj in range(len(objects)):
                print('  ' + os.path.join(args.sci_path, spec1dfiles[ii]) + ' ' + meta_tbl['name'][jj], file=f2)

        print('flux end', file=f1)
        f1.close()

        print('coadd1d end', file=f2)
        f2.close()


