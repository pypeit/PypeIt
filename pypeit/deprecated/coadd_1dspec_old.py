#!/usr/bin/env python
"""
Wrapper to the linetools XSpecGUI
"""
import argparse



def parser(options=None):
    parser = argparse.ArgumentParser(description='Script to coadd a set of spec1D files and 1 or more slits and 1 or more objects. Current defaults use Optimal + Fluxed extraction. [v1.1]')
    parser.add_argument("infile", type=str, help="Input file (YAML)")
    parser.add_argument("--debug", default=False, action='store_true', help="Turn debugging on")
    parser.add_argument("--show", default=False, action='store_true', help="Show the coadded spectra")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args, unit_test=False, path=''):
    """ Runs the XSpecGui on an input file
    path : str, optional
      Mainly for running the unit test
    """

    import glob
    import yaml
    import IPython

    import numpy as np
    from numpy import isnan

    from astropy.io import fits

    from pypeit import msgs
    from pypeit.core import coadd
    from pypeit import specobjs
    from pypeit.spectrographs import util

    # Load the input file
    with open(args.infile, 'r') as infile:
        coadd_dict = yaml.load(infile)

    # Spectrograph
    spectrograph = util.load_spectrograph(coadd_dict.pop('spectrograph'))

    # Grab object names in the spectra
    filelist = coadd_dict.pop('filenames')
    # Allow for wildcards
    files = []
    for ifl in filelist:
        if '*' in ifl:
            files += glob.glob(path+ifl)
        else:
            files += [path+ifl]
    # Load spectra
    if len(files) == 0:
        msgs.error("No files match your input list")
    else:
        msgs.info("Coadding {:d} data frames".format(len(files)))
        # figure out whether it is Echelle or Longslit
        header0 = fits.getheader(files[0],0)
        pypeline = header0['PYPELINE']
        # also need norder for Echelle data
        if pypeline == 'Echelle':
            #ext_final = fits.getheader(files[0], -1)
            #norder = ext_final['ECHORDER'] + 1
            ext_first = fits.getheader(files[0], 1)
            ext_final = fits.getheader(files[0], -1)
            norder = abs(ext_final['ECHORDER'] - ext_first['ECHORDER']) + 1
            order_vec = np.arange(np.fmax(ext_first['ECHORDER'],ext_final['ECHORDER']),
                                  np.fmin(ext_first['ECHORDER'],ext_final['ECHORDER'])-1,-1)
    fdict = {}
    for ifile in files:
        # Open file
        hdulist = fits.open(ifile)
        # Grab objects
        objects = [hdu.name for hdu in hdulist][1:]
        fdict[ifile] = objects

    # Global parameters?
    if 'global' in coadd_dict.keys():
        gparam = coadd_dict.pop('global')
    else:
        gparam = {}
    if args.debug:
        gparam['debug'] = True
    if args.show:
        gparam['show'] = True
    sv_gparam = gparam.copy()
    # Extraction
    if 'extract' in coadd_dict.keys():
        ex_value = coadd_dict.pop('extract')
    else:
        ex_value = 'OPT'
    msgs.info("Using {:s} extraction".format(ex_value))
    # Fluxed data?
    if 'flux' in coadd_dict.keys():
        flux_value = coadd_dict.pop('flux')
    else:
        flux_value = True
    # sensfunc weight
    if 'sensfile' in coadd_dict.keys():
        sensfile = coadd_dict.pop('sensfile')
    else:
        sensfile = None

    # Loop on sources
    for key in coadd_dict.keys():
        # Re-init gparam
        gparam = sv_gparam.copy()
        iobj = coadd_dict[key]['object']
        # Check iobj input
        if isinstance(iobj, list):
            if len(iobj) != len(files):
                raise IOError("Input list of object names must have same length as files")
        #
        outfile = coadd_dict[key]['outfile']

        # Scale
        if 'scale' in coadd_dict[key]:
            scale_dict = coadd_dict[key]['scale']
        else:
            scale_dict = None

        # Generate local keywords
        try:
            local_kwargs = coadd_dict[key]['local']
        except KeyError:
            local_kwargs = {}
        else:
            for lkey in local_kwargs:
                gparam[lkey] = local_kwargs[lkey]

        if unit_test:
            return gparam, ex_value, flux_value, iobj, outfile, files, local_kwargs


        # Loop on spec1d files
        gdfiles = []
        extensions = []
        gdobj = []

        for fkey in fdict:
            # Input as str or list
            if not isinstance(iobj, list) == 1:  # Simple single object
                use_obj = iobj
            else:
                ind = files.index(fkey)
                use_obj = iobj[ind]

            if pypeline == 'Echelle':
                gdfiles.append(fkey)
                gdobj += [use_obj]
            else:
                # Find object indices
                # FW: mtch_obj_to_objects will return None when no matching and raise TypeError: cannot unpack non-iterable NoneType object
                try:
                    mtch_obj, idx = specobjs.mtch_obj_to_objects(use_obj, fdict[fkey], **local_kwargs)
                except TypeError:
                    mtch_obj = None
                if mtch_obj is None:
                    msgs.info("No object {:s} in file {:s}".format(iobj, fkey))
                elif len(mtch_obj) == 1:
                    #Check if optimal extraction is present in all objects.
                    # If not, warn the user and set ex_value to 'box'.
                    hdulist = fits.open(fkey)
                    try: #In case the optimal extraction array is a NaN array
                        if flux_value is True: # If we have a fluxed spectrum, look for flam
                            obj_opt = hdulist[mtch_obj[0]].data['OPT_FLAM']
                        else: # If not, look for counts
                            obj_opt = hdulist[mtch_obj[0]].data['OPT_COUNTS']
                        if any(isnan(obj_opt)):
                            msgs.warn("Object {:s} in file {:s} has a NaN array for optimal extraction. Boxcar will be used instead.".format(mtch_obj[0],fkey))
                            ex_value = 'box'
                    except KeyError: #In case the array is absent altogether.
                        msgs.warn("Object {:s} in file {:s} doesn't have an optimal extraction. Boxcar will be used instead.".format(mtch_obj[0],fkey))
                        try:
                            if flux_value is True: # If we have a fluxed spectrum, look for flam
                                hdulist[mtch_obj[0]].data['BOX_FLAM']
                            else: # If not, look for counts
                                hdulist[mtch_obj[0]].data['BOX_COUNTS']
                        except KeyError:
                            #In case the boxcar extract is also absent
                            msgs.error("Object {:s} in file {:s} doesn't have a boxcar extraction either. Co-addition cannot be performed".format(mtch_obj[0],fkey))
                        ex_value = 'box'
                    gdfiles.append(fkey)
                    gdobj += mtch_obj
                    extensions.append(idx[0]+1)
                else:
                    raise ValueError("Multiple matches to object {:s} in file {:s}".format(iobj, fkey))

        # Load spectra
        if len(gdfiles) == 0:
            msgs.error("No files match your input criteria")

        # QA file name
        exten = outfile.split('.')[-1]  # Allow for hdf or fits or whatever
        qafile = outfile.replace(exten, 'pdf')

        ## The following part will be replaced with the new coadd code
        if pypeline == 'Echelle':

            # Check whether the scale_dict is in the right shape.
            if 'orderscale' in gparam.keys():
                orderscale_value = gparam['orderscale']
            else:
                orderscale_value = 'median'

            if (scale_dict is not None) and (orderscale_value=='photometry'):
                if len(scale_dict) != norder:
                    raise IOError("You need to specifiy the photometric information for every order.")

            wave_stack, flux_stack, ivar_stack, mask_stack = coadd.ech_combspec(
                gdfiles, gdobj, sensfile=sensfile, ex_value=ex_value, flux_value=flux_value, phot_scale_dicts=scale_dict,
                outfile=outfile, qafile=qafile, **gparam)

        else:
            wave_stack, flux_stack, ivar_stack, mask_stack = coadd.multi_combspec(
                gdfiles, gdobj, ex_value=ex_value, flux_value=flux_value, phot_scale_dicts=scale_dict,
                outfile=outfile, qafile=qafile, **gparam)
