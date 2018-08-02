*********
Code Flow
*********

Overview
========

Setup
=====

Flow
----

Below is the code flow for the :ref:`pypeit-setup` script.  The
following are nearly all function names or object methods.
The module name is typically the first item, e.g. arparse.init
is a method in arparse.py.  Here goes::

   ├── pypeit_setup
   |  ├── pyputils.make_pypeit_file(pypeit_file, spectrograph, dfnames)
   |  ├── run_pypeit.parser(options)
   |  ├── pypeit.PypeIt(args)
   |  |  ├── load_input(pypeit_file)
   |  |  |  ++ generates pyp_dict
   |  |  ├── arparse.get_argflag_class()
   |  |  |  ++ generates argf Class
   |  |  ├── argf.init_param()
   |  |  ├── plines = argf.load_lines(parlines)
   |  |  ├── argf.set_paramlist(plines)
   |  |  ├── arparse.get_spect_class()
   |  |  |  ++ generates spect Class
   |  |  ├── spect.load_file(base=True)  # default
   |  |  ├── spect.set_paramlist(lines)
   |  |  ├── spect.load_file()  # instrument specific
   |  |  ├── spect.set_paramlist(lines)
   |  |  ├── spect.set_paramlist(plines)   # using pyp_dict
   |  |  ├── spect.load_lines(spclines)      # command line
   |  |  ├── argf.save()
   |  |  ├── spect.save()
   |  |  ├── arparse.init(argf, spect)  # Saved into settings.
   |  |  ├── # fitsdict created
   |  |  ├── fitsdict = arload.load_headers(datlines)
   |  |  |  ├── # Checks on FITS files
   |  |  |  ├── settings.spect['check'][ch]
   |  |  ├── # Flip dispersion direction (if needed)
   |  ├── ARMLSD  # This is formally below PypeIt, but I want to reduce the indentation here
   |  |  ├── armbase.SetupScience(fitsdict)
   |  |  |  ├── filesort = arsort.sort_data(fitsdict)
   |  |  |  |  ├── find_standard_file()
   |  |  |  |  ++ Generates filesort dict
   |  |  |  ├── arsort.match_science(fitsdict, filesort)
   |  |  |  |  ++ Written to settings.spect[ftag]['index']
   |  |  |  ├── arsciexp.ScienceExposure(i, fitsdict)
   |  |  |  |  ++ Generates sciexp list of ScienceExposure objects
   |  |  |  ++ setup_dict generated
   |  |  |  ├── arsort.instr_setup(sciexp, kk+1, fitsdict, setup_dict)
   |  |  |  |  ++ Generates setupIDs
   |  |  |  ++ group_dict generated
   |  |  |  ├── arsort.write_sorted(group_dict, setup_dict)
   |  |  |  ├── arsort.write_setup(setup_dict)


Items
-----

Items created and carried around::

    filesort
    fitsdict
    settings.spect
    settings.argf
    setup_dict
    group_dict
    sciexp

