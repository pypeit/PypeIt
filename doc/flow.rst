*********
Code Flow
*********

Overview
========

Setup
=====

Flow
----

Here is the code flow for the :ref:`pypit-setup` script::

   ├── pypit_setup
   |  ├── pyputils.make_pypit_file()
   |  ├── run_pypit.parser()
   |  ├── pypit.PYPIT()
   |  |  ├── load_input()
   |  |  |  ++ generates pyp_dict
   |  |  ├── arparse.get_argflag_class()
   |  |  |  ++ generates argf Class
   |  |  ├── argf.init_param()
   |  |  ├── argf.set_param()
   |  |  ├── plines = argf.load_lines()
   |  |  ├── argf.set_paramlist(plines)
   |  |  ├── arparse.get_spect_class()
   |  |  |  ++ generates spect Class
   |  |  ├── spect.load_file(base=True)  # default
   |  |  ├── spect.set_paramlist()
   |  |  ├── spect.load_file()  # instrument specific
   |  |  ├── spect.set_paramlist()
   |  |  ├── spect.set_paramlist()   # using pyp_dict
   |  |  ├── spect.load_lines()      # command line
   |  |  ├── spect.save()
   |  |  ├── argf.save()
   |  |  ├── arparse.init(argf, spect)  # Saved into settings.
   |  |  ├── # fitsdict created
   |  |  ├── fitsdict = arload.load_headers(datlines)
   |  |  |  ├── # Checks on FITS files
   |  |  |  ├── settings.spect['check'][ch]
   |  |  ├── # Flip dispersion direction (if needed)
   |  ├── ARMLSD  # This is formally below PYPIT, but I want to reduce the indentation here
   |  |  ├── armbase.SetupScience(fitsdict)
   |  |  |  ├── arsort.sort_data()
   |  |  |  |  ├── find_standard_file()
   |  |  |  ├── arsort.match_science()
   |  |  |  |  ++ Written to settings.spect[ftag]['index']
   |  |  |  ├── arsciexp.ScienceExposure()
   |  |  |  |  ++ Generates sciexp list of ScienceExposure objects
   |  |  |  ++ setup_dict generated
   |  |  |  ├── arsort.instr_setup()
   |  |  |  |  ++ Generates setupIDs
   |  |  |  ++ group_dict generated
   |  |  |  ├── arsort.write_sorted(group_dict, setup_dict)
   |  |  |  ├── arsort.write_setup(setup_dict)


Items
-----

Items created and carried around::

    fitsdict
    settings.spect
    settings.argf
    setup_dict
    group_dict
    sciexp

