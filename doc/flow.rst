*********
Code Flow
*********

Overview
========

Setup
=====

Here is the code flow for the
`pypit-setup`_ script::

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

