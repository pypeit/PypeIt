.. highlight:: rest

*****
Setup
*****

This document will describe how to organize your
raw files and the settings files for running PYPIT.


Raw Files
=========


Settings File
=============
To begin reducing your data, you'll need a settings file
with the extension ".pypit" that will provide the information
necessary for PYPIT to know how to proceed with reducing your
particular set of data.

This section will instruct you on how to build a .pypit
setting file from scratch. For reference, there are
existing settings files in PYPIT development suite.
The PYPIT development suite is recommended for download in
:doc:`installing`, and the relevant settings files are located
in::

    ~/PYPIT-development-suite/pypit_files/

However, let's look at how to make one from scratch.

Settings File, line by line
+++++++++++++++++++++++++++
Create a .pypit file. Name it anything you want, but for example,
it's useful to have: the instrument name, the grating or grism used,
the dichroic, etc. For example, we could call our settings file
'lris_blue_long_600_4000_d560.pypit', for our data as collected
on LRIS's blue arm, in long slit mode, using the 600/4000 grism
and d560 dichroic.

You can make any comments in your .pypit settings file with a
pound sign::

    # This is a comment line

The first thing you'll need to include is any changes on the
default settings. Most of these, except the spectrograph, aren't
truly necessary in the most basic settings file (there are existing
set defaults), but will be useful for running PYPIT::

    # Change the default settings
    run ncpus 1                     # number of CPUs to use; can also negative integers,
                                    so -1 means all but one CPU
    run spectrograph lris_blue      # the spectrograph (+arm, if necessary) this set of data is from;
                                    see README for list of available instruments
    output verbosity 2                   # level of screen output; 0 = no output, 1 = low level of output;
                                    2 = output everything
    output overwrite True              # overwrite any existing output files?
    output sorted lris_blue_long_600_4000_d560     # name of output files

Next, tell PYPIT where your raw data lives!::

    # Read in the data
    data read
     /Users/path/to/your/raw/data/*.fits
    data end

Then, give PYPIT some information about your raw data. For
example, PYPIT only accepts calibration files if they were
created within a time window of the science frame of interest.
You can set your own time window here. PYPIT also requires a
certain number of each type of calibration file to be matched
with the science frame, and here you can set what you want the
minimum to be::

    spect read
     #fits calwin 1000.     # calibration window; default window is 12 hrs;
                            here it is changed to 1000. hrs
     pixelflat number 1       # number of pixel flats needed for data reduction
     bias number 3          # number of bias frames; note that in this case,
                            PYPIT will combine the 3 biases into a master bias
     arc number 1           # number of arcs
     trace number 1         # number of trace frames
    spect end

Settings File, as a whole
+++++++++++++++++++++++++
With that, the most basic settings file looks something like this::

    # Change the default settings
    run ncpus 1
    run spectrograph lris_blue
    output verbosity 2
    output overwrite True
    output sorted lris_blue_long_600_4000_d560

    # Read in the data
    data read
     /Users/path/to/your/raw/data/*.fits
    data end

    spect read
     #fits calwin 1000.

     pixelflat number 1
     bias number 3
     arc number 1
     trace number 1
    spect end

You can now run PYPIT with this .pypit settings file! See how in
:doc:`running`.

Additional parameters for the Settings File
+++++++++++++++++++++++++++++++++++++++++++
In addition to the basic settings file as shown above, there
are other parameters that you can tell PYPIT to run by::

    # Reduce


    spect read
     # not needed if everything runs smoothly. Check your .lst file and make sure that each frame was identified correctly (that each file is properly identified as a bias, arc, slitflat, standard, science). If any file was misidentified, you can force the file type to be something different below (note that you can also identify your various calibration and science files below if you don't want to deal with the .lst file):

     #set bias     b150910_2036.fits.gz
     #set bias     b150910_2037.fits.gz
     #set bias     b150910_2038.fits.gz
     #set pixelflat  b150910_2051.fits.gz
     #set trace    b150910_2051.fits.gz
     #set standard b150910_2083.fits.gz
     ################################
    spect end