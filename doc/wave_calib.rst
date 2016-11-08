.. _wavecalib:

.. highlight:: rest

**********************
Wavelength Calibration
**********************

.. index:: wave_calib

Basic Algorithms
================

These notes will describe the algorithms used to perform
wavelength calibration with PYPIT.

Line Lists
==========

Without exception, arc line wavelengths are taken from
the `NIST database <http://physics.nist.gov/PhysRefData`_,
in vacuum. These data are stored as ASCII tables in data/arc_lines/NIST.
Here are the available lamps:

======  ==========  =============
Lamp    Range (A)   Last updated
======  ==========  =============
ArI     3000-10000  21 April 2016
CdI     3000-10000  21 April 2016
CuI     3000-10000  13 June 2016
HeI     2900-12000  2 May 2016
HgI     3000-10000  21 April 2016
KrI     4000-12000  14 April 2016
NeI     3000-10000  21 April 2016
XeI     4000-12000  14 April 2016
ZnI     2900-8000   2 May 2016
======  ==========  =============

By-Hand Calibration
===================

If the automatic algorithm is failing (heaven forbid; and you should
probably raise an Issue on PYPIT if you are sure it isn't your fault),
you can input a set of pixel, wavelength values as a crutch in
your .pypit setup file.  Here is the recommended approach:

#. Run PYPIT with --debug_arc on. This will force the code to stop inside ararc.py
#. Print the pixel values to the screen

   *  (Pdb) tcent

#. Plot the arc spectrum.

   *  (Pdb) plt.plot(yprep)
   *  (Pdb) plt.show()

#. Compare that spectrum with a known one and ID a few lines.  Write down.  Better be using vacuum wavelengths
#. Add pixel values and wavelengths to your .pypit file, e.g.

   * arc calibrate IDpixels 872.062,902.7719,1931.0048,2452.620,3365.25658,3887.125
   * arc calibrate IDwaves 3248.4769,3274.905,4159.763,4610.656,5402.0634,5854.110

Validation
==========

See the iPython Notebook under test_suite for a comparison of the
wavelength solution for PYPIT vs. LowRedux.

Adding a new grating to existing instrument
===========================================

This section describes how to add a new
wavelength solution for a new instrument and/or
grating.

Open ararc.py to add information to an already
existing instrument.

In the method setup_param, add the new disperser to the
list, under the appropriate instrument. In setup_param,
all the defaults for all instruments and gratings are listed
first::

    # Defaults
    arcparam = dict(llist='',
        disp=0.,           # Ang/unbinned pixel
        b1=0.,               # Pixel fit term (binning independent)
        b2=0.,               # Pixel fit term
        wvmnx=[2900.,12000.],# Guess at wavelength range
        disp_toler=0.1,      # 10% tolerance
        match_toler=3.,      # Matcing tolerance (pixels)
        func='legendre',     # Function for fitting
        n_first=1,           # Order of polynomial for first fit
        n_final=4,           # Order of polynomial for final fit
        nsig_rej=2.,         # Number of sigma for rejection
        nsig_rej_final=3.0,  # Number of sigma for rejection (final fit)
        Nstrong=13)          # Number of lines for auto-analysis

Find the instrument you'd like to add a grating to. For
example, to add the 1200/5000 grating to KAST red, the
section of interest would be::

    elif sname=='kast_red':
        lamps = ['HgI','NeI','ArI']
        #arcparam['llist'] = slf._argflag['run']['pypitdir'] + 'data/arc_lines/kast_red.lst'

And the following lines should be added::

        elif disperser == '1200/5000':
            arcparam['disp']=1.17 # This information is on the instrument's website
            arcparam['b1']= 1./arcparam['disp']/slf._msarc[det-1].shape[0]
            arcparam['wvmnx'][1] = 5000. # This is a guess at max wavelength covered
            arcparam['n_first']=2 # Should be able to lock on

Now in armlsd.py, put a stop after wavelength calibration
to check that the arc lines were correctly identified for
this new disperser. To do this, in method ARMLSD, find::

                # Extract arc and identify lines
                wv_calib = ararc.simple_calib(slf, det)
                slf.SetFrame(slf._wvcalib, wv_calib, det)
                slf._qa.close()
                debugger.set_trace()

Note that the last two lines were added so that the QA
plots can be correctly closed, and the process stopped.
Run PYPIT, and check in the QA plots that the arc lines
identified by PYPIT are consistent with a pre-existing
arc line mapping, and you're done!

Flexure Correction
==================

By default, the code will calculate a flexure shift based on the
extracted sky spectrum (boxcar). See :doc:`flexure` for
further details.


