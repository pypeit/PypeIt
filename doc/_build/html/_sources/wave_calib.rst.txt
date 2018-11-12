.. _wavecalib:

.. highlight:: rest

**********************
Wavelength Calibration
**********************

.. index:: wave_calib

Basic Algorithms
================

These notes will describe the algorithms used to perform
wavelength calibration in 1D (i.e. down the slit/order)
with PypeIt.   The basic steps are:

 1. Extract 1D arc spectra down the center of each slit/order
 2. Load the parameters guiding wavelength calibration
 3. Generate the 1D wavelength fits

For the primary step (#3), the preferred approach is a
new pattern-searching algorithm currently packaged in the PypeIt
repository named `arclines`.  It is designed to estimate
the dispersion and wavelength coverage of the spectrum with
limited inputs and then automatically identify the known
arc lines.

The code is guided by the WaveCalib class, partially described
by this `WaveCalib.ipynb <https://github.com/pypeit/pypeit/blob/master/doc/nb/WaveCalib.ipynb>`_
Notebook.


Line Lists
==========

Without exception, arc line wavelengths are taken from
the `NIST database <http://physics.nist.gov/PhysRefData`_,
*in vacuum*. These data are stored as ASCII tables in the
`arclines` repository. Here are the available lamps:

======  ==========  =============
Lamp    Range (A)   Last updated
======  ==========  =============
ArI     3000-10000  21 April 2016
CdI     3000-10000  21 April 2016
CuI     3000-10000  13 June 2016
HeI     2900-12000  2 May 2016
HgI     3000-10000  May 2018
KrI     4000-12000  May 2018
NeI     3000-10000  May 2018
XeI     4000-12000  May 2018
ZnI     2900-8000   2 May 2016
======  ==========  =============

By-Hand Calibration
===================

If the automatic algorithm is failing (heaven forbid; and you should
probably raise an Issue on PypeIt if you are sure it isn't your fault),
you can input a set of pixel, wavelength values as a crutch in
your .pypeit setup file.  Here is the recommended approach:

#. Run PypeIt with --debug_arc on. This will force the code to stop inside ararc.py
#. Print the pixel values to the screen

   *  (Pdb) tcent

#. Plot the arc spectrum.

   *  (Pdb) plt.plot(yprep)
   *  (Pdb) plt.show()

#. Compare that spectrum with a known one and ID a few lines.  Write down.  Better be using vacuum wavelengths
#. Add pixel values and wavelengths to your .pypeit file, e.g.

   * arc calibrate IDpixels 872.062,902.7719,1931.0048,2452.620,3365.25658,3887.125
   * arc calibrate IDwaves 3248.4769,3274.905,4159.763,4610.656,5402.0634,5854.110


Flexure Correction
==================

By default, the code will calculate a flexure shift based on the
extracted sky spectrum (boxcar). See :doc:`flexure` for
further details.

Wavelength Frame
================

PypeIt offers several frames of reference that can used for the
wavelength scale. The first choice is whether you would like the
data to be calibrated to air or vacuum wavelengths. This option
is controlled by the argument::

    reduce calibrate wavelength air

where the default value is to calibrate to vacuum. You can also
specify 'pixel', which will save the pixel values instead of the
wavelength values (i.e. a wavelength calibration will not be
performed).  The calibration follows the Ciddor schema
(Ciddor 1996, Applied Optics 62, 958).


You can also choose if you want the wavelength scale corrected
to the heliocentric (Sun-centered), barycentric (Solar system
barycentre), or topocentric (telescope centered). None is also
an option, but this defaults to topocentric. This option
is governed by the command::

    reduce calibrate refframe barycentric

where the default value is a heliocentric wavelength scale.
More details are provided in :doc:`heliocorr`.


Developers
==========

Adding a new grating to existing instrument
-------------------------------------------

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
        #arcparam['llist'] = slf._argflag['run']['pypeitdir'] + 'data/arc_lines/kast_red.lst'

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
Run PypeIt, and check in the QA plots that the arc lines
identified by PypeIt are consistent with a pre-existing
arc line mapping, and you're done!


Validation
==========

See the iPython Notebook under test_suite for a comparison of the
wavelength solution for PypeIt vs. LowRedux.
