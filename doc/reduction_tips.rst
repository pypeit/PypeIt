==============
Reduction Tips
==============

Overview
========

This document describes commonly adjusted
:doc:`pypeit_par` related to the sky-subtraction
and extraction of the spectra.

It assumes:

1. You have already generated a set of :doc:`out_spec2D` files.

2. You have inspected the spectra wish to improve aspects of them.


Object Finding
==============

Refer to :doc:`object_finding` for full details on the algorithm.
This process is guided by the :ref:`pypeit_par:FindObjPar Keywords`.

The most common to modify is **sig_thresh** which sets the
search for any source with *peak* flux in excess of **sig_thresh**
times the RMS.  The default is 10 and you may wish to
reduce this parameter, e.g.::

    [reduce]
      [[findobj]]
          sig_thresh = 3.

This will search for any source with peak flux 3-sigma above the
estimated RMS in the smashed slit profile.

