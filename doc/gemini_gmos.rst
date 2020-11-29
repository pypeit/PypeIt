***********
Gemini GMOS
***********


Overview
========

This file summarizes several instrument specific
settings that are related to the Gemini/GMOS spectrograph.


Nod and Shuffle
===============

For the time being, we have enabled reductions of data
taken in Nod+Shuffle mode by simply replicating the calibrations.
That is, we do *not* subtract the nodded images but reduce
it as if it were a separate slit.

Arcs
====

The CuAr lamps are pretty faint in the blue which lead
to some "unique" challenges.  At present we have
lowered the default tracethresh parameter to 10, i.e.::

    par['calibrations']['tilts']['tracethresh'] = 10.  # Deals with faint CuAr lines

It is possible you will want to increase this, but unlikely.
