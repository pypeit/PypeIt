***********
Gemini GMOS
***********


Overview
========

This file summarizes several instrument specific
settings that are related to the Gemini/GMOS spectrograph.


Arcs
====

The CuAr lamps are pretty faint in the blue which lead
to some "unique" challenges.  At present we have
lowered the default tracethresh parameter to 10, i.e.::

    par['calibrations']['tilts']['tracethresh'] = 10.  # Deals with faint CuAr lines

It is possible you will want to increase this, but unlikely.
