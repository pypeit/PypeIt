.. highlight:: rest

**************
Object Tracing
**************

This document describes how the code traces
each object found within a slit/order.

Parameters
==========

The following table describes parameters related to tracing the
objects down the slit/order.  These parameters also follow a
prefix of `trace object`.

============== ==============================  ==================================================
Parameter      Options                         Description
============== ==============================  ==================================================
order          int; default=2                  Order of the polynomial function to be used to fit
                                               the object trace in each slit.
function       polynomial,legendre, chebyshev  Function to be used to trace the object in each
============== ==============================  ==================================================
