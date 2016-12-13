.. highlight:: rest

****************
Bias Subtraction
****************


Overview
========

The code can perform bias subtraction using input bias frames
or by analyzing and subtracting off an estimate from the overscan
region(s).  The default for ARMLSD is to use bias
frames and require that several be provided.
A future implementation will combine the two approaches.

Methods
=======

Subtract Bias Frame
-------------------

This method combines the set of input bias frames and
subtracts the resulting MasterBias from all other frames.
This is the default in ARMLSD when bias frames are not
provided.  It can be explicitly enforced by adding::

    bias useframe bias

to the PYPIT reduction file.

Subtract Overscan
-----------------

Analyze the overscan region and subtract a replicated
version of the result from the full image.  This
method is applied by adding::

    bias useframe overscan

to the PYPIT reduction file.

Overscan Algorithms
+++++++++++++++++++

SavGol, Polynomial
