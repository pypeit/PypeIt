.. highlight:: rest

****************
Bias Subtraction
****************


Overview
========

The code can perform bias subtraction using input bias frames
or by analyzing and subtracting off an estimate from the overscan
region(s).  The latter is strongly recommended at the present time
and is the default for armlsd.

A future implementation will be to combine the two approaches.

Overscan Algorithms
===================

SavGol, Polynomial
