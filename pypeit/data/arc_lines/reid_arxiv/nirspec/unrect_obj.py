#!/usr/bin/env python

from basics import *
from DataData import *
from SimpleQueue import *
import re
import glob

stars = glob.glob("*tdo.fits")

for star in stars:
    xstr = "unrect -xdist -ydist -wdist %s" % (star)
    verbex(xstr,dummy=0)

