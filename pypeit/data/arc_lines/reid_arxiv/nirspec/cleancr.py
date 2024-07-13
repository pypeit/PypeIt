#!/usr/bin/env python

import os
import re
L = open("Use","r").readlines()
D = [l.split() for l in L]
specs = [l.split()[0] for l in L]


for spec in specs:
    inn = "%stdfocn.fits" % spec
    outcn = "%stdfocnc.fits" % spec
    outcrn = "%stdfocnc_cr.fits" % spec
    com = "~/src/dcr/dcr %s %s %s" % (inn,outcn,outcrn)
    print com
    os.system(com)

