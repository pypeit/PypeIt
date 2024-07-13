#!/usr/bin/env python

from basics import *
from DataData import *
from SimpleQueue import *
import re

def rectdiff(spec,back):

    outs = re.sub(r"td\.fits",r"tdo.fits",spec)
    ins = re.sub(r"td\.fits",r"t.fits",spec)
    outb = re.sub(r"td\.fits",r"tdo.fits",back)
    mods = re.sub(r"td\.fits",r"tdm.fits",spec)
    modb = re.sub(r"td\.fits",r"tdm.fits",back)

    com = "getrect -xdist -ydist -x 1 -y 4 -nx 8 -dy 7 -b 3 -nsub2 4 -niter 1  "
    c = com + ins
    verbex(c,dummy=0)
    copyit = "copyrect -xdist %s %s %s " % (ins,spec,back)
    verbex(copyit,dummy=0)

    com = "skyrect -xdist -ydist -b 1 -kx 3 -ky 3 -dkx 0.5 -skyint 1024 -check 7 -pct 50 -skycut 5 -spass -5 -rn 25 -g 2 "
    c = com + spec
    verbex(c,dummy=0)
    c = com + back
    verbex(c,dummy=0)



    diff = "efits %s %s \"i1 - i2\" %s" % (back, modb, outb)
    verbex(diff,dummy=0)

    diff = "efits %s %s \"i1 - i2\" %s" % (spec, mods, outs)
    verbex(diff,dummy=0)
        
    diffn = re.sub(r"td\.fits",r"tdod.fits",spec)
    diff = "efits %s %s \"i1 - i2\" %s" % (outs, outb, diffn)
    verbex(diff,dummy=0)



pairlist = sys.argv[1]

L = open(pairlist,"r").readlines()
specs = [l.split()[0] for l in L]

for nspec in arange(len(specs)-1):

    nback = nspec + 1
    back = specs[nback]
    spec = specs[nspec]

    rectdiff(spec,back)
    
nback = len(specs)-2
nspec = len(specs)-1

back = specs[nback]
spec = specs[nspec]
rectdiff(spec,back)
    
