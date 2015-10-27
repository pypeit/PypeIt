#!/usr/bin/env python
#
'''
This needs to be run in the TEST_SUITES directory.
Usually found in a Dropbox.
'''
import sys, os, os.path
import pdb
import subprocess
import warnings 

#sys.path.append(os.getenv('PYPIT')+'/src/')
#import pypit, arload

# Point to all sub-folders 
walk = os.walk('./')
instruments = next(walk)[1]

pwd = os.getcwd()

# Loop on instruments
for instr in instruments:
    #if instr in ['Kast_blue','Kast_red']: # For testing
    #    continue
    # Setups
    setups = next(os.walk(instr))[1]
    for setup in setups:
        # Look for redux file in PYPIT
        redfile = os.getenv('PYPIT')+'/test_suite/'+instr.lower()+'_'+setup.lower()+'.red'
        if not os.path.exists(redfile):
            warnings.warn('No redux file: {:s}'.format(redfile))
            warnings.warn('Not testing..')
            continue
        # Edit data directory 
        with open(redfile, 'r') as infile:
            lines = infile.readlines()
        for kk,iline in enumerate(lines):
            if 'data read' in iline:
                dpth = lines[kk+1]
                i0 = dpth.rfind('/')
                newdpth = ' '+pwd+'/'+instr+'/'+setup+dpth[i0:]
                lines[kk+1] = newdpth
        # Generate folder as need be
        wdir = os.getenv('TST_PYPIT')+'/'+instr
        if not os.path.exists(wdir):
            os.makedirs(wdir)
        # Write to TST_PYPIT
        outfile = wdir+'/'+instr.lower()+'_'+setup.lower()+'.red'
        with open(outfile, 'w') as ofile:
            for iline in lines:
                ofile.writelines(iline)
        # Run       
        logfile = wdir+'/'+instr.lower()+'_'+setup.lower()+'.log'
        print('Running pypit on {:s} --- '.format(outfile))
        with open(logfile,'w') as f:
            subprocess.call(['python', os.getenv('PYPIT')+'/src/pypit.py', outfile], stderr=f, cwd=wdir)#, shell=True)
        print('Done running pypit on {:s} --- '.format(outfile))
        subprocess.call(['tail', logfile])
        # Need some merit of success..


# cd to Test Suite folder
#os.chdir(os.getenv('TST_PYPIT'))

