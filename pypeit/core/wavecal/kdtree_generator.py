"""This script is used to generate the KD Tree that is needed for
the kdtree pattern matching wavelength calibration algorithm. At
present, this method is only used for calibrating ThAr lamps.

You should not run this script unless you know what you're doing,
since you could mess up the ThAr patterns that are used in the
wavelength calibration routine. This script should not be called
from within PypeIt - it should be run as a standalone script, and
it's only purpose is to generate a KD Tree with the desired patterns.
"""

# NOTE: No longer used.  Use KD tree in scikit-learn:
#   https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.KDTree.html
# See benchmarks here:
#   https://jakevdp.github.io/blog/2013/04/29/benchmarking-nearest-neighbor-searches-in-python/

import os

from pypeit.core.wavecal import waveio
from astropy.table import vstack
#import numba as nb
from scipy.spatial import cKDTree
import numpy as np
import pickle

from pypeit import data

def trigon(linelist, numsrch, maxlin):
    """ Generate a series of trigon patterns, given an input list of detections or lines from a linelist

    linelist : ndarray
      list of wavelength calibration lines (must be sorted by ascending wavelength)
    numsrch : int
      Number of consecutive detected lines used to generate a pattern. For
      example, if numsrch is 4, there are four lines (called 1 2 3 4). The following
      patterns will be generated (assuming line #1 is the left anchor):
      1 2 3  (in this case line #3 is the right anchor)
      1 2 4  (in this case line #4 is the right anchor)
      1 3 4  (in this case line #4 is the right anchor)
    maxlin : float
      Value (in pixels in the case of detections or Angstroms in the case of a linelist)
      over which the wavelength solution can be considered linear.
    """

    nptn = 3  # Number of lines used to create a pattern

    sz_l = linelist.shape[0]

    # Count the number of patterns that will be created
    cnt = 0
    for l in range(0, sz_l - nptn + 1):
        nup = (l + nptn - 1) + numsrch
        if nup > sz_l: nup = sz_l
        for ll in range(l + nptn - 1, nup):
            if (linelist[ll] - linelist[l]) > maxlin: continue
            for x in range(l + 1, ll):
                cnt += 1

    index = np.zeros((cnt, nptn), dtype=np.uint64)
    pattern = np.zeros((cnt, nptn - 2),dtype=float)

    # Generate the patterns
    cnt = 0
    for l in range(0, sz_l - nptn + 1):
        nup = (l + nptn - 1) + numsrch
        if nup > sz_l: nup = sz_l
        for ll in range(l + nptn - 1, nup):
            if (linelist[ll] - linelist[l]) > maxlin: continue
            # Create a pattern with these two endpoints
            for x in range(l + 1, ll):
                index[cnt, 0] = l
                index[cnt, 1] = x
                index[cnt, 2] = ll
                pattern[cnt, 0] = (linelist[x] - linelist[l]) / (linelist[ll] - linelist[l])
                cnt += 1
    return pattern, index



def tetragon(linelist, numsrch, maxlin):
    """ Generate a series of tetragon patterns, given an input list of detections or lines from a linelist

    linelist : ndarray
      list of wavelength calibration lines (must be sorted by ascending wavelength)
    numsrch : int
      Number of consecutive detected lines used to generate a pattern. For
      example, if numsrch is 5, there are four lines (called 1 2 3 4 5). The following
      patterns will be generated (assuming line #1 is the left anchor):
      1 2 3 4 (in this case line #4 is the right anchor)
      1 2 3 5 (in this case line #5 is the right anchor)
      1 2 4 5 (in this case line #5 is the right anchor)
      1 3 4 5 (in this case line #5 is the right anchor)
    maxlin : float
      Value (in pixels in the case of detections or Angstroms in the case of a linelist)
      over which the wavelength solution can be considered linear.
    """

    nptn = 4  # Number of lines used to create a pattern

    sz_l = linelist.shape[0]

    # Count the number of patterns that will be created
    cnt = 0
    for l in range(0, sz_l - nptn + 1):
        nup = (l + nptn - 1) + numsrch
        if nup > sz_l: nup = sz_l
        for ll in range(l + nptn - 1, nup):
            if (linelist[ll] - linelist[l]) > maxlin: continue
            for x in range(l + 1, ll - 2):
                for xx in range(x + 1, ll):
                    cnt += 1

    index = np.zeros((cnt, nptn), dtype=np.uint64)
    pattern = np.zeros((cnt, nptn - 2),dtype=float)

    # Generate the patterns
    cnt = 0
    for l in range(0, sz_l - nptn + 1):
        nup = (l + nptn - 1) + numsrch
        if nup > sz_l: nup = sz_l
        for ll in range(l + nptn - 1, nup):
            if (linelist[ll] - linelist[l]) > maxlin: continue
            # Create a pattern with these two endpoints
            for x in range(l + 1, ll - 2):
                for xx in range(x + 1, ll):
                    index[cnt, 0] = l
                    index[cnt, 1] = x
                    index[cnt, 2] = xx
                    index[cnt, 3] = ll
                    pattern[cnt, 0] = (linelist[x] - linelist[l]) / (linelist[ll] - linelist[l])
                    pattern[cnt, 1] = (linelist[xx] - linelist[l]) / (linelist[ll] - linelist[l])
                    cnt += 1
    return pattern, index


def pentagon(linelist, numsrch, maxlin):
    """
    see trigon and tetragon for an example docstring
    """
    nptn = 5  # Number of lines used to create a pattern

    sz_l = linelist.shape[0]

    # Count the number of patterns that will be created
    cnt = 0
    for l in range(0, sz_l - nptn + 1):
        nup = (l + nptn - 1) + numsrch
        if nup > sz_l: nup = sz_l
        for ll in range(l + nptn - 1, nup):
            if (linelist[ll] - linelist[l]) > maxlin: continue
            for x in range(l + 1, ll - 3):
                for xx in range(x + 1, ll - 2):
                    for xxx in range(xx + 1, ll - 1):
                        cnt += 1

    index = np.zeros((cnt, nptn), dtype=np.uint64)
    pattern = np.zeros((cnt, nptn - 2),dtype=float)

    # Generate the patterns
    cnt = 0
    for l in range(0, sz_l - nptn + 1):
        nup = (l + nptn - 1) + numsrch
        if nup > sz_l: nup = sz_l
        for ll in range(l + nptn - 1, nup):
            if (linelist[ll] - linelist[l]) > maxlin: continue
            # Create a pattern with these two endpoints
            for x in range(l + 1, ll - 3):
                for xx in range(x + 1, ll - 2):
                    for xxx in range(xx + 1, ll - 1):
                        index[cnt, 0] = l
                        index[cnt, 1] = x
                        index[cnt, 2] = xx
                        index[cnt, 3] = xxx
                        index[cnt, 4] = ll
                        pattern[cnt, 0] = (linelist[x] - linelist[l]) / (linelist[ll] - linelist[l])
                        pattern[cnt, 1] = (linelist[xx] - linelist[l]) / (linelist[ll] - linelist[l])
                        pattern[cnt, 2] = (linelist[xxx] - linelist[l]) / (linelist[ll] - linelist[l])
                        cnt += 1
    return pattern, index


def hexagon(linelist, numsrch, maxlin):
    """
    see trigon and tetragon for an example docstring
    """

    # Number of lines used to create a pattern
    nptn = 6

    sz_l = linelist.shape[0]
    # Count the number of patterns that will be created
    cnt = 0
    for l in range(0, sz_l - nptn + 1):
        nup = (l + nptn - 1) + numsrch
        if nup > sz_l: nup = sz_l
        for ll in range(l + nptn - 1, nup):
            if (linelist[ll] - linelist[l]) > maxlin: continue
            for x in range(l + 1, ll - 4):
                for xx in range(x + 1, ll - 3):
                    for xxx in range(xx + 1, ll - 2):
                        for xxxx in range(xxx + 1, ll - 1):
                            cnt += 1

    index = np.zeros((cnt, nptn),dtype=np.uint64)
    pattern = np.zeros((cnt, nptn - 2),dtype=float)

    # Generate the patterns
    cnt = 0
    for l in range(0, sz_l - nptn + 1):
        nup = (l + nptn - 1) + numsrch
        if nup > sz_l: nup = sz_l
        for ll in range(l + nptn - 1, nup):
            if (linelist[ll] - linelist[l]) > maxlin: continue
            # Create a pattern with these two endpoints
            for x in range(l + 1, ll - 4):
                for xx in range(x + 1, ll - 3):
                    for xxx in range(xx + 1, ll - 2):
                        for xxxx in range(xxx + 1, ll - 1):
                            index[cnt, 0] = l
                            index[cnt, 1] = x
                            index[cnt, 2] = xx
                            index[cnt, 3] = xxx
                            index[cnt, 4] = xxxx
                            index[cnt, 5] = ll
                            pattern[cnt, 0] = (linelist[x] - linelist[l]) / (linelist[ll] - linelist[l])
                            pattern[cnt, 1] = (linelist[xx] - linelist[l]) / (linelist[ll] - linelist[l])
                            pattern[cnt, 2] = (linelist[xxx] - linelist[l]) / (linelist[ll] - linelist[l])
                            pattern[cnt, 3] = (linelist[xxxx] - linelist[l]) / (linelist[ll] - linelist[l])
                            cnt += 1
    return pattern, index


def main(polygon, numsearch=8, maxlinear=100.0, use_unknowns=True, leafsize=30, verbose=False,
         ret_treeindx=False, outname=None, ):
    """Driving method for generating the KD Tree

    Parameters
    ----------
    polygon : int
      Number of sides to the polygon used in pattern matching
    numsearch : int
      Number of adjacent lines to use when deriving patterns
    maxlinear : float
      Over how many Angstroms is the solution deemed to be linear
    use_unknowns : bool
      Include unknown lines in the wavelength calibration (these may arise from lines other than Th I/II and Ar I/II)
    leafsize : int
      The leaf size of the tree
    """

    # Load the ThAr linelist
    line_lists_all = waveio.load_line_lists(['ThAr'])
    line_lists = line_lists_all[np.where(line_lists_all['ion'] != 'UNKNWN')]
    unknwns = line_lists_all[np.where(line_lists_all['ion'] == 'UNKNWN')]
    if use_unknowns:
        tot_list = vstack([line_lists, unknwns])
    else:
        tot_list = line_lists
    wvdata = np.array(tot_list['wave'].data)  # Removes mask if any
    wvdata.sort()

    # NIST_lines = (line_lists_all['NIST'] > 0) & (np.char.find(line_lists_all['Source'].data, 'MURPHY') >= 0)
    # wvdata = line_lists_all['wave'].data[NIST_lines]
    # wvdata.sort()

    if polygon == 3:
        if verbose: print("Generating patterns for a trigon")
        pattern, index = trigon(wvdata, numsearch, maxlinear)
    elif polygon == 4:
        if verbose: print("Generating patterns for a tetragon")
        pattern, index = tetragon(wvdata, numsearch, maxlinear)
    elif polygon == 5:
        if verbose: print("Generating patterns for a pentagon")
        pattern, index = pentagon(wvdata, numsearch, maxlinear)
    elif polygon == 6:
        if verbose: print("Generating patterns for a hexagon")
        pattern, index = hexagon(wvdata, numsearch, maxlinear)
    else:
        if verbose: print("Patterns can only be generated with 3 <= polygon <= 6")
        return None

    if outname is None:
        outname = data.get_linelist_filepath(f'ThAr_patterns_poly{polygon}_search{numsearch}.kdtree')
    outindx = outname.replace('.kdtree', '.index')
    print("Generating Tree")
    tree = cKDTree(pattern, leafsize=leafsize)
    print("Saving Tree")
    pickle.dump(tree, open(outname, 'wb'))
    print("Written KD Tree file:\n{0:s}".format(outname))
    np.save(outindx, index)
    print("Written index file:\n{0:s}".format(outindx))
    #_ = pickle.load(open(outname, 'rb'))
    #print("loaded successfully")
    if ret_treeindx:
        return tree, index

# Test
if __name__ == '__main__':
    """Set the number of sides to the polygon. Some examples:
    =========================================================    
    A trigon (polygon=3) contains
      (1) a starting point (s),
      (2) an end point (e), and
      (3) something in between (b)

                    |
                    |      |
    |               |      |
    |               |      |
    s               b      e

    Then, the value (b-s)/(e-s) is in the same coordinate system
    for both detlines and linelist.
    =========================================================    
    A tetragon (polygon=4) contains
      (1) a left line (l),
      (2) a right line (r), and
      (3) two lines in between (a, b)


                    |
              |     |      |
    |         |     |      |
    |         |     |      |
    l         a     b      r

    Then, the values (a-ll)/(r-ll) and (b-ll)/(r-ll) are in the same
    coordinate system for both detlines and linelist.
    """
    polygon = 4
    numsearch = 10
    main(polygon, numsearch=numsearch, verbose=True)
