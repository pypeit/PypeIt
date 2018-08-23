"""This script is used to generate the KD Tree that is needed for
the kdtree pattern matching wavelength calibration algorithm. At
present, this method is only used for calibrating ThAr lamps.

You should not run this script unless you know what you're doing,
since you could mess up the ThAr patterns that are used in the
wavelength calibration routine. This script should not be called
from within PypeIt - it should be run as a standalone script, and
it's only purpose is to generate a KD Tree with the desired patterns.
"""

from pypeit.core.wavecal import waveio
from astropy.table import vstack
import numba as nb
from scipy.spatial import cKDTree
import numpy as np
import pickle


@nb.jit(nopython=True, cache=True)
def trigon(linelist, numsrch, maxlin):
    """
    linelist : list of wavelength calibration lines (must be sorted by ascending wavelength)
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
            for x in range(l + 1, ll - 1):
                cnt += 1

    index = np.zeros((cnt, nptn), dtype=nb.types.uint64)
    pattern = np.zeros((cnt, nptn - 2), dtype=nb.types.float64)

    # Generate the patterns
    cnt = 0
    for l in range(0, sz_l - nptn + 1):
        nup = (l + nptn - 1) + numsrch
        if nup > sz_l: nup = sz_l
        for ll in range(l + nptn - 1, nup):
            if (linelist[ll] - linelist[l]) > maxlin: continue
            # Create a pattern with these two endpoints
            for x in range(l + 1, ll - 1):
                index[cnt, 0] = l
                index[cnt, 1] = x
                index[cnt, 2] = ll
                pattern[cnt, 0] = (linelist[x] - linelist[l]) / (linelist[ll] - linelist[l])
                cnt += 1
    return pattern, index


@nb.jit(nopython=True, cache=True)
def tetragon(linelist, numsrch, maxlin):
    """
    linelist : list of wavelength calibration lines (must be sorted by ascending wavelength)
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
                for xx in range(x + 1, ll - 1):
                    cnt += 1

    index = np.zeros((cnt, nptn), dtype=nb.types.uint64)
    pattern = np.zeros((cnt, nptn - 2), dtype=nb.types.float64)

    # Generate the patterns
    cnt = 0
    for l in range(0, sz_l - nptn + 1):
        nup = (l + nptn - 1) + numsrch
        if nup > sz_l: nup = sz_l
        for ll in range(l + nptn - 1, nup):
            if (linelist[ll] - linelist[l]) > maxlin: continue
            # Create a pattern with these two endpoints
            for x in range(l + 1, ll - 2):
                for xx in range(x + 1, ll - 1):
                    index[cnt, 0] = l
                    index[cnt, 1] = x
                    index[cnt, 2] = xx
                    index[cnt, 3] = ll
                    pattern[cnt, 0] = (linelist[x] - linelist[l]) / (linelist[ll] - linelist[l])
                    pattern[cnt, 1] = (linelist[xx] - linelist[l]) / (linelist[ll] - linelist[l])
                    cnt += 1
    return pattern, index


@nb.jit(nopython=True, cache=True)
def pentagon(linelist, numsrch, maxlin):
    """
    linelist : list of wavelength calibration lines (must be sorted by ascending wavelength)
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

    index = np.zeros((cnt, nptn), dtype=nb.types.uint64)
    pattern = np.zeros((cnt, nptn - 2), dtype=nb.types.float64)

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


@nb.jit(nopython=True, cache=True)
def hexagon(linelist, numsrch, maxlin):
    """
    linelist : list of wavelength calibration lines (must be sorted by ascending wavelength)
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

    index = np.zeros((cnt, nptn), dtype=nb.types.uint64)
    pattern = np.zeros((cnt, nptn - 2), dtype=nb.types.float64)

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


def main(polygon, numsearch=8, maxlinear=15.0, use_unknowns=True, leafsize=30):
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

    if polygon == 3:
        print("Generating patterns for a trigon")
        pattern, index = trigon(wvdata, numsearch, maxlinear)
    elif polygon == 4:
        print("Generating patterns for a tetragon")
        pattern, index = tetragon(wvdata, numsearch, maxlinear)
    elif polygon == 5:
        print("Generating patterns for a pentagon")
        pattern, index = pentagon(wvdata, numsearch, maxlinear)
    elif polygon == 6:
        print("Generating patterns for a hexagon")
        pattern, index = hexagon(wvdata, numsearch, maxlinear)
    else:
        print("Patterns can only be generated with 3 <= polygon <= 6")
        return None

    outname = '../../data/arc_lines/lists/ThAr_patterns_poly{0:d}_search{1:d}.kdtree'.format(polygon, numsearch)

    print("Generating Tree")
    tree = cKDTree(pattern, leafsize=leafsize)
    print("Saving Tree")
    pickle.dump(tree, open(outname, 'wb'))
    print("Written KD Tree file:\n{0:s}".format(outname))
    #_ = pickle.load(open(outname, 'rb'))
    #print("loaded successfully")


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
    numsearch = 8
    main(polygon, numsearch=numsearch)
