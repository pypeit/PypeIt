from IPython import embed
import numpy as np

from pypeit.core.fitting import robust_fit
from pypeit.core import trace

# Left and right edges uses for testing
# Pulled from keck_hires_hs1700+6416_h45ah_red_b2_ech_0.00_xd_-0.00_1x2
def ex_lefts():
    return np.array([ 387.48466467,  419.74944269,  452.7626946 ,  486.59039802,
                      521.24024038,  556.68944113,  592.94414902,  630.05184547,
                      668.02586614,  706.90703733,  746.73026992,  787.53311534,
                      829.35374693,  872.23568857,  916.22050729,  961.34921725,
                     1007.74046978, 1105.14334103, 1155.39584898, 1207.07871502,
                     1260.25865643, 1315.01036114, 1371.41261315, 1429.54436242,
                     1489.49852871, 1551.36901195, 1615.26125598, 1681.28568972,
                     1749.55792592, 1820.21637138, 1893.40015359, 1969.26214779,
                     2048.00900207, 2130.08866216, 2215.01472161, 2303.39330558,
                     2395.47431655, 2491.51937196, 2591.80187799, 2696.6263941 ,
                     2806.40894244, 2921.55324348])

def ex_rights():
    return np.array([ 417.3750106 ,  450.0306145 ,  483.27774126,  517.18974649,
                      551.84420416,  587.25718106,  623.53074822,  660.61932348,
                      698.58047385,  737.44928937,  777.25633863,  818.0477269 ,
                      859.85375605,  902.71292489,  946.67158948,  991.79414868,
                     1038.49793985, 1135.55065652, 1185.77954643, 1237.44393944,
                     1290.60859833, 1345.34978922, 1401.7396414 , 1459.86011192,
                     1519.8027442 , 1581.66637517, 1645.55134019, 1711.58372343,
                     1779.84232834, 1850.50167339, 1923.68888709, 1999.56313505,
                     2078.21066646, 2160.49703097, 2245.45457385, 2333.85552002,
                     2425.93486802, 2521.95453935, 2622.25346894, 2727.15466796,
                     2837.05077696, 2952.32806799])

def test_find_missing_orders():
    left = ex_lefts()
    right = ex_rights()

    gap = left[1:] - right[:-1]
    width = right - left
    cen = (right + left)/2

    width_fit = robust_fit(cen, width, 3, function='legendre',
                           lower=3., upper=3., maxiter=5, sticky=True)
    gap_fit = robust_fit(cen[:-1], gap, 3, function='legendre',
                         lower=3., upper=3., maxiter=5, sticky=True)

    order_cen, order_missing = trace.find_missing_orders(cen, width_fit, gap_fit)

    assert order_cen.size > cen.size, 'Should be missing orders'
    assert np.sum(order_missing) == 1, 'Should find 1 missing order'
    assert np.where(order_missing)[0][0] == 17, 'Missing order incorrect'


def test_extrapolate_orders():
    left = ex_lefts()
    right = ex_rights()

    gap = left[1:] - right[:-1]
    width = right - left
    cen = (right + left)/2

    width_fit = robust_fit(cen, width, 3, function='legendre',
                           lower=3., upper=3., maxiter=5, sticky=True)
    gap_fit = robust_fit(cen[:-1], gap, 3, function='legendre',
                         lower=3., upper=3., maxiter=5, sticky=True)

    min_spat = 200.
    max_spat = 3100.
    lower_order_cen, upper_order_cen \
            = trace.extrapolate_orders(cen, width_fit, gap_fit, min_spat, max_spat)
    
    assert lower_order_cen.size == 6, 'Incorrect number of lower extrapolations'
    assert upper_order_cen.size == 1, 'Incorrect number of upper extrapolations'
    assert np.all(lower_order_cen > min_spat), 'Extrapolations outside specified range'
    assert np.all(upper_order_cen < max_spat), 'Extrapolations outside specified range'



