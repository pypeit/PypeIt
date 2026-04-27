
from IPython import embed

import numpy as np
from pypeit.core import slitdesign_matching

def test_match_positions():

    # Random number generator
    rng = np.random.default_rng(99)

    # Fake nominal positions
    nominal = np.linspace(0, 1, 5)
    d = np.mean(np.diff(nominal))

    #----------
    # Test 1: Same number of measurements

    # Vectors used to rearrange and then recover a set of "measured" positions
    isrt = np.arange(nominal.size)
    rng.shuffle(isrt)
    srt = np.argsort(isrt)

    # Generate fake data with no matching ambiguity
    measured = rng.uniform(low=-d/3, high=d/3, size=nominal.size) + nominal[isrt]

    match = slitdesign_matching.match_positions_1D(measured, nominal)

    assert np.array_equal(match, srt), 'Bad match'

    #----------
    # Test 2: Add ambiguity by appending another measurement that is a better
    # match than the existing one
    measured = np.append(measured, [nominal[2]])

    match = slitdesign_matching.match_positions_1D(measured, nominal)
    assert match[2] == nominal.size, 'Should use new measurement'
    assert np.sum(np.absolute(nominal - measured[match])) \
                < np.sum(np.absolute(nominal - measured[srt])), 'Should get a better match'

    #----------
    # Test 3: Fake a missing measurement
    measured = measured[1:5]
    match = slitdesign_matching.match_positions_1D(measured, nominal)

    assert match.size == nominal.size, 'Match vector should be the same size as nominal'
    assert match[isrt[0]] == -1, 'Should not match missing element'

    #----------
    # Test 4: Missing a measurement and effectively a repeat of an existing one
    measured = np.append(measured, [nominal[1]])
    match = slitdesign_matching.match_positions_1D(measured, nominal)

    assert match.size == nominal.size, 'Match vector should be the same size as nominal'
    assert np.all(match > -1), 'All should match'
    assert np.sum(np.absolute(nominal - measured[match]) > d/2) == 1, 'Should be one large outlier'

    #----------
    # Test 5: Same as test 4, but include a tolerance that should remove the outlier
    match = slitdesign_matching.match_positions_1D(measured, nominal, tol=d/2)

    assert match.size == nominal.size, 'Match vector should be the same size as nominal'
    assert match[isrt[0]] == -1, 'Should not match missing element'


