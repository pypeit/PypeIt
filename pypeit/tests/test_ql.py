
from IPython import embed

from pypeit.metadata import PypeItMetaData
from pypeit.scripts import ql

def test_merge():

    # There are two existing setups, 'A', 'B', and the new reductions use the
    # same setup as 'B'.
    existing_setups = ['A', 'B']
    calib_match = {'A': dict(setup='B')}
    frame_setup = ['A', 'A', 'A']
    new_frame_setup = ql.merge_setups(existing_setups, calib_match, frame_setup)
    assert new_frame_setup == ['B', 'B', 'B'], 'Bad setup replacement'

    # There is one existing setup, 'A', but the new reductions are from a
    # different setup.  This assigns that new setup to 'B'.
    existing_setups = ['A']
    calib_match = {'A': None}
    frame_setup = ['A', 'A', 'A']
    new_frame_setup = ql.merge_setups(existing_setups, calib_match, frame_setup)
    assert new_frame_setup == ['B', 'B', 'B'], 'Bad setup replacement'

    # Somehow a previous set of calibration set was deleted, such that the setup
    # identifier is irregular.  This creates the new setup identifier in sequence.
    existing_setups = ['A', 'C']
    calib_match = {'A': None}
    frame_setup = ['A', 'A', 'A']
    new_frame_setup = ql.merge_setups(existing_setups, calib_match, frame_setup)
    assert new_frame_setup == ['B', 'B', 'B'], 'Bad setup replacement'

    # There are two existing setups, 'A', 'B', and two new setups to be reduced.
    # The new 'A' setup matches the existing 'A' setup, but the new 'B' setup
    # does not match either of the existing ones.  This matches the relevant
    # frames to the existing setup and adds the new setup identifier.
    existing_setups = ['A', 'B']
    calib_match = {'A': dict(setup='A'), 'B': None}
    frame_setup = ['A', 'A', 'A', 'B', 'B', 'B']
    new_frame_setup = ql.merge_setups(existing_setups, calib_match, frame_setup)
    assert new_frame_setup == ['A', 'A', 'A', 'C', 'C', 'C'], 'Bad setup replacement'

    # Adding multiple new setups.
    existing_setups = ['A', 'B']
    calib_match = {'A': None, 'B': None}
    frame_setup = ['A', 'A', 'A', 'B', 'B', 'B']
    new_frame_setup = ql.merge_setups(existing_setups, calib_match, frame_setup)
    assert new_frame_setup == ['C', 'C', 'C', 'D', 'D', 'D'], 'Bad setup replacement'

    # Make sure that the function can handle the identifiers with two
    # characters.
    existing_setups = ['AA']
    calib_match = {'A': dict(setup='AA')}
    frame_setup = ['A', 'A', 'A']
    new_frame_setup = ql.merge_setups(existing_setups, calib_match, frame_setup)
    assert new_frame_setup == ['AA', 'AA', 'AA'], 'Bad setup replacement'



