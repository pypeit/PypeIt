"""
Module to run tests on PypeItPar classes
"""
import os

from IPython import embed

import pytest

from pypeit.par import newparset
from pypeit.par import newpypeitpar
from pypeit.par import pypeitpar
from pypeit.par import parset
from pypeit.par import util
from pypeit.spectrographs.util import load_spectrograph


def test_eval_tuple():
    # one tuple
    t = ['(1', '2', '3)']
    assert util.eval_tuple(t) == [(1,2,3)], 'Bad tuple evaluation'
    # two tuples
    t = ['(1', '2)', '(3', '4)']
    assert util.eval_tuple(t) == [(1,2),(3,4)], 'Bad tuple evaluation'
    # a tuple with a mix of int and float types
    t = ['(1', '2.3)', '(3', '4)']
    assert util.eval_tuple(t) == [(1,2.3),(3,4)], 'Bad tuple evaluation'
    # a tuple with a string
    t = ['(1', '2.3)', '(test', '4)']
    assert util.eval_tuple(t) == [(1, 2.3), ('test', 4)], 'Bad tuple evaluation'


def test_eval_list():
    # one list
    t = ['[1', '2', '3]']
    assert util.eval_list(t) == [[1,2,3]], 'Bad list evaluation'
    # two lists
    t = ['[1', '2]', '[3', '4]']
    assert util.eval_list(t) == [[1,2],[3,4]], 'Bad list evaluation'
    # a list with a mix of int and float types
    t = ['[1', '2.3]', '[3', '4]']
    assert util.eval_list(t) == [[1,2.3],[3,4]], 'Bad list evaluation'
    # a list with a string
    t = ['[1', '2.3]', '[test', '4]']
    assert util.eval_list(t) == [[1, 2.3], ['test', 4]], 'Bad list evaluation'


def test_framegroup():
    pypeitpar.FrameGroupPar()

def test_framegroup_types():
    t = pypeitpar.FrameGroupPar.valid_frame_types()
    assert 'bias' in t, 'Expected to find \'bias\' in list of valid frame types'

def test_processimages():
    pypeitpar.ProcessImagesPar()

def test_flatfield():
    pypeitpar.FlatFieldPar()

def test_flexure():
    pypeitpar.FlexurePar()

def test_coadd1d():
    pypeitpar.Coadd1DPar()

def test_coadd2d():
    pypeitpar.Coadd2DPar()

def test_cube():
    pypeitpar.CubePar()

def test_fluxcalibrate():
    pypeitpar.FluxCalibratePar()

def test_sensfunc():
    pypeitpar.SensFuncPar()

def test_sensfuncuvis():
    pypeitpar.SensfuncUVISPar()

def test_telluric():
    pypeitpar.TelluricPar()

# TODO: Valid spectrographs are not longer read by pypeit.pypeitpar; it causes
# a circular import.
#def test_spectrographs():
#    s = pypeitpar.ReduxPar.valid_spectrographs()
#    assert 'keck_lris_blue' in s, 'Expected to find keck_lris_blue as a spectrograph!'

def test_redux():
    pypeitpar.ReduxPar()

def test_wavelengthsolution():
    pypeitpar.WavelengthSolutionPar()

def test_edgetrace():
    pypeitpar.EdgeTracePar()

def test_wavetilts():
    pypeitpar.WaveTiltsPar()

def test_reduce():
    pypeitpar.ReducePar()

def test_findobj():
    pypeitpar.FindObjPar()

def test_skysub():
    pypeitpar.SkySubPar()

def test_extraction():
    pypeitpar.ExtractionPar()

def test_calibrations():
    pypeitpar.CalibrationsPar()

def test_pypeit():
    pypeitpar.PypeItPar()

def test_fromcfgfile():
    pypeitpar.PypeItPar.from_cfg_file()

def test_writecfg():
    default_file = 'default.cfg'
    if os.path.isfile(default_file):
        os.remove(default_file)
    pypeitpar.PypeItPar.from_cfg_file().to_config(cfg_file=default_file)
    assert os.path.isfile(default_file), 'No file written!'
    os.remove(default_file)

def test_readcfg():
    default_file = 'default.cfg'
    if not os.path.isfile(default_file):
        pypeitpar.PypeItPar.from_cfg_file().to_config(cfg_file=default_file)
    pypeitpar.PypeItPar.from_cfg_file('default.cfg')
    os.remove(default_file)

def print_diff_header(name1, name2):
    print("----------------------------------------------------------------")
    print(f"Comparing {name1} to {name2}")
    print("----------------------------------------------------------------")

def diff_pars(par1, name1, par2, name2):
    """Print human readable diff output when comparing pars"""
    
    printed_header = False
    result = True
    first_par_keys = set(par1.keys())
    second_par_keys = set(par2.keys())
    missing_from_second = first_par_keys - second_par_keys
    if len(missing_from_second) > 0:
        if not printed_header:
            print_diff_header(name1, name2)
            printed_header=True
        print(f"{name1} keys missing from {name2}: {','.join(missing_from_second)}")
        result=False

    missing_from_first = second_par_keys - first_par_keys
    if len(missing_from_first) > 0:
        if not printed_header:
            print_diff_header(name1, name2)
            printed_header=True
        print(f"{name2} keys missing from {name1}: {','.join(missing_from_second)}")
        result=False

    for key in first_par_keys & second_par_keys:
        first_value = par1[key]
        second_value = par2[key]
        if isinstance(first_value, parset.ParSet) and isinstance(second_value, parset.ParSet):
            if not diff_pars(first_value, f"{name1}->{key}", second_value, f"{name2}->{key}"):
                result=False
        else:
            if not first_value == second_value:
                result = False
                if not printed_header:
                    print_diff_header(name1, name2)
                    printed_header=True
                else:
                    print("\n")
                print(f"{name1}->{key} differs from {name2}->{key}")
                print(f"{name1}->{key} is type {type(first_value)}, value {first_value}")
                print(f"{name2}->{key} is type {type(second_value)}, value {second_value}")
    
    return result

def test_readwritecfg():
    # Test writing and re-reading a file results in the same
    # values
    default_file = 'default.cfg'
    if os.path.isfile(default_file):
        os.remove(default_file)
    default_par = pypeitpar.PypeItPar()
    default_par.to_config(cfg_file=default_file)
    assert os.path.isfile(default_file), 'No file written!'
    read_par = pypeitpar.PypeItPar.from_cfg_file(default_file)

    # Comparison that produces nice output when they differ
    assert diff_pars(default_par, "Default", read_par, "Read")
    os.remove(default_file) 

def test_mergecfg():
    # Create a file with the defaults
    user_file = 'user_adjust.cfg'
    if os.path.isfile(user_file):
        os.remove(user_file)
    p = pypeitpar.PypeItPar.from_cfg_file()

    # Make some modifications
    p['rdx']['spectrograph'] = 'keck_lris_blue'
    p['calibrations']['biasframe']['useframe'] = 'overscan'

    # Write the modified config
    p.to_config(cfg_file=user_file)
    assert os.path.isfile(user_file), 'User file was not written!'

    # Read it in as a config to merge with the defaults
    p = pypeitpar.PypeItPar.from_cfg_file(merge_with=(user_file))

    # Check the values are correctly read in
    assert p['rdx']['spectrograph'] == 'keck_lris_blue', 'Test spectrograph is incorrect!'
    assert p['calibrations']['biasframe']['useframe'] == 'overscan', \
                'Test biasframe:useframe is incorrect!'

    # Clean-up
    os.remove(user_file)

def test_sync():
    p = pypeitpar.PypeItPar()
    proc = pypeitpar.ProcessImagesPar()
    proc['combine'] = 'median'
#    proc['cr_sigrej'] = 20.5
    p.sync_processing(proc)
    assert p['scienceframe']['process']['combine'] == 'median'
    assert p['calibrations']['biasframe']['process']['combine'] == 'median'
    # Sigma rejection of cosmic rays for arc frames is already turned
    # off by default
#    assert p['calibrations']['arcframe']['process']['cr_sigrej'] < 0
#    assert p['calibrations']['traceframe']['process']['cr_sigrej'] == 20.5

def test_telescope():
    pypeitpar.TelescopePar()

def test_fail_badpar():
    p = load_spectrograph('gemini_gnirs_echelle').default_pypeit_par()

    # Faults because there's no junk parameter
    cfg_lines = ['[calibrations]', '[[biasframe]]', '[[[process]]]', 'junk = True']
    with pytest.raises(ValueError):
        _p = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=p.to_config(), 
                                                merge_with=cfg_lines) # Once as list
    
def test_fail_badlevel():
    p = load_spectrograph('gemini_gnirs_echelle').default_pypeit_par()

    # Faults because process isn't at the right level (i.e., there's no
    # process parameter for CalibrationsPar)
    cfg_lines = ['[calibrations]', '[[biasframe]]', '[[process]]', 'cr_reject = True']
    with pytest.raises(ValueError):
        _p = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=p.to_config(), 
                                                merge_with=(cfg_lines,))  #Once as tuple


def test_lists():
    # Initialise the parset
    p = load_spectrograph('keck_kcwi').default_pypeit_par()

    # Test with a single element list
    p['calibrations']['alignment']['locations'] = [0.5]
    _p = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=p.to_config())  # Once as tuple
    assert(isinstance(_p['calibrations']['alignment']['locations'], list))
    assert(len(_p['calibrations']['alignment']['locations']) == 1)
    assert (_p['calibrations']['alignment']['locations'][0] == 0.5)

    # Test with a multi-element list
    p['calibrations']['alignment']['locations'] = [0.0, 1.0]
    _p = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=p.to_config())  # Once as tuple
    assert(isinstance(_p['calibrations']['alignment']['locations'], list))
    assert(len(_p['calibrations']['alignment']['locations']) == 2)
    assert (_p['calibrations']['alignment']['locations'][0] == 0.0)
    assert (_p['calibrations']['alignment']['locations'][1] == 1.0)

    # Test something that should fail
    with pytest.raises(TypeError):
        p['calibrations']['alignment']['locations'] = 0.0
        _p = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=p.to_config())  # Once as tuple


# def test_framegrp_to_hdr():
#     biaspar = pypeitpar.FrameGroupPar(frametype='bias')
# 
#     embed()
#     exit()
# 
# 
# test_framegrp_to_hdr()

#def test_new_parset():
#    #p = pypeitpar.NewScatteredLightPar()
#
#    p = pypeitpar.NewProcessImagesPar()
#    embed()
#    exit()

# test_new_parset()


def compare_old_new(opar, npar):
    # opar = ocls()
    # npar = ncls()

    okeys = sorted(opar.keys())
    nkeys = sorted(npar.keys())

    comp = f'{opar.__class__.__name__}/{npar.__class__.__name__}'

    # Check that the keys are the same
    diff_keys = set(okeys) - set(nkeys)
    if len(diff_keys) != 0:
        print(f'{comp}: Keys are different: {diff_keys}')
#    assert len(diff_keys) == 0, f'{comp}: Keys are different: {diff_keys}'

    common_keys = set(okeys) & set(nkeys)

    # Check that the types are the same
    odtype = {key:opar.dtype[key] for key in common_keys}
    ndtype = {key:npar.parameters[key]['dtype'] for key in common_keys}
    if odtype != ndtype:
        print(f'{comp}: dtype dictionaries are different')
        diff_keys = [key for key in common_keys if odtype[key] != ndtype[key]]
        print(f'Different dtypes are for : {diff_keys}')
#    assert opar.dtype == ndtype, f'{comp}: dtype dictionaries are different'

    # Check that the defaults are the same
    # Ignore defaults that are parsets
    ndefault = {key:npar.parameters[key]['default'] for key in common_keys
                if not isinstance(npar.parameters[key]['default'], newparset.ParSet)}
    odefault = {key:opar.default[key] for key in common_keys
                if not isinstance(opar.default[key], parset.ParSet)}
    diff_keys = set(odefault.keys()) - set(ndefault.keys())
    assert len(diff_keys) == 0, \
            f'{comp}: Keys after excluding ParSet defaults are different: {diff_keys}'
    if odefault != ndefault:
        print(f'{comp}: default dictionaries are different')
        diff_keys = [key for key in common_keys if odefault[key] != ndefault[key]]
        print(f'Different default are for : {diff_keys}')
#    assert opar.default == ndefault, f'{comp}: default dictionaries are different'

    # Check that the options are the same
    noptions = {key:npar.parameters[key]['options'] for key in common_keys}
    ooptions = {key:opar.options[key] for key in common_keys}
    if ooptions != noptions:
        print(f'{comp}: option dictionaries are different')
        diff_keys = [key for key in common_keys if ooptions[key] != noptions[key]]
        print(f'Different options are for : {diff_keys}')
#    assert opar.options == noptions, f'{comp}: option dictionaries are different'

    # Check that the descriptions are the same
    ndescr = {key:npar.parameters[key]['descr'] for key in common_keys}
    odescr = {key:opar.descr[key] for key in common_keys}
    if odescr != ndescr:
        print(f'{comp}: description dictionaries are different')
        diff_keys = [key for key in common_keys if odescr[key] != ndescr[key]]
        print(f'Different descriptions are for : {diff_keys}')
#    assert opar.descr == ndescr, f'{comp} descr dictionaries are different'

def test_old_new_parsets():
    print('    ')
    print('    ')
    compare_old_new(pypeitpar.ScatteredLightPar(), newpypeitpar.ScatteredLightPar())
    print('    ')
    print('    ')
    compare_old_new(pypeitpar.ProcessImagesPar(), newpypeitpar.ProcessImagesPar())
    print('    ')
    print('    ')
    compare_old_new(pypeitpar.FrameGroupPar(), newpypeitpar.FrameGroupPar())
    print('    ')
    print('    ')
    par = pypeitpar.PypeItPar()
    compare_old_new(par['calibrations']['biasframe'], newpypeitpar.BiasFramePar())
    print('    ')
    print('    ')
    compare_old_new(par['calibrations']['darkframe'], newpypeitpar.DarkFramePar())
    print('    ')
    print('    ')
    compare_old_new(par['calibrations']['scattlightframe'], newpypeitpar.ScatteredLightFramePar())
    print('    ')
    print('    ')
    compare_old_new(par['calibrations']['pixelflatframe'], newpypeitpar.PixelFlatFramePar())
    print('    ')
    print('    ')
    compare_old_new(par['calibrations']['illumflatframe'], newpypeitpar.IllumFlatFramePar())
    print('    ')
    print('    ')
    compare_old_new(par['calibrations']['lampoffflatsframe'], newpypeitpar.LampOffFlatsFramePar())
    print('    ')
    print('    ')
    compare_old_new(par['calibrations']['slitless_pixflatframe'], newpypeitpar.SlitlessPixFlatFramePar())
    print('    ')
    print('    ')
    compare_old_new(par['calibrations']['pinholeframe'], newpypeitpar.PinholeFramePar())
    print('    ')
    print('    ')
    compare_old_new(par['calibrations']['alignframe'], newpypeitpar.AlignFramePar())
    print('    ')
    print('    ')
    compare_old_new(par['calibrations']['arcframe'], newpypeitpar.ArcFramePar())
    print('    ')
    print('    ')
    compare_old_new(par['calibrations']['tiltframe'], newpypeitpar.TiltFramePar())
    print('    ')
    print('    ')
    compare_old_new(par['calibrations']['traceframe'], newpypeitpar.TraceFramePar())
    print('    ')
    print('    ')
    compare_old_new(par['calibrations']['standardframe'], newpypeitpar.StandardFramePar())
    print('    ')
    print('    ')
    compare_old_new(par['calibrations']['skyframe'], newpypeitpar.SkyFramePar())
    print('    ')
    print('    ')
    compare_old_new(par['scienceframe'], newpypeitpar.ScienceFramePar())
    print('    ')
    print('    ')
    compare_old_new(pypeitpar.FlatFieldPar(), newpypeitpar.FlatFieldPar())
    print('    ')
    print('    ')
    compare_old_new(pypeitpar.FlexurePar(), newpypeitpar.FlexurePar())
    print('    ')
    print('    ')
    compare_old_new(pypeitpar.AlignPar(), newpypeitpar.AlignPar())
    print('    ')
    print('    ')
    compare_old_new(pypeitpar.Coadd1DPar(), newpypeitpar.Coadd1DPar())
    print('    ')
    print('    ')
    compare_old_new(pypeitpar.Coadd2DPar(), newpypeitpar.Coadd2DPar())
    print('    ')
    print('    ')
    compare_old_new(pypeitpar.CubePar(), newpypeitpar.CubePar())
    print('    ')
    print('    ')
    compare_old_new(pypeitpar.FluxCalibratePar(), newpypeitpar.FluxCalibratePar())
    print('    ')
    print('    ')
    compare_old_new(pypeitpar.SensfuncUVISPar(), newpypeitpar.SensfuncUVISPar())
    print('    ')
    print('    ')
    compare_old_new(pypeitpar.TelluricPar(), newpypeitpar.TelluricPar())
    print('    ')
    print('    ')
    compare_old_new(pypeitpar.SensFuncPar(), newpypeitpar.SensFuncPar())
    print('    ')
    print('    ')
    compare_old_new(pypeitpar.SlitMaskPar(), newpypeitpar.SlitMaskPar())
    print('    ')
    print('    ')
    compare_old_new(pypeitpar.ReduxPar(), newpypeitpar.ReduxPar())
    print('    ')
    print('    ')
    compare_old_new(pypeitpar.WavelengthSolutionPar(), newpypeitpar.WavelengthSolutionPar())
    print('    ')
    print('    ')
    compare_old_new(pypeitpar.EdgeTracePar(), newpypeitpar.EdgeTracePar())
    print('    ')
    print('    ')
    compare_old_new(pypeitpar.WaveTiltsPar(), newpypeitpar.WaveTiltsPar())
    print('    ')
    print('    ')
    compare_old_new(pypeitpar.FindObjPar(), newpypeitpar.FindObjPar())
    print('    ')
    print('    ')
    compare_old_new(pypeitpar.SkySubPar(), newpypeitpar.SkySubPar())
    print('    ')
    print('    ')
    compare_old_new(pypeitpar.ExtractionPar(), newpypeitpar.ExtractionPar())
    print('    ')
    print('    ')
    compare_old_new(pypeitpar.Collate1DPar(), newpypeitpar.Collate1DPar())
    print('    ')
    print('    ')
    compare_old_new(pypeitpar.TelescopePar(), newpypeitpar.TelescopePar())
    print('    ')
    print('    ')
    compare_old_new(pypeitpar.ReducePar(), newpypeitpar.ReducePar())
    print('    ')
    print('    ')
    compare_old_new(pypeitpar.CalibrationsPar(), newpypeitpar.CalibrationsPar())
    print('    ')
    print('    ')
    compare_old_new(pypeitpar.PypeItPar(), newpypeitpar.PypeItPar())
    print('    ')
    print('    ')

#    embed()
#    exit()

test_old_new_parsets()

