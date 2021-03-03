
from IPython import embed

from pypeit import sensfunc_new

def test_telluricdata_init():

    td = telluric.TelluricData(norders=1, nspec=2048, n_obj_par=3)
    assert td.model['WAVE'].shape == (1,2048), 'bad wavelength array shape'
    assert td.model['TELL_THETA'].shape == (1,7), 'bad or changed telluric parameter shape'
    assert td.model['OBJ_THETA'].shape == (1,3), 'bad or changed object parameter shape'
    

def test_telluricdata_io():

    td = telluric.TelluricData(norders=1, nspec=2048, n_obj_par=3)

    

if __name__ == '__main__':
    test_telluricdata_io()

