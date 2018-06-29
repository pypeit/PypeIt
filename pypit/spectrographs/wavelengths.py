""" Spectrograph specific info/methods for wavelengths
"""
import numpy as np

from astropy.io import fits

from pypit import msgs
from pypit.spectrographs import deimos
from pypit.spectrographs import lris
from pypit import arparse

from pypit import ardebug as debugger


def setup_param(spectrograph, msarc_shape, fitstbl, arc_idx,
                calibrate_lamps=None):
    """ Setup for arc analysis

    Parameters
    ----------
    spectrograph : str
    msarc_shape : tuple
    fitstbl : Table
      Contains relevant information from fits header files
    arc_idx : int
      Index of the arc frame in the fitstbl
    calibrate_lamps : str, optional
       List of lamps used

    """
    # Defaults
    arcparam = dict(llist='',
                    disp=0.,             # Ang/unbinned pixel
                    b1=0.,               # Pixel fit term (binning independent)
                    b2=0.,               # Pixel fit term
                    lamps=[],            # Line lamps on
                    wv_cen=0.,           # Estimate of central wavelength
                    wvmnx=[2900.,12000.],# Guess at wavelength range
                    disp_toler=0.1,      # 10% tolerance
                    match_toler=3.,      # Matching tolerance (pixels)
                    min_ampl=300.,       # Minimum amplitude
                    func='legendre',     # Function for fitting
                    n_first=1,           # Order of polynomial for first fit
                    n_final=4,           # Order of polynomial for final fit
                    nsig_rej=2.,         # Number of sigma for rejection
                    nsig_rej_final=3.0,  # Number of sigma for rejection (final fit)
                    Nstrong=13)          # Number of lines for auto-analysis

    modify_dict = None
    # Instrument/disperser specific
    disperser = fitstbl["dispname"][arc_idx]
    binspatial, binspectral = arparse.parse_binning(fitstbl['binning'][arc_idx])
    if spectrograph == 'shane_kast_blue':
        # Could have the following depend on lamps that were turned on
        lamps = ['CdI','HgI','HeI']
        if disperser == '600/4310':
            arcparam['disp']=1.02
            arcparam['b1']=6.88935788e-04
            arcparam['b2']=-2.38634231e-08
            arcparam['wvmnx'][1] = 6000.
            arcparam['wv_cen'] = 4250.
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))
    elif spectrograph=='shane_kast_red':
        lamps = ['NeI','HgI','HeI','ArI']
        if disperser == '600/7500':
            arcparam['disp']=1.30
            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0] / binspectral
            arcparam['wvmnx'][0] = 5000.
            arcparam['n_first']=2 # Should be able to lock on
        elif disperser == '1200/5000':
            arcparam['disp']=0.63
            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0] / binspectral
            arcparam['wvmnx'][0] = 5000.
            arcparam['n_first']=2 # Should be able to lock on
            arcparam['wv_cen'] = 6600.
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))
    elif spectrograph=='shane_kast_red_ret':
        lamps = ['NeI','HgI','HeI','ArI']
        if disperser == '600/7500':
            arcparam['disp']=2.35
            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0] / binspectral
            arcparam['wvmnx'][0] = 5000.
            arcparam['n_first']=2 # Should be able to lock on
        elif disperser == '1200/5000':
            arcparam['disp']=1.17
            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0] / binspectral
            arcparam['wvmnx'][0] = 5000.
            arcparam['n_first']=2 # Should be able to lock on
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))
    elif spectrograph=='keck_lris_blue':
        lamps = ['NeI', 'ArI', 'CdI', 'KrI', 'XeI', 'ZnI','CdI','HgI']
        if disperser == '600/4000':
            arcparam['n_first']=2 # Too much curvature for 1st order
            arcparam['disp']=0.63 # Ang per pixel (unbinned)
            arcparam['b1']= 4.54698031e-04
            arcparam['b2']= -6.86414978e-09
            arcparam['wvmnx'][1] = 6000.
            arcparam['wv_cen'] = 4000.
        elif disperser == '400/3400':
            arcparam['n_first']=2 # Too much curvature for 1st order
            arcparam['disp']=1.02
            arcparam['b1']= 2.72694493e-04
            arcparam['b2']= -5.30717321e-09
            arcparam['wvmnx'][1] = 6000.
        elif disperser == '300/5000':
            arcparam['n_first'] = 2
            arcparam['wv_cen'] = 4500.
            arcparam['disp'] = 1.43
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))
    elif spectrograph=='keck_lris_red':
        arcparam['wv_cen'] = fitstbl['wavecen'][arc_idx]
        lamps = ['ArI','NeI','HgI','KrI','XeI']  # Should set according to the lamps that were on
        if disperser == '600/7500':
            arcparam['n_first']=3 # Too much curvature for 1st order
            arcparam['disp']=0.80 # Ang per pixel (unbinned)
            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0] / binspectral
            arcparam['wvmnx'][1] = 11000.
        elif disperser == '600/10000':
            arcparam['n_first']=2 # Too much curvature for 1st order
            arcparam['disp']=0.80 # Ang per pixel (unbinned)
            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0] / binspectral
            arcparam['wvmnx'][1] = 12000.
        elif disperser == '400/8500':
            arcparam['n_first']=2 # Too much curvature for 1st order
            arcparam['disp']=1.19 # Ang per pixel (unbinned)
            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0] / binspectral
            arcparam['wvmnx'][1] = 11000.
            arcparam['min_ampl'] = 3000.  # Lines tend to be very strong
            arcparam['nsig_rej_final'] = 5.
        elif disperser == '900/5500':
            arcparam['n_first']=2 # Too much curvature for 1st order
            arcparam['disp']=0.53 # Ang per pixel (unbinned)
            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0] / binspectral
            arcparam['wvmnx'][1] = 7000.
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))
    elif spectrograph=='keck_deimos':
        arcparam['wv_cen'] = fitstbl['dispangle'][arc_idx]
        # TODO -- Should set according to the lamps that were on
        lamps = ['ArI','NeI','KrI','XeI']
        if disperser == '830G': # Blaze 8640
            arcparam['n_first']=2 # Too much curvature for 1st order
            arcparam['disp']=0.47 # Ang per pixel (unbinned)
            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0]
            arcparam['wvmnx'][0] = 550.
            arcparam['wvmnx'][1] = 11000.
            arcparam['min_ampl'] = 3000.  # Lines tend to be very strong
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))
    elif spectrograph=='wht_isis_blue':
        modify_dict = dict(NeI={'min_wave': 3000.,'min_intensity': 299,
                                'min_Aki': 0.},ArI={'min_intensity': 399.})
        lamps=['CuI','NeI','ArI']
        if fitstbl["dichroic"][arc_idx].strip() == '5300':
            arcparam['wvmnx'][1] = 6000.
        else:
            msgs.error('Not ready for this dichroic {:s}!'.format(disperser))
        if disperser == 'R300B':
            arcparam['n_first']=1  #
            arcparam['disp']=0.80  # Ang per pixel (unbinned)
            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0]
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))
    elif spectrograph == 'tng_dolores':
        lamps = ['NeI', 'HgI']
        if disperser == 'LR-R':
            arcparam['n_first'] = 2  # Too much curvature for 1st order
            arcparam['disp'] = 2.61  # Ang per pixel (unbinned)
            arcparam['disp_toler'] = 0.1  # Ang per pixel (unbinned)
            arcparam['wvmnx'][0] = 4470.0
            arcparam['wvmnx'][1] = 10073.0
            arcparam['wv_cen'] = 7400.
            arcparam['b1'] = 1. / arcparam['disp'] / msarc_shape[0] / binspectral
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))
    else:
        msgs.error('ararc.setup_param: Not ready for this instrument {:s}!'.format(spectrograph))
    # Load linelist
    #if settings.argflag['arc']['calibrate']['lamps'] is not None:
    if calibrate_lamps is not None:
        arcparam['lamps'] = calibrate_lamps
    else:
        arcparam['lamps'] = lamps
    slmps = lamps[0]
    for lamp in lamps[1:]:
        slmps=slmps+','+lamp
    msgs.info('Loading line list using {:s} lamps'.format(slmps))
    #    arcparam['llist'] = ararclines.load_arcline_list(slf, idx, lamps, disperser,
    arcparam['llist'] = ararclines.load_arcline_list(lamps, disperser, spectrograph,
                                                     wvmnx=arcparam['wvmnx'],
                                                     modify_parse_dict=modify_dict)
    # Binning
    arcparam['disp'] *= binspectral

    # Return
    return arcparam
