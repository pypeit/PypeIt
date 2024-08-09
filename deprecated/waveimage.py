"""
Module for guiding construction of the Wavelength Image

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import inspect

import numpy as np
import os

from pypeit import msgs
from pypeit import utils
from pypeit import datamodel

from IPython import embed


class WaveImage(datamodel.DataContainer):
    version = '1.0.0'

    # I/O
    output_to_disk = None #('WVTILTS_IMAGE', 'WVTILTS_FULLMASK', 'WVTILTS_DETECTOR_CONTAINER')
    hdu_prefix = None

    # Master fun
    master_type = 'Wave'
    master_file_format = 'fits'

    datamodel = {
        'image':  dict(otype=np.ndarray, atype=np.floating, desc='2D Wavelength image'),
        'PYP_SPEC': dict(otype=str, desc='PypeIt spectrograph name'),
    }

    def __init__(self, image, PYP_SPEC=None):
        # Parse
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        d = dict([(k,values[k]) for k in args[1:]])
        # Setup the DataContainer
        datamodel.DataContainer.__init__(self, d=d)



class BuildWaveImage(object):
    """
    Class to generate the Wavelength Image

    Args:
        slits (:class:`pypeit.edgetrace.SlitTraceSet`):
            Object holding the slit edge locations
        tilts (np.ndarray or None):
            Tilt image
        wv_calib (dict or None): wavelength solution dictionary
            Parameters are read from wv_calib['par']
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            The `Spectrograph` instance that sets the
            instrument used to take the observations.  Used to set
            :attr:`spectrograph`.
        det (int or None):

    Attributes:
        image (np.ndarray): Wavelength image
        steps (list): List of the processing steps performed

    """
    master_type = 'Wave'

#    @classmethod
#    def from_master_file(cls, master_file):
#        """
#
#        Args:
#            master_file (str):
#
#        Returns:
#            waveimage.WaveImage:
#
#        """
#        # Spectrograph
#        spectrograph, extras = masterframe.items_from_master_file(master_file)
#        head0 = extras[0]
#        # Master info
#        master_dir = head0['MSTRDIR']
#        master_key = head0['MSTRKEY']
#        # Instantiate
#        slf = cls(None, None, None, spectrograph, None, master_dir=master_dir,
#                  master_key=master_key, reuse_masters=True)
#        slf.image = slf.load(ifile=master_file)
#        # Return
#        return slf

    # TODO: Is maskslits ever anything besides slits.mask? (e.g., see calibrations.py call)
    def __init__(self, slits, tilts, wv_calib, spectrograph, det):

        # MasterFrame
        #masterframe.MasterFrame.__init__(self, self.master_type, master_dir=master_dir,
        #                                 master_key=master_key, reuse_masters=reuse_masters)
        # Required parameters
        self.spectrograph = spectrograph
        self.det = det

        # TODO: Do we need to assign slits to self?
        self.slits = slits
        self.tilts = tilts
        self.wv_calib = wv_calib
        if self.slits is None:
            self.slitmask = None
            self.slit_spat_pos = None
        else:
            # NOTE: This uses the pad defined by EdgeTraceSetPar
            self.slitmask = self.slits.slit_img()
            # This selects the coordinates for the tweaked edges if
            # they exist, original otherwise.
            self.slit_spat_pos = self.slits.spatial_coordinates()

        # For echelle order, primarily
        # TODO: only echelle is ever used.  Do we need to keep the whole
        # thing?
        self.par = wv_calib['par'] if wv_calib is not None else None

        # Main output
        self.image = None
        self.steps = []

    def build_wave(self):
        """
        Main algorithm to build the wavelength image

        Returns:
            `numpy.ndarray`_: The wavelength image.

        """
        # Loop on slits
        ok_slits = np.where(np.invert(self.slits.mask))[0]
        self.image = np.zeros_like(self.tilts)
        nspec = self.slitmask.shape[0]

        # Error checking on the wv_calib
        #if (nspec-1) != int(self.wv_calib[str(0)]['fmax']):
        #    msgs.error('Your wavelength fits used inconsistent normalization. Something is wrong!')

        # If this is echelle print out a status message and do some error checking
        if self.par['echelle']:
            msgs.info('Evaluating 2-d wavelength solution for echelle....')
            if len(self.wv_calib['fit2d']['orders']) != len(ok_slits):
                msgs.error('wv_calib and ok_slits do not line up. Something is very wrong!')

        # Unpack some 2-d fit parameters if this is echelle
        for slit in ok_slits:
            thismask = (self.slitmask == slit)
            if self.par['echelle']:
                # TODO: Put this in `SlitTraceSet`?
                order, indx = self.spectrograph.slit2order(self.slit_spat_pos[slit])
                # evaluate solution
                self.image[thismask] = utils.func_val(self.wv_calib['fit2d']['coeffs'],
                                                       self.tilts[thismask],
                                                       self.wv_calib['fit2d']['func2d'],
                                                       x2=np.ones_like(self.tilts[thismask])*order,
                                                       minx=self.wv_calib['fit2d']['min_spec'],
                                                       maxx=self.wv_calib['fit2d']['max_spec'],
                                                       minx2=self.wv_calib['fit2d']['min_order'],
                                                       maxx2=self.wv_calib['fit2d']['max_order'])
                self.image[thismask] /= order
            else:
                iwv_calib = self.wv_calib[str(slit)]
                self.image[thismask] = utils.func_val(iwv_calib['fitc'], self.tilts[thismask],
                                                       iwv_calib['function'],
                                                       minx=iwv_calib['fmin'],
                                                       maxx=iwv_calib['fmax'])

        # Return
        self.steps.append(inspect.stack()[0][3])
        return WaveImage(self.image, PYP_SPEC=self.spectrograph.spectrograph)

    def __repr__(self):
        # Generate sets string
        txt = '<{:s}: >'.format(self.__class__.__name__)
        return txt

