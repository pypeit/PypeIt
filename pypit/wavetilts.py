# Module for guiding Arc/Sky line tracing
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import os

from importlib import reload

from pypit import msgs
from pypit import ardebug as debugger
from pypit import ararc
from pypit.core import artracewave
from pypit import arutils
from pypit import masterframe
from pypit import ginga


# For out of PYPIT running
if msgs._debug is None:
    debug = debugger.init()
    debug['develop'] = True
    msgs.reset(debug=debug, verbosity=2)

frametype = 'tilts'

# Place these here or elsewhere?
#  Wherever they be, they need to be defined, described, etc.
default_settings = dict(tilts={'idsonly': False,
                               'trthrsh': 1000.,
                               'order': 1,
                               'function': 'legendre',
                               'method': 'spca',
                               'params': [1,1,0],
                               }
                        )
settings_spect = dict(det01={'saturation': 60000., 'nonlinear': 0.9})
#settings_spect[dnum]['saturation']*settings_spect[dnum]['nonlinear']

#  See save_master() for the data model for output


class WaveTilts(masterframe.MasterFrame):
    """Class to guide slit/order tracing

    Parameters
    ----------
    msarc : ndarray
      Arc image

    Attributes
    ----------
    frametype : str
      Hard-coded to 'tilts'

    """
    def __init__(self, msarc, settings=None, det=None, setup='',
                 lordloc=None, rordloc=None, pixlocn=None, pixcen=None):

        # Required parameters (but can be None)
        self.msarc = msarc
        self.lordloc = lordloc
        self.rordloc = rordloc
        self.pixlocn = pixlocn
        self.pixcen = pixcen

        # Optional parameters
        self.det = det
        if settings is None:
            self.settings = default_settings.copy()
        else:
            self.settings = settings
            if 'tilts' not in self.settings:
                self.settings.update(default_settings.copy())

        # Attributes
        self.frametype = frametype

        # Main outputs

        # Key Internals

        # MasterFrame
        masterframe.MasterFrame.__init__(self, self.frametype, setup, self.settings)


    def show(self, attr, slit=None, display='ginga'):
        """
        Display an image or spectrum in TraceSlits

        Parameters
        ----------
        attr : str
          'edges' -- Show the mstrace image and the edges
          'edgearr' -- Show the edgearr image
          'siglev' -- Show the Sobolev image
        display : str (optional)
          'ginga' -- Display to an RC Ginga
        """
        if attr == 'sslit':
            debugger.chk_arc_tilts(self.msarc, self.trcdict,
                                   sedges=(self.lordloc[:,slit], self.rordloc[:,slit]))
            msgs.info("Green = ok line;  red=rejected")

    def master(self):
        """ Mainly for PYPIT running

        Parameters
        ----------

        Returns
        -------
        loaded : bool

        """
        # Load master frame?
        loaded = False
        if self._masters_load_chk():
            loaded = self.load_master()
        # Return
        return loaded

    def _analyze_tilt_traces(self, slit):
        reload(artracewave)
        self.badlines, self.maskrows = artracewave.analyze_spec_lines(self.msarc, slit,
                                                            self.trcdict,
                                                            self.pixcen,
                                                            self.settings)
        return self.badlines, self.maskrows

    def _extract_arcs(self):
        # Extract an arc down each slit/order
        self.arccen, self.arc_maskslit, _ = ararc.get_censpec(self.lordloc, self.rordloc,
                                                              self.pixlocn, self.msarc, self.det,
                                                              gen_satmask=False)
        self.satmask = np.zeros_like(self.msarc)
        return self.arccen, self.arc_maskslit

    def _trace_tilts(self, slit):
        #reload(artracewave)
        # Determine the tilts for this slit
        self.trcdict = artracewave.trace_tilt(self.pixcen, self.rordloc, self.lordloc, self.det,
                                         self.msarc, slit, settings_spect, self.settings,
                                         censpec=self.arccen[:, slit], nsmth=3,
                                              trthrsh=self.settings['tilts']['trthrsh'])
        # Return
        return self.trcdict

    def run(self, maskslits):
        """ Main driver for tracing arc lines

        Parameters
        ----------

        Returns
        -------
        """
        #if settings.argflag['trace']['slits']['tilts']['method'].lower() == "zero":
        # If the user sets no tilts, return here
        if self.settings['method'].lower() == "zero":
            # Assuming there is no spectral tilt
            self.tilts = np.outer(np.linspace(0.0, 1.0, self.msarc.shape[0]), np.ones(self.msarc.shape[1]))
            return self.tilts, None, None

        # Extract the arc spectra for all slits
        self.arccen, self.arc_maskslit = self._extract_arc()

        # Setup
        fitxy = [self.settings['order'], 1]

        # maskslit
        if maskslits is not None:
            mask = maskslits & (self.arc_maskslit==1)
        else:
            mask = self.arc_maskslit
        # TODO -- NEED TO PASS THIS BACK!?
        #slf._maskslits[det-1] = mask
        gdslits = np.where(mask == 0)[0]

        # Loop on all slits
        for slit in gdslits:
            pass

        # Final tilts image
        final_tilts = np.zeros_like(self.msarc)

    def _qa(self, use_slitid=True):
        """
        QA
          Wrapper to artraceslits.slit_trace_qa()

        Returns
        -------

        """
        #artraceslits.slit_trace_qa(self.mstrace, self.lcen,
        #                           self.rcen, self.extrapord, self.setup,
        #                           desc="Trace of the slit edges D{:02d}".format(self.det),
        #                           use_slitid=use_slitid)


    def __repr__(self):
        # Generate sets string
        txt = '<{:s}: '.format(self.__class__.__name__)
        if len(self.steps) > 0:
            txt+= ' steps: ['
            for step in self.steps:
                txt += '{:s}, '.format(step)
            txt = txt[:-2]+']'  # Trim the trailing comma
        txt += '>'
        return txt

