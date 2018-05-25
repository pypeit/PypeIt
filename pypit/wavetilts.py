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
                               'order': 2,
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
        self.all_trcdict = {}
        self.polytilts = None

        # MasterFrame
        masterframe.MasterFrame.__init__(self, self.frametype, setup, self.settings)



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
        self.badlines, self.maskrows, self.tcoeff, self.all_tilts, self.model2, = artracewave.analyze_spec_lines(
            self.msarc, slit, self.all_trcdict[str(slit)], self.pixcen, self.settings)
        if self.badlines > 0:
            msgs.warn("There were {0:d} additional arc lines that should have been traced".format(self.badlines) +
                      msgs.newline() + "(perhaps lines were saturated?). Check the spectral tilt solution")
        return self.badlines

    def _extract_arcs(self):
        # Extract an arc down each slit/order
        self.arccen, self.arc_maskslit, _ = ararc.get_censpec(self.lordloc, self.rordloc,
                                                              self.pixlocn, self.msarc, self.det,
                                                              gen_satmask=False)
        self.satmask = np.zeros_like(self.msarc)
        return self.arccen, self.arc_maskslit

    def _prepare_polytilts(self, skip_QA=False, show_QA=False):
        reload(artracewave)
        self.polytilts, self.outpar = artracewave.prepare_polytilts(
            self.msarc, self.maskrows, self.tcoeff, self.all_tilts, self.settings,
            setup=self.setup, skip_QA=skip_QA, show_QA=show_QA)
        #
        return self.polytilts

    def _trace_tilts(self, slit):
        #reload(artracewave)
        # Determine the tilts for this slit
        trcdict = artracewave.trace_tilt(self.pixcen, self.rordloc, self.lordloc, self.det,
                                         self.msarc, slit, settings_spect, self.settings,
                                         censpec=self.arccen[:, slit], nsmth=3,
                                              trthrsh=self.settings['tilts']['trthrsh'])
        # Load up
        self.all_trcdict[str(slit)] = trcdict.copy()
        # Return
        return trcdict

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

    def _qa(self, slit):
        """
        QA
          Wrapper to artraceslits.slit_trace_qa()

        Returns
        -------

        """
        reload(artracewave)
        self.tiltsplot, self.ztilto, self.xdat = artracewave.prep_tilts_qa(
            self.msarc, self.all_tilts, self.tilts, self.all_trcdict[str(slit)]['arcdet'],
            self.pixcen, slit)

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
            debugger.chk_arc_tilts(self.msarc, self.all_trcdict[str(slit)],
                                   sedges=(self.lordloc[:,slit], self.rordloc[:,slit]))
            msgs.info("Green = ok line;  red=rejected")
        elif attr == 'model2':
            tmp = self.all_trcdict[str(slit)].copy()
            tmp['xtfit'] = [np.arange(self.model2.shape[0])]*self.model2.shape[1]
            tmp['ytfit'] = [self.model2[:,ii] for ii in range(self.model2.shape[1])]
            debugger.chk_arc_tilts(self.msarc, tmp, sedges=(self.lordloc[:,slit], self.rordloc[:,slit]))
        elif attr == 'polytilt_img':
            if self.polytilts is not None:
                debugger.show_image(self.polytilts)
        elif attr == ['tilts','polytilts']:
            if slit is None:
                msgs.warn("Need to input a slit")
                return
            tmp = self.all_trcdict[str(slit)].copy()
            if attr == 'tilts':
                tilts = self.tilts
            else:
                tilts = self.polytilts
            self.tiltsplot, self.ztilto, self.xdat = artracewave.prep_tilts_qa(
                self.msarc, self.all_tilts, tilts, tmp['arcdet'],
                self.pixcen, slit)
            tmp['xtfit'] = [self.xdat[:,ii] for ii in range(self.xdat.shape[1])]
            tmp['ytfit'] = [self.ztilto[:,ii] for ii in range(self.ztilto.shape[1])]
            debugger.chk_arc_tilts(self.msarc, tmp,
                                   sedges=(self.lordloc[:,slit], self.rordloc[:,slit]))


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

