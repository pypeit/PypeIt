# Module for guiding Arc/Sky line tracing
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import os

from importlib import reload

from astropy.io import fits

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
                               'function': 'legendre',       # Function for arc line fits
                               'yorder': 2,
                               'poly_2D': False,           # Use 2D polynomial for polytilts
                               'poly_2Dfunc': 'polynomial',  # Function for 2D fit
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
        if self.lordloc is not None:
            self.nslit = self.lordloc.shape[1]
        else:
            self.nslit = 0

        # Main outputs
        self.final_tilts = None

        # Key Internals
        self.all_trcdict = [None]*self.nslit
        self.polytilts = None
        self.tilts = None

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
        self.badlines, self.maskrows, self.tcoeff, self.all_tilts = artracewave.analyze_spec_lines(
            self.msarc, slit, self.all_trcdict[slit], self.pixcen, self.settings)
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
        reload(arutils)
        self.polytilts, self.outpar = artracewave.prepare_polytilts(
            self.msarc, self.maskrows, self.tcoeff, self.all_tilts, self.settings,
            setup=self.setup, skip_QA=skip_QA, show_QA=show_QA)
        #
        return self.polytilts

    def _tilts_spca(self, slit):
        reload(artracewave)
        self.tilts = artracewave.tilts_spca(self.msarc, self.polytilts, self.pixcen, slit,
            self.all_trcdict[slit]['arcdet'], self.all_trcdict[slit]['aduse'],
            self.rordloc, self.lordloc)

    def _tilts_spline(self, slit):
        reload(artracewave)
        self.tilts = artracewave.tilts_spline(self.all_tilts,
                                              self.all_trcdict[slit]['arcdet'], self.all_trcdict[slit]['aduse'],
                                              self.polytilts, self.msarc)

    def _trace_tilts(self, slit):
        #reload(artracewave)
        # Determine the tilts for this slit
        trcdict = artracewave.trace_tilt(self.pixcen, self.rordloc, self.lordloc, self.det,
                                         self.msarc, slit, settings_spect, self.settings,
                                         censpec=self.arccen[:, slit], nsmth=3,
                                              trthrsh=self.settings['tilts']['trthrsh'])
        # Load up
        self.all_trcdict[slit] = trcdict.copy()
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
        self.arccen, self.arc_maskslit = self._extract_arcs()

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
            self.msarc, self.all_tilts, self.tilts, self.all_trcdict[slit]['arcdet'],
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
        reload(ginga)
        if attr == 'sslit':
            ginga.chk_arc_tilts(self.msarc, self.all_trcdict[slit],
                                   sedges=(self.lordloc[:,slit], self.rordloc[:,slit]))
            msgs.info("Green = ok line;  red=rejected")
        elif attr == 'model':
            tmp = self.all_trcdict[slit].copy()
            tmp['xtfit'] = self.all_trcdict[slit]['xmodel']
            tmp['ytfit'] = self.all_trcdict[slit]['ymodel']
            ginga.chk_arc_tilts(self.msarc, tmp, sedges=(self.lordloc[:,slit], self.rordloc[:,slit]))
            msgs.info("Ignore the color scheme")
        elif attr == 'polytilt_img':
            if self.polytilts is not None:
                ginga.show_image(self.polytilts)
        elif attr in ['polytilts', 'tilts']:
            tmp = self.all_trcdict[slit].copy()
            tmp['xtfit'] = []
            tmp['ytfit'] = []
            # arcdet is only the approximately nearest pixel (not even necessarily)
            for arcdet in self.all_trcdict[slit]['arcdet']:
                tmp['xtfit'].append(np.arange(self.msarc.shape[1]))
                if attr == 'tilts':  # Need to offset here as we finally fit to the lines not the row
                    yval = (2*self.tilts[arcdet, self.pixcen[arcdet, slit]]-self.tilts[arcdet,:])*(self.msarc.shape[0]-1)
                    # Offset onto the centroid arc line
                    yval += self.polytilts[arcdet, self.pixcen[arcdet, slit]]*(self.msarc.shape[0]-1) - arcdet
                else:
                    yval = self.polytilts[arcdet,:] * (self.msarc.shape[0]-1)
                # Save
                tmp['ytfit'].append(yval)
            # Show
            ginga.chk_arc_tilts(self.msarc, tmp,
                                sedges=(self.lordloc[:,slit], self.rordloc[:,slit]))

    def save_master(self, outfile=None, use_tilts_as_final=False):
        """

        Parameters
        ----------
        outfile
        use_tilts_as_final

        Returns
        -------

        """
        if outfile is None:
            outfile = self.ms_name
        #
        if use_tilts_as_final:
            msgs.warn("Using tilts as final.  Better know what you are doing!")
            self.final_tilts = self.tilts
        #
        if self.final_tilts is None:
            msgs.warn("final_tilts not yet created.  Make it!")
            return
        #
        hdu0 = fits.PrimaryHDU(self.final_tilts)
        hdul = [hdu0]

        for slit in range(self.nslit):
            # fweight and model
            xtfits = self.all_trcdict[slit]['xtfit']  # For convenience
            xszs = [len(xtfit) for xtfit in xtfits]
            maxx = np.max(xszs)
            fwm_img = np.zeros((maxx, len(xtfits), 4))
            # Fill fweight and model
            model_cnt = 0
            for kk, xtfit in enumerate(xtfits):
                fwm_img[0:xszs[kk], kk, 0] = xtfit
                fwm_img[0:xszs[kk], kk, 1] = self.all_trcdict[slit]['ytfit'][kk]
                #
                if self.all_trcdict[slit]['aduse'][kk]:
                    fwm_img[0:xszs[kk], kk, 2] = self.all_trcdict[slit]['xmodel'][model_cnt]
                    fwm_img[0:xszs[kk], kk, 3] = self.all_trcdict[slit]['ymodel'][model_cnt]
                    model_cnt += 1
            hdu1 = fits.ImageHDU(fwm_img)
            hdu1.name = 'FWM{:03d}'.format(slit)
            hdul.append(hdu1)
        # Finish
        hdulist = fits.HDUList(hdul)
        hdulist.writeto(outfile, clobber=True)


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

