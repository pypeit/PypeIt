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
                               'yorder': 4,
                               'func2D': 'legendre',  # Function for 2D fit
                               'method': 'spca',  # defunct
                               'params': [1,1,0],  # defunct
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
                 lordloc=None, rordloc=None, pixlocn=None, pixcen=None, slitpix=None):

        # Required parameters (but can be None)
        self.msarc = msarc
        self.lordloc = lordloc
        self.rordloc = rordloc
        self.pixlocn = pixlocn
        self.pixcen = pixcen
        self.slitpix = slitpix

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
        self.steps = []

        # Main outputs
        self.final_tilts = None

        # Key Internals
        self.all_trcdict = [None]*self.nslit
        self.polytilts = None
        self.tilts = None

        # MasterFrame
        masterframe.MasterFrame.__init__(self, self.frametype, setup, self.settings)


    def _analyze_lines(self, slit):
        reload(artracewave)
        self.badlines, self.all_tilts = artracewave.analyze_lines(
            self.msarc, self.all_trcdict[slit], slit, self.pixcen, self.settings)
        if self.badlines > 0:
            msgs.warn("There were {0:d} additional arc lines that should have been traced".format(self.badlines) +
                      msgs.newline() + "(perhaps lines were saturated?). Check the spectral tilt solution")
        # Step
        self.steps.append(inspect.stack()[0][3])
        return self.badlines

    def _extract_arcs(self):
        # Extract an arc down each slit/order
        self.arccen, self.arc_maskslit, _ = ararc.get_censpec(self.lordloc, self.rordloc,
                                                              self.pixlocn, self.msarc, self.det,
                                                              gen_satmask=False)
        self.satmask = np.zeros_like(self.msarc)
        # Step
        self.steps.append(inspect.stack()[0][3])
        return self.arccen, self.arc_maskslit

    def _fit_tilts(self, slit, show_QA=False):
        reload(artracewave)
        reload(arutils)
        self.tilts, self.outpar = artracewave.fit_tilts(self.msarc, slit, self.all_tilts,
                                                        self.settings, setup=self.setup, show_QA=show_QA)
        # Step
        self.steps.append(inspect.stack()[0][3])
        return self.tilts

    '''
    def _tilts_jxp_spca(self, slit):
        reload(artracewave)
        self.tilts = artracewave.tilts_spca(self.msarc, self.polytilts, self.pixcen, slit,
                                            self.all_trcdict[slit]['arcdet'], self.all_trcdict[slit]['aduse'],
                                            self.rordloc, self.lordloc, self.all_tilts)

    def _tilts_spca(self, slit):
        reload(artracewave)
        self.tilts = artracewave.tilts_spca(self.msarc, self.polytilts, self.pixcen, slit,
            self.all_trcdict[slit]['arcdet'], self.all_trcdict[slit]['aduse'],
            self.rordloc, self.lordloc)

    def _tilts_jxp_spline(self, slit):
        reload(artracewave)
        self.tilts = artracewave.tilts_spline(self.all_tilts,
                                              self.all_trcdict[slit]['arcdet'], self.all_trcdict[slit]['aduse'],
                                              self.polytilts, self.msarc, use_mtilt=True)
    def _tilts_spline(self, slit):
        reload(artracewave)
        self.tilts = artracewave.tilts_spline(self.all_tilts,
                                              self.all_trcdict[slit]['arcdet'], self.all_trcdict[slit]['aduse'],
                                              self.polytilts, self.msarc)
    '''

    def _trace_tilts(self, slit):
        reload(artracewave)
        # Determine the tilts for this slit
        trcdict = artracewave.trace_tilt(self.pixcen, self.rordloc, self.lordloc, self.det,
                                         self.msarc, slit, settings_spect, self.settings,
                                         censpec=self.arccen[:, slit], nsmth=3,
                                              trthrsh=self.settings['tilts']['trthrsh'])
        # Load up
        self.all_trcdict[slit] = trcdict.copy()
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return trcdict

    def run(self, maskslits=None):
        """ Main driver for tracing arc lines

        Code flow:
           1.  Extract an arc spectrum down the center of each slit/order
           2.  Loop on slits/orders
             i.   Trace the arc lines (fweight is the default)
             ii.  Fit the individual arc lines
             iii.  2D Fit to the offset from pixcen
             iv. Save

        Parameters
        ----------

        Returns
        -------
        """
        # If the user sets no tilts, return here
        if self.settings['tilts']['method'].lower() == "zero":
            # Assuming there is no spectral tilt
            self.final_tilts = np.outer(np.linspace(0.0, 1.0, self.msarc.shape[0]), np.ones(self.msarc.shape[1]))
            return self.final_tilts, None, None

        if maskslits is None:
            maskslits = np.zeros(self.nslit, dtype=bool)

        # Extract the arc spectra for all slits
        self.arccen, self.arc_maskslit = self._extract_arcs()

        # maskslit
        mask = maskslits & (self.arc_maskslit==1)
        gdslits = np.where(mask == 0)[0]

        # Final tilts image
        self.final_tilts = np.zeros_like(self.msarc)
        # Loop on all slits
        for slit in gdslits:
            # Trace
            _ = self._trace_tilts(slit)

            # Model line-by-line
            _ = self._analyze_lines(slit)

            # 2D model of the tilts
            #   Includes QA
            self.tilts = self._fit_tilts(slit)

            # Save to final image
            word = self.slitpix == slit+1
            self.final_tilts[word] = self.tilts[word]

        return self.final_tilts, maskslits

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
            ginga.chk_arc_tilts(self.msarc, tmp, sedges=(self.lordloc[:,slit], self.rordloc[:,slit]), all_green=True)
        elif attr == 'polytilt_img':
            if self.polytilts is not None:
                ginga.show_image(self.polytilts)
        elif attr in ['polytilts', 'tilts']:
            tmp = self.all_trcdict[slit].copy()
            tmp['xtfit'] = []
            tmp['ytfit'] = []

            # arcdet is only the approximately nearest pixel (not even necessarily)
            for idx in np.where(self.all_trcdict[slit]['aduse'])[0]:
                tmp['xtfit'].append(np.arange(self.msarc.shape[1]))
                xgd = self.all_trcdict[slit]['xtfit'][idx][self.all_trcdict[slit]['xtfit'][idx].size//2]
                ycen = self.all_tilts[1][int(xgd),idx]
                # This is not exact.  Could make a correction.  Probably is close enough
                yval = ycen + self.tilts[int(ycen),:]
                tmp['ytfit'].append(yval)
            '''
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
            '''
            # Show
            ginga.chk_arc_tilts(self.msarc, tmp,
                                sedges=(self.lordloc[:,slit], self.rordloc[:,slit]))

    def save_master(self, outfile=None):
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

