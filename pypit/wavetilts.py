# Module for guiding Arc/Sky line tracing
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import os

#from importlib import reload

from astropy.io import fits

from pypit import msgs
from pypit import ardebug as debugger
from pypit.core import ararc
from pypit.core import artracewave
from pypit.core import armasters
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
def default_settings():
    default_settings = dict(tilts={'idsonly': False,
                               'tracethresh': 1000.,
                               'order': 2,
                               'function': 'legendre',       # Function for arc line fits
                               'yorder': 4,
                               'func2D': 'legendre',  # Function for 2D fit
                               'method': 'spca',  # defunct
                               'params': [1,1,0],  # defunct
                               }
                        )
    return default_settings
#  See save_master() for the data model for output


class WaveTilts(masterframe.MasterFrame):
    """Class to guide slit/order tracing

    Parameters
    ----------
    msarc : ndarray
      Arc image
    tslits_dict : dict
      Input from TraceSlits
    settings_det : dict
      Detector settings -- Needed for arc line saturation
    det : int
      Detector index
    settings : dict
      Tilts settings

    Attributes
    ----------
    frametype : str
      Hard-coded to 'tilts'
    steps : list
    mask : ndarray, bool
      True = Ignore this slit
    all_trcdict : list of dict
      All trace dict's
    tilts : ndarray
      Tilts for a single slit/order
    all_ttilts : list of tuples
      Tuple of tilts ndarray's
    final_tilts : ndarray
      Final tilts image

    """
    def __init__(self, msarc, settings=None, det=None, settings_det=None, setup='',
                 tslits_dict=None, pixlocn=None):

        # Parameters (but can be None)
        self.msarc = msarc
        self.tslits_dict = tslits_dict
        self.pixlocn = pixlocn
        self.settings_det = settings_det

        # Optional parameters
        self.det = det
        if settings is None:
            self.settings = default_settings()
        else:
            self.settings = settings
            if 'tilts' not in self.settings:
                self.settings.update(default_settings())

        # Attributes
        self.frametype = frametype
        if self.tslits_dict is not None:
            self.nslit = self.tslits_dict['lcen'].shape[1]
        else:
            self.nslit = 0
        self.steps = []

        # Main outputs
        self.final_tilts = None

        # Key Internals
        self.mask = None
        self.all_trcdict = [None]*self.nslit
        self.tilts = None
        self.all_ttilts = [None]*self.nslit
        self.final_tilts = None

        # MasterFrame
        masterframe.MasterFrame.__init__(self, self.frametype, setup, self.settings)

    @classmethod
    def from_master_files(cls, setup, mdir='./'):
        """
        Build the class from Master frames

        Parameters
        ----------
        setup : str
        mdir : str, optional

        Returns
        -------
        slf

        """
        # Arc
        msarc_file = armasters.core_master_name('arc', setup, mdir)
        msarc, _, _ = armasters._load(msarc_file)

        # Instantiate
        slf = cls(msarc)

        # Tilts
        mstilts_file = armasters.core_master_name('tilts', setup, mdir)
        hdul = fits.open(mstilts_file)
        slf.final_tilts = hdul[0].data
        slf.tilts = slf.final_tilts

        # Dict
        slf.all_trcdict = []
        islit = 0
        for hdu in hdul[1:]:
            if hdu.name == 'FWM{:03d}'.format(islit):
                # Setup
                fwm_img = hdu.data
                narc = fwm_img.shape[1]
                trcdict = dict(xtfit=[], ytfit=[], xmodel=[], ymodel=[], ycen=[], aduse=np.zeros(narc, dtype=bool))
                # Fill  (the -1 are for ycen which is packed in at the end)
                for iarc in range(narc):
                    trcdict['xtfit'].append(fwm_img[:-1,iarc,0])
                    trcdict['ytfit'].append(fwm_img[:-1,iarc,1])
                    trcdict['ycen'].append(fwm_img[-1,iarc,1])  # Many of these are junk
                    if np.any(fwm_img[:-1,iarc,2] > 0):
                        trcdict['xmodel'].append(fwm_img[:-1,iarc,2])
                        trcdict['ymodel'].append(fwm_img[:-1,iarc,3])
                        trcdict['aduse'][iarc] = True
                #
                slf.all_trcdict.append(trcdict.copy())
            else:
                slf.all_trcdict.append(None)
            islit += 1
        # FInish
        return slf


    def _analyze_lines(self, slit):
        """
        Analyze the tilts of the arc lines in a given slit/order

        Wrapper to artracewave.analyze_lines()

        Parameters
        ----------
        slit : int

        Returns
        -------
        self.badlines

        """
        self.badlines, self.all_ttilts[slit] = artracewave.analyze_lines(
            self.msarc, self.all_trcdict[slit], slit, self.tslits_dict['pixcen'], self.settings)
        if self.badlines > 0:
            msgs.warn("There were {0:d} additional arc lines that should have been traced".format(self.badlines) +
                      msgs.newline() + "(perhaps lines were saturated?). Check the spectral tilt solution")
        # Step
        self.steps.append(inspect.stack()[0][3])
        return self.badlines

    def _extract_arcs(self, gen_satmask=False):
        """
        Extract the arcs down each slit/order

        Wrapper to ararc.get_censpec()

        Returns
        -------
        self.arccen
        self.arc_maskslit

        """
        # Extract an arc down each slit/order
        self.arccen, self.arc_maskslit, _ = ararc.get_censpec(self.tslits_dict['lcen'], self.tslits_dict['rcen'],
                                                              self.pixlocn, self.msarc, self.det,
                                                              gen_satmask=gen_satmask)
        self.satmask = np.zeros_like(self.msarc)
        # Step
        self.steps.append(inspect.stack()[0][3])
        return self.arccen, self.arc_maskslit

    def _fit_tilts(self, slit, show_QA=False, doqa=True):
        """

        Parameters
        ----------
        slit : int
        show_QA : bool, optional
          Show the QA plot (e.g. in a Notebook)
        doqa : bool, optional
          Perform the QA

        Returns
        -------
        self.tilts : ndarray

        """
        self.tilts, self.outpar = artracewave.fit_tilts(self.msarc, slit, self.all_ttilts[slit],
                                                        self.settings, setup=self.setup,
                                                        show_QA=show_QA, doqa=doqa)
        # Step
        self.steps.append(inspect.stack()[0][3])
        return self.tilts

    def _trace_tilts(self, slit, wv_calib=None):
        """

        Parameters
        ----------
        slit : int
        wv_calib : dict, optional
          Used only for avoiding ghosts


        Returns
        -------
        trcdict : dict
          Filled in self.all_trcdict[]

        """
        # Determine the tilts for this slit
        tracethresh_in = self.settings['tilts']['tracethresh']
        if isinstance(tracethresh_in,(float, int)):
            tracethresh = tracethresh_in
        elif isinstance(tracethresh_in, (list, np.ndarray)):
            tracethresh = tracethresh_in[slit]
        else:
            raise ValueError('Invalid input for parameter tracethresh')

        trcdict = artracewave.trace_tilt(self.tslits_dict['pixcen'], self.tslits_dict['lcen'],
                                         self.tslits_dict['rcen'], self.det,
                                         self.msarc, slit, self.settings_det, self.settings,
                                         censpec=self.arccen[:, slit], nsmth=3,
                                         tracethresh=tracethresh,
                                         wv_calib=wv_calib)
        # Load up
        self.all_trcdict[slit] = trcdict.copy()
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return trcdict

    def run(self, maskslits=None, doqa=True, wv_calib=None, gen_satmask=False):
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
        maskslits : ndarray (bool), optional
        doqa : bool
        wv_calib : dict
        gen_satmask : bool, optional
          Generate a saturation mask?

        Returns
        -------
        self.final_tilts
        maskslits
        """
        # If the user sets no tilts, return here
        if self.settings['tilts']['method'].lower() == "zero":
            # Assuming there is no spectral tilt
            self.final_tilts = np.outer(np.linspace(0.0, 1.0, self.msarc.shape[0]), np.ones(self.msarc.shape[1]))
            return self.final_tilts, None, None

        if maskslits is None:
            maskslits = np.zeros(self.nslit, dtype=bool)

        # Extract the arc spectra for all slits
        self.arccen, self.arc_maskslit = self._extract_arcs(gen_satmask=gen_satmask)

        # maskslit
        self.mask = maskslits & (self.arc_maskslit==1)
        gdslits = np.where(self.mask == 0)[0]

        # Final tilts image
        self.final_tilts = np.zeros_like(self.msarc)
        # Loop on all slits
        for slit in gdslits:
            # Trace
            _ = self._trace_tilts(slit, wv_calib=wv_calib)

            # Model line-by-line
            _ = self._analyze_lines(slit)

            # 2D model of the tilts
            #   Includes QA
            self.tilts = self._fit_tilts(slit, doqa=doqa)

            # Save to final image
            word = self.tslits_dict['slitpix'] == slit+1
            self.final_tilts[word] = self.tilts[word]

        return self.final_tilts, maskslits

    def _qa(self, slit):
        """
        QA
          Wrapper to artraceslits.slit_trace_qa()

        Parameters
        ----------
        slit : int

        Returns
        -------

        """
        self.tiltsplot, self.ztilto, self.xdat = artracewave.prep_tilts_qa(
            self.msarc, self.all_ttilts[slit], self.tilts, self.all_trcdict[slit]['arcdet'],
            self.pixcen, slit)

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
            # Bad slit?
            if self.mask[slit]:
                continue
            # fweight and model
            xtfits = self.all_trcdict[slit]['xtfit']  # For convenience
            xszs = [len(xtfit) if xtfit is not None else 0 for xtfit in xtfits]
            maxx = np.max(xszs)
            # Add 1 to pack in ycen
            fwm_img = np.zeros((maxx+1, len(xtfits), 4)) - 9999999.9
            # Fill fweight and model
            model_cnt = 0
            for kk, xtfit in enumerate(xtfits):
                if xtfit is None:
                    continue
                #
                fwm_img[0:xszs[kk], kk, 0] = xtfit
                fwm_img[0:xszs[kk], kk, 1] = self.all_trcdict[slit]['ytfit'][kk]
                #
                if self.all_trcdict[slit]['aduse'][kk]:
                    szmod = self.all_trcdict[slit]['xmodel'][model_cnt].size # Slits on edge can be smaller
                    fwm_img[0:szmod, kk, 2] = self.all_trcdict[slit]['xmodel'][model_cnt]
                    fwm_img[0:szmod, kk, 3] = self.all_trcdict[slit]['ymodel'][model_cnt]
                    model_cnt += 1
                    # ycen
                    xgd = self.all_trcdict[slit]['xtfit'][kk][self.all_trcdict[slit]['xtfit'][kk].size//2]
                    ycen = self.all_ttilts[slit][1][int(xgd),kk]
                    fwm_img[-1, kk, 1] = ycen
            hdu1 = fits.ImageHDU(fwm_img)
            hdu1.name = 'FWM{:03d}'.format(slit)
            hdul.append(hdu1)
        # Finish
        hdulist = fits.HDUList(hdul)
        hdulist.writeto(outfile, clobber=True)

    def show(self, attr, slit=None, display='ginga', cname=None):
        """
        Display an image or spectrum in TraceSlits

        Parameters
        ----------
        attr : str
          'fweight'  -- Show the msarc image and the tilts traced by fweight
          'model'    -- Show the msarc image and the poylynomial model fits to the individual arc lines that
                        were traced by fweight.
          'arcmodel -- This illustrates the global final 2-d model fit to the indivdiaul models of each traced fweight arc line
                       tilts evaluated at the location of the specific arclines that wered use for the fit.
          'final_tilts' -- Show the final 2-d tilt model for all the slits that were fit.
        slit : int, optional
                    -- The slit to plot. This needs to be an integer between 1 and nslit
        display : str (optional)
          'ginga' -- Display to an RC Ginga
        """
        # ToDO I don't see why we are not looping over all slits for all of this. Why should we restrict to an individual fit?
        if (self.tslits_dict['lcen'] is not None) and (slit is not None):
            sedges=(self.tslits_dict['lcen'][:,slit], self.tslits_dict['rcen'][:,slit])
        else:
            sedges = None
        if attr == 'fweight':
            if slit is None:
                msgs.error("Need to provide the slit with this option")
            ginga.chk_arc_tilts(self.msarc, self.all_trcdict[slit],
                                sedges=sedges)
            msgs.info("Green = ok line;  red=not used")
        elif attr == 'model':
            if slit is None:
                msgs.error("Need to provide the slit with this option")
            tmp = self.all_trcdict[slit-1].copy()
            tmp['xtfit'] = self.all_trcdict[slit-1]['xmodel']
            tmp['ytfit'] = self.all_trcdict[slit-1]['ymodel']
            ginga.chk_arc_tilts(self.msarc, tmp, sedges=sedges, all_green=True)
        elif attr in ['arcmodel']:
            if slit is None:
                msgs.error("Need to provide the slit with this option")
            tmp = self.all_trcdict[slit].copy()
            tmp['xtfit'] = []
            tmp['ytfit'] = []

            ynorm = np.outer(np.linspace(0., 1., self.msarc.shape[0]), np.ones(self.msarc.shape[1]))
            polytilts = (ynorm-self.tilts)*(self.msarc.shape[0]-1)
            slitpix = self.tslits_dict['slitpix']
            # arcdet is only the approximately nearest pixel (not even necessarily)
            for idx in np.where(self.all_trcdict[slit-1]['aduse'])[0]:
                xnow = np.arange(self.msarc.shape[1])
                if self.all_ttilts is not None:  # None if read from disk
                    xgd = self.all_trcdict[slit-1]['xtfit'][idx][self.all_trcdict[slit]['xtfit'][idx].size//2]
                    ycen = self.all_ttilts[slit-1][1][int(xgd),idx]
                else:
                    ycen = self.all_trcdict[slit-1]['ycen'][idx]
                ynow = ycen + polytilts[int(ycen),:]
                # Only plot the xnow, ynow values that are on this slit
                onslit = (slitpix[int(np.rint(xnow)),int(np.rint(ynow))]) == slit
                tmp['xtfit'].append(xnow[onslit])
                tmp['ytfit'].append(ynow[onslit])

            # Show
            msgs.warn("Display via tilts is not exact")  # Could make a correction.  Probably is close enough
            ginga.chk_arc_tilts(self.msarc, tmp, sedges=sedges, all_green=True, cname=cname)
        elif attr == 'final_tilts':
            if self.final_tilts is not None:
                ginga.show_image(self.final_tilts)
        else:
            msgs.error('Unrecognized attribute')

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

'''
def get_wv_tilts(det, setup, tilt_settings, settings_det, tslits_dict,
                 pixlocn, msarc, wv_calib, maskslits):
    """
    Load/Generate the tilts image

    Parameters
    ----------
    det : int
      Required for processing
    setup : str
      Required for MasterFrame loading
    tilt_settings : dict
      Tilt specific settings
    settings_det : dict
      Detector settings
    tslits_dict : dict
      Slits dict; required for processing
    pixlocn : ndarray
      Required for processing
    msarc : ndarray
      Required for processing
    wv_calib : dict
      1D wavelength fits
    maskslits : ndarray (bool)
      Indicates which slits are masked

    Returns
    -------
    mstilts : ndarray
      Tilt image
    wv_maskslits : ndarray (bool)
      Indicates slits that were masked (skipped)
    waveTilts : WaveTilts object
    """
    # Instantiate
    waveTilts = WaveTilts(msarc, settings=tilt_settings,
                                    det=det, setup=setup,
                                    tslits_dict=tslits_dict, settings_det=settings_det,
                                    pixlocn=pixlocn)
    # Master
    mstilts = waveTilts.master()
    if mstilts is None:
        mstilts, wt_maskslits = waveTilts.run(maskslits=maskslits,
                                              wv_calib=wv_calib)
        waveTilts.save_master()
    else:
        wt_maskslits = np.zeros_like(maskslits, dtype=bool)
    # Return
    return mstilts, wt_maskslits, waveTilts
'''
