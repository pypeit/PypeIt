# Module for guiding Arc/Sky line tracing
from __future__ import absolute_import, division, print_function

import os
import inspect
import numpy as np

#from importlib import reload

from astropy.io import fits

from pypeit import msgs
from pypeit import masterframe
from pypeit import ginga
from pypeit.core import arc
from pypeit.core import tracewave
from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph
from pypeit import spectrographs
import copy

from pypeit import debugger


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
    
    # Frametype is a class attribute
    frametype = 'tilts'

    def __init__(self, msarc, tslits_dict, spectrograph=None, par=None, wavepar = None, det=None, setup=None, master_dir=None,
                 mode=None, redux_path=None, bpm=None):

        # TODO: (KBW) Why was setup='' in this argument list and
        # setup=None in all the others?  Is it because of the
        # from_master_files() classmethod below?  Changed it to match
        # the rest of the MasterFrame children.

        # Instantiate the spectograph
        # TODO: (KBW) Do we need this?  It's only used to get the
        # non-linear counts and the name of the master directory

        self.spectrograph = load_spectrograph(spectrograph)
        self.par = pypeitpar.WaveTiltsPar() if par is None else par
        self.wavepar = pypeitpar.WavelengthSolutionPar() if wavepar is None else wavepar

        # MasterFrame
        masterframe.MasterFrame.__init__(self, self.frametype, setup,
                                         master_dir=master_dir, mode=mode)


        # Parameters (but can be None)
        self.msarc = msarc
        self.tslits_dict = tslits_dict
        self.bpm = bpm
        self.inmask = (self.bpm == 0) if self.bpm is not None else None

        # Optional parameters
        self.det = det
        self.redux_path = redux_path

        # Attributes
        if self.tslits_dict is not None:
            self.nslit = self.tslits_dict['lcen'].shape[1]
            self.slitcen = self.tslits_dict['slitcen']
        else:
            self.nslit = 0
            self.slitcen = None

        self.steps = []
        self.slitmask = None

        # Key Internals
        self.mask = None
        self.all_trace_dict = [None]*self.nslit
        self.tilts = None
        # 2D fits are stored as a dictionary rather than list because we will jsonify the dict
        self.all_fit_dict = [None]*self.nslit

        # Main outputs
        self.final_tilts = None
        self.fit_dict = None
        self.trace_dict = None

    # This method does not appear finished
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

        # Instantiate
        slf = cls(None, setup=setup)
        msarc_file = masterframe.master_name('arc', setup, mdir)
        # Arc
        msarc, _, _ = slf.load_master(msarc_file)
        slf.msarc = msarc


        # Tilts
        mstilts_file = masterframe.master_name('tilts', setup, mdir)
        hdul = fits.open(mstilts_file)
        slf.final_tilts = hdul[0].data
        slf.tilts = slf.final_tilts
        slf.coeffs = slf.hdu[1].data

        # Dict
        slf.all_trcdict = []
        islit = 0
        for hdu in hdul[2:]:
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


    def _extract_arcs(self):
        """
        Extract the arcs down each slit/order

        Wrapper to arc.get_censpec()

        Returns
        -------
        self.arccen
        self.arc_maskslit

        """
        # Extract an arc down each slit/order
        slitmask = self.spectrograph.slitmask(self.tslits_dict) if self.slitmask is None else self.slitmask

        self.arccen, self.arc_maskslit = arc.get_censpec(self.tslits_dict['lcen'], self.tslits_dict['rcen'],
                                                         slitmask, self.msarc, inmask = self.inmask)
        # Step
        self.steps.append(inspect.stack()[0][3])
        return self.arccen, self.arc_maskslit

    def _find_lines(self, arcspec, slit_cen, slit, debug_lines=False):

        # Find good lines for the tilts
        nonlinear_counts = self.spectrograph.detector[self.det-1]['saturation']*self.spectrograph.detector[self.det-1]['nonlinear']

        if self.par['idsonly'] is not None:
            # Put in some hook here for getting the lines out of the wave calib for i.e. LRIS ghosts.
            only_these_lines = None
            pass
        else:
            only_these_lines = None

        tracethresh = self._parse_param(self.par, 'tracethresh', slit)
        lines_spec, lines_spat = tracewave.tilts_find_lines(
            arcspec, slit_cen, tracethresh=tracethresh, sigdetect=self.par['sigdetect'],
            nfwhm_neigh=self.par['nfwhm_neigh'],only_these_lines=only_these_lines, fwhm=self.wavepar['fwhm'],
            nonlinear_counts=nonlinear_counts, debug_lines = debug_lines)

        self.steps.append(inspect.stack()[0][3])
        return lines_spec, lines_spat



    def _fit_tilts(self, trc_tilt_dict, slit_cen, spat_order, spec_order, slit, show_QA=False, doqa=True, debug=False):
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
        coeffs

        """

        # Now perform a fit to the tilts
        tilt_fit_dict = tracewave.fit_tilts(
            trc_tilt_dict, spat_order=spat_order, spec_order=spec_order,maxdev=self.par['maxdev2d'],
            sigrej=self.par['sigrej2d'],func2d=self.par['func2d'],doqa=doqa,setup=self.setup,slit=slit, show_QA=show_QA,
            out_dir=self.redux_path, debug=debug)

        # Evaluate the fit
        tilts = tracewave.fit2tilts((tilt_fit_dict['nspec'], tilt_fit_dict['nspat']),slit_cen,tilt_fit_dict['coeff2'], tilt_fit_dict['func'])

        # Step
        self.all_fit_dict[slit] = copy.deepcopy(tilt_fit_dict)
        self.steps.append(inspect.stack()[0][3])
        return tilts, tilt_fit_dict['coeff2']

    def _trace_tilts(self, arcimg, lines_spec, lines_spat, thismask, slit):
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

        trace_dict = tracewave.trace_tilts(
            arcimg, lines_spec, lines_spat, thismask, inmask=self.inmask, fwhm=self.wavepar['fwhm'],
            spat_order=self.par['spat_order'], maxdev_tracefit=self.par['maxdev_tracefit'],
            sigrej_trace=self.par['sigrej_trace'])

        # Load up
        self.all_trace_dict[slit] = copy.deepcopy(trace_dict)
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return trace_dict


    def run(self, maskslits=None, doqa=True, show_QA=True, debug=True):
        """ Main driver for tracing arc lines

            Code flow:
               1.  Extract an arc spectrum down the center of each slit/order
               2.  Loop on slits/orders
                 i.   Trace and fit the arc lines (This is done twice, once with trace_crude as the tracing crutch, then
                      again with a PCA model fit as the crutch)
                 iii.  2D Fit to the offset from slitcen
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

        if maskslits is None:
            maskslits = np.zeros(self.nslit, dtype=bool)

        self.slitmask = self.spectrograph.slitmask(self.tslits_dict)

        # Extract the arc spectra for all slits
        self.arccen, self.arc_maskslit = self._extract_arcs()

        # maskslit
        self.mask = maskslits & (self.arc_maskslit==1)
        gdslits = np.where(self.mask == 0)[0]

        # Final tilts image
        self.final_tilts = np.zeros_like(self.msarc)
        max_spat_dim = (np.asarray(self.par['spat_order']) + 1).max()
        max_spec_dim = (np.asarray(self.par['spec_order']) + 1).max()
        self.coeffs = np.zeros((max_spat_dim, max_spec_dim,self.nslit))
        self.spat_order = np.zeros(self.nslit, dtype=int)
        self.spec_order = np.zeros(self.nslit, dtype=int)

        # Loop on all slits
        for slit in gdslits:
            # Identify lines for tracing tilts
            self.lines_spec, self.lines_spat = self._find_lines(self.arccen[:,slit], self.slitcen[:,slit], slit)

            thismask = self.slitmask == slit
            # Trace
            self.trace_dict = self._trace_tilts(self.msarc, self.lines_spec, self.lines_spat, thismask, slit)

            self.spat_order[slit] = self._parse_param(self.par, 'spat_order', slit)
            self.spec_order[slit] = self._parse_param(self.par, 'spec_order', slit)
            # 2D model of the tilts, includes construction of QA

            self.tilts, coeff_out = self._fit_tilts(self.trace_dict, self.slitcen[:,slit], self.spat_order[slit],
                                                    self.spec_order[slit], slit,doqa=doqa, show_QA = show_QA, debug=debug)
            self.coeffs[0:self.spat_order[slit]+1, 0:self.spec_order[slit]+1 , slit] = coeff_out
            # Save to final image
            self.final_tilts[thismask] = self.tilts[thismask]

        self.tilts_dict = {'tilts':self.final_tilts, 'coeffs':self.coeffs, 'slitcen': self.slitcen, 'func2d':self.par['func2d'],
                           'nslit': self.nslit, 'spat_order': self.spat_order, 'spec_order': self.spec_order}
        return self.tilts_dict, maskslits

    def load_master(self, filename, exten = 0, force = False):


        # Does the master file exist?
        if not os.path.isfile(filename):
            msgs.warn("No Master frame found of type {:s}: {:s}".format(self.frametype, filename))
            if force:
                msgs.error("Crashing out because reduce-masters-force=True:" + msgs.newline() + filename)
            return None
        else:
            msgs.info("Loading a pre-existing master calibration frame of type: {:}".format(self.frametype) + " from filename: {:}".format(filename))
            hdu = fits.open(filename)
            head0 = hdu[0].header
            tilts = hdu[0].data
            head1 = hdu[1].header
            coeffs = hdu[1].data
            slitcen = hdu[2].data
            spat_order = hdu[3].data
            spec_order = hdu[4].data
            tilts_dict = {'tilts':tilts,'coeffs':coeffs,'slitcen':slitcen,'func2d': head1['FUNC2D'], 'nslit': head1['NSLIT'],
                          'spat_order':spat_order, 'spec_order':spec_order}
            return tilts_dict

    # JFH THis routine does not follow the current master protocol of taking a data argument. There is no reason to
    # save all this other information here
    def save_master(self, tilts_dict, outfile=None, steps=None, overwrite=True):
        """

        Parameters
        ----------
        outfile
        use_tilts_as_final

        Returns
        -------
        """


        _outfile = self.ms_name if outfile is None else outfile
        # Additional keywords for the Header
        keywds = None if steps is None else dict(steps=','.join(steps))
        # Check for existing
        if os.path.exists(_outfile) and (not overwrite):
            msgs.warn("This file already exists.  Use overwrite=True to overwrite it")
            return
        #
        msgs.info("Saving master {0:s} frame as:".format(self.frametype) + msgs.newline() + _outfile)
        hdu0 = fits.PrimaryHDU(tilts_dict['tilts'])
        hdul = [hdu0]
        hdu_coeff = fits.ImageHDU(tilts_dict['coeffs'])
        hdu_coeff.header['FUNC2D'] = tilts_dict['func2d']
        hdu_coeff.header['NSLIT'] = tilts_dict['nslit']
        hdul.append(hdu_coeff)
        hdu_slitcen = fits.ImageHDU(tilts_dict['slitcen'])
        hdul.append(hdu_slitcen)
        hdu_spat_order = fits.ImageHDU(tilts_dict['spat_order'])
        hdul.append(hdu_spat_order)
        hdu_spec_order = fits.ImageHDU(tilts_dict['spec_order'])
        hdul.append(hdu_spec_order)
        # Finish
        hdulist = fits.HDUList(hdul)
        hdulist.writeto(_outfile, clobber=True)

    def _parse_param(self, par, key, slit):

        # Find good lines for the tilts
        param_in = par[key]
        if isinstance(param_in, (float, int)):
            param = param_in
        elif isinstance(param_in, (list, np.ndarray)):
            param = param_in[slit]
        else:
            raise ValueError('Invalid input for parameter {:s}'.format(key))

        return param

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

        viewer, ch = ginga.show_image(self.arcimg*(self.slitmask == slit), chname='Tilts')
        ginga.show_tilts(
            viewer, ch, self.trace_dict, sedges=(self.tslits_dict['lcen'][:,slit],self.tslits_dict['rcen'][:,slit]),
            points = True, clear_canvas=True)

        # TODO Need to update the show function!

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
                onslit = (slitmask[int(np.rint(xnow)),int(np.rint(ynow))]) == slit
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
        """

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



#    def _qa(self, slit):
#        """
#        QA
#          Wrapper to traceslits.slit_trace_qa()
#
#        Parameters
#        ----------
#        slit : int
#
#        Returns
#        -------
#
#        """
#        self.tiltsplot, self.ztilto, self.xdat = tracewave.prep_tilts_qa(
#            self.msarc, self.all_ttilts[slit], self.tilts, self.all_trcdict[slit]['arcdet'],
#            self.pixcen, slit)

