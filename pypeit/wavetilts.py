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
from pypeit.core import tracewave, pixels
from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph
import copy



class WaveTilts(masterframe.MasterFrame):
    """Class to guide slit/order tracing

    Args:
        msarc (ndarray): Arc image
        tslits_dict (dict): dict from TraceSlits class (e.g. slitpix)
        par (:class:`pypeit.par.pypeitpar.WaveTiltsPar`):
            The parameters used to fuss with the tilts
        wavepar (:class:`pypeit.par.pypeitpar.WaveSolutionPar`):
            The parameters used for the wavelength solution
        det (int): Detector index

    Attributes:
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

    def __init__(self, msarc, tslits_dict, spectrograph, par, wavepar, det=1,
                 master_key=None, master_dir=None, reuse_masters=False, redux_path=None, bpm=None):

        self.spectrograph = spectrograph
        self.par = par # pypeitpar.WaveTiltsPar() if par is None else par
        self.wavepar = wavepar # pypeitpar.WavelengthSolutionPar() if wavepar is None else wavepar
        # MasterFrame
        masterframe.MasterFrame.__init__(self, self.frametype, master_key,
                                         master_dir, reuse_masters=reuse_masters)

        # Parameters (but can be None)
        self.msarc = msarc
        self.tslits_dict = tslits_dict
        self.bpm = bpm
        # Optional parameters
        self.det = det
        self.redux_path = redux_path
        #
        if self.spectrograph is not None:
            self.nonlinear_counts = self.spectrograph.detector[self.det-1]['saturation']*self.spectrograph.detector[self.det-1]['nonlinear']
        else:
            self.nonlinear_counts=1e10
        # Set the slitmask and slit boundary related attributes that the code needs for execution. This also deals with
        # arcimages that have a different binning then the trace images used to defined the slits
        if self.tslits_dict is not None and self.msarc is not None:
            self.slitmask_science = pixels.tslits2mask(self.tslits_dict)
            inmask = (self.bpm == 0) if self.bpm is not None else np.ones_like(self.slitmask_science, dtype=bool)
            self.shape_science = self.slitmask_science.shape
            self.shape_arc = self.msarc.shape
            self.nslits = self.tslits_dict['slit_left'].shape[1]
            self.slit_left = arc.resize_slits2arc(self.shape_arc, self.shape_science, self.tslits_dict['slit_left'])
            self.slit_righ = arc.resize_slits2arc(self.shape_arc, self.shape_science, self.tslits_dict['slit_righ'])
            self.slitcen   = arc.resize_slits2arc(self.shape_arc, self.shape_science, self.tslits_dict['slitcen'])
            self.slitmask  = arc.resize_mask2arc(self.shape_arc, self.slitmask_science)
            self.inmask = (arc.resize_mask2arc(self.shape_arc, inmask)) & (self.msarc < self.nonlinear_counts)
        else:
            self.slitmask_science = None
            self.shape_science = None
            self.shape_arc = None
            self.nslits = 0
            self.slit_left = None
            self.slit_righ = None
            self.slitcen = None
            self.slitmask = None
            self.inmask = None

        # Key Internals
        self.mask = None
        self.all_trace_dict = [None]*self.nslits
        self.tilts = None
        # 2D fits are stored as a dictionary rather than list because we will jsonify the dict
        self.all_fit_dict = [None]*self.nslits
        self.steps = []
        # Main outputs
        self.final_tilts = None
        self.fit_dict = None
        self.trace_dict = None


    def extract_arcs(self, slitcen, slitmask, msarc, inmask):
        """
        Extract the arcs down each slit/order

        Wrapper to arc.get_censpec()

        Args:
            slitcen (ndarray): Image for tracing
            slitmask (ndarray):
            msarc (ndarray):
            inmask (ndarray):

        Returns:
            ndarray, ndarray:  Extracted arcs

        """
        arccen, arc_maskslit = arc.get_censpec(slitcen, slitmask, msarc,
                                               inmask=inmask, nonlinear_counts=self.nonlinear_counts)
        # Step
        self.steps.append(inspect.stack()[0][3])
        return arccen, arc_maskslit

    def find_lines(self, arcspec, slit_cen, slit, debug=False):
        """
        Find the lines for tracing

        Wrapper to tracewave.tilts_find_lines()

        Args:
            arcspec:
            slit_cen:
            slit (int):
            debug:

        Returns:
            ndarray, ndarray:  Spectral, spatial positions of lines to trace

        """


        if self.par['idsonly'] is not None:
            # Put in some hook here for getting the lines out of the wave calib for i.e. LRIS ghosts.
            only_these_lines = None
            pass
        else:
            only_these_lines = None

        tracethresh = self._parse_param(self.par, 'tracethresh', slit)
        lines_spec, lines_spat = tracewave.tilts_find_lines(
            arcspec, slit_cen, tracethresh=tracethresh, sig_neigh=self.par['sig_neigh'],
            nfwhm_neigh=self.par['nfwhm_neigh'],only_these_lines=only_these_lines, fwhm=self.wavepar['fwhm'],
            nonlinear_counts=self.nonlinear_counts, debug_peaks=False, debug_lines = debug)

        self.steps.append(inspect.stack()[0][3])
        return lines_spec, lines_spat



    def fit_tilts(self, trc_tilt_dict, thismask, slit_cen, spat_order, spec_order, slit, show_QA=False, doqa=True, debug=False):
        """
        Fit the tilts

        Args:
            trc_tilt_dict (dict): Contains information from tilt tracing
            slit_cen (ndarray): (nspec,) Central trace for this slit
            spat_order (int): Order of the 2d polynomial fit for the spatial direction
            spec_order (int): Order of the 2d polytnomial fit for the spectral direction
            slit (int): integer index for the slit in question

        Optional Args:
            show_QA: bool, default = False
                show the QA instead of writing it out to the outfile
            doqa: bool, default = True
                Construct the QA plot
            debug: bool, default = False
                Show additional plots useful for debugging.

        Returns:
           (tilts, coeffs)
            tilts: ndarray (nspec, nspat)
               tilts image
            coeff: ndarray (spat_order + 1, spec_order+1)
               Array containing the coefficients for the 2d legendre polynomial fit
        """



        # Now perform a fit to the tilts
        tilt_fit_dict, trc_tilt_dict_out = tracewave.fit_tilts(
            trc_tilt_dict, thismask, slit_cen, spat_order=spat_order, spec_order=spec_order,maxdev=self.par['maxdev2d'],
            sigrej=self.par['sigrej2d'],func2d=self.par['func2d'],doqa=doqa,master_key=self.master_key,slit=slit, show_QA=show_QA,
            out_dir=self.redux_path, debug=debug)

        # Evaluate the fit
        #tilts = tracewave.fit2tilts((tilt_fit_dict['nspec'], tilt_fit_dict['nspat']),slit_cen,tilt_fit_dict['coeff2'], tilt_fit_dict['func'])

        # Populate the fit dict, and update the all_trace_dict
        self.all_fit_dict[slit] = copy.deepcopy(tilt_fit_dict)
        self.all_trace_dict[slit] = copy.deepcopy(trc_tilt_dict_out)

        self.steps.append(inspect.stack()[0][3])
        return tilt_fit_dict['coeff2']

    def trace_tilts(self, arcimg, lines_spec, lines_spat, thismask, slit_cen):
        """
        Trace the tilts

        Args:
            arcimg (ndarray): (nspec, nspat) Arc image
            lines_spec (ndarray): (nlines) Array containing the spectral pixel location of each line found for this slit
            lines_spat (ndarray): (nlines)
               Array containing the spatial pixel location of each line, which is the slitcen evaluate at the spectral position position
               of the line stored in lines_spec
            thismask (ndarray): (nspec, nspat), type=bool
               Image indicating which pixels lie on the slit in equation. True = on the slit. False = not on slit
            slit (int): Integer index indicating the slit in question


        Returns:
            dict: Dictionary containing informatin on the traced tilts required to fit the filts.

        """

        trace_dict = tracewave.trace_tilts(
            arcimg, lines_spec, lines_spat, thismask, slit_cen, inmask=self.inmask, fwhm=self.wavepar['fwhm'],
            spat_order=self.par['spat_order'], maxdev_tracefit=self.par['maxdev_tracefit'],
            sigrej_trace=self.par['sigrej_trace'])

        # Load up
        #self.all_trace_dict[slit] = copy.deepcopy(trace_dict)
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return trace_dict


    def run(self, maskslits=None, doqa=True, debug=False, show=False):
        """ Main driver for tracing arc lines

            Code flow:
               1.  Extract an arc spectrum down the center of each slit/order
               2.  Loop on slits/orders
                 i.   Trace and fit the arc lines (This is done twice, once with trace_crude as the tracing crutch, then
                      again with a PCA model fit as the crutch)
                 iii.  2D Fit to the offset from slitcen
                 iv. Save

        Keyword Args:
            maskslits (ndarray of bool, optional):
            doqa (bool):
            show (bool):

        Returns:
            dict, ndarray:  Tilts dict and maskslits array

        """

        if maskslits is None:
            maskslits = np.zeros(self.nslits, dtype=bool)

        # Extract the arc spectra for all slits
        self.arccen, self.arc_maskslit = self.extract_arcs(self.slitcen, self.slitmask, self.msarc, self.inmask)

        # maskslit
        self.mask = maskslits & (self.arc_maskslit==1)
        gdslits = np.where(self.mask == 0)[0]

        # Final tilts image
        self.final_tilts = np.zeros(self.shape_science,dtype=float)
        max_spat_dim = (np.asarray(self.par['spat_order']) + 1).max()
        max_spec_dim = (np.asarray(self.par['spec_order']) + 1).max()
        self.coeffs = np.zeros((max_spec_dim, max_spat_dim,self.nslits))
        self.spat_order = np.zeros(self.nslits, dtype=int)
        self.spec_order = np.zeros(self.nslits, dtype=int)

        # TODO sort out show methods for debugging
        #if show:
        #    viewer,ch = ginga.show_image(self.msarc*(self.slitmask > -1),chname='tilts')

        # Loop on all slits
        for slit in gdslits:
            msgs.info('Computing tilts for slit {:d}/{:d}'.format(slit,self.nslits-1))
            # Identify lines for tracing tilts
            self.lines_spec, self.lines_spat = self.find_lines(self.arccen[:,slit], self.slitcen[:,slit], slit, debug=debug)

            thismask = self.slitmask == slit
            # Trace
            self.trace_dict = self.trace_tilts(self.msarc, self.lines_spec, self.lines_spat, thismask, self.slitcen[:,slit])
            #if show:
            #    ginga.show_tilts(viewer, ch, self.trace_dict)

            self.spat_order[slit] = self._parse_param(self.par, 'spat_order', slit)
            self.spec_order[slit] = self._parse_param(self.par, 'spec_order', slit)
            # 2D model of the tilts, includes construction of QA
            coeff_out = self.fit_tilts(self.trace_dict, thismask, self.slitcen[:,slit], self.spat_order[slit],
                                       self.spec_order[slit], slit,doqa=doqa, show_QA = show, debug=show)
            self.coeffs[0:self.spec_order[slit]+1, 0:self.spat_order[slit]+1 , slit] = coeff_out
            # Tilts are created with the size of the original slitmask, which corresonds to the same binning
            # as the science images, trace images, and pixelflats etc.
            self.tilts = tracewave.fit2tilts(self.slitmask_science.shape, coeff_out, self.par['func2d'])
            # Save to final image
            thismask_science = self.slitmask_science == slit
            self.final_tilts[thismask_science] = self.tilts[thismask_science]

        self.tilts_dict = {'tilts':self.final_tilts, 'coeffs':self.coeffs, 'slitcen': self.slitcen, 'func2d':self.par['func2d'],
                           'nslit': self.nslits, 'spat_order': self.spat_order, 'spec_order': self.spec_order}
        return self.tilts_dict, maskslits

    def load_master(self, filename, exten=0, force=False):
        """
        Load the master frame

        Args:
            filename (str):  Master file
            exten (int, optinal):
            force (bool, optional):

        Returns:
            dict:  Tilts dict

        """


        # Does the master file exist?
        if not os.path.isfile(filename):
            msgs.warn("No Master frame found of type {:s}: {:s}".format(self.frametype, filename))
            if force:
                msgs.error("Crashing out because reduce-masters-force=True:" + msgs.newline() + filename)
            return None, None
        else:
            msgs.info("Loading a pre-existing master calibration frame of type: {:}".format(self.frametype) + " from filename: {:}".format(filename))
            hdu = fits.open(filename)
            head0 = hdu[0].header
            tilts = hdu[0].data
            coeffs = hdu[1].data
            slitcen = hdu[2].data
            spat_order = hdu[3].data
            spec_order = hdu[4].data
            tilts_dict = {'tilts':tilts,'coeffs':coeffs,'slitcen':slitcen,'func2d': head0['FUNC2D'], 'nslit': head0['NSLIT'],
                          'spat_order':spat_order, 'spec_order':spec_order}
            return tilts_dict, head0

    def save_master(self, tilts_dict, outfile=None, steps=None, overwrite=True):
        """
        Save the tilts dict to a MasterFile

        Args:
            tilts_dict (dict):  To be saved
            outfile (str):

        Returns:

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
        hdu0.name='TILTS'
        hdu0.header['FUNC2D'] = tilts_dict['func2d']
        hdu0.header['NSLIT'] =  tilts_dict['nslit']
        hdul = [hdu0]
        hdu_coeff = fits.ImageHDU(tilts_dict['coeffs'])
        hdu_coeff.name='COEFFS'
        hdul.append(hdu_coeff)
        hdu_slitcen = fits.ImageHDU(tilts_dict['slitcen'])
        hdu_slitcen.name = 'SLITCEN'
        hdul.append(hdu_slitcen)
        hdu_spat_order = fits.ImageHDU(tilts_dict['spat_order'])
        hdu_spat_order.name = 'SPAT_ORDER'
        hdul.append(hdu_spat_order)
        hdu_spec_order = fits.ImageHDU(tilts_dict['spec_order'])
        hdu_spec_order.name = 'SPEC_ORDER'
        hdul.append(hdu_spec_order)
        # Finish
        hdulist = fits.HDUList(hdul)
        hdulist.writeto(_outfile, clobber=True)

    def _parse_param(self, par, key, slit):
        """
        Grab a parameter for a given slit

        Args:
            par (ParSet):
            key (str):
            slit (int):

        Returns:
            object:  Value of the parameter

        """

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
            viewer, ch, self.trace_dict, sedges=(self.tslits_dict['slit_left'][:,slit],self.tslits_dict['slit_righ'][:,slit]),
            points = True, clear_canvas=True)

        # TODO Need to update the show function!

        """
        # ToDO I don't see why we are not looping over all slits for all of this. Why should we restrict to an individual fit?
        if (self.tslits_dict['slit_left'] is not None) and (slit is not None):
            sedges=(self.tslits_dict['slit_left'][:,slit], self.tslits_dict['slit_righ'][:,slit])
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



def load_tilts(filename):
    """
    Utility function which enables one to load the tilts from a master file in one line of code without
    instantiating the class.

    Args:
        filename (str): Master file name

    Returns:
        dict:  The trace slits dict

    """

    waveTilts = WaveTilts(None, None, None, None, None)
    tilts_dict, _ = waveTilts.load_master(filename)
    return tilts_dict['tilts']

