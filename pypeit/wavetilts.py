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

    def __init__(self, msarc, spectrograph=None, par=None, det=None, setup=None, master_dir=None,
                 mode=None, tslits_dict=None, redux_path=None, bpm=None):

        # TODO: (KBW) Why was setup='' in this argument list and
        # setup=None in all the others?  Is it because of the
        # from_master_files() classmethod below?  Changed it to match
        # the rest of the MasterFrame children.

        # Instantiate the spectograph
        # TODO: (KBW) Do we need this?  It's only used to get the
        # non-linear counts and the name of the master directory

        self.spectrograph = load_spectrograph(spectrograph)

        # MasterFrame
        masterframe.MasterFrame.__init__(self, self.frametype, setup,
                                         master_dir=master_dir, mode=mode)

        self.par = pypeitpar.WaveTiltsPar() if par is None else par

        # Parameters (but can be None)
        self.msarc = msarc
        if bpm is None:
            self.bpm = np.zeros_like(msarc)
        else:
            self.bpm = bpm
        self.tslits_dict = tslits_dict

        # Optional parameters
        self.det = det
        self.redux_path = redux_path

        # Attributes
        if self.tslits_dict is not None:
            self.nslit = self.tslits_dict['lcen'].shape[1]
        else:
            self.nslit = 0
        self.steps = []
        self.slitmask = None

        # Key Internals
        self.mask = None
        self.all_trcdict = [None]*self.nslit
        self.tilts = None
        self.all_ttilts = [None]*self.nslit

        # Main outputs
        self.final_tilts = None
        self.coeffs = None
        self.tilts_dict = None

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


    def _analyze_lines(self, slit):
        """
        Analyze the tilts of the arc lines in a given slit/order

        Wrapper to tracewave.analyze_lines()

        Parameters
        ----------
        slit : int

        Returns
        -------
        self.badlines

        """
        self.badlines, self.all_ttilts[slit] \
                = tilts.analyze_lines(self.msarc, self.all_trcdict[slit], slit,
                                            self.tslits_dict['pixcen'], order=self.par['order'],
                                            function=self.par['function'])
        if self.badlines > 0:
            msgs.warn('There were {0:d} additional arc lines that '.format(self.badlines) +
                      'should have been traced' + msgs.newline() + '(perhaps lines were '
                      'saturated?). Check the spectral tilt solution')
        # Step
        self.steps.append(inspect.stack()[0][3])
        return self.badlines

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
        inmask = (self.bpm == 0) if self.bpm is not None else None
        self.arccen, self.arc_maskslit = arc.get_censpec(self.tslits_dict['lcen'], self.tslits_dict['rcen'],
                                                         self.slitmask, self.msarc, inmask = inmask)
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
        coeffs

        """
        self.tilts, coeffs, self.outpar = tilts.fit_tilts(self.msarc, slit, self.all_ttilts[slit],
                                                        order=self.par['order'],
                                                        yorder=self.par['yorder'],
                                                        func2D=self.par['func2D'],
                                                        setup=self.setup, show_QA=show_QA,
                                                        doqa=doqa, out_dir=self.redux_path)
        # Step
        self.steps.append(inspect.stack()[0][3])
        return self.tilts, coeffs

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
        tracethresh_in = self.par['tracethresh']
        if isinstance(tracethresh_in,(float, int)):
            tracethresh = tracethresh_in
        elif isinstance(tracethresh_in, (list, np.ndarray)):
            tracethresh = tracethresh_in[slit]
        else:
            raise ValueError('Invalid input for parameter tracethresh')

        nonlinear_counts = self.spectrograph.detector[self.det-1]['saturation'] \
                                * self.spectrograph.detector[self.det-1]['nonlinear']

        # JFH Code block starts here
        ########
        # def get_tilts(npca = 2, fwhm=4.0, ncoeff=5, maxdev_tracefit=0.1,percentile_reject=0.10, max_badpix_frac=0.20, tcrude_maxerr=1.0,
        # tcrude_maxshift=3.0, tcrude_maxshift0=3.0, tcrude_nave=5,)

        from pypeit.core import pixels
        from pypeit.core import extract

        show_tilts = True

        nspat = self.msarc.shape[1]
        nspec = self.msarc.shape[0]
        arcimg = self.msarc
        arc_spec = self.arccen[:, slit]
        slit_left = self.tslits_dict['lcen'][:,slit].copy()
        slit_righ = self.tslits_dict['rcen'][:,slit].copy()
        inmask = (self.bpm == False)
        # Center of the slit
        slit_cen = (slit_left + slit_righ)/2.0

        slitmask = self.spectrograph.slitmask(self.tslits_dict)
        thismask = slitmask == slit

        # Tilt specific Optional parameters
        tracethresh = 20.0 # significance threshold for an arc line to be traced
        sigdetect = 5.0 # This is the significance threshold for finding neighboring lines. The number of these neighboring lines
        #  determines how many nsig > tracethresh lines that may be removed because they are too close.
        npix_neigh = 7.0
        only_these_lines = None # These are lines from the wavelength solution that are known to be good. If set code only uses these
        # identified lines to avoid problems with ghosts (is this still necessary with new tracing code?)
        ncoeff = 3 # Order of legendre polynomial fits to the tilts
        spec_order = 2
        maxdev_2dfit = 1.0 #
        debug = True
        maxdev_tracefit = 1.0 # maximum deviation
        sigrej_trace = 3.0 # From each line we compute a median absolute deviation of the trace from the polynomial fit. We then
        # analyze the distribution of mads for all the lines, and reject sigrej_trace outliers from that distribution.
        max_badpix_frac = 0.20 # Maximum fraction of total pixels masked by the trace_gweight algorithm (because the residuals are too large)
        # Trace Crude parameters
        tcrude_maxerr = 1.0 #
        tcrude_maxshift = 3.0
        tcrude_maxshift0 = 3.0,
        tcrude_nave = 5
        show_tracefits = False # show the trace fits

        # Optional Parameters for tilts_find_lines
        fwhm = 4.0 # expected spectral fwhm of the arc lines
        fit_frac_fwhm = 1.25
        cont_frac_fwhm = 1.0
        max_frac_fwhm = 2.5
        cont_samp = 30
        niter_cont = 3
        debug_lines = False

        lines_spec, lines_spat = tracewave.tilts_find_lines(
            arc_spec, slit_cen, tracethresh=tracethresh, sigdetect=sigdetect, npix_neigh=npix_neigh,
            only_these_lines=only_these_lines, fwhm=fwhm, nonlinear_counts=nonlinear_counts, fit_frac_fwhm=fit_frac_fwhm,
            cont_frac_fwhm=cont_frac_fwhm, max_frac_fwhm=max_frac_fwhm, cont_samp=cont_samp, niter_cont=niter_cont,
            debug = debug_lines)

        slit_width = int(np.ceil((slit_righ - slit_left).max()))

        trc_tilt_dict0 = tracewave.trace_tilts(arcimg, lines_spec, lines_spat, slit_width, thismask, inmask=inmask,
                                               tilts_guess=None,fwhm=fwhm,ncoeff=ncoeff, maxdev_tracefit=maxdev_tracefit,
                                               sigrej_trace=sigrej_trace, max_badpix_frac=max_badpix_frac,
                                               tcrude_maxerr=tcrude_maxerr, tcrude_maxshift=tcrude_maxshift,
                                               tcrude_maxshift0=tcrude_maxshift0,
                                               tcrude_nave=tcrude_nave,show_fits=show_tracefits)
        # Now evaluate the model of the tilts for all of our lines
        # Testing!!!!
        # Now perform a fit to the tilts
        tilt_fit_dict0 = tracewave.fit_tilts(trc_tilt_dict0, spat_order=3, spec_order=2, maxdev=maxdev_2dfit, debug=True,
                                            doqa=True,setup='test',slit=0, show_QA=False, out_dir='./')
        # Now evaluate the model of the tilts for all of our lines
        piximg0 = tracewave.fit2piximg(tilt_fit_dict0)


        # Do a PCA fit, which rejects some outliers
        iuse = trc_tilt_dict0['use_tilt']
        nuse =np.sum(iuse)
        msgs.info('PCA modeling {:d} good tilts'.format(nuse))
        coeff_npoly_pca = None
        #pca_explained_var = None
        # TODO Should we truncate this PCA by hand, or just let it explain variance
        pca_explained_var = 99.8
        npca = 2
        # PCA fit good orders, and predict bad orders
        pca_fit, poly_fit_dict, pca_mean, pca_vectors = extract.pca_trace(
            trc_tilt_dict0['tilts_sub_fit'], predict=np.invert(iuse), npca = npca, coeff_npoly=coeff_npoly_pca,
            order_vec=lines_spec,xinit_mean=lines_spec, minv=0.0, maxv=float(trc_tilt_dict0['nsub'] - 1), debug=True)

        # Now trace again with the PCA predictions as the starting crutches
        trc_tilt_dict1 = tracewave.trace_tilts(arcimg, lines_spec, lines_spat, slit_width, thismask, inmask=inmask,
                                               tilts_guess=pca_fit, fwhm=fwhm,ncoeff=ncoeff, maxdev_tracefit=maxdev_tracefit,
                                               sigrej_trace=sigrej_trace, max_badpix_frac=max_badpix_frac,
                                               show_fits=show_tracefits)


        # Now perform a fit to the tilts
        tilt_fit_dict1 = tracewave.fit_tilts(trc_tilt_dict1, spat_order=ncoeff, spec_order=spec_order, maxdev=maxdev_2dfit, debug=True,
                                            doqa=True,setup='test',slit=0, show_QA=False, out_dir='./')
        # Now evaluate the model of the tilts for all of our lines
        piximg1 = tracewave.fit2piximg(tilt_fit_dict1)
        # Since the y independent variable is the tilts in the way we do the 2d fit, and soem tilts are spurios, it does
        # no good to evaluate the global fit at these spurious tilts to get the new tracing crutches. The following is a
        # hack to determine this from the piximg. Perhaps there is a better way.
        spec_img = np.outer(np.arange(nspec), np.ones(nspat))
        spec_vec = np.arange(nspec)
        nlines=len(lines_spec)
        interp_flag = np.ones(nlines,dtype=bool)
        tilts_crutch = np.zeros((nspat, nlines))
        spat_min = trc_tilt_dict1['spat_min']
        spat_max = trc_tilt_dict1['spat_max']
        piximg1[np.invert(thismask)] = 1e10
        # Is there a faster more clever way to do this?
        for iline in range(nlines):
            min_spat = np.fmax(spat_min[iline], 0)
            max_spat = np.fmin(spat_max[iline], nspat - 1)
            min_spec = int(np.fmax(np.round(lines_spec[iline]) - np.round(0.01*nspec),0))
            max_spec = int(np.fmin(np.round(lines_spec[iline]) + np.round(0.01*nspec),nspec-1))
            piximg_sub = piximg1[min_spec:max_spec,:]
            for ispat in range(min_spat,max_spat):
                if np.any(np.diff(piximg_sub[:,ispat] < 0)):
                    # If we ever encounter an unsorted piximg_sub the logic below makes no sense so just use the
                    # previous polynomial fit as the crutch and exit this loop
                    tilts_crutch[ispat,iline] = trc_tilt_dict1['tilts_fit'][:,iline]
                    msgs.warn('piximg is not monotonic. Your tilts are probably bad')
                    break
                else:
                    tilts_crutch[ispat,iline] = np.interp(lines_spec[iline],piximg_sub[:,ispat],spec_vec[min_spec:max_spec])

        trc_tilt_dict1['tilts_crutch'] = tilts_crutch
        # Now trace again with the PCA predictions as the starting crutches
        trc_tilt_dict2 = tracewave.trace_tilts(arcimg, lines_spec, lines_spat, slit_width, thismask, inmask=inmask,
                                               tilts_guess=tilts_crutch, fwhm=fwhm,ncoeff=ncoeff, maxdev_tracefit=maxdev_tracefit,
                                               sigrej_trace=sigrej_trace, max_badpix_frac=max_badpix_frac,
                                               show_fits=show_tracefits)

        # Now perform a second fit to the tilts
        tilt_fit_dict2 = tracewave.fit_tilts(trc_tilt_dict2, spat_order=ncoeff, spec_order=spec_order, maxdev=maxdev_2dfit, debug=True,
                                            doqa=True,setup='test',slit=0, show_QA=False, out_dir='./')
        # Now evaluate the model of the tilts for all of our lines
        piximg2 = tracewave.fit2piximg(tilt_fit_dict2)

        from IPython import embed
        embed()

        # Now trace again with the piximg model as the starting crutches


        if show_tilts:
            viewer, ch = ginga.show_image(arcimg * thismask, chname='Tilts')
            # ginga.show_tilts(viewer, ch, tilts,tilts_spat, tilts_mask, tilts_err, sedges = (slit_left, slit_righ))
            ginga.show_tilts(viewer, ch, trc_tilt_dict1, crutch=False, sedges=(slit_left, slit_righ), points = True, clear_canvas=True)

#        from matplotlib import pyplot as plt
#        tilt_mask = trc_tilt_dict1['tilts_mask']
#        plt.plot(trc_tilt_dict1['tilts_spat'][tilt_mask], tilts_crutch[tilt_mask], 'ko', markersize=2.0)
#        plt.plot(trc_tilt_dict1['tilts_spat'][tilt_mask], trc_tilt_dict1['tilts_fit'][tilt_mask], 'ro', mfc='none', markersize=2.0)


        # Now do a fit


        """
        for iline in range(nlines):
            min_spat = np.fmax(spat_min[iline], 0)
            max_spat = np.fmin(spat_max[iline], nspat - 1)
            sub_img = (piximg*thismask)[:, min_spat:max_spat]
            spec_img_sub = spec_img[  :, min_spat:max_spat]
            ispec_min = np.argmin(np.abs(sub_img - lines_spec[iline]), axis=0)
            tilts_crutch[min_spat:max_spat,iline] = spec_img_sub[ispec_min,np.arange(sub_img.shape[1])]
        """




        # Load up
        self.all_trcdict[slit] = trc_tilt_dict.copy()
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
        if self.par['method'].lower() == "zero":
            # Assuming there is no spectral tilt
            self.final_tilts = np.outer(np.linspace(0.0, 1.0, self.msarc.shape[0]), np.ones(self.msarc.shape[1]))
            return self.final_tilts, None, None

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
        self.coeffs = np.zeros((self.par['order'] + 2,self.par['yorder'] +1,self.nslit))
        # Loop on all slits
        for slit in gdslits:
            # Trace
            _ = self._trace_tilts(slit, wv_calib=wv_calib)

            # Model line-by-line
            _ = self._analyze_lines(slit)

            # 2D model of the tilts
            #   Includes QA
            self.tilts, self.coeffs[:,:,slit] = self._fit_tilts(slit, doqa=doqa)

            # Save to final image
            word = self.slitmask == slit
            self.final_tilts[word] = self.tilts[word]

        self.tilts_dict = {'tilts':self.final_tilts, 'coeffs':self.coeffs, 'func2D':self.par['func2D']}
        return self.tilts_dict, maskslits

    def _qa(self, slit):
        """
        QA
          Wrapper to traceslits.slit_trace_qa()

        Parameters
        ----------
        slit : int

        Returns
        -------

        """
        self.tiltsplot, self.ztilto, self.xdat = tracewave.prep_tilts_qa(
            self.msarc, self.all_ttilts[slit], self.tilts, self.all_trcdict[slit]['arcdet'],
            self.pixcen, slit)

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
            tilts_dict = {'tilts':tilts,'coeffs':coeffs,'func2D': head1['FUNC2D']} # This is the tilts_dict
            return tilts_dict #, head0, [filename]

    # JFH THis routine does not follow the current master protocol of taking a data argument. There is no reason to
    # save all this other information here
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
        hdu_coeff = fits.ImageHDU(self.coeffs)
        hdu_coeff.header['FUNC2D'] = self.par['func2D']
        hdul.append(hdu_coeff)

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

