""" Module for flexure routines

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import inspect
import pathlib

from astropy.io import ascii
from IPython import embed
import matplotlib
from matplotlib import pyplot as plt
import numpy as np

from pypeit import dataPaths
from pypeit import log
from pypeit import specobjs
from pypeit.core import fitting
from pypeit.datamodel import DataContainer
from pypeit.images.detector_container import DetectorContainer
from pypeit.images.mosaic import Mosaic


def sky_em_residuals(
    wave:np.ndarray, flux:np.ndarray, ivar:np.ndarray, sky_waves:np.ndarray, plot=False, noff=5.,
    nfit_min=20
):
    """
    Calculate residuals and other metrics for a set of input sky emission lines.

    Args:
        wave (`numpy.ndarray`_):
            Wavelengths (in air!)
        flux (`numpy.ndarray`_):
            Fluxes
        ivar (`numpy.ndarray`_):
            Inverse variance
        sky_waves (`numpy.ndarray`_):
            Skyline wavelengths (in air!)
        plot (bool, optional):
            If true, plot the residuals
        noff (int, optional):
            Range in Ang to analyze labout emission line. Defaults to 5.
        nfit_min (int, optional):
            Minimum number of pixels required to do a fit. Defaults to 20.

    Returns:
        tuple of `numpy.ndarray`_ -- sky line wavelength of good lines, wavelength offset,
            error in wavelength offset, sky line width,
            error in sky line width
    """


    dwave = []
    diff = []
    diff_err  = []
    los = []
    los_err= []

    good_ivar = ivar > 0

    # Loop on known sky lines
    for line in sky_waves: 
        wline = [line-noff,line+noff] 
        mw    = (wave > wline[0]) & (wave < wline[1]) & good_ivar
        
        # Reuire minimum number
        if np.sum(mw) <= nfit_min:
            continue

        p=[0,0,0,0]
        # Guess
        p0 = list(fitting.guess_gauss(wave[mw], flux[mw]))
        # Fit
        try:
            p, pcov = fitting.fit_gauss(wave[mw], flux[mw], w_out=1./np.sqrt(ivar[mw]),
                                        guesses=p0, nparam=4)
        except RuntimeError as e:
            log.warning('First attempt at Gaussian fit failed, ending with RuntimeError.  Original '
                      f'exception: {e.args[0]}  Assuming this is because it hit the maximum '
                      'number of function evaluations.  Trying again with a maximum of 10000.')
            # Try again with larger limit on the number of function evaluations
            p, pcov = fitting.fit_gauss(wave[mw], flux[mw], w_out=1./np.sqrt(ivar[mw]),
                                        guesses=p0, nparam=4, maxfev=10000)

        perr = np.sqrt(np.diag(pcov))
        #except:
        #    p=p0
        #    p[2] = -99
        #    perr=p0

        # Continue
        d = p[2] - line

        # For debugging
        if plot:
            gfit = fitting.gauss_4deg(wave[mw],*p)
            plt.figure(figsize=(8,3)) 
            plt.plot(wave[mw],gfit,'g')
            plt.plot(wave[mw],flux[mw])
            plt.title('{} {:0.2f} diff= {:0.3f}'.format(line,p[3],d))
            plt.show()

        # Check
        if not np.isfinite(perr[2]):
            perr[2] = 1000.
        # Save
        dwave = np.append(dwave,line)
        diff = np.append(diff,d)
        diff_err = np.append(diff_err,perr[2])
        los = np.append(los,p[3])
        los_err = np.append(los_err,perr[3])

    # Cut on quality
    m=(diff_err < 0.1) & (diff_err > 0.0)
    # Return
    return dwave[m], diff[m], diff_err[m], los[m], los_err[m]


# TODO -- Consider separating the methods from the DataContainer as per calibrations
class MultiSlitFlexure(DataContainer):
    """
    Class to perform multi-detector flexure analysis.

    Based on code written by Marla Geha for DEIMOS.

    The datamodel attributes are:

    .. include:: ../include/class_datamodel_multislitflexure.rst
    """

    # Set the version of this class
    version = '1.1.0'

    datamodel = {'s1dfile': dict(otype=str, descr='spec1d filename'), 
                 'PYP_SPEC': dict(otype=str, descr='PypeIt spectrograph name'),
                 'ndet': dict(otype=int, descr='Number of detectors per spectrum'),
                 'nslits': dict(otype=int, descr='Number of slits'),
                 'is_msc': dict(otype=np.ndarray, atype=(int, np.integer),
                                descr='Flag that the "det" is the mosaic ID (ndet, nslits)'),
                 'det': dict(otype=np.ndarray, atype=(int, np.integer),
                             descr='Integer identifiers for the detector or mosaic (ndet, nslits)'),
                 'SN': dict(otype=np.ndarray, atype=np.floating, descr='S/N (ndet, nslits)'),
                 'slitid': dict(otype=np.ndarray, atype=np.floating, descr='Slit ID (nslits)'),
                 'mn_wv': dict(otype=np.ndarray, atype=np.floating,
                               descr='Mininum wavelength of the slit [Ang] (nslits)'),
                 'indiv_fit_slope': dict(otype=np.ndarray, atype=np.floating,
                                         descr='Fits to each slit individually (nslits)'),
                 'indiv_fit_b': dict(otype=np.ndarray, atype=np.floating,
                                     descr='Same as above but for b (nslits)'),
                 'indiv_fit_los': dict(otype=np.ndarray, atype=np.floating,
                                       descr='Same as above but for line width (nslits)'),
                 'fit_slope': dict(otype=np.ndarray, atype=np.floating,
                                   descr='Fitted slope (nslits)'),
                 'fit_b': dict(otype=np.ndarray, atype=np.floating,
                               descr='Fitted b value(nslits)'),
                 'fit_los': dict(otype=np.ndarray, atype=np.floating,
                                 descr='Fitted line width(nslits)'),
                 'resid_sky': dict(otype=np.ndarray, atype=np.floating,
                                   descr='Residuals of flexure model on sky lines (nslits)'),
                 'objra': dict(otype=np.ndarray, atype=np.floating, descr='Object RA (nslits)'),
                 'objdec': dict(otype=np.ndarray, atype=np.floating, descr='Object DEC (nslits)'),
                 'maskdef_id': dict(otype=np.ndarray, atype=np.integer, descr='Mask ID (nslits)'),
                 'rms_arc': dict(otype=np.ndarray, atype=np.floating,
                                 descr='RMS of fit (ndet, nslits)')}

    internals = ['flex_par',        # Parameters (FlexurePar)
                 'spectrograph',    # spectrograph
                 'specobjs',        # SpecObjs object
                 'sobj_idx',        # (ndet, nslits); Index to specobjs (tuple of arrays)
                 'sky_table',       # Sky line table
                 # 2D models
                 'pmodel_m',
                 'pmodel_b',
                 'pmodel_l'
                ]

    def __init__(self, s1dfile=None, PYP_SPEC=None, nslits=None, det=None, 
                 SN=None, slitid=None, mn_wv=None, fit_slope=None, fit_b=None,
                 fit_los=None, objra=None, objdec=None, maskdef_id=None, rms_arc=None, 
                 resid_sky=None, indiv_fit_slope=None, indiv_fit_b=None,
                 indiv_fit_los=None):

        # Setup the DataContainer
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        _d = {k: values[k] for k in args[1:]}
        # Init
        super().__init__(d=_d)

        # Load up specobjs
        # NOTE: specobjs here *is* the pypeit module
        self.specobjs = specobjs.SpecObjs.from_fitsfile(self.s1dfile, chk_version=False) 
        #  Sky lines -- This one is ASCII, so don't use load_sky_spectrum()
        sky_file = 'sky_single_mg.dat'
        self.sky_table = ascii.read(dataPaths.sky_spec.get_file_path(sky_file))

    # NOTE: If you make changes to how this object is bundled into the output
    # datamodel, make sure you update the documentation in
    # doc/calibrations/flexure.rst!
    def _bundle(self):
        """
        Override the base class method simply to set the HDU extension name.
        """
        return super()._bundle(ext='FLEXURE')

    def init(self, spectrograph, par):
        """ Initialize this and that about the slits, par, spectrograph
        e.g. RA, DEC, S/N

        Args:
            spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
                The spectrograph instance that sets the instrument used to take
                the observations.  Used to set :attr:`spectrograph`.
            par (:class:`~pypeit.par.pypeitpar.FlexurePar`):
                The parameters used for the flexure processing
        """
        # Internals
        self.spectrograph = spectrograph
        self.flex_par = par
        # Set
        self.PYP_SPEC = self.spectrograph.name
        self.sobj_idx = self.spectrograph.spec1d_match_spectra(self.specobjs)
        #
        self.nslits = len(self.sobj_idx[0])
        self.ndet = len(self.sobj_idx)
        
        # Fill in 1D
        self['slitid'] = self.specobjs[self.sobj_idx[0]]['SLITID'].astype(float)
        self['objra'] = self.specobjs[self.sobj_idx[0]]['RA']
        self['objdec'] = self.specobjs[self.sobj_idx[0]]['DEC']
        #self['slitname'] = self.specobjs[self.sobj_idx[0]]['MASKDEF_OBJNAME']
        self['maskdef_id'] = self.specobjs[self.sobj_idx[0]]['MASKDEF_ID']

        # Compile the list of detector *names* once
        DETs = self.specobjs.DET
        # Find which ones are actually mosaics
        is_msc = np.array([Mosaic.name_prefix in d for d in DETs]).astype(np.uint16)
        # Use the relevant parser to get the integer identifier
        det_msc_num = np.array([Mosaic.parse_name(d) if m else DetectorContainer.parse_name(d) 
                                    for d,m in zip(DETs, is_msc)])
        # Then assign the attributes
        self.is_msc = np.vstack(tuple(is_msc[self.sobj_idx[det]] for det in range(self.ndet)))
        self.det = np.vstack(tuple(det_msc_num[self.sobj_idx[det]] for det in range(self.ndet)))

        # S/N and mn_wv from the spectra
        self['SN'] = np.zeros((self.ndet, self.nslits), dtype=float)
        self['mn_wv'] = np.zeros((self.ndet, self.nslits), dtype=float)
        for det in range(self.ndet):
            self['SN'][det] = [sobj.S2N for sobj in self.specobjs[self.sobj_idx[det]]]
            self['mn_wv'][det] = [sobj.mnx_wave[0] for sobj in self.specobjs[self.sobj_idx[det]]]

    def fit_mask_surfaces(self):
        """
        Fit 2D model to linear flexure models from each slit as a function of
        RA, DEC.
        """
        # Cut on S/N
        good_SN = self['SN'] > self.flex_par['multi_min_SN']
        good_slit = np.sum(good_SN, axis=0) == self.ndet

        # Basic stats
        mu = np.median(self['indiv_fit_slope'][good_slit])
        sd = np.std(self['indiv_fit_slope'][good_slit])
        mu2 = np.median(self['indiv_fit_b'][good_slit])
        sd2 = np.std(self['indiv_fit_b'][good_slit])

        # Cut down to +/- 2sigma
        mgood = (np.abs(self['indiv_fit_slope']-mu) < 2.*sd) \
                    & ( np.abs(self['indiv_fit_b']-mu2) < 2.*sd2) & good_slit

        # Fit me (without additional rejection)
        # TODO -- Allow for x,y position instead of RA, DEC
        self.pmodel_m = fitting.robust_fit(self['objra'][mgood],
                                           self['indiv_fit_slope'][mgood], (2,2),
                                           function='polynomial2d',
                                           x2=self['objdec'][mgood])
        self.pmodel_b = fitting.robust_fit(self['objra'][mgood],
                                           self['indiv_fit_b'][mgood], (2,2),
                                           function='polynomial2d',
                                           x2=self['objdec'][mgood])
        self.pmodel_l = fitting.robust_fit(self['objra'][mgood],
                                           self['indiv_fit_los'][mgood], (2,2),
                                           function='polynomial2d',
                                           x2=self['objdec'][mgood])

    def measure_sky_lines(self):
        """Main method to analyze the sky lines for all the slits
        """

        # Init
        for key in ['indiv_fit_slope', 'indiv_fit_b', 'indiv_fit_los']:
            self[key] = np.zeros(self.nslits)

        # Loop on slits
        for i in np.arange(0,self.nslits,1):
            if (i % 10) == 0:
                log.info("Working on slit {} of {}".format(i, self.nslits))

            if not np.all(self['SN'][:,i] > 1.):
                continue

            # Loop on detectors
            sky_lines, sky_diffs, sky_ediffs, sky_loss = [], [], [], []
            for det in range(self.ndet):
                sobj = self.specobjs[self.sobj_idx[det][i]]

                # Measure em
                # The following will break if only boxcar...
                # TODO -- Allow for boxcar
                sky_line, sky_diff, sky_ediff, los, _ = sky_em_residuals(
                    sobj['OPT_WAVE'], 
                    sobj['OPT_COUNTS_SKY'], 
                    sobj['OPT_COUNTS_IVAR'],
                    self.sky_table['Wave'])

                # Hold em
                sky_lines.append(sky_line)
                sky_diffs.append(sky_diff)
                sky_ediffs.append(sky_ediff)
                sky_loss.append(los)

            # Concatenate
            sky_lines = np.concatenate(sky_lines)
            sky_diffs = np.concatenate(sky_diffs)
            sky_ediffs = np.concatenate(sky_ediffs)
            sky_loss = np.concatenate(sky_loss)
            
            # FIT SINGLE SLIT SKY LINES WITH A LINE           
            linear_fit = fitting.robust_fit(sky_lines,
                                            sky_diffs,
                                            weights=1./sky_ediffs**2,  
                                            function='polynomial', 
                                            order=1,
                                            maxrej=1,  # Might increase
                                            lower=3., upper=3.)
            # Save 
            self['indiv_fit_b'][i]     = linear_fit.fitc[0]
            self['indiv_fit_slope'][i] = linear_fit.fitc[1]
            self['indiv_fit_los'][i]   = np.median(sky_loss)

    def update_fit(self):
        """Update fits for each slit based on 2D model
        """
        # Do it
        self['fit_slope'] = self.pmodel_m.eval(self['objra'],x2=self['objdec'])
        self['fit_b']     = self.pmodel_b.eval(self['objra'],x2=self['objdec'])
        self['fit_los']   = self.pmodel_l.eval(self['objra'],x2=self['objdec'])

        # CALCULATE RESIDUALS FROM FIT
        #   Only for QA (I think)
        resid_sky = []
        for i in range(self.nslits):

            # Require sufficient S/N in reddest detector
            if self['SN'][-1,i] > 0:
                # Load up the full spectrum
                tmp_wave, all_flux, all_sky, all_ivar = np.ndarray(0), \
                    np.ndarray(0), np.ndarray(0), np.ndarray(0)
                # TODO -- Allow for Boxcar
                for det in range(self.ndet):
                    sobj = self.specobjs[self.sobj_idx[det][i]]
                    tmp_wave = np.concatenate((tmp_wave, sobj.OPT_WAVE))
                    all_flux = np.concatenate((all_flux, sobj.OPT_COUNTS))
                    all_sky = np.concatenate((all_sky, sobj.OPT_COUNTS_SKY))
                    all_ivar = np.concatenate((all_ivar, sobj.OPT_COUNTS_IVAR))
                
                # Massage
                fitwave  = self['fit_slope'][i]*tmp_wave + self['fit_b'][i]
                all_wave = tmp_wave - fitwave

                # TRIM ENDS
                all_wave=all_wave[5:-15]
                all_flux=all_flux[5:-15]
                all_ivar=all_ivar[5:-15]
                all_sky=all_sky[5:-15]

                # REMOVE CRAZY 500-SIGMA VALUES
                cmask = (all_sky > np.percentile(all_sky,0.1)) & (all_sky < np.percentile(all_sky,99.9))

                m=np.median(all_sky[cmask])
                s=np.std(all_sky[cmask])
                mm = (all_sky > 500.*s + m) | (all_sky < m-50.*s)
                all_sky[mm] = m
                all_ivar[mm] = 1e6
                if (np.sum(mm) > 10):
                    log.warning('Removing more than 10 pixels of data')
                
                _,diff,diff_err,_,_ = sky_em_residuals(all_wave, all_sky, all_ivar,
                                                       self.sky_table['Wave'])
                m = np.isfinite(diff)
                sky_mean = np.average(np.abs(diff[m]), weights = 1./diff_err[m]**2)
                resid_sky = np.append(resid_sky,sky_mean)

            else:
                resid_sky = np.append(resid_sky,-1)

        self['resid_sky'] = resid_sky

    def qa_plots(self, plot_dir:str, root:str):
        """Generate QA plots

        Args:
            plot_dir (str): Top-lvel folder for QA
                QA/ is generated beneath this, as needed
            root (str): Root for output files
        """

        # Generate QA folder as need be
        qa_dir = pathlib.Path(plot_dir) / 'QA'
        if not qa_dir.is_dir():
            qa_dir.mkdir(parents=True)
        
        '''
        # Slopes
        pdf2 = matplotlib.backends.backend_pdf.PdfPages(os.path.join(qa_dir, 'flex_slits_'+root+'.pdf'))
        plt.rcParams.update({'figure.max_open_warning': 0})
        for i in np.arange(0,self.nslits,1):

            if not np.all(self['SN'][:,i] > 0.):
                continue


            # SKY LINES FIRST
            r_sky_line, r_sky_diff,r_sky_ediff,r_los,r_elos = sky_em_residuals(hdu[r].data['OPT_WAVE'], \
                                                    hdu[r].data['OPT_COUNTS_SKY'],\
                                                    hdu[r].data['OPT_COUNTS_IVAR'])

            b_sky_line, b_sky_diff,b_sky_ediff,b_los,b_elos = sky_em_residuals(hdu[b].data['OPT_WAVE'], \
                                                    hdu[b].data['OPT_COUNTS_SKY'],\
                                                    hdu[b].data['OPT_COUNTS_IVAR'])

            fig, (ax1,ax2) = plt.subplots(1, 2,figsize=(20,4))
            ax1.plot(r_sky_line,r_sky_diff,'ro',alpha=0.8,label='Red chip: Sky Emission')
            ax1.plot(b_sky_line,b_sky_diff,'bo',alpha=0.8,label='Blue chip: Sky Emission')
            ax1.errorbar(b_sky_line,b_sky_diff,yerr=b_sky_ediff,fmt='none',ecolor='b',alpha=0.5)
            ax1.errorbar(r_sky_line,r_sky_diff,yerr=r_sky_ediff,fmt='none',ecolor='r',alpha=0.5)
            ax1.text(6320,0,'{}'.format(b),fontsize=11)
            ax1.text(8500,0,'{}'.format(r),fontsize=11)
            ax1.set_ylim(-0.45,0.45)

            x=np.arange(6000,9000,1)
            l1 = slits['fit_slope'][i]*x + slits['fit_b'][i]
            l2 = fslits['fit_slope'][i]*x + fslits['fit_b'][i]
            ax1.plot(x,l1,'-')
            ax1.plot(x,l2,'--')
            ax1.axhline(linewidth=1, color='grey',alpha=0.5)
            ax1.set_ylabel('Wavelength offset (AA)')
            ax1.set_xlabel('Wavelength (AA)')
            ax1.set_xlim(6300,9100)
            t = 'Sky Line Fits , resid = {:0.4f} AA, arc = {:0.2f}'.format(slits['resid_sky'][i],0.32*slits['rms_arc_r'][i])
            ax1.set_title(t)

            sky_diff  = np.concatenate((r_sky_diff,b_sky_diff),axis=None)
            sky_lines = np.concatenate((r_sky_line,b_sky_line),axis=None)
            sky_ediff = np.concatenate((r_sky_ediff,b_sky_ediff),axis=None)
            sky_los   = np.concatenate((r_los,b_los),axis=None)


            ax2.plot(r_sky_line,r_los,'ro',alpha=0.8,label='Red chip: Sky Emission')
            ax2.plot(b_sky_line,b_los,'bo',alpha=0.8,label='Blue chip: Sky Emission')
            ax2.errorbar(r_sky_line,r_los,yerr=r_elos,fmt='none',ecolor='r',alpha=0.5)
            ax2.errorbar(b_sky_line,b_los,yerr=b_elos,fmt='none',ecolor='b',alpha=0.5)
            ax2.axhline(fslits['fit_los'][i],linewidth=1, color='grey',alpha=0.5)

            ax2.set_title('Line widths')
            ax2.set_xlabel('Wavelength (AA)')
            ax2.set_ylim(0.3,0.8)
            ax2.set_xlim(6300,9100)

            pdf2.savefig()
        pdf2.close()
        plt.close('all')
        '''

        #########################################################################
        # CREATE FULL MASK FITS
        pdf = matplotlib.backends.backend_pdf.PdfPages(
            plot_dir+'QA/flex_mask_'+root+'.pdf')
        xslit = self['objra']
        yslit = self['objdec']
        t=2.

        mu =  np.median(self['indiv_fit_slope'])
        sd =  np.std(self['indiv_fit_slope'])
        mu2 =  np.median(self['indiv_fit_b'])
        sd2 =  np.std(self['indiv_fit_b'])
        mu3 =  np.median(self['indiv_fit_los'])
        sd3 =  np.std(self['indiv_fit_los'])

        # PLOT FITTED VALUES
        fig, (ax1,ax2,ax3) = plt.subplots(1, 3,figsize=(22,5))
    
        mm1=-0.00005
        mm2=0.00005
        print(mu-t*sd,mu+t*sd)
        ax1.scatter(xslit,yslit,c=self['indiv_fit_slope'],
                    cmap="cool",vmin = mm1,vmax=mm2 )# mu-t*sd,vmax=mu+t*sd)
        ax1.set_ylabel('Dec [deg]')
        ax1.set_xlabel('RA [deg]')
        ax1.set_title('Wave MEASURE: line slope')
        #cax, _ = matplotlib.colorbar.make_axes(ax1)
        #normalize = matplotlib.colors.Normalize(vmin = mu-t*sd,vmax=mu+t*sd)
        #cbar = matplotlib.colorbar.ColorbarBase(cax, cmap='cool',norm=normalize)


        ax2.scatter(xslit,yslit,c=self['indiv_fit_b'],cmap="summer",
                    vmin = mu2-t*sd2,vmax=mu2+t*sd2)
        ax2.set_ylabel('Dec [deg]')
        ax2.set_xlabel('RA [deg]')
        ax2.set_title('Wave MEASURE: line intercept')
        cax, _ = matplotlib.colorbar.make_axes(ax2)
        normalize = matplotlib.colors.Normalize(vmin = mu2-t*sd2,vmax=mu2+t*sd2)
        #cbar = matplotlib.colorbar.ColorbarBase(cax, cmap='summer',norm=normalize)


        ax3.scatter(xslit,yslit,c=self['indiv_fit_los'],cmap="cool",vmin = mu3-t*sd3,vmax=mu3+t*sd3)
        ax3.set_ylabel('Dec [deg]')
        ax3.set_xlabel('RA [deg]')
        ax3.set_title('Wave MEASURE: line width')
        cax, _ = matplotlib.colorbar.make_axes(ax3)
        normalize = matplotlib.colors.Normalize(vmin = mu3-t*sd3,vmax=mu3+t*sd3)
        #cbar = matplotlib.colorbar.ColorbarBase(cax, cmap='cool',norm=normalize)

        pdf.savefig()
        
        #######################
        # PLOT MEASURED VALUES
        fig, (ax1,ax2,ax3) = plt.subplots(1, 3,figsize=(22,5))
    
        ax1.scatter(xslit,yslit,c=self['fit_slope'],
                    cmap="cool",vmin = mu-t*sd,vmax=mu+t*sd)

        ax1.set_ylabel('Dec [deg]')
        ax1.set_xlabel('RA [deg]')
        ax1.set_title('Wave fit: line slope')
        cax, _ = matplotlib.colorbar.make_axes(ax1)
        normalize = matplotlib.colors.Normalize(vmin = mu-t*sd,vmax=mu+t*sd)
        #cbar = matplotlib.colorbar.ColorbarBase(cax, cmap='cool',norm=normalize)


        ax2.scatter(xslit,yslit,c=self['fit_b'],
                    cmap="summer",vmin = mu2-t*sd2,vmax=mu2+t*sd2)
        ax2.set_ylabel('Dec [deg]')
        ax2.set_xlabel('RA [deg]')
        ax2.set_title('Wave fit: line intercept')
        cax, _ = matplotlib.colorbar.make_axes(ax2)
        normalize = matplotlib.colors.Normalize(vmin = mu2-t*sd2,vmax=mu2+t*sd2)
        #cbar = matplotlib.colorbar.ColorbarBase(cax, cmap='summer',norm=normalize)


        ax3.scatter(xslit,yslit,c=self['fit_los'],
                    cmap="cool",vmin = mu3-t*sd3,vmax=mu3+t*sd3)
        ax3.set_ylabel('Dec [deg]')
        ax3.set_xlabel('RA [deg]')
        ax3.set_title('Wave fit: line width')
        cax, _ = matplotlib.colorbar.make_axes(ax3)
        normalize = matplotlib.colors.Normalize(vmin = mu3-t*sd3,vmax=mu3+t*sd3)
        #cbar = matplotlib.colorbar.ColorbarBase(cax, cmap='cool',norm=normalize)

        
        pdf.close()
