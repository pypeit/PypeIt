# Module for the ScienceImage class
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import os

from importlib import reload
import datetime

from astropy.time import Time

from pypit import msgs
from pypit import processimages
from pypit import arspecobj
from pypit.core import arprocimg
from pypit.core import arskysub
from pypit.core import arextract
from pypit import artrace
from pypit import ginga

from pypit import ardebug as debugger

# For out of PYPIT running
if msgs._debug is None:
    debug = debugger.init()
    debug['develop'] = True
    msgs.reset(debug=debug, verbosity=2)

# Does not need to be global, but I prefer it
frametype = 'science'


class ScienceImage(processimages.ProcessImages):
    """
    This class will generat the pixel-level FlatField
      The master() method returns the image

    Parameters
    ----------
    file_list : list
      List of raw files to produce the flat field
    spectrograph : str
    settings : dict-like
    msbias : ndarray or str or None
    tslits_dict : dict
      dict from TraceSlits class (e.g. slitpix)
    tilts : ndarray
      tilts from WaveTilts class
    det : int
    setup : str

    Attributes
    ----------
    frametype : str
      Set to 'pixelflat'
    mspixelflat : ndarray
      Stacked image
    mspixelflatnrm : ndarray
      Normalized flat
    extrap_slit
    msblaze : ndarray
      Blaze function fit to normalize
    blazeext :
    slit_profiles : ndarray
      Slit profile(s)
    self.ntckx : int
      Number of knots in the spatial dimension
    self.ntcky : int
      Number of knots in the spectral dimension

    """
    # Keep order same as processimages (or else!)
    def __init__(self, file_list=[], spectrograph=None, settings=None,
                 tslits_dict=None, tilts=None, det=None, setup=None, datasec_img=None,
                 bpm=None, maskslits=None, pixlocn=None, objtype='science',
                 fitstbl=None, scidx=0, mswave=None):

        # Parameters unique to this Object
        self.det = det
        self.setup = setup
        self.tslits_dict = tslits_dict
        self.tilts = tilts
        self.maskslits = maskslits
        self.fitstbl = fitstbl
        self.pixlocn = pixlocn
        self.objtype = objtype
        self.scidx = scidx

        # Start us up
        processimages.ProcessImages.__init__(self, file_list, spectrograph=spectrograph,
                                             settings=settings, det=det, datasec_img=datasec_img,
                                             bpm=bpm)

        # Attributes (set after init)
        self.frametype = frametype

        # Key outputs/internals
        self.sciframe = None
        self.rawvarframe = None
        self.modelvarframe = None
        self.global_sky = None
        self.tracelist = []


        # Settings
        # The copy allows up to update settings with user settings without changing the original
        if settings is None:
            # Defaults
            self.settings = processimages.default_settings()
        else:
            self.settings = settings.copy()

        # Child-specific Internals
        #    See ProcessImages
        self.crmask = None

    def init_time_names(self, camera, timeunit='mjd'):
        tbname = None
        try:
            if "T" in self.fitstbl['date'][self.scidx]:
                tbname = self.fitstbl['date'][self.scidx]
        except IndexError:
            debugger.set_trace()
        else:
            if tbname is None:
                if timeunit == "mjd":
                    # Not ideal, but convert MJD into a date+time
                    timval = Time(self.fitstbl['time'][self.scidx] / 24.0, scale='tt', format='mjd')
                    tbname = timval.isot
                else:
                    # Really not ideal... just append date and time
                    tbname = self.fitstbl['date'][self.scidx] + "T" + str(self.fitstbl['time'][self.scidx])
        tval = Time(tbname, format='isot')#'%Y-%m-%dT%H:%M:%S.%f')
        dtime = datetime.datetime.strptime(tval.value, '%Y-%m-%dT%H:%M:%S.%f')
        # Finish
        self._inst_name = camera
        self._target_name = self.fitstbl['target'][self.scidx].replace(" ", "")
        self._basename = self._target_name+'_'+self._inst_name+'_'+ \
                         datetime.datetime.strftime(dtime, '%Y%b%dT') + \
                         tbname.split("T")[1].replace(':','')
        # Save Time object
        self._time = tval
        # Return time
        return self._time, self._basename

    def _build_specobj(self):
        self.specobjs = arspecobj.init_exp(self.tslits_dict['lcen'],
                                           self.tslits_dict['rcen'],
                                           self.sciframe.shape,
                                           self.maskslits,
                                           self.det, self.scidx, self.fitstbl,
                                           self.tracelist, self.settings,
                                           objtype=self.objtype)
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.specobjs

    def boxcar(self, mswave):
        msgs.info("Performing boxcar extraction")
        self.skycorr_box = arextract.boxcar(self.specobjs, self.sciframe,
                                            self.modelvarframe, self.bpm,
                                            self.global_sky, self.crmask,
                                            self.tracelist, mswave,
                                            self.maskslits, self.tslits_dict['slitpix'])
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.skycorr_box

    def original_optimal(self, mswave):
        msgs.info("Attempting optimal extraction with model profile")
        # Profile
        arextract.obj_profiles(self.det, self.specobjs, self.sciframe-self.global_sky-self.skycorr_box,
                               self.modelvarframe, self.crmask, self.tracelist, self.tilts,
                               self.maskslits, self.tslits_dict['slitpix'], doqa=False)
        # Extract
        self.obj_model = arextract.optimal_extract(self.specobjs,
                                              self.sciframe-self.global_sky-self.skycorr_box,
                                              self.modelvarframe, self.crmask,
                                              self.tracelist, self.tilts,
                                              mswave, self.maskslits, self.tslits_dict['slitpix'])
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.obj_model

    def extraction(self, mswave):
        reload(arspecobj)
        reload(arextract)

        # Init specobjs
        #  Nested -- self.specobjs[slit][object]
        self.specobjs = self._build_specobj()

        # Boxcar -- Fills specobj.boxcar in place
        self.skycorr_box = self.boxcar(mswave)
        self.finalsky = self.global_sky+self.skycorr_box

        if self.objtype != 'standard':
            # Optimal (original recipe)
            self.obj_model = self.original_optimal(mswave)
            #
            msgs.info("Update model variance image (and trace?) and repeat")
            _ = self._build_modelvar(skyframe=self.finalsky, objframe=self.obj_model)
            self.obj_model = self.original_optimal(mswave)

            # Final variance image
            self.finalvar = self._build_modelvar(skyframe=self.finalsky, objframe=self.obj_model)
        else:
            self.finalvar = self._grab_varframe()

        # Step
        self.steps.append(inspect.stack()[0][3])

        # Return
        return self.specobjs, self.finalvar, self.finalsky

    def find_objects(self, doqa=False):
        reload(artrace)

        varframe = self._grab_varframe()

        nslit = len(self.maskslits)
        gdslits = np.where(~self.maskslits)[0]
        self.tracelist = []

        # Loop on good slits
        for slit in range(nslit):
            if slit not in gdslits:
                self.tracelist.append({})
                continue
            tlist = artrace.trace_objects_in_slit(self.det, slit, self.tslits_dict,
                                                  self.sciframe, self.global_sky,
                                                  varframe, self.crmask, self.settings)
            # Append
            self.tracelist += tlist

        # QA?
        if doqa: # and (not msgs._debug['no_qa']):
            obj_trace_qa(slf, sciframe, det, tracelist, root="object_trace", normalize=False)

        self.nobj = np.sum([tdict['nobj'] for tdict in self.tracelist if 'nobj' in tdict.keys()])
        return self.tracelist, self.nobj


    def global_skysub(self, settings_skysub, use_tracemask=False):
        reload(arskysub)

        self.global_sky = np.zeros_like(self.sciframe)
        gdslits = np.where(~self.maskslits)[0]

        varframe = self._grab_varframe()

        if use_tracemask:
            tracemask = self._build_tracemask()
        else:
            tracemask = None

        for slit in gdslits:
            msgs.info("Working on slit: {:d}".format(slit))
            # TODO -- Replace this try/except when a more stable b-spline is used..
            try:
                slit_bgframe = arskysub.bg_subtraction_slit(self.tslits_dict, self.pixlocn,
                                                            slit, self.tilts, self.sciframe,
                                                            varframe, self.bpm,
                                                            self.crmask, settings_skysub,
                                                            tracemask=tracemask)
            except ValueError:  # Should have been bspline..
                msgs.warn("B-spline sky subtraction failed.  Slit {:d} will no longer be processed..".format(slit))
                #msgs.warn("Continue if you wish..")
                debugger.set_trace()
                self.maskslits[slit] = True
            else:
                self.global_sky += slit_bgframe

        # Build model variance
        msgs.info("Building model variance from the Sky frame")
        self.modelvarframe = self._build_modelvar()
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.global_sky, self.modelvarframe

    def _build_modelvar(self, skyframe=None, objframe=None):
        if skyframe is None:
            skyframe = self.global_sky
        self.modelvarframe = arprocimg.variance_frame(
            self.datasec_img, self.det, self.sciframe, skyframe=skyframe,
            settings_det=self.settings['detector'], objframe=objframe)
        return self.modelvarframe

    def _process(self, bias_subtract, pixel_flat, apply_gain=True, dnoise=0.):
        """ Process the image
        Wrapper to ProcessImages.process()

        Needed in part to set self.sciframe, although I could kludge it another way..
        """
        self.sciframe = self.process(bias_subtract=bias_subtract, apply_gain=apply_gain, pixel_flat=pixel_flat)

        # Variance image
        self.rawvarframe = self.build_rawvarframe(dnoise=dnoise)

        # CR mask
        self.crmask = self.build_crmask()

        return self.sciframe, self.rawvarframe, self.crmask

    def _grab_varframe(self):
        #
        if self.modelvarframe is None:
            msgs.info("Using raw variance frame for sky sub")
            varframe = self.rawvarframe
        else:
            msgs.info("Using model variance frame for sky sub")
            varframe = self.modelvarframe
        return varframe

    def _build_tracemask(self):
        # Create a trace mask of the object
        self.trcmask = np.zeros_like(self.sciframe)
        for sl in range(len(self.tracelist)):
            if 'nobj' in self.tracelist[sl].keys():
                if self.tracelist[sl]['nobj'] > 0:
                    self.trcmask += self.tracelist[sl]['object'].sum(axis=2)
        self.trcmask[np.where(self.trcmask > 0.0)] = 1.0
        return self.trcmask


    def show(self, attr, image=None, display='ginga'):
        """
        Show one of the internal images
          Should probably put some of these in ProcessImages

        Parameters
        ----------
        attr : str
          bgframe -- Show the combined flat image, unnormalized
          norm -- Show the combined normalized flat image
        display : str, optional

        Returns
        -------

        """
        if attr == 'global':
            if self.global_sky is not None:
                ginga.show_image(self.global_sky, chname='Global')
        elif attr == 'sci':
            if self.sciframe is not None:
                ginga.show_image(self.sciframe)
        elif attr == 'rawvar':
            if self.rawvarframe is not None:
                ginga.show_image(self.rawvarframe)
        elif attr == 'modelvar':
            if self.modelvarframe is not None:
                ginga.show_image(self.modelvarframe)
        elif attr == 'crmasked':
            ginga.show_image(self.sciframe*(1-self.crmask), chname='CR')
        elif attr == 'skysub':
            ginga.show_image(self.sciframe-self.global_sky, chname='SkySub')
        elif attr == 'image':
            ginga.show_image(image)
        else:
            msgs.warn("Not an option for show")
