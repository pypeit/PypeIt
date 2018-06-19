# Module for the ScienceImage class
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np

#from importlib import reload
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
    This class will organize and run actions related to
    a Science or Standard star exposure

    Parameters
    ----------
    file_list : list
      List of raw files to produce the flat field
    spectrograph : str
    settings : dict-like
    tslits_dict : dict
      dict from TraceSlits class (e.g. slitpix)
    tilts : ndarray
      tilts from WaveTilts class
      used for sky subtraction and object finding
    det : int
    setup : str
    datasec_img : ndarray
      Identifies pixels to amplifiers
    bpm : ndarray
      Bad pixel mask
    maskslits : ndarray (bool)
      Specifies masked out slits
    pixlocn : ndarray
    objtype : str
      'science'
      'standard'
    fitstbl : Table
      Header info
    scidx : int
      Row in the fitstbl corresponding to the exposure

    Attributes
    ----------
    frametype : str
      Set to 'science'
    sciframe : ndarray
      Processed 2D frame
    rawvarframe : ndarray
      Variance generated without a sky (or object) model
    modelvarframe : ndarray
      Variance generated with a sky model
    finalvar : ndarray
      Final variance frame
    global_sky : ndarray
      Sky model across the slit/order
    skycorr_box : ndarray
      Local corrections to the sky model
    final_sky : ndarray
      Final sky model; may include 'local' corrections
    obj_model : ndarray
      Model of the object flux
    trcmask : ndarray
      Masks of objects for sky subtraction
    tracelist : list
      List of traces for objects in slits
    inst_name : str
      Short name of the spectrograph, e.g. KASTb
    target_name : str
      Parsed from the Header
    basename : str
      Combination of camera, target, and time
      e.g. J1217p3905_KASTb_2015May20T045733.56
    time : Time
      time object
    specobjs : list
      List of specobjs


    """
    # Keep order same as processimages (or else!)
    def __init__(self, file_list=[], spectrograph=None, settings=None,
                 tslits_dict=None, tilts=None, det=None, setup=None, datasec_img=None,
                 bpm=None, maskslits=None, pixlocn=None, objtype='science',
                 fitstbl=None, scidx=0):

        # Parameters unique to this Object
        self.det = det
        self.setup = setup
        self.tslits_dict = tslits_dict
        self.tilts = tilts
        self.maskslits = maskslits
        self.pixlocn = pixlocn
        self.objtype = objtype
        self.fitstbl = fitstbl
        self.scidx = scidx

        # Start us up
        processimages.ProcessImages.__init__(self, file_list, spectrograph=spectrograph,
                                             settings=settings, det=det,
                                             datasec_img=datasec_img,
                                             bpm=bpm)

        # Attributes (set after init)
        self.frametype = frametype

        # Key outputs/internals
        self.sciframe = None
        self.rawvarframe = None
        self.modelvarframe = None
        self.obj_model = None
        self.global_sky = None
        self.skycorr_box = None
        self.finalvar = None
        self.finalsky = None
        self.trcmask = None
        self.tracelist = []

        self.time = None
        self.inst_name = None
        self.target_name = None
        self.basename = None

        self.specobjs = []


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
        """
        Setup the basename (for file output mainly)
        and time objects (for heliocentric)

        Parameters
        ----------
        camera : str
          Taken from settings['mosaic']['camera']
        timeunit : str
          mjd

        Returns
        -------
        self.time : Time
        self.basename : str

        """
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
        # Time
        tval = Time(tbname, format='isot')#'%Y-%m-%dT%H:%M:%S.%f')
        dtime = datetime.datetime.strptime(tval.value, '%Y-%m-%dT%H:%M:%S.%f')
        self.time = tval
        # Basename
        self.inst_name = camera
        self.target_name = self.fitstbl['target'][self.scidx].replace(" ", "")
        self.basename = self.target_name+'_'+self.inst_name+'_'+ \
                         datetime.datetime.strftime(dtime, '%Y%b%dT') + \
                         tbname.split("T")[1].replace(':','')
        # Return
        return self.time, self.basename

    def _build_specobj(self):
        """
        Initialize the specobjs for all slits
        Key input is self.tracelist

        Wrapper to arspecobj.init_exp

        Returns
        -------
        self.specobjs : list

        """
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

    def _build_modelvar(self, skyframe=None, objframe=None):
        """
        Generate a model variance image using the sky model
        and (optional) object model

        Wrapper to arprocimg.variance_frame

        Parameters
        ----------
        skyframe : ndarray
          Sky model
        objframe : ndarray
          Object model

        Returns
        -------
        self.modelvarframe : ndarray
          Model variance image

        """
        if skyframe is None:
            skyframe = self.global_sky
        self.modelvarframe = arprocimg.variance_frame(
            self.datasec_img, self.det, self.sciframe, skyframe=skyframe,
            settings_det=self.settings['detector'], objframe=objframe)
        return self.modelvarframe

    def boxcar(self, mswave):
        """
        Perform boxcar extraction

        Wrapper to arextract.boxcar

        Parameters
        ----------
        mswave : ndarray
          Wavelength image

        Returns
        -------
        self.skycorr_box : ndarray
          Local corrections to the sky model

        Extractions are ingested within self.specobjs

        """
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
        """
        Perform optimal extraction using the 'original' PYPIT algorithm

        Wrapper to arextract.obj_profiles and arextract.optimal_extract

        Parameters
        ----------
        mswave : ndarray
          Wavelength image

        Returns
        -------
        self.obj_model : ndarray
          Model of the object flux; used for improving the variance estimate

        """
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
        """
        Perform the extraction
          Boxcar and optimal (as desired)

        Code flow:
          1. Instantiate the specobjs list -- This should be done when finding objects
          2. Boxcar extraction
            i. Update sky model
          3. Optimal extraction
            i. Iterative
            ii. With an update to the variance image
          4. One last update to the variance image

        Parameters
        ----------
        mswave : ndarray
          Wavelength image

        Returns
        -------
        self.specobjs
        self.finalvar
        self.finalsky

        """
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
            self.finalvar = self.rawvarframe

        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.specobjs, self.finalvar, self.finalsky

    def find_objects(self, doqa=False):
        """
        Find objects in the slits
        This is currently setup only for ARMS

        Wrapper to artrace.trace_objects_in_slit

        Parameters
        ----------
        doqa : bool, optional
          Generate QA?

        Returns
        -------
        self.tracelist : list
          List of tracedict objects genreated by object finding (to be Refactored)
        self.nobj : int
          Total number of objects in all slits

        """
        # Grab the 'best' variance frame
        varframe = self.rawvarframe

        # Prepare to loop on slits
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

        # Finish
        self.nobj = np.sum([tdict['nobj'] for tdict in self.tracelist if 'nobj' in tdict.keys()])
        if doqa:  # QA?
            obj_trace_qa(slf, sciframe, det, tracelist, root="object_trace", normalize=False)
        # Steps
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.tracelist, self.nobj

    def global_skysub(self, settings_skysub, use_tracemask=False):
        """
        Perform global sky subtraction, slit by slit

        Wrapper to arskysub.bg_subtraction_slit

        Parameters
        ----------
        settings_skysub : dict
          Guides sky subtraction algorithm(s)
        use_tracemask : bool, optional
          Mask objects (requires they were found previously)

        Returns
        -------
        self.global_sky : ndarray
        self.modelvarframe : ndarray
          Variance image using the sky model

        """
        # Prep
        self.global_sky = np.zeros_like(self.sciframe)
        gdslits = np.where(~self.maskslits)[0]
        # Grab the currently 'best' variance frame
        if self.modelvarframe is not None:
            # JFH does not quite approve of this..
            varframe = self.modelvarframe
        else:
            varframe = self.rawvarframe

        # Mask objects?
        if use_tracemask:
            tracemask = self._build_tracemask()
        else:
            tracemask = None

        # Loop on slits
        for slit in gdslits:
            msgs.info("Working on slit: {:d}".format(slit))
            slit_bgframe = arskysub.bg_subtraction_slit(
                slit+1, self.tslits_dict['slitpix'], self.tslits_dict['edge_mask'],
                self.sciframe, varframe, self.tilts,
                bpm=self.bpm, crmask=self.crmask, tracemask=tracemask)
            #slit_bgframe = arskysub.bg_subtraction_slit(self.tslits_dict, self.pixlocn,
            #                                            slit, self.tilts, self.sciframe,
            #                                            varframe, self.bpm,
            #                                            self.crmask, settings_skysub,
            #                                            tracemask=tracemask)
            #except ValueError:  # Should have been bspline..
            #    msgs.warn("B-spline sky subtraction failed.  Slit {:d} will no longer be processed..".format(slit))
            #    #msgs.warn("Continue if you wish..")
            #    debugger.set_trace()
            #    self.maskslits[slit] = True
            #else:
            self.global_sky += slit_bgframe

        # Build model variance
        msgs.info("Building model variance from the Sky frame")
        self.modelvarframe = self._build_modelvar()
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.global_sky, self.modelvarframe

    def _process(self, bias_subtract, pixel_flat, apply_gain=True, dnoise=0.):
        """ Process the image

        Wrapper to ProcessImages.process()

        Needed in part to set self.sciframe, although I could kludge it another way..

        Returns
        -------
        self.sciframe
        self.rawvarframe
        self.crmask

        """
        # Process
        self.sciframe = self.process(bias_subtract=bias_subtract, apply_gain=apply_gain, pixel_flat=pixel_flat)

        # Construct raw variance image
        self.rawvarframe = self.build_rawvarframe(dnoise=dnoise)

        # Build CR mask
        self.crmask = self.build_crmask()

        return self.sciframe, self.rawvarframe, self.crmask

    '''
    def _grab_varframe(self):
        """
        Simple algorithm to grab the currently 'best' variance image

        Order of priority --
          finalvar
          modelvarframe
          rawvarframe

        Returns
        -------
        varframe : ndarray

        """
        # Make the choice
        if self.finalvar is not None:
            varframe = self.finalvar
        elif self.modelvarframe is not None:
            varframe = self.modelvarframe
        else:
            varframe = self.rawvarframe
        # Return it
        return varframe
    '''

    def _build_tracemask(self):
        """
        Create a tracemask for the object
        Used in sky subtraction

        Returns
        -------
        self.trcmask

        """
        self.trcmask = np.zeros_like(self.sciframe)
        # Create a trace mask of the object
        for sl in range(len(self.tracelist)):
            if 'nobj' in self.tracelist[sl].keys():
                if self.tracelist[sl]['nobj'] > 0:
                    self.trcmask += self.tracelist[sl]['object'].sum(axis=2)
        self.trcmask[np.where(self.trcmask > 0.0)] = 1.0
        # Return
        return self.trcmask

    def show(self, attr, image=None, display='ginga'):
        """
        Show one of the internal images
          Should probably put some of these in ProcessImages

        Parameters
        ----------
        attr : str
          global -- Sky model (global)
          sci -- Processed science image
          rawvar -- Raw variance image
          modelvar -- Model variance image
          crmasked -- Science image with CRs set to 0
          skysub -- Science image with global sky subtracted
          image -- Input image
        display : str, optional
        image : ndarray, optional
          User supplied image to display

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

    def __repr__(self):
        txt = '<{:s}: nimg={:d}'.format(self.__class__.__name__,
                                        self.nfiles)
        if len(self.steps) > 0:
            txt+= ' steps: ['
            for step in self.steps:
                txt += '{:s}, '.format(step)
            txt = txt[:-2]+']'  # Trim the trailing comma
        txt += '>'
        return txt

