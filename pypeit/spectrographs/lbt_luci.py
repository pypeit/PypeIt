""" Module for LBT/LUCI specific codes
"""
import numpy as np


from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph



class LBTLUCISpectrograph(spectrograph.Spectrograph):
    """
    Class to handle LBT/LUCI specific code
    """
    def __init__(self):
        # Get it started
        super(LBTLUCISpectrograph, self).__init__()
        self.spectrograph = 'lbt_luci'
        self.telescope = telescopes.LBTTelescopePar()
        self.timeunit = 'isot'

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for  LBT/LUCI reductions.

        OLD CODE from LBT MODS
        """
        par = pypeitpar.PypeItPar()
        # Scienceimage default parameters
        par['reduce'] = pypeitpar.ReducePar()
        # Always flux calibrate, starting with default parameters
        par['fluxcalib'] = pypeitpar.FluxCalibratePar()
        # Always correct for flexure, starting with default parameters
        par['flexure'] = pypeitpar.FlexurePar()
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['pixelflatframe']['exprng'] = [0, None]
        par['calibrations']['traceframe']['exprng'] = [0, None]
        par['calibrations']['arcframe']['exprng'] = [None, 60]
        par['calibrations']['standardframe']['exprng'] = [1, 200]
        par['scienceframe']['exprng'] = [200, None]
        return par

    def init_meta(self):
        """
        Generate the meta data dict
        Note that the children can add to this

        Returns:
            self.meta: dict (generated in place)

        """
        meta = {}
        # Required (core)
        meta['ra'] = dict(ext=0, card='OBJRA')
        meta['dec'] = dict(ext=0, card='OBJDEC')
        meta['target'] = dict(ext=0, card='OBJECT')
        meta['decker'] = dict(ext=0, card='MASKID')
        meta['binning'] = dict(ext=0, card=None, default='1,1')
        meta['filter1'] = dict(ext=0, card='FILTERS')
        meta['idname'] = dict(card=None, compound=True)
        meta['mjd'] = dict(ext=0, card='MJD-OBS')
        meta['exptime'] = dict(ext=0, card='EXPTIME')
        meta['airmass'] = dict(ext=0, card='AIRMASS')
        meta['dispname'] = dict(ext=0, card='GRATNAME')


        self.meta = meta

    def compound_meta(self, headarr, meta_key):
        """

        Args:
            headarr: list
            meta_key: str

        Returns:
            value

        """

        # Populate the idname based on the header information of LUCI
        # This is an implicit way of pre-typing without adding too many
        # variables to the self.meta.
        if meta_key == 'idname':
            targetname = (headarr[0].get('OBJECT'))
            dispname = (headarr[0].get('GRATNAME'))
            calib_unit = (headarr[0].get('CALIB'))
            filter1 = (headarr[0].get('FILTER1'))
            filter2 = (headarr[0].get('FILTER2'))
            lamp1 = (headarr[0].get('LAMP1'))
            lamp2 = (headarr[0].get('LAMP2'))
            lamp3 = (headarr[0].get('LAMP3'))
            lamp4 = (headarr[0].get('LAMP4'))
            lamp5 = (headarr[0].get('LAMP5'))
            lamp6 = (headarr[0].get('LAMP6'))

            # object frame -> will be typed as science
            # This currently includes sky flats, science and standard images
            # We will guess standards using the beginning of their names.
            if ((dispname != 'Mirror') and
                (calib_unit == False) and
                (lamp1 == False) and
                (lamp2 == False) and
                (lamp3 == False) and
                (lamp4 == False) and
                (lamp5 == False) and
                (lamp6 == False)):
                if (targetname[:3] == 'HIP' or
                    targetname[:2] == 'HD' or
                    targetname[:5] == 'Feige'):
                    return 'standard'
                else:
                    return 'object'
            # flat frames -> will be typed as pixelflat, trace
            elif ((calib_unit == True) and
                  ((lamp4 == True) or
                   (lamp5 == True) or
                   (lamp6 == True))):
                return 'flat'
            # arcs -> will be typed as arc, tilt
            elif ((dispname != 'Mirror') and
                  (calib_unit == True) and
                  ((lamp1 == True) or
                   (lamp2 == True) or
                   (lamp3 == True))):
                return 'arc'
            # pixelflat off -> will be typed as bias
            elif ((dispname != 'Mirror') and
                (calib_unit == True) and
                (lamp1 == False) and
                (lamp2 == False) and
                (lamp3 == False) and
                (lamp4 == False) and
                (lamp5 == False) and
                (lamp6 == False) and
                (filter1 != 'blind') and
                (filter2 != 'blind')):
                return 'flat_off'
            # darks -> will not be typed currently
            elif ((filter1 == 'blind') or
                  (filter2 == 'blind')):
                return 'dark'

        else:
            msgs.error("Not ready for this compound meta")

    # Uses parent metadata keys

    def configuration_keys(self):
        return ['decker', 'dispname']

    def pypeit_file_keys(self):
        pypeit_keys = super(LBTLUCISpectrograph, self).pypeit_file_keys()
        pypeit_keys += ['calib', 'comb_id', 'bkg_id', 'idname']
        return pypeit_keys

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        # ATTENTION: Standards have to be added manually for LUCI because
        # there is not unique flag that allows to distinguish between targets
        # and standards
        if ftype in ['science']:
            return good_exp & (fitstbl['idname'] == 'object')
        if ftype in ['standard']:
            return good_exp & (fitstbl['idname'] == 'standard')
        if ftype == 'bias':
            # for NIR data we type off lamp flats as biases
            return good_exp & (fitstbl['idname'] == 'flat_off')
        if ftype in ['pixelflat', 'trace']:
            # Flats and trace frames are typed together
            return good_exp & (fitstbl['idname'] == 'flat')
        if ftype in ['dark']:
            # NOT Typing dark frames
            # return np.zeros(len(fitstbl), dtype=bool)
            # for testing dark typing uncommen the following line and comment
            # out the previous line
            return good_exp & (fitstbl['idname'] == 'dark')
        if ftype in ['arc', 'tilt']:
            return (good_exp & ((fitstbl['idname'] == 'object') |
                    (fitstbl['idname'] == 'arc')))



        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)


class LBTLUCI1Spectrograph(LBTLUCISpectrograph):
    """
    Child to handle LBT/LUCI1 specific code
    """
    def __init__(self):
        # Get it started
        super(LBTLUCI1Spectrograph, self).__init__()
        self.spectrograph = 'lbt_luci1'
        self.camera = 'LUCI1'
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         = 0,
                            specaxis        = 1,
                            specflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.25,
                    # Dark current nominally is < 360 electrons per hours
                    # but the dark subtraction will effectively bring this to 0
                            darkcurr        = 0.0,
                    # Saturation is 55000, but will be set to dummy value for
                    # now
                            saturation      = 1e+8,
                    # NIR detectors are non-linear even in lower percentages
                    # of the full well, thus for precision measurements one
                    # should take into account a general non-linearity
                    # correction.
                            nonlinear       = 0.80,
                    # In fact there are 32 amplifiers, which gain and ronoise
                    # are extremely similar to each other, thus it will be
                    # mimicked as 1
                            numamplifiers   = 1,
                    # The readout noise for LUCI are different for
                    # different readout modes. The LIR mode values will be
                    # commented in and the MER values will be uncommented:
                            gain= 2.0,
                            # ronoise= 10.3,
                            ronoise         = 4.61,
                            datasec='[5:2044,5:2044]',
                    # For Luci the first 4 pixels on each side can
                    # technically be used for as a biassec. This is not
                    # included here.
                            oscansec= '[5:2044,1:4]',
                            suffix          = '_luci1'
                            )]
        self.numhead = 1


    def default_pypeit_par(self):
        """
         Set default parameters for Keck/MOSFIRE
         """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'lbt_luci1'

        # for key in par['calibrations'].keys():
        #     print(key)
        #     par['calibrations'][key]['process']['overscan'] = 'none'

        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths'][
            'rms_threshold'] = 0.20  # 0.20  # Might be grating dependent..
        par['calibrations']['wavelengths']['sigdetect'] = 5.0
        par['calibrations']['wavelengths']['fwhm'] = 5.0
        par['calibrations']['wavelengths']['n_final'] = 4
        par['calibrations']['wavelengths']['lamps'] = ['OH_NIRES']
        par['calibrations']['wavelengths']['nonlinear_counts'] = \
        self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['method'] = 'holy-grail'
        # Reidentification parameters
        par['calibrations']['slitedges']['edge_thresh'] = 300.
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Flats
        par['calibrations']['flatfield']['illumflatten'] = True

        # Extraction
        # Model full slit currently turned on
        par['reduce']['extraction']['model_full_slit'] = False
        # Tailored profile nsigma parameter for the standard, trying 100 (30
        # was standard
        par['reduce']['extraction']['std_prof_nsigma'] = 100.
        # Do not perform global sky subtraction for standard stars
        par['reduce']['skysub']['global_sky_std'] = True
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['extraction']['sn_gauss'] = 4.0

        # Flexure
        par['flexure']['method'] = 'skip'

        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] = 'nothing'
        # par['scienceframe']['process']['satpix'] = 'reject'

        par['scienceframe']['process']['overscan'] = 'none'
        # par['standardframe']['process']['overscan'] = 'none'

        return par

    def check_headers(self, headers):
        """
        Check headers match expectations for an LBT LUCI1 exposure.

        See also
        :func:`pypeit.spectrographs.spectrograph.Spectrograph.check_headers`.

        Args:
            headers (list):
                A list of headers read from a fits file
        """
        expected_values = { '0.INSTRUME': 'LUCI1',
                            '0.NAXIS': 2 }
        super(LBTLUCI1Spectrograph, self).check_headers(headers,
                                                              expected_values=expected_values)

class LBTLUCI2Spectrograph(LBTLUCISpectrograph):
    """
    Child to handle LBT/LUCI2 specific code
    """
    def __init__(self):
        # Get it started
        super(LBTLUCI2Spectrograph, self).__init__()
        self.spectrograph = 'lbt_luci2'
        self.camera = 'LUCI2'
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         = 0,
                            specaxis        = 1,
                            specflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.25,
                            darkcurr        = 0.0,
                            # Saturation is 55000, but will be set to dummy value for
                            # now
                            saturation=1e+8,
                            nonlinear       = 0.80,
                            numamplifiers   = 1,
                            gain            = 2.0,
                            # ronoise         = 10.0,
                            ronoise         = 4.47,
                            datasec='[5:2044,5:2044]',
                            oscansec='[5:2044,1:4]',
                            suffix          = '_luci2'
                            )]
        self.numhead = 1


    def default_pypeit_par(self):
        """
         Set default parameters for Keck/MOSFIRE
         """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'lbt_luci2'

        # for key in par['calibrations'].keys():
        #     print(key)
        #     par['calibrations'][key]['process']['overscan'] = 'none'

        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths'][
            'rms_threshold'] = 0.20  # 0.20  # Might be grating dependent..
        par['calibrations']['wavelengths']['sigdetect'] = 5.0
        par['calibrations']['wavelengths']['fwhm'] = 5.0
        par['calibrations']['wavelengths']['n_final'] = 4
        par['calibrations']['wavelengths']['lamps'] = ['OH_NIRES']
        par['calibrations']['wavelengths']['nonlinear_counts'] = \
            self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['method'] = 'holy-grail'


        par['calibrations']['slitedges']['edge_thresh'] = 300
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        par['calibrations']['slitedges']['fit_order'] = 8

        # Flats
        par['calibrations']['flatfield']['illumflatten'] = True
        # par['calibration']['flatfield']['tweak_slits'] = False

        # Extraction
        # Model full slit currently turned on
        par['reduce']['extraction']['model_full_slit'] = True
        # Tailored profile nsigma parameter for the standard
        par['reduce']['extraction']['std_prof_nsigma'] = 100.
        # Do not perform global sky subtraction for standard stars
        par['reduce']['skysub']['global_sky_std'] = False
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['extraction']['sn_gauss'] = 4.0

        # Flexure
        par['flexure']['method'] = 'skip'

        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] = 'nothing'
        # par['scienceframe']['process']['satpix'] = 'reject'

        par['scienceframe']['process']['overscan'] = 'none'
        # par['standardframe']['process']['overscan'] = 'none'

        return par



    def check_headers(self, headers):
        """
        Check headers match expectations for an LBT LUCI1 exposure.

        See also
        :func:`pypeit.spectrographs.spectrograph.Spectrograph.check_headers`.

        Args:
            headers (list):
                A list of headers read from a fits file
        """
        expected_values = { '0.INSTRUME': 'LUCI2',
                            '0.NAXIS': 2 }
        super(LBTLUCI1Spectrograph, self).check_headers(headers,
                                                              expected_values=expected_values)
