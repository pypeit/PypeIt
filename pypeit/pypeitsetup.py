"""
Class for organizing PypeIt setup
"""
import os
import inspect
import datetime
import numpy as np

from astropy.table import hstack, Table

from pypeit import msgs
from pypeit.metadata import PypeItMetaData

from pypeit.par import PypeItPar
from pypeit.par.util import parse_pypeit_file, make_pypeit_file
from pypeit.spectrographs.util import load_spectrograph

from pypeit import debugger

from spit import classifier
from spit import preprocess as spit_p
from spit import labels as spit_l

from tensorflow import keras
from PIL import Image

os.environ['SPIT_MODEL_PATH'] = 'C:/Users/kjk11/Desktop/'


class PypeItSetup(object):
    """
    Prepare for a pypeit run.

    .. todo::
        - This is now mostly a wrapper for PypeItMetaData.  Should we
          remove this class, or merge PypeItSetup and PypeItMetaData.

    The main deliverables are the set of parameters used for pypeit's
    algorithms (:attr:`par`), an :obj:`astropy.table.Table` with the
    details of the files to be reduced (:attr:`fitstbl`), and a
    dictionary with the list of instrument setups.

    Args:
        file_list (list):
            A list of strings with the full path to each file to be
            reduced.
        frametype (:obj:`dict`, optional):
            A dictionary that associates the name of the file (just the
            fits file name without the full path) to a specific frame
            type (e.g., arc, bias, etc.).  If None, this is determined
            by the :func:`get_frame_types` method.
        usrdata (:obj:`astropy.table.Table`, optional):
            A user provided set of data used to supplement or overwrite
            metadata read from the file headers.  The table must have a
            `filename` column that is used to match to the metadata
            table generated within PypeIt.
        setups (:obj:`list`, optional):
            A list of setups that each file can be associated with.  If
            None, all files are expected to be for a single setup.
        cfg_lines (:obj:`list`, optional):
            A list of strings that provide a set of user-defined
            parameters for executing pypeit.  These are the lines of a
            configuration file.  See the documentation for the
            `configobj`_ package.  One of the user-level inputs should
            be the spectrograph that provided the data to be reduced.
            One can get the list of spectrographs currently served by
            running::

                from pypeit.spectrographs.util import valid_spectrographs
                print(valid_spectrographs())

            To use all the default parameters when reducing data from a
            given spectrogaph, you can set `cfg_lines = None`, but you
            then *must* provide the `spectrograph_name`.
        spectrograph_name (:obj:`str`, optional):
            If not providing a list of configuration lines
            (`cfg_lines`), this sets the spectrograph.  The spectrograph
            defined in `cfg_lines` takes precedent over anything
            provided by this argument.
        pypeit_file (:obj:`str`, optional):
            The name of the pypeit file used to instantiate the
            reduction.  This can be None, and will lead to default names
            for output files (TODO: Give list).  Setting `pypeit_file`
            here *only sets the name of the file*.  To instantiate a
            `PypitSetup` object directly from a pypeit file (i.e. by
            reading the file), use the :func:`from_pypeit_file` method;
            i.e.::

                setup = PypitSetup.from_pypeit_file('myfile.pypeit')

    Attributes:
        file_list (list):
            See description of class argument.
        frametype (dict):
            See description of class argument.
        setups (list):
            See description of class argument.
        pypeit_file (str):
            See description of class argument.
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            An instance of the `Spectograph` class used throughout the
            reduction procedures.
        par (:class:`pypeit.par.pypeitpar.PypitPar`):
            An instance of the `PypitPar` class that provides the
            parameters to all the algorthms that pypeit uses to reduce
            the data.
        fitstbl (:class:`pypeit.metadata.PypeItMetaData`):
            A `Table` that provides the salient metadata for the fits
            files to be reduced.
        setup_dict (dict):
            The dictionary with the list of instrument setups.
        steps (list):
            The steps run to provide the pypeit setup.

    .. _configobj: http://configobj.readthedocs.io/en/latest/
    """

    def __init__(self, file_list, path=None, frametype=None, usrdata=None, setups=None,
                 cfg_lines=None, spectrograph_name=None, pypeit_file=None):

        # The provided list of files cannot be None
        if file_list is None or len(file_list) == 0:
            msgs.error('Must provide a list of files to be reduced!')

        # Save input
        self.file_list = file_list
        self.path = os.getcwd() if path is None else path
        self.frametype = frametype
        self.usrdata = usrdata
        self.setups = setups
        self.pypeit_file = pypeit_file
        self.user_cfg = cfg_lines

        # Determine the spectrograph name
        _spectrograph_name = spectrograph_name if cfg_lines is None \
            else PypeItPar.from_cfg_lines(merge_with=cfg_lines)['rdx']['spectrograph']

        # Cannot proceed without spectrograph name
        if _spectrograph_name is None:
            msgs.error(
                'Must provide spectrograph name directly or using configuration lines.')

        # Instantiate the spectrograph
        self.spectrograph = load_spectrograph(_spectrograph_name)

        # Get the spectrograph specific configuration to be merged with
        # the user modifications.
        spectrograph_cfg_lines = self.spectrograph.default_pypeit_par().to_config()

        # Instantiate the pypeit parameters.  The user input
        # configuration (cfg_lines) can be None.
        self.par = PypeItPar.from_cfg_lines(
            cfg_lines=spectrograph_cfg_lines, merge_with=cfg_lines)

        # Prepare internals for execution
        self.fitstbl = None
        self.setup_dict = None
        self.steps = []

    @classmethod
    def from_pypeit_file(cls, filename):
        """
        Instantiate the :class:`PypeitSetup` object using a pypeit file.

        Args:
            filename (str):
                Name of the pypeit file to read.  Pypit files have a
                specific set of valid formats. A description can be
                found `here`_ (include doc link).

        Returns:
            :class:`PypeitSetup`: The instance of the class.
        """
        cfg_lines, data_files, frametype, usrdata, setups = parse_pypeit_file(
            filename)
        return cls(data_files, frametype=frametype, usrdata=usrdata, setups=setups,
                   cfg_lines=cfg_lines, pypeit_file=filename)

    @classmethod
    def from_file_root(cls, root, spectrograph, extension='.fits', output_path=None):
        """
        Instantiate the :class:`PypeItSetup` object by providing a file
        root.

        This is based on first writing a vanilla PypeIt file for the
        provided spectrograph and extension to a file in the provided
        path.

        Args:
            root (:obj:`str`):
                The root path to all the files for PypeIt to reduce.
                This should be everything up to the wild-card before the
                file extension to use to find the relevant files.  The
                root itself can have wild cards to read through multiple
                directories.
            spectrograph (:obj:`str`):
                The PypeIt name of the spectrograph used to take the
                observations.  This should be one of the available
                options in
                :func:`pypeit.spectrographs.valid_spectrographs`.
            extension (:obj:`str`, optional):
                The extension common to all the fits files to reduce.
                Default is '.fits', meaning anything with `root*.fits*`
                will be be included.
            output_path (:obj:`str`, optional):
                Path to use for the output.  If None, the default is
                './setup_files'.  If the path doesn't yet exist, it is
                created.

        Returns:
            :class:`PypitSetup`: The instance of the class.
        """
        # Set the output directory
        outdir = os.path.join(
            os.getcwd(), 'setup_files') if output_path is None else output_path
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        # Set the output file name
        date = str(datetime.date.today().strftime('%Y-%b-%d'))
        pypeit_file = os.path.join(
            outdir, '{0}_{1}.pypeit'.format(spectrograph, date))
        msgs.info(
            'A vanilla pypeit file will be written to: {0}'.format(pypeit_file))

        # Generate the pypeit file
        cls.vanilla_pypeit_file(
            pypeit_file, root, spectrograph, extension=extension)

        # Now setup PypeIt using that file
        return cls.from_pypeit_file(pypeit_file)

    @staticmethod
    def vanilla_pypeit_file(pypeit_file, root, spectrograph, extension='.fits'):
        """
        Write a vanilla PypeIt file.

        Args:
            pypeit_file (str):
              Name of PypeIt file
            root (str):
            spectrograph (str):
              Name of spectrograph
            extension (str, optional):
              File extension

        Returns:

        """
        # Generate
        dfname = os.path.join(root, '*{0}*'.format(extension)) \
            if os.path.isdir(root) else '{0}*{1}*'.format(root, extension)
        # configuration lines
        cfg_lines = ['[rdx]']
        cfg_lines += ['    spectrograph = {0}'.format(spectrograph)]
#        cfg_lines += ['    sortroot = {0}'.format(root)] # check if this needs to be uncommented
        make_pypeit_file(pypeit_file, spectrograph, [
                         dfname], cfg_lines=cfg_lines, setup_mode=True)

    @property
    def nfiles(self):
        """The number of files to reduce."""
        if self.fitstbl is None:
            msgs.warn('No fits files have been read!')
        return 0 if self.fitstbl is None else len(self.fitstbl)

    def __repr__(self):
        return '<{:s}: nfiles={:d}>'.format(self.__class__.__name__, self.nfiles)

    def build_fitstbl(self, strict=True):
        """
        Construct the table with metadata for the frames to reduce.

        Largely a wrapper for :func:`pypeit.core.load.create_fitstbl`.

        Args:
            strict (:obj:`bool`, optional):
                Function will fault if :func:`fits.getheader` fails to
                read the headers of any of the files in
                :attr:`file_list`.  Set to False to only report a
                warning and continue.

        Returns:
            :obj:`astropy.table.Table`: Table with the metadata for each
            fits file to reduce.  Note this is different from
            :attr:`fitstbl` which is a :obj:`PypeItMetaData` object
        """
        # Build and sort the table
        self.fitstbl = PypeItMetaData(self.spectrograph, par=self.par, file_list=self.file_list,  # CHANGE THIS BACK TO files WHEN YOU ARE DONE WITH CODING EVERYTHING ELSE
                                      usrdata=self.usrdata, strict=strict)
        # Sort by the time
        if 'time' in self.fitstbl.keys():
            self.fitstbl.sort('time')

        # Add this to the completed steps
        self.steps.append(inspect.stack()[0][3])

        # Return the table
        return self.fitstbl.table

    def get_frame_types(self, flag_unknown=False, use_header_id=False):
        """
        Include the frame types in the metadata table.

        This is mainly a wrapper for
        :func:`PypeItMetaData.get_frame_types`.

        .. warning::

            Because this merges the frame types with the existing
            :attr:`fitstbl` this should only be run once.

        Args:
            flag_unknown (:obj:`bool`, optional):
                Allow for frames to have unknown types instead of
                crashing.  This should be True for initial setup and
                False otherwise.  Passed to get_frame_types()
            use_header_id (bool, optional):
                Passed to get_frame_types()

        """
        # Use PypeItMetaData methods to get the frame types
        _ = self.fitstbl.get_frame_types(flag_unknown=flag_unknown, user=self.frametype,
                                         useIDname=use_header_id)
        # Include finished processing step
        self.steps.append(inspect.stack()[0][3])

    def load_metadata(self, fits_file):
        """
        Load the fitstbl from disk (a binary FITS table)

        Args:
            fits_file (str):  Name of PypeItMetaData file

        Returns:
            obj:`PypeItMetaData`: The so-called fitstbl

        """
        self.fitstbl = PypeItMetaData(
            self.spectrograph, self.par, data=Table.read(fits_file))
        msgs.info("Loaded fitstbl from {:s}".format(fits_file))
        return self.fitstbl.table

    def write_metadata(self, sort_dir=None, ofile=None):
        """
        Write the :class:`astropy.table.Table` object in :attr:`fitstbl`
        to a file.

        If an output file is provided, the file is used.  If that file
        name inclues '.fits', the output will be a fits file; otherwise
        the format is ascii.fixed_width.

        If no output file, the default is an ascii table with an '.lst'
        extension.  The root name of the file is either the spectrograph
        name or the root name of the pypeit file, if the latter exists.

        If a `sort_dir` is provided, the directory of the nominal output
        file is changed to this directory.

        Args:
            sort_dir (:obj:`str`, optional):
                The full root of the name for the metadata table
                ('.lst') file.
            ofile (:obj:`str, optional):
                The name of the file to write.  See description above.
        """
        if ofile is None:
            ofile = self.spectrograph.spectrograph + '.lst' if self.pypeit_file is None \
                else self.pypeit_file.replace('.pypeit', '.lst')
            if sort_dir is not None:
                ofile = os.path.join(sort_dir, os.path.split(ofile)[1])

        format = None if '.fits' in ofile else 'ascii.fixed_width'
        self.fitstbl.write(ofile,
                           #columns=None if format is None else self.spectrograph.pypeit_file_keys(),
                           format=format, overwrite=True)

    def run(self, setup_only=False, calibration_check=False,
            use_header_id=False, sort_dir=None, write_bkg_pairs=False):
        """
        Once instantiated, this is the main method used to construct the
        object.

        The code flow is as follows::
            - Build the fitstbl from an input file_list (optional)
            - Type the files (bias, arc, etc.)
            - Match calibration files to the science files
            - Generate the setup_dict
                - Write group info to disk
                - Write calib info to disk (if main run)

        It is expected that a user will run this three times if they're
        being careful.  Once with `setup_only=True` to confirm the
        images are properly typed and grouped together for calibration.
        A second time with `calibration_check=True` to confirm the
        appropriate calibrations frames are available.  And a third time
        to do the actual setup before proceeding with the reductions.

        Args:
            setup_only (:obj:`bool`, optional):
                Only this setup will be performed.  Pypit is expected to
                execute in a way that ends after this class is fully
                instantiated such that the user can inspect the results
                before proceeding.  This has the effect of providing
                more output describing the success of the setup and how
                to proceed, and provides warnings (instead of errors)
                for issues that may cause the reduction itself to fail.
            calibration_check (obj:`bool`, optional):
                Only check that the calibration frames are appropriately
                setup and exist on disk.  Pypit is expected to execute
                in a way that ends after this class is fully
                instantiated such that the user can inspect the results
                before proceeding. 
            use_header_id (:obj:`bool`, optional):
                Allow setup to use the frame types drawn from single,
                instrument-specific header keywords set to `idname` in
                the metadata table (:attr:`fitstbl`).
            sort_dir (:obj:`str`, optional):
                The directory to put the '.sorted' file.

        Returns:
            :class:`pypeit.par.pypeitpar.PypeItPar`,
            :class:`pypeit.spectrographs.spectrograph.Spectrograph`,
            :class:`astropy.table.Table`: Returns the attributes
            :attr:`par`, :attr:`spectrograph`, :attr:`fitstbl`
            If running with `setup_only` or
            `calibrations_check`, these are all returned as `None`
            values.
        """
        # Kludge
        pypeit_file = '' if self.pypeit_file is None else self.pypeit_file

        # Build fitstbl
        if self.fitstbl is None:
            self.build_fitstbl(strict=not setup_only)  # , bkg_pairs=bkg_pairs)

        # File typing
        a = self.get_frame_types(flag_unknown=setup_only or calibration_check,
                                 use_header_id=use_header_id)

        # frame types aren't being printed, you have to figure out where you want prediction to be stored, and if you want filepath to spit to be parameter for prediction method

        # spit classification
        spit_frame_types, spit_framebits = self.spit_classify()

        # compare spit to rules-based
        for index, framebit in enumerate(self.fitstbl['framebit']):
                # check if they're different
            if spit_framebits[index] != framebit:
                # if so, notify user of difference, and set to spit's values
                # add at which files they differ so the user can check
                msgs.warn("SPIT and rules-based framebits differ at %s" %
                          (self.fitstbl['filename'][0]))
                self.fitstbl['framebit'][index] = spit_framebits[index]
                self.fitstbl['frametype'][index] = spit_frame_types[index]

        # Determine the configurations and assign each frame to the
        # specified configuration
        ignore_frames = ['bias', 'dark']
        cfgs = self.fitstbl.unique_configurations(ignore_frames=ignore_frames)
        self.fitstbl.set_configurations(cfgs, ignore_frames=ignore_frames)

        # Assign frames to calibration groups
        self.fitstbl.set_calibration_groups(global_frames=['bias', 'dark'])

        # Set default comb_id
        self.fitstbl.set_combination_groups()

        # Assign science IDs based on the calibrations groups (to be
        # deprecated)
        self.fitstbl['failures'] = False                    # TODO: placeholder

        if setup_only:
            # Collate all matching files and write .sorted Table (on pypeit_setup only)
            sorted_file = self.spectrograph.spectrograph + '.sorted' \
                if pypeit_file is None or len(pypeit_file) == 0 \
                else pypeit_file.replace('.pypeit', '.sorted')
            if sort_dir is not None:
                sorted_file = os.path.join(
                    sort_dir, os.path.split(sorted_file)[1])
            self.fitstbl.write_sorted(
                sorted_file, write_bkg_pairs=write_bkg_pairs)
            msgs.info("Wrote sorted file data to {:s}".format(sorted_file))

        else:
            # Write the calib file
            calib_file = self.spectrograph.spectrograph + '.calib' \
                if pypeit_file is None or len(pypeit_file) == 0 \
                else pypeit_file.replace('.pypeit', '.calib')
            if sort_dir is not None:
                calib_file = os.path.join(
                    sort_dir, os.path.split(calib_file)[1])
            self.fitstbl.write_calib(calib_file)

        # Finish (depends on PypeIt run mode)
        # TODO: Do we need this functionality?
        if calibration_check:
            msgs.info("Inspect the .calib file: {:s}".format(calib_file))
            msgs.info("*********************************************************")
            msgs.info("Calibration check complete and successful!")
#            msgs.info("Set 'run calcheck False' to continue with data reduction")
            msgs.info("*********************************************************")
            # TODO: Why should this not return the same as when setup_only is True

        if setup_only:
            for idx in np.where(self.fitstbl['failures'])[0]:
                msgs.warn("No Arc found: Skipping object {:s} with file {:s}".format(
                    self.fitstbl['target'][idx], self.fitstbl['filename'][idx]))
            msgs.info("Setup is complete.")
            msgs.info("Inspect the .sorted file")
            return None, None, None

        return self.par, self.spectrograph, self.fitstbl

    def spit_classify(self, spit_path=os.getenv('SPIT_MODEL_PATH')):
        """
        Get frame types of all SPIT predictions.

        Parameters:
            spit_path:
                Path to the best spit model. Environmental variable by default.
                ***User must have model downloaded, and must manually set path.***
        Returns:
            spit_frame_types:
                1-d list containing string frame type for each file
            spit_framebits:
                1-d list containing frame bit integers for each file
        """

        # load best model, SPIT_MODEL_PATH must be set by user elsewhere until best_model.h5 is in spit or pypeit repo
        best_model = keras.models.load_model(spit_path+'best_model.h5')

        # initialize predictions list
        predictions = []

        # loop through all files in fitstbl['filename']
        for directory, file in zip(self.fitstbl['directory'], self.fitstbl['filename']):
            # get image data from fits file
            raw_frame, hdu = self.spectrograph.load_raw_frame(
                os.path.join(directory, file))

            # perform pre-processing on the image
            preproc_dict = spit_p.original_preproc_dict()
            img_np = spit_p.process_image(raw_frame, preproc_dict)
            
            img_np = np.uint8(img_np)
            # classify with spit
            pred_list, pred = self.predict_one_image(img_np, best_model)

            # add to list of predictions
            predictions.append(pred)

        # get frame types and framebits from spit predictions
        spit_frame_types, spit_framebits = self.predictions_to_framebits(
            predictions)

        # return predictions
        return spit_frame_types, spit_framebits

    def predict_one_image(self, image, model):
        """
        Predict one image's label with a given neural network.

        Parameters:
            image:
                2-d numpy array containing image data
            model: 
                keras model, likely densenet
        Returns:
            prediction[0]: 
                1-d array containing probability of each label being the true label
            maxIndex[0][0]: 
                integer containing the index at which the highest probability (predicted label) is contained
        """
        # get shape of the image
        ishape = list(image.shape)

        # reshape to what is necessary for classifier
        image = image.reshape([1]+ishape+[1])

        # get predictions for labels
        prediction = model.predict(image)

        # get max value of prediction
        maxIndex = np.where(prediction[0] == np.amax(prediction[0]))

        # return correct values within respective return values
        return prediction[0], maxIndex[0][0]

    def predictions_to_framebits(self, predictions):
        """
        Get framebits and frame types of all SPIT predictions.

        Parameters:
            predictions: 
                1-d list containing integers representing index of label
        Returns:
            spit_frame_types:
                1-d list containing string frame type for each file
                type_bits:
                        1-d list containing frame bit integers for each file
        """
        # convert integers to string labels using label_dict and frame_types, store in numpy array with each index having strings of frametypes
        spit_labels = self.get_spit_labels()
        spit_frame_types = self.get_spit_frame_types(spit_labels, predictions)

        # use some of the code from metadata.py lines 1267-1295
        type_bits = np.zeros(
            len(self.fitstbl), dtype=self.fitstbl.type_bitmask.minimum_dtype())

        # create dictionary to be used for setting framebit
        user = {}
        for index, file in enumerate(self.fitstbl['filename']):
            user[file] = spit_frame_types[index]

        # Use the user-defined frame types from the input dictionary, could eventually integrate this and other methods in metadata.py
        if len(user.keys()) != len(self.fitstbl):
            raise ValueError(
                'The user-provided dictionary does not match table length.')
        msgs.info('Converting SPIT frame types.')
        for ifile, ftypes in user.items():
            indx = self.fitstbl['filename'] == ifile
            # turn on relevant framebits
            type_bits[indx] = self.fitstbl.type_bitmask.turn_on(
                type_bits[indx], flag=ftypes.split(','))

        return spit_frame_types, type_bits

    def get_spit_labels(self):
        """
        Get SPIT labels in framematch format.

        Parameters:
                None
        Returns:
            classes:
                1-d list containing cleaned up SPIT labels for frame typing.
        """

        # initialize label dict and classes array, could optimize this
        label_dict = spit_l.kast_label_dict()
        classes = []

        for label in label_dict:
            if '_label' in label:
                label = label.replace('_label', '')
            classes.append(label)
        return classes

    def get_spit_frame_types(self, spit_labels, predictions):
        """
        Get frame types of all SPIT predictions.

        Parameters:
            spit_labels:
                Cleaned up labels from SPIT's label dict (strings)
            predictions: 
                1-d list containing integers representing index of label
        Returns:
            spit_frame_types:
                1-d list containing string frame type for each file
        """
        # initialize empty list
        spit_frame_types = []

        # get all frame types for each label
        for p in predictions:
            spit_frame_type = spit_labels[p]  # a string

            # map spit labels to proper frames for specific cases
            if 'flat' in spit_frame_type:
            	spit_frame_type = spit_frame_type.replace('flat','pixelflat,trace')
            elif 'arc' in spit_frame_type:
            	spit_frame_type += ',tilt'
            spit_frame_types.append(spit_frame_type)
        return spit_frame_types
