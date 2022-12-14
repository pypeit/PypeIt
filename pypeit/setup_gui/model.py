import os
import re
import traceback
import enum

import numpy as np

from qtpy.QtCore import QAbstractTableModel, QAbstractItemModel, QModelIndex, Qt, Signal, QObject, QThread

from configobj import ConfigObj

from pypeit import msgs, spectrographs
from pypeit.pypeitsetup import PypeItSetup
from pypeit.inputfiles import PypeItFile

class ModelState(enum.Enum):
    NEW = enum.auto()        
    UNCHANGED = enum.auto()
    CHANGED = enum.auto()


class OpCanceledError(Exception):
    def __init__(self):
        super().__init__()

class LogWatcher():
    def __init__(self, log_file):
        if log_file is not None:
            self._log = open(os.fspath(log_file), "w")
        else:
            self._log = None
        self._watches = dict()

    def write(self, message):
        if self._log is not None:
            self._log.write(message)
            self._log.flush()
        for watch in self._watches.items():
            match = watch[1][0].search(message)
            if match is not None:
                watch[1][1](watch[0], match)


    def close(self):
        if self._log is not None:
            self._log.close()

    def watch(self, name, compiled_re, callback):
        self._watches[name] = (compiled_re, callback)

    def unwatch(self, name):
        del self._watches[name]

class SpectrographProxy():
    
    @staticmethod
    def available_spectrographs():
        return spectrographs.available_spectrographs


class PypeItMetadataProxy(QAbstractTableModel):
    def __init__(self, config="ObsLog"):
        super().__init__()
        self._config = config
        # Setup the default columns for an empty table
        self.reset()

    def reset(self):
        super().beginResetModel()
        self._colnames = ['filename', 'frametype', 'ra', 'dec', 'target', 'dispname', 'decker', 'binning', 'mjd', 'airmass', 'exptime']
        self._metadata = None
        self.sortColumn = 8 # "mjd"
        self.sortOrder = Qt.AscendingOrder
        super().endResetModel()


    def rowCount(self, parent_index):
        """Overridden method from QAbstractTableModel. Returns number of rows in metadta"""
        if (parent_index.isValid() or # Per Qt docs for a table model
            self._metadata is None):
            return 0
        else:            
            return len(self._metadata.table[self._displayIdx])

    def columnCount(self, parent_index):
        if parent_index.isValid():
            # Per Qt docs for a table model
            return 0
        else:
            return len(self._colnames)

    def data(self, index, role):

        if role == Qt.TextAlignmentRole:
            return Qt.AlignRight

        elif role != Qt.DisplayRole or self._metadata is None or index.row() > len(self._metadata.table) or index.column() > len(self._metadata.table.colnames):
            return None # Nothing to display at this index
            
        else:
            # The columns being displayed are a subset of the metadata based on PypeItMetadata.set_pypeit_cols,
            # So we use the column name instead of the index
            colname = self._colnames[index.column()]
            value = self._metadata[colname][self._displayIdx][self._sortIdx][index.row()]
            # Round floating point values to look better
            if isinstance(value, np.float64) or isinstance(value, np.float32):
                value = round(value, 3)
            return str(value)

    def sort(self, column, order=Qt.AscendingOrder):
        if self._metadata is not None and column < len(self._colnames):
            self.sortColumn = column
            self.sortOrder = order
            colname = self._colnames[self.sortColumn]
            self._sortIdx = self._metadata.table[self._displayIdx].argsort(colname, reverse=(order==Qt.DescendingOrder))
            self.dataChanged.emit(self.index(0, 0), self.index(len(self._displayIdx)-1, len(self._colnames)-1))

    def headerData(self, section, orientation, role):
        if role == Qt.DisplayRole:
            if orientation == Qt.Orientation.Horizontal and section < len(self._colnames):
                return self._colnames[section]
            else:
                return str(section + 1)
        elif role == Qt.InitialSortOrderRole and section == self.sortColumn:
            return self.sortOrder
        else:
            return None

    def setMetadata(self, metadata):
        """
        Sets a new PypeItMetadata object for the proxy model.
        """
        # Tell views a model reset is happening
        super().beginResetModel()

        # Preseve the original sorting column before clearing the column names
        if self.sortColumn < len(self._colnames):
            origianl_sort_colname = self._colnames[self.sortColumn]
        else:
            origianl_sort_colname = "mjd"

        self._metadata = None
        self._colnames = []
        
        # Now set the new column names
        self._colnames = metadata.set_pypeit_cols(write_bkg_pairs=True)

        # Finally set the new data, signalling views about the new rows.
        self._metadata = metadata
        if self._config != "ObsLog":
            self._displayIdx = metadata.table['setup'] == self._config
        else:
            self._displayIdx = [True for x in metadata.table]

        # Try to set sorting parameters to what they were before.
        if origianl_sort_colname in self._colnames:
            self.sortColumn = self._colnames.index(origianl_sort_colname)
        
        # If that isn't possible, try to reset to the default "mjd" column
        elif "mjd" in self._colnames:
            self.sortColumn = self._colnames.index("mjd")
            self.sortOrder = Qt.AscendingOrder

        # Otherwise just sort by the first column
        else:
            self.sortColumn = 0
            self.sortOrder = Qt.AscendingOrder
        
        self._sortIdx = self._metadata.table[self._displayIdx].argsort(self._colnames[self.sortColumn], reverse=(self.sortOrder == Qt.DescendingOrder))
        super().endResetModel()

class _UserConfigTreeNode:
    """
    Internal class used for the internal pointer in PypeItParamsProxy
    """
    def __init__(self, node, key=None, parent=None):
        
        self.parent = parent
        self.key = key
        if isinstance(node, dict):
            self.children = [_UserConfigTreeNode(key=k, node=node[k], parent=self) for k in node.keys()]
            self.value=None
        else:
            self.value=node
            self.children = []
    

class PypeItParamsProxy(QAbstractItemModel):
    """
    A Proxy model that maps a PypeIt PypeItPar to a QAbstractItemModel suitable for a Qt QTreeView

    Args:
    par (pypeit.par.parset.PypeItPar): The PypeIt parameters being displayed as a tree.

    """
    def __init__(self, pypeit_setup):
        super().__init__(None)

        self.par = pypeit_setup.par
        self._userConfigTree = _UserConfigTreeNode(ConfigObj(pypeit_setup.user_cfg))

    def rowCount(self, parent):
        """
        Returns the number of items under parent. Overridden from QAbstractItemModel.
        
        Returns: (int) - The number of items under parent. Can be 0 if parent has no child items.

        """
        #msgs.info(f"rowCount: {repr(parent)}")
        if not parent.isValid():
            # If parent is invalid, it's the root of the tree
            #msgs.info("rowCount invalid")
            node = self._userConfigTree
        else:
            # Otherwise, if the parent is an index created by the index() method, it
            # points to the parent node in its internalPointer()
            #msgs.info("rowCount valid")
            node = parent.internalPointer()

        return len(node.children)

    def index(self, row, column, parent):
        """
        Creates a QModelIndex that points to an item in the ParSet tree. Overridden from QAbstractItemModel.
        
        Args:
            row    (int): The row of the item underneath parent to create an index for.

            column (int): The column of the item. Typically 0 is used for the name of the ParSet or 
                          value, and 1 is used for actual parameter values (1 is never used for a ParSet).

            parent (QModelIndex): An index for the parent. If this is invalid, it refers to the root of the tree.

        Returns:
            A new index object to point to the item.
        
        """
        #msgs.info(f"index: {row} {column} {repr(parent)}")
        if not parent.isValid():
            # Use the root of the config tree as the parent
            parent_node = self._userConfigTree
        else:
            # The parent is valid, we're creating an index to one of its children
            # Get the parent from the internal data created by this method earlier.
            parent_node = parent.internalPointer()

        # Find the child using the passed in row
        child_node = parent_node.children[row]

        # Column 1 is for actual config values, not nested parameter sets
        # (Or Column 1 is for leaves of the tree if that makes more sense)
        # Return an invalid index to indicate there's nothing to display here
        if column == 1 and child_node.value is None:
            return QModelIndex()
        else:
            # Otherwise create the index.
            return super().createIndex(row, column, child_node)

    def data(self, index, role):
        """Returns data for a given index. Overridden from QAbstractItemModel.

        Args:
        index (QModelIndex): The index of the item, as returned by index()

        role (Qt.ItemDataRole): The role of the data being requested. This method supports
                                Qt.DisplayRole (for displaying text).
        
        
        Returns: (str or None): A string if there's something to display at this index for a DisplayRole, or
                                None if there's nothing to display or an unsupported role
        """
        #msgs.info(f"data: {repr(index)}, {repr(role)}")
        if role == Qt.DisplayRole:
            if index.column() == 0:
                # Display the name (aka key) of the item
                # This is the reason for including the key in the UserConfigTreeNode.
                return index.internalPointer().key
            else:
                # Display the value, which will be None if the index isn't
                # pointing to a leaf node
                value = index.internalPointer().value
                if value is not None:
                    return str(value)

        return None

    def headerData(self, section, orientation, role):
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            if section == 0:
                return "Setting"
            elif section == 1:
                return "Value"
        return super().headerData(section, orientation, role)


    def columnCount(self, parent):
        """Return the number of columns, which is always 2 because we use column 0 for the ParSet tree and
        column 1 for the configuration values. Overridden from QAbstractItemModel.

        Args:
        parent (QModelIndex): Parent to return the column count for. (Not used)

        Return: (int) The number of columns. Always 2.        
        """
        #msgs.info(f"columnCount: {repr(parent)}")
        return 2

    def parent(self, index):
        """Return the parent of an item in the model. Overridden from QAbstractItemModel.

        Args:
        index (QModelIndex): The index of the item, as returned by the index() method.

        Return: (QModelIndex): An index for the items parent. Maybe be invalid if item is at the top level.
        """
        # This method is why the parent node is included in the _UserConfigTreeNode data.
        #msgs.info(f"parent: {repr(index)}")
        node = index.internalPointer()
        if node.parent is None:
            # Root node, there is no parent
            return QModelIndex()
        else:
            # Need to know the row # for the parent, which is in the grandparent
            grandparent = node.parent.parent
            if grandparent is None:
                # The parent is the root node, it's row is 0
                row = 0
            else:
                parent_and_siblings = [x.key for x in grandparent.children ]
                row = parent_and_siblings.index(node.parent.key)
            # Column is always 0, because 1 is reserved for leaf nodes which can't be parents
            return super().createIndex(row, 0, node.parent)

class SetupConfigurationModel(QObject):
    """
    A model representing a single PypeIt configuration, consisting of file metadata, a common observing setup for all those files,
    and PypeIt parameters. This is the information needed for a single .pypeit file.

    Args: 
    pypeit_setup (pypeit.pypeitsetup.PypeItSetup): The PypeItSetup object this configuration is a part of.
    config_name (str): The name of this configuration.
    config_dict (dict): The metadata for this configurations, as returned by PypeItMetadata.unique_configurations.
    state (ModelState): The state of the model.  "CHANGED" if it was built by running setup or "UNCHANGED" if it came from loading
                        an existing .pypeit file.

    Signals:
    stateChanged (str): Sent when the state of the configuration changes. The name of the configuration is sent as the signal parameter.
    """
    stateChanged = Signal(str)

    def __init__(self, pypeit_setup, config_name, config_dict, state):
        super().__init__()
        self._pypeit_setup = pypeit_setup
        self.name = config_name
        self.config_dict = config_dict
        self.metadata_model = PypeItMetadataProxy(config_name)
        self.metadata_model.setMetadata(pypeit_setup.fitstbl)
        self.params_model = PypeItParamsProxy(pypeit_setup)
        self.state = state

    def save(self, outdir):
        try:
            self._pypeit_setup.fitstbl.write_pypeit(outdir, cfg_lines=self._pypeit_setup.user_cfg,
                                                    configs = [self.name])
        except Exception as e:
            msgs.warn(f"Failed saving setup {self.name} to {outdir}.")
            msgs.warn(traceback.format_exc())
            # Raise an exception that will look nice when displayed to the GUI
            raise RuntimeError(f"Failed saving setup {self.name} to {outdir}.\nException: {e}")
        self.state = ModelState.UNCHANGED
        self.stateChanged.emit(self.name)
        

class PypeItSetupProxy(QObject):
    operation_progress = Signal(str)
    operation_complete = Signal(str)
    raw_data_dir_changed = Signal(str)
    configs_deleted = Signal(list)
    configs_added = Signal(list)
    spectrograph_changed = Signal(str)


    def __init__(self, logname):
        super().__init__()
        self._spectrograph = None
        self._pypeit_setup = None
        self.metadata_model = PypeItMetadataProxy()
        self.raw_data_directory = None
        self.configs = dict()

    def setup_logging(self, logname, verbosity):
        self._log_watcher = LogWatcher(logname)
        msgs.reset(verbosity=verbosity, log_object=self._log_watcher)


    @property
    def available_spectrographs(self):
        return spectrographs.available_spectrographs


    @property
    def state(self):
        if len(self.configs.values()) == 0:
            return ModelState.NEW
        else:
            if any([config.state==ModelState.CHANGED for config in self.configs.values()]):
                return ModelState.CHANGED
            else:
                return ModelState.UNCHANGED

    def set_spectrograph(self, new_spec):
        msgs.info(f"Spectrograph is now {new_spec}")
        self._spectrograph = new_spec
        self.spectrograph_changed.emit(new_spec)

    @property
    def spectrograph(self):
        return self._spectrograph

    def set_raw_data_directory(self, new_directory):
        """
        Change the model's current raw data directory.

        Args:
        new_directory (str): The new directory containing raw data.

        Signals:
        raw_data_dir_chasnged(new_directory)
        
        """
        msgs.info(f"Setting raw directory: {new_directory}, current spec is {self._spectrograph}")
        spectrograph = self._spectrograph
        self.reset()
        self.raw_data_directory = new_directory
        self._spectrograph = spectrograph
        self._pypeit_setup = PypeItSetup.from_file_root(new_directory, self._spectrograph) 
        # Set the default configuration params for a new pypeit file
        self._pypeit_setup.user_cfg = ['[rdx]', f'spectrograph = {self._spectrograph}']
        self.raw_data_dir_changed.emit(self.raw_data_directory)
        self.spectrograph_changed.emit(self._spectrograph)

    def get_num_raw_files(self):
        return 0 if self._pypeit_setup is None else len(self._pypeit_setup.file_list)

    def reset(self):
        msgs.info(f"Resetting to empty setup.")
        self._pypeit_setup = None
        self.metadata_model.reset()        
        self.raw_data_directory = None
        self.raw_data_dir_changed.emit(None)
        self._spectrograph = None
        self.spectrograph_changed.emit(None)
        self._setConfigurations({})

    def generate_obslog(self):

        try:
            if self._pypeit_setup is None:
                raise ValueError("No PypeItSetup object. set_raw_data_directory should be called before calling generate_obslog.")

            added_metadata_re = re.compile("Adding metadata for (\S.*)$")
            self._log_watcher.watch("added_metadata", added_metadata_re, self._addedMetadata)

            # These were taken from the default parameters in pypeit_obslog
            self._pypeit_setup.run(setup_only=True,
                                write_files=False, 
                                groupings=True,
                                clean_config=False)
            self._log_watcher.unwatch("added_metadata")
            self.metadata_model.setMetadata(self._pypeit_setup.fitstbl)
            self._setConfigurations(self._pypeit_setup.fitstbl.unique_configurations())
            self.operation_complete.emit("")
        except OpCanceledError:
            # The operation was canceled, reset pypeit_setup and return
            msgs.info("OpCanceled")
            self._pypeit_setup = None
            self.operation_complete.emit("CANCEL")
        except Exception as e:
            msgs.info("Exception")
            msgs.info(traceback.format_exc())
            # Any other exception is an error reading the metadata
            self._pypeit_setup = None
            self.operation_complete.emit("Could not read metadata. Are you sure the spectrograph is set correctly?\nCheck the logs for more information.")

    def save_all(self, outdir):
        for config in self.configs.values():
            config.save(outdir)


    def save_config(self, config_name, outdir):
        config = self.configs[config_name]
        config.save(outdir)

    def open_pypeit_file(self, pypeit_file):
        
        # We can't just create a PypeItSetup using from_pypeit_file because we 
        # need the paths for raw_data_directory.
        pf = PypeItFile.from_file(pypeit_file)

        if len(pf.file_paths) > 0:
            # TODO: Support multiple file paths
            self.raw_data_directory = pf.file_paths[0]
        else:
            raise ValueError(f"PypeIt input file {pypeit_file} is missing a path entry.")

        self._pypeit_setup = PypeItSetup(pf.filenames,
                                         usrdata=pf.data,
                                         setups=[pf.setup_name],
                                         cfg_lines=pf.cfg_lines,
                                         pypeit_file=pypeit_file,
                                         setup_dict=pf.setup)

        # Setup our proxy models and notify the view widgets
        self.raw_data_dir_changed.emit(self.raw_data_directory)
        self._pypeit_setup.build_fitstbl()        
        self._pypeit_setup.fitstbl.set_configurations(fill=pf.setup_name)
        self._pypeit_setup.fitstbl.get_frame_types(user=dict(zip(pf.data['filename'], pf.data['frametype'])))
        self.metadata_model.setMetadata(self._pypeit_setup.fitstbl)        
        self._spectrograph = self._pypeit_setup.spectrograph.name
        self.spectrograph_changed.emit(self._spectrograph)
        self._setConfigurations(self._pypeit_setup.fitstbl.unique_configurations(), state=ModelState.UNCHANGED)

    def _addedMetadata(self, name, match):
        if QThread.currentThread().isInterruptionRequested():
            raise OpCanceledError()

        self.operation_progress.emit(match.group(1))

    def _setConfigurations(self, unique_configs, state=ModelState.CHANGED):
        """
        Private method to reset configurations to a new set. Any previous configurations are deleted.

        Args:
        unique_config (dict): A the new configurations to set.

        state (ModelState): The state of the new model. Defaults to CHANGED, meaning it has not been saved.

        Return: None
        """
        msgs.info(f"Unique Configs {unique_configs}")

        # Delete previous configurations
        self.configs_deleted.emit(list(self.configs.keys()))

        config_names = list(unique_configs.keys())
        self.configs = dict()
        for config_name in config_names:
            new_config_model = SetupConfigurationModel(self._pypeit_setup, config_name, unique_configs[config_name], state=state)
            self.configs[config_name] = new_config_model

        msgs.info(f"Self configs: {self.configs}")
        self.configs_added.emit(list(self.configs.values()))