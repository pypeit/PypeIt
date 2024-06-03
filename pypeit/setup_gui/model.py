"""
The model portion of the PypeIt Setup GUI.  Responsible for holding the data of the Setup GUI, notifying the view of any changes
to that data, and translating PypeIt datastructures to a form usable by Qt.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import os
from collections import deque
import copy
import traceback
import enum
import glob
import numpy as np
import astropy.table
import io
import typing
from pathlib import Path
from qtpy.QtCore import QAbstractTableModel, QAbstractItemModel, QAbstractListModel, QModelIndex, Qt, Signal, QObject, QThread, QStringListModel
import qtpy
from configobj import ConfigObj

from pypeit import msgs, spectrographs
from pypeit.spectrographs import available_spectrographs
from pypeit.pypeitsetup import PypeItSetup
from pypeit.metadata import PypeItMetaData
from pypeit.inputfiles import PypeItFile
import pypeit.core

class ModelState(enum.Enum):
    """The state values for a model object."""

    NEW = enum.auto()
    """The model is in a fresh, uninitialized state."""

    UNCHANGED = enum.auto()
    """The model contains data that has not been changed since being read from or saved to a file."""

    CHANGED = enum.auto()
    """The model contains data that has not been saved to a file."""




class LogBuffer(io.TextIOBase):
    """Imitation file object that is passed to the PypeIt msgs logging system. It maintains a buffer
    of log messages that the user can view through the GUI. It is also used to monitor progress of 
    background operations, by registering regular expressions to watch the log for.

    Args:
        log_file (str): The log file to receive the log messages. If this is None the log is not
                        written to a file.
        verbosity (int): The verbosity of log messages to pass on. 0 = No logging. 1 = INFO,
                         BUG, WARNING, and ERROR only, 2 = All.
        max_len (int,Optional):   The maximum number of log lines to buffer. Defaults to 1000.
    """

    def __init__(self, log_file, verbosity, max_len=1000):
        super().__init__()

        if log_file is not None:
            self._log = open(os.fspath(log_file), "w")
        else:
            self._log = None
        self._verbosity = verbosity
        self._watches = dict()
        self.maxlen = max_len
        self._buffer = deque([], maxlen=max_len)

    def write(self, message):
        """Simulates the write method of a file object to monitors log messages for for matching messages.
        The messages are then sent to a log file (if one has been configured).
        
        Args:
            message (str): The log message being written to the log.
        """

        # Store the message
        self._buffer.append(message)

        # Notify clients for specific watched log messages
        for watch in self._watches.items():
            re = watch[1][0]
            if re is None:
                watch[1][1](message)
            else:
                match = re.search(message)
                if match is not None:
                    watch[1][1](watch[0], match)

        # Enforce verbosity
        if self._verbosity == 0:
            return

        # Write to the log file if one was given
        if self._log is not None:
            self._log.write(message)
            self._log.flush()


    def close(self):
        """Closes the log file (if any)."""
        if self._log is not None:
            self._log.close()

    def watch(self, name, compiled_re, callback):
        """Registers a regular expression to watch the log for.
        If a log message that matches the expression is logged, the client callback 
        is called.

        Args:
            name (str):                                      Name to register the regular expression under. This can be
                                                             passed to :meth:`unwatch` to stop monitoring for that expression.
            callback (:class:`collections.abc.Callable`):    A function or object that will be called when a matching log
                                                             message is found. It is called with two arguments: the name
                                                             used to register the expression and the Match object regurned by
                                                             the regular expression.
            compiled_re (:class:`typing.Pattern`, Optional): A compiled Python regular expression to match log messages
                                                             against. If this is not given, the caller is notified of all log messages.
        """
        self._watches[name] = (compiled_re, callback)

    def unwatch(self, name):
        """Stops monitoring the log for a previously registered regular expression.

        Args:
            name (str): The name passed to :meth:`watch` to register the regular expression.
        """
        if name in self._watches:
            del self._watches[name]

    def __iter__(self):
        """Allow iteration through the log buffer.
        
        Returns:
            (:obj:`collections.abc.Iterator`): An iterator over the lines in  the log buffer."""
        return self._buffer.__iter__()

    def __len__(self):
        """Return the number of lines in the buffer.
        
        Returns:
            int
        """
        return len(self._buffer)

    def __bool__(self):
        """Return a true status indicating we're ready to receive data"""
        return True

class PypeItMetadataUniquePathsProxy(QAbstractListModel):
    """A Proxy model filtering the content of a PypeItMetadataModel object to only show the
    unique paths within it to Qt views.
    
    Args:
        metadata_model (PypeItMetaData): The model being filtered.
    """
    def __init__(self, metadata_model):
        super().__init__()
        self.metadata = metadata_model.metadata
        self._setUniqueIndex()
        metadata_model.modelReset.connect(self._setUniqueIndex)
        metadata_model.rowsInserted.connect(self._setUniqueIndex)
        metadata_model.rowsRemoved.connect(self._setUniqueIndex)
        metadata_model.dataChanged.connect(self._setUniqueIndex)

    def _setUniqueIndex(self, *args, **kwargs):
        """Sets the Numpy index array for the unique paths within the metadata."""
        self.beginResetModel()
        
        if self.metadata is not None:
            items, self._unique_index = np.unique(self.metadata['directory'],return_index=True)
        else:
            self._unique_index = []

        self.endResetModel()
    
    def rowCount(self, parent_index=QModelIndex()):
        """Returns the number of unique paths. Inherited from QAbstractItemModel."""
        return len(self._unique_index)
    
    def data(self, index, role):
        """Returns the path for a given QModelIndex.
        
        Args:
            index (QModelIndex):   The QModelIndex for the row to get data for.
            role (Qt.DisplayRole): The role to return data for. This method only supports "DisplayRole".
        """
        if role == Qt.DisplayRole:
            if index.isValid() and self.metadata is not None and index.row() < len(self._unique_index):
                return str(self.metadata['directory'][self._unique_index][index.row()])

        return None    

    def getPaths(self):
        """Returns all the paths known to this model.
        
        Return: 
            list of str: The list of paths, or an empty list of the model is empty.
        """
        return [] if self.metadata is None else list(self.metadata['directory'][self._unique_index])

class PypeItMetadataModel(QAbstractTableModel):
    """
    Provides a Qt model interface for a PypeItMetaData object.  This proxy implements
    a QAbstractItemModel interface to present file metadata to Qt views. 

    It also supports editing. 

    Args:
        metadata: The PypeItMetaData object being wrapped. If this is None, the
                  model is in a "NEW" state.
    """
    def __init__(self, metadata : typing.Union[PypeItMetaData, None]):
        super().__init__()

        self.metadata = metadata
        self.editable_columns=['calib', 'comb_id', 'bkg_id', 'frametype']

        self.colnames = []

        # How large a character is in a unicode dtype. Used to determine string maximum length
        temp_dtype = np.dtype("U1")
        self._dtype_char_size = temp_dtype.itemsize

        self.reset()
        
    @property
    def spectrograph(self) -> spectrographs.spectrograph:
        """The spectrograph object for the metadata. None if no metadata has been set."""
        return None if self.metadata is None else self.metadata.spectrograph

    def getColumnNumFromName(self, colname):
        """Return the column number that has the passed in name.
        
        Args:
            colname (str): The name of the column.

        Return:
            int: 
                The column number (0 based) of the named column. -1 if the column
                name is not in this model.
        """
        if colname in self.colnames:
            return self.colnames.index(colname)
        else:
            return -1

    def getColumnNameFromNum(self, index):
        """Return the name of the column with the given number.
        
        Args:
            index (str): The 0 based column number.

        Return:
            str: 
                The name of the column. If the passed in index is out of bounds,
                an IndexError value is raised.
        """

        return self.colnames[index.column()]

    def getAllFrameTypes(self):
        """Return the allowable values for the frametype column.
        
        Return:
            list of str: List of names of the allowable frame types.
        """
        return pypeit.core.framematch.FrameTypeBitMask().keys()

    def rowCount(self, parent_index=QModelIndex()):
        """Returns number of rows under a parent. Overridden method from QAbstractItemModel.
        
        Args:
            parent_index (QModelIndex): The index of the parent. Not applicable to a table model
                                        as there's no parent/child heirarchy.

        Returns:
            int: The number of rows in the table.
        """
        if (parent_index.isValid() or # Per Qt docs for a table model
            self.metadata is None):
            return 0
        else:
            return len(self.metadata)

    def columnCount(self, parent_index=QModelIndex()):
        """Returns number of columns in under a parent. Overridden method from QAbstractItemModel. 

        Args:
            parent_index (QModelIndex): The index of the parent. Not applicable to a table model
                                        as there's no parent/child heirarchy.

        Returns:
            int: The number of columns in the table.
        """
        if parent_index.isValid():
            # Per Qt docs for a table model
            return 0
        else:
            return len(self.colnames)

    def data(self, index, role):
        """Returns data for a given role at a given index. Overridden method from QAbstractItemModel. 
        
        Args:
            index (QModelIndex): The index in the table to return data form.
            role (Qt.DisplayRole): The role to return data for. This method supports the "TextAlignmentRole"
                                   for returning alignment information for a table cell, the "EditRole" 
                                   for data to edit, and the "DisplayRole" for displaying the data within a table cell.

        Return:
            Object: The requested data, or None if there is no applicable data.
        """

        if role == Qt.TextAlignmentRole:
            return Qt.AlignLeft

        elif self.metadata is not None and (role == Qt.EditRole or role==Qt.DisplayRole):
            # The columns being displayed are a subset of the metadata,
            # So we use the column name instead of the column number
            colname = self.colnames[index.column()]
            value = self.metadata[colname][index.row()]

            if role == Qt.DisplayRole:
                # Convert to a string for displaying

                # Check for commented out filenames and remove the # characters
                if colname == 'filename' and value.lstrip().startswith("#"):
                    value = value.strip("# ")

                # Round floating point values to look better
                elif isinstance(value, np.float64) or isinstance(value, np.float32):
                    value = round(value, 3)
                return str(value)
            else:
                return value
        # No data for any other cases
        return None

    def setData(self, index, value, role=Qt.EditRole):
        """Modifies data for a given role at a given index. Overridden method from QAbstractItemModel. 
        
        Args:
            index (QModelIndex):   The index in the table to return data form.
            value (Any):           The value to set in the metadata. 
            role (Qt.DisplayRole): The role to set data for. Only the "EditRole" is supported.
                                   Defaults to Qt.EditRole.  

        Return:
            True if the edit succeeded, false if it failed.
        """
        if role==Qt.EditRole and self.metadata is not None:
            colname = self.colnames[index.column()]
            if colname in self.editable_columns:
                max_length = self.getStringColumnSize(colname)
                if max_length is not None:
                    # This is a string type, does the new data fit
                    value = str(value)
                    if len(value) > max_length:
                        # We need to increase the string size to make the new value fit
                        self.resizeStringColumn(colname, len(value))
                
                try:
                    self.metadata[colname][index.row()] = value
                except ValueError as e:
                    msgs.warn(f"Failed to set {colname} row {index.row()} to '{value}'. ValueError: {e}")

                self.dataChanged.emit(index,index,[Qt.DisplayRole, Qt.EditRole])
                return True

        return False

    def flags(self, index):
        """Returns item flags for a given index. Overridden method from QAbstractItemModel
        
        Args:
            index (QModelIndex): The index to get flags for


        Return:
            flags (Qt.ItemFlag): 
                The bitmask flags for the index. Currently only Qt.ItemFlag.ItemIsEditable
                is supported.
        """
        base_flags = super().flags(index)
        if self.colnames[index.column()] in self.editable_columns:
            base_flags |= Qt.ItemFlag.ItemIsEditable
        return base_flags

    def headerData(self, section, orientation, role):
        """Display header data for the table. For columns we give a column name, for rows we return
        a number. Overridden method from QAbstractItemModel.

        Args:
            section (int): The "section" the header is used for. For table models this is the 
                           zero based row or column number.
            orientation (Qt.Orientation): Whether the horizontal (column header) or vertical (row header)
                                          is being requested.
            role (Qt.Role): The role header data is being requested for. This method supports
                            Qt.InitialSortOrderRole and Qt.DisplayRole.

        Returns:
            str: The name for the column/row, or None if there is no applicable name for the given section/orientation/role.
        """
        if role == Qt.DisplayRole:
            # Display the actual name of the header row/column
            if orientation == Qt.Orientation.Horizontal and section < len(self.colnames):
                # Columns have propper names
                return self.colnames[section]
            else:
                return " "
        else:
            # A non-applicable role or a sort order request for a column that we're not sorted by.
            return None

    def getSetup(self, row : int) -> dict:
        """Return the configuration/setup values for the given row of metadata.
        
        Args:
            row: The row within the model to return the setup values for
        Return:
            A dictionary of configuration key/value pairs. An empty dictionary is returned
            if there is no metadata or the spectrograph has no configuration keys
        """
        if self.metadata is None or row >= len(self.metadata):
            return {}
        return self.metadata.get_configuration(row)



    def getDefaultColumns(self):
        """Return the default columns to display to the user. This can vary based
        on the current spectrograph.
        
        Return:
        A list of column names, in order to display.
        """
        if self.metadata is None:
            return ['filename', 'frametype', 'ra', 'dec', 'target', 'dispname', 'decker', 'binning', 'mjd', 'airmass', 'exptime']
        return self.metadata.set_pypeit_cols(write_bkg_pairs=True)

    def getStringColumnSize(self, colname: str) -> typing.Union[int,None]:
        """
        Return the maximum size of a string column.

        Args:
            colname: The name of the column
        Return:
            The maximum size of the column, or None if it is not a string column.
        """
        dt = self.metadata[colname].dtype
        if dt.kind == "U":
            # Divide the size by the size for a single character column
            return int(dt.itemsize / self._dtype_char_size)
        else:
            # Not a string
            return None

    def resizeStringColumn(self, colname : str, new_size : int)->None:
        """Resize a string column to fit a new size.
        
        Args:
            colname: The name of the string column.
        """
        # Create a new column with the new size.
        # We do this before setting the column, so that an exception will leave the
        # metadata in the previous state
        new_column = astropy.table.Column(self.metadata[colname],dtype=f'<U{new_size}')

        self.metadata[colname] = new_column

    def reset(self):
        """
        Reset the proxy assuming the metadata object has completely changed
        """
        # Tell views a model reset is happening
        self.beginResetModel()

        # Reset column names
        self.colnames = self.getDefaultColumns()

        self.endResetModel()

    def setMetadata(self, metadata):
        """Sets the PypeItMetaData object being wrapped.
        
        Args:
            metadata (:obj:`pypeit.metadata.PypeItMetaData`): The metadata being wrapped.
        """
        self.metadata=metadata
        self.reset()

    def removeMetadataRows(self, rows):
        """Removes rows from the metadata
        
        Args:
            row (Sequence of int): A sequence of integer row indexes of the rows to remove.

        """
        if self.metadata is None:
            return None

        # Remove the row, making sure views are notified
        # We do this in reverse sorted order so that indices don't change before deletes
        for row in sorted(rows,reverse=True):
            # We could try to group these rows into ranges, but
            # it doesn't seem worth it
            self.beginRemoveRows(QModelIndex(), row, row)
            msgs.info(f"Removing metadata row {row}")
            self.metadata.remove_rows([row])
            self.endRemoveRows()


    def pasteFrom(self, other_metadata_model):
        if self.metadata is None:
            # Pasting into an empty model, just copy from the other metadata
            self.beginResetModel()
            self.colnames = other_metadata_model.colnames
            self.spectrograph = other_metadata_model.spectrograph
            table_copy = other_metadata_model.metadata.table.copy(True)
            self.metadata = PypeItMetaData(self.spectrograph, other_metadata_model.metadata.par, data=table_copy)
            self.endResetModel()
        else:
            # Make sure the two metadata models are compatible
            num_rows = 0 if other_metadata_model.metadata is None else len(other_metadata_model.metadata.table)
            if num_rows == 0:
                # Nothing to paste
                return
            elif self.spectrograph.name != other_metadata_model.spectrograph.name:
                raise RuntimeError(f"Can't paste file metadata between different spectrographs.")
            elif self.metadata.table.colnames != other_metadata_model.metadata.table.colnames:
                # Note we compare the underlying table's columns, since tables loaded from a file
                # will have fewer columns than ones generated by running setup
                raise RuntimeError(f"Can't paste file metadata with different columns.")
            else:
                # Paste the data. First, make sure our text columns are big enough for the corresponding
                # column in the source metadata.
                for colname in self.metadata.table.colnames:
                    our_max_size = self.getStringColumnSize(colname)
                    if our_max_size is not None:
                        # A String column
                        other_max_size = other_metadata_model.getStringColumnSize(colname)
                        if other_max_size is None:
                            # It isn't a string, something's wrong
                            raise(f"Can't paste metadata, incompatible types for column '{colname}'")
                        if other_max_size > our_max_size:
                            self.resizeStringColumn(colname, other_max_size)
                # 
                # Second update the rows that have changed
                first_updated_row = None
                last_updated_row = None
                indx = np.ones(self.rowCount(),dtype=bool)
                current_files = self.metadata.frame_paths(indx)
                indx = np.ones(other_metadata_model.rowCount(),dtype=bool)
                other_files = other_metadata_model.metadata.frame_paths(indx)
                
                indices_to_add = [] # Keep track  of new rows to add

                for i in range(len(other_files)):
                    if other_files[i] in current_files:
                        # Duplicate row, don't add a new row, rather overwrite
                        indx = current_files.index(other_files[i])
                        self.metadata.table[indx] = other_metadata_model.metadata.table[i]
                        if first_updated_row is None:
                            first_udpated_row = indx
                        last_updated_row = indx
                    else:
                        indices_to_add.append(i)

                # Signal clients for the range of updated rows. This may include extra rows
                # but that's ok
                if first_updated_row is not None:
                    self.dataChanged.emit(self.index(first_udpated_row, 0),self.index(last_updated_row,self.columnCount()-1))
                
                # Third: add new rows
                if len(indices_to_add) > 0:
                    self.beginInsertRows(QModelIndex(), self.rowCount(), self.rowCount()+len(indices_to_add)-1)
                    for i in indices_to_add:
                        self.metadata.table.add_row(other_metadata_model.metadata.table[i])
                    self.endInsertRows()
                

    def createCopyForConfig(self, config_name):
        """Create a new copy of this metadata model containing only rows for a given config.
        
        Args:
            config_name (str): Name of one of the unique configurations in the metadata.

        Return:
            PypeItMetadataModel: A deep copy of the meatdata matching the config_name
        """
        msgs.info(f"Creating new metadata for config {config_name}")

        config_rows  = [ config_name in setup for setup in self.metadata.table['setup'] ]
        return self.createCopyForRows(config_rows)
    
    def createCopyForRows(self, rows):
        """Createa a new copy of this metadata model containing only a given set of rows.
        
        Args:
            rows (list): A list of the indexes of the rows to include in the copy.

        Return:
            PypeItMetadataModel: A deep copy of the meatdata in the given rows.
        """
        # If this is an empty model, return another empty model
        if self.metadata is None:
            return PypeItMetadataModel(None)
        else: # Copy the selected rows
            table_copy = self.metadata.table[rows].copy(True)                
            copy = PypeItMetadataModel(metadata=PypeItMetaData(self.spectrograph, self.metadata.par, data=table_copy))
            return copy
    
    def getPathsModel(self):
        """Returns a model for displaying the paths in this metadata.
        
        Return:
            PypeItMetadataUniquePathsProxy: 
                A proxy model that only contains the unique 'directory' 
                entries of this metadata.
        """
        return PypeItMetadataUniquePathsProxy(self)

    def commentMetadataRows(self, rows):
        """Comment out metadata rows. Rows that are already commented out are not affected.
        
        Args:
            rows (list of int): The indices of the rows to comment out.
            
        Return:
            None
        """

        if self.metadata is not None:

            # Recreate the entire column of filenames to avoid issues with astropy tables having a fixed length string dtype
            # We are careful to not comment out a row that has already been commented out.
            new_filenames = ["# " + filename if i in rows and not filename.lstrip().startswith("#") else filename for i, filename in enumerate(self.metadata['filename'])]
            self.metadata['filename'] = new_filenames

            # We need a upper left, lower right for the data changed event
            min_row = min(rows)
            max_row = max(rows)
            col = self.colnames.index('filename')
            self.dataChanged.emit(self.index(min_row, col), self.index(max_row, col))

    def uncommentMetadataRows(self, rows):
        """Uncomment out metadata rows. Rows that are not commented out are not affected.
        
        Args:
            rows (list of int): The indices of the rows to comment out.
            
        Return:
            None
        """

        if self.metadata is not None:

            # Recreate the entire column of filenames to avoid issues with astropy tables having a fixed length string dtype
            new_filenames = [filename.strip("# ") if i in rows else filename for i, filename in enumerate(self.metadata['filename'])]
            self.metadata['filename'] = new_filenames

            # We need a upper left, lower right for the data changed event
            min_row = min(rows)
            max_row = max(rows)
            col = self.colnames.index('filename')
            self.dataChanged.emit(self.index(min_row, col), self.index(max_row, col))

    def isCommentedOut(self, index : QModelIndex) -> bool:
        """Return whether a particular row is commented out.
        
        Args:
            index: The index of the row to check.
        
        Return:
            bool:
                True if the row is commented out, False otherwise.
        """
        return self.metadata is not None and self.metadata[index.row()]['filename'].lstrip().startswith("#")

class _UserConfigTreeNode:
    """
    Internal class used to represent the tree structure of PypeIt parameters.
    It is used as the internal pointer in PypeItParamsProxy.

    Take this example config::

        [rdx]
            spectrograph = shane_kast_blue
        [calibrations]
            bpm_usebias = true
            [[biasframe]]
                frametype = bias

    The resulting data structure is::

        Node 1: 
            parent = None
            key = None
            value = None
            children = [Node 2, Node 4]
        Node 2:
            parent = None
            key = "rdx"
            value = None
            children = [Node 3]
        Node 3:
            parent = Node 2
            key = "spectrograph"
            value = "shane_kast_blue"
            children = []
        Node 4:
            parent = None
            key = "calibrations"
            value = None
            children = [Node 5, Node 6]
        Node 5:
            parent = Node 4
            key = "bpm_usebias"
            value = True
            children = []
        Node 6:
            parent = Node 4
            key = biasframe
            value = None
            children = [Node 7]
        Node 7:
            parent = Node 6
            key = frametype
            value = bias
            children = []

    Of note is that the "root" in the data model used by Qt is an invalid entry, i.e. it isn't shown in 
    the view as a node. This matches the PypeIt parameters format in that there can be multiple
    top level entries (e.g. [rdx] and [calibrations] above). This is why nodes 2 and 4 have parent = None 
    even though Node 1 is their parent in the data structure. It prevents them from reporting a parent node to Qt and 
    causing issues.

    Args:
        node (dict or value):  The data for a node within a tree. If it's a dict, this
                               is a parent node with no value but each entry in the dict
                               is a child. Otherwise it's treated as a value 
                               for the given key.
        key (str): Name of the value or section (e.g. "rdx" or "spectrograph"). If None, this is the
                   root entry which has no equivalent in Qt.
        parent (_UserConfigTreeNode): The parent of this entry. None if the entry is the root entry or
                                      a top level entry.
        
    """
    def __init__(self, node, key=None, parent=None):
        
        self.parent = parent
        self.key = key
        if isinstance(node, dict):
            self.children = [_UserConfigTreeNode(key=k, node=node[k], parent=self if key is not None else None) for k in node.keys()]
            self.value=None
        else:
            self.value=node
            self.children = []

    def getConfigLines(self, level=0):
        """Return the user's configuration lines starting at a given level.
        
        Args:
            level (int,optional): The level to provide configuration lines for, with level 0 being the top level.
            
        Return: 
            list[str]: 
                The configuration lines, starting at the given level and descending down 
                into all of the lower levels in the tree.
        """
        indent = " " * (level-1)
        if len(self.children) == 0:
            return [f"{indent}{self.key} = {self.value}"]
        else:
            if self.key is not None:                
                lines = [ f"{indent}{'['*level}{self.key}{']'*level}"]
            else:
                # Only the root node should have a None key
                lines = []

            for child in self.children:
                lines += child.getConfigLines(level+1)

            return lines

class PypeItParamsProxy(QAbstractItemModel):
    """
    A Proxy model that maps a PypeIt :class:`PypeItPar` to a QAbstractItemModel suitable for a Qt QTreeView.
    It uses the _UserConfigTreeNode recursive data structure to represent the parameters.

    Args:
        pypeit_setup (PypeItSetup): The PypeItSetup object containing the parameters to represent.

    """
    def __init__(self, par, cfg_lines=None):
        super().__init__(None)

        # TODO is this needed? Currently self.par isn't used but maybe it could be used
        # to show default values to clients? Or help info?
        self.par = par 
        if cfg_lines is None:
            # Create default config lines
            cfg_lines = ["[rdx]",f"    spectrograph = {par['rdx']['spectrograph']}"]

        self.cfg_lines = cfg_lines
        self._userConfigTree = _UserConfigTreeNode(ConfigObj(cfg_lines))

    def getConfigLines(self):
        return self.cfg_lines

    def rowCount(self, parent=QModelIndex()):
        """
        Returns the number of items under a parent node. Overridden from QAbstractItemModel.
        
        Returns:
            int: The number of items under parent. Can be 0 if parent has no child items.

        """
        if not parent.isValid():
            # If parent is invalid, it's the top level of the tree
            node = self._userConfigTree
        else:
            # Otherwise, if the parent is an index created by the index() method, it
            # points to the parent node in its internalPointer()
            if parent.column() == 1:
                # Column 1 does not have children
                return 0

            #msgs.info("rowCount valid")
            node = parent.internalPointer()

        return len(node.children)

    def index(self, row, column, parent=QModelIndex()):
        """
        Creates a QModelIndex that points to an item in the PypeItPar parameter tree. Overridden from QAbstractItemModel.
        
        Args:
            row    (int): The row of the item underneath parent to create an index for.

            column (int): The column of the item. Typically 0 is used for the name of the :class:`PypeItPar` or 
                          value, and 1 is used for actual parameter values (1 is never used for a :class:`PypeItPar`).

            parent (QModelIndex): An index for the parent. If this is invalid, it refers to the root of the tree.

        Returns:
            QModelIndex: A new index object to point to the item.
        
        """
        if not parent.isValid():
            # Use the root of the config tree as the parent
            parent_node = self._userConfigTree
        else:
            if parent.column() == 1:
                # Column one does not have children
                return QModelIndex()

            # The parent is valid, we're creating an index to one of its children
            # Get the parent from the internal data created by this method earlier.
            parent_node = parent.internalPointer()


        # Find the child using the passed in row
        child_node = parent_node.children[row]

        # Create the index, using child_node as the intenralPointer
        return super().createIndex(row, column, child_node)

    def data(self, index, role=Qt.DisplayRole):
        """Returns data for a given index. Overridden from QAbstractItemModel.

        Args:
            index (QModelIndex): The index of the item, as returned by index()

            role (Qt.ItemDataRole): The role of the data being requested. This method supports
                                    Qt.DisplayRole (for displaying text).
        
        
        Returns: 
            str: 
                A string if there's something to display at the given index and role, or
                None if there's nothing to display.
        """
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
        """Return data for the header row of the :class:`PypeItPar` tree.  For the horizontal
        header we call section 0 aka column 0 "Setting" and section 1 aka column 1 "Value".

        Args:
            section (int): The section aka column of the header.
            orientation(Qt.Orientation): The orientation of the header. This model only supports
                                         the horizontal header row.
            role (Qt.DisplayRole): The display role of the data to return. Only Qt.DisplayRole is supported
                                   for showing the text labels of the header.

        Returns:
            str: The name of a column or None if headerData isn't applicable for the
                 section/orientation/role passed in.
        
        """
        if role == Qt.DisplayRole and orientation == Qt.Horizontal:
            if section == 0:
                return "Setting"
            elif section == 1:
                return "Value"
        return super().headerData(section, orientation, role)


    def columnCount(self, parent):
        """Return the number of columns, which is always 2 because we use column 0 for the name of the
        parameter or :class:`PypeItPar` and column 1 for the configuration values. 
        Overridden from QAbstractItemModel.

        Args:
        parent (QModelIndex): Parent to return the column count for. (Not used)

        Return:
            int: The number of columns. Always 2.        
        """
        return 2

    def parent(self, index):
        """Return the parent of an item in the model. Overridden from QAbstractItemModel.

        Args:
            index (QModelIndex): The index of the item, as returned by the index() method.

        Returns: 
            QModelIndex: An index for the items parent. Maybe be invalid if the item is at the top level.
        """
        # This method is why the parent node is included in the _UserConfigTreeNode data.
        node = index.internalPointer()
        if node is None or node.parent is None:
            # Root node or top level node, there is no parent
            return QModelIndex()
        else:
            # Need to know the row # for the parent, which is in the grandparent
            grandparent = node.parent.parent
            if grandparent is None:
                # The grandparent is the root node
                grandparent = self._userConfigTree

            parent_and_siblings = [x.key for x in grandparent.children ]
            row = parent_and_siblings.index(node.parent.key)
            # Column is always 0, because 1 is reserved for leaf nodes which can't be parents
            return super().createIndex(row, 0, node.parent)


class PypeItFileModel(QObject):
    """
    A model representing the contents of a single .pypeit file. This involves a spectrograph, configuration values, 
    file metadata, and PypeIt parameters.

    Args: 
        pypeit_setup (PypeItSetup): The PypeItSetup object this configuration is a part of.
        config_name (str):  The name of this configuration.
        config_dict (dict): The metadata for this configuration, as returned by PypeItMetadata.unique_configurations.
        state (ModelState): The state of the model.  "CHANGED" if it was built by running setup or "UNCHANGED" if it came from loading
                            an existing .pypeit file.
    """

    stateChanged = Signal(str, ModelState)
    """Signal(str): Signal sent when the state of the file model changes. The config name of the file is sent as the signal parameter."""

    def __init__(self, metadata_model, state, name_stem, config_dict=None, params_model=None):
        super().__init__()
        self.state = state
        self.params_model=params_model
        self.metadata_model = None
        self.setMetadataModel(name_stem, config_dict, metadata_model)
        self.save_location = None

    def setMetadataModel(self, name_stem, config_dict, metadata_model):
        """Sets the metadata model containing the list of raw data files and associated metadata.
        
        Args:
            name_stem (str):      The name_stem of the metadata which will be used in the name of the pypeit file. This is typially the
                                  setup name (e.g. "A").
            config_dict (dict):   The metadata for this configuration, as returned by PypeItMetadata.unique_configurations.
            metadata_model (PypeItMetadataModel): The model object to set.            
        """
        if self.metadata_model is not None:
            self.metadata_model.dataChanged.disconnect(self._update_state)
            self.metadata_model.modelReset.disconnect(self._update_state)
            self.metadata_model.rowsInserted.disconnect(self._update_state)
            self.metadata_model.rowsRemoved.disconnect(self._update_state)            

        self.metadata_model = metadata_model
        self.name_stem = name_stem
        self.paths_model = metadata_model.getPathsModel()

        self._spectrograph = metadata_model.spectrograph

        if self.state == ModelState.NEW:
            self.config_values = {}
        else:            

            # Build a list of the configuration key/value pairs
            self.config_values = config_dict

        # Monitor the model for changes
        self.metadata_model.dataChanged.connect(self._update_state, Qt.ConnectionType.DirectConnection)
        self.metadata_model.modelReset.connect(self._update_state, Qt.ConnectionType.DirectConnection)
        self.metadata_model.rowsInserted.connect(self._update_state, Qt.ConnectionType.DirectConnection)
        self.metadata_model.rowsRemoved.connect(self._update_state, Qt.ConnectionType.DirectConnection)

    def _update_state(self):
        """Signal handler that detects changes to the metadata model and updates
        our model state to ModelState.CHANGED
        """
        msgs.info("Updating state")
        if self.state == ModelState.NEW:
            # Only move to "CHANGED" if there are rows in the metadata.
            if self.metadata_model.rowCount() > 0:
                self.state = ModelState.CHANGED
                # Update our configuration values
                self.config_values = self.metadata_model.getSetup(0)
        else:
            self.state = ModelState.CHANGED

        self.stateChanged.emit(self.name_stem,self.state)

    @property
    def filename(self):
        """str: The full path name of the pypeit file.
        """
        # TODO figure out when to use default name, vs a set name, vs "New Name.pypeit" or whatever.        
        # If name is key in parent model, how do we change it???
        return os.path.join("" if self.save_location is None else self.save_location, f"{self.spec_name}_{self.name_stem}.pypeit")

    @property
    def spec_name(self):
        """str: The name of the spectrograph for this pypeit file."""
        return self._spectrograph.name

    def save(self):
        """Save a .pypeit file from this model. The file name is chosen per PypeIt standards.
        """
        try:
            if self.metadata_model.metadata is None:
                metadata_table = None
                setup_dict = {}
            else:
                metadata_table = self.metadata_model.metadata.table[self.metadata_model.colnames]
                # Sorting is consistent with PypeItMetadata.write_pypeit
                metadata_table.sort(['frametype','filename'])

                # A pypeit file can only have one configuration. We choose the first one, if the user
                # has somehow created a file with multiple configurations, we hope they know what they're doing
                configs = self.metadata_model.metadata.unique_configurations()
                config_to_save = list(configs.keys())[0]
                setup_dict = {f'Setup {self.name_stem}':configs[config_to_save]}
    
            pf = PypeItFile(self.params_model.getConfigLines(),self.metadata_model.getPathsModel().getPaths(), metadata_table, setup_dict,vet=False,preserve_comments=True)    

            msgs.info(f"Saving filename: {self.filename}")
            if self.save_location is not None:
                os.makedirs(self.save_location,exist_ok=True)
            pf.write(self.filename) 
            
        except Exception as e:
            msgs.warn(f"Failed saving setup {self.name_stem} to {self.save_location}.")
            msgs.warn(traceback.format_exc())
            # Raise an exception that will look nice when displayed to the GUI
            raise RuntimeError(f"Failed saving setup {self.name_stem} to {self.save_location}.\nException: {e}")
        self.state = ModelState.UNCHANGED
        self.stateChanged.emit(self.name_stem, self.state)
        

class PypeItObsLogModel(QObject):
    """Model representing an obslog to the GUI.
    """

    spectrograph_changed = Signal(str)
    """Signal(str): Signal sent when the spectrgraph for the setup has changed. Sends the name of the new spectrograph."""

    def __init__(self):
        super().__init__()
        self._spectrograph = None
        self.metadata_model = PypeItMetadataModel(None)
        self.paths_model = QStringListModel()
        self.raw_data_files = []
        self.default_extension = ".fits"

    @property
    def state(self):
        """ModelState: The state of the obslog model. Either NEW or UNCHANGED."""
        if self.metadata_model.metadata is None:
            msgs.info("Obslog state is NEW")
            return ModelState.NEW
        else:
            msgs.info("Obslog state is UNCHANGED")
            return ModelState.UNCHANGED

    def set_spectrograph(self, new_spec):
        """Set the current spectrograph.

        Args:
            new_spec (str): The name of the new spectrograph.
        """
        msgs.info(f"Spectrograph is now {new_spec}")
        self._spectrograph = spectrographs.util.load_spectrograph(new_spec)
        if self.metadata_model.spectrograph is not None and self.metadata_model.spectrograph.name != new_spec:
            self.metadata_model.reset()
        self.spectrograph_changed.emit(self._spectrograph.name)

    def setMetadata(self, metadata):
        """Sets the obslog to reflect a new PypeItMetadata object."""
        self.metadata_model.setMetadata(metadata)
        self._spectrograph = self.metadata_model.spectrograph
        self.paths_model.setStringList(np.unique(metadata.table['directory']))
        self.spectrograph_changed.emit(self._spectrograph.name)

    @property
    def spec_name(self):
        """str: The name of the current spectrograph, or None if not set. """
        return None if self._spectrograph is None else self._spectrograph.name

    @property
    def spectrograph(self):
        """:obj:`pypeit.spectrographs.spectrograph`: The currently selected spectrograph, or None if not set."""
        return self._spectrograph

    @property
    def raw_data_directories(self):
        """list of str: The list directories containing raw data for the PypeItSetup object."""
        return self.paths_model.stringList()

    def add_raw_data_directory(self, new_directory):        
        """
        Adds a new directory to the model's list of directories.

        Args:
            new_directory (str): The new directory containing raw data.        
        """
        msgs.info(f"Adding raw directory: {new_directory}, current spec is {self._spectrograph}")
        if new_directory not in self.paths_model.stringList():
            row_number = self.paths_model.rowCount()
            self.paths_model.insertRows(row_number, 1)
            self.paths_model.setData(self.paths_model.index(row_number,0),new_directory)

    def scan_raw_data_directories(self):
        """
        Scans all of the raw data directories for raw data files.

        Returns:
            list[str]: The raw data files found.
        """
        allowed_extensions = self._spectrograph.allowed_extensions
        if allowed_extensions is None or len(allowed_extensions) == 0:
            # Most spectrographs don't set the allowed extensions, just use the
            # default from the command line. Append a "*" to match compressed files
            allowed_extensions = [self.default_extension + "*"]

        raw_data_files = []
        for directory in self.paths_model.stringList():
            msgs.info(f"Scanning directory: {directory}")
            for extension in allowed_extensions:
                #  The command line may set a root, which isn't a directory but a prefix
                if not os.path.isdir(directory):
                    glob_pattern = directory + "*" + extension
                else:
                    glob_pattern = os.path.join(directory, "*" + extension)

                msgs.info(f"Searching for raw data files with {glob_pattern}")
                raw_data_files += glob.glob(glob_pattern)

        return raw_data_files


    def reset(self):
        """Reset the model to an empty state."""

        msgs.info(f"Resetting to empty setup.")
        self.raw_data_files = []
        self.raw_data_dirs = []
        self.metadata_model.setMetadata(None)     
        self.paths_model.setStringList([])
        self._spectrograph = None
        self.spectrograph_changed.emit(None)

class PypeItSetupGUIModel(QObject):
    """
    Maintains the state of the overall PypeItSetupGUI.
    """
    filesDeleted = Signal(list)
    """Signal(list): Signal sent when PypeIt files within the GUI have been closed. Sends the list of
                     removed files."""

    filesAdded = Signal(list)
    """Signal(list): Signal sent when PypeIt files within the GUI have been created or opened. Sends the list of the
                     new files."""

    stateChanged = Signal()
    """Signal(): Signal sent when the state attribute of the model has changed."""

    def __init__(self):
        super().__init__()
        self._pypeit_setup = None
        self.pypeit_files = dict()
        self.log_buffer = None
        self.obslog_model = PypeItObsLogModel()
        self._clipboard = PypeItMetadataModel(None)

    def setup_logging(self, logname, verbosity):
        """
        Setup the PypeIt logging mechanism and a log buffer
        for monitoring the progress of operations and
        viewing the log.

        Args:
            logname (str): The filename to log to.
            verbosity (int): The verbosity level to log at.
        """
        self.log_buffer = LogBuffer(logname,verbosity)
        msgs.reset(verbosity=verbosity, log=self.log_buffer, log_to_stderr=False)
        msgs.info(f"QT Version: {qtpy.QT_VERSION}")
        msgs.info(f"PySide version: {qtpy.PYSIDE_VERSION}")
        msgs.info(f"PyQt version: {qtpy.PYQT_VERSION}")
        msgs.info(f"QtPy API_NAME: {qtpy.API_NAME}")

    @property
    def state(self):
        """ModelState: The state of the model.
        """
        if len(self.pypeit_files) == 0:
            # No pypeit files, just use the obslog state
            msgs.info(f"GUI state is {self.obslog_model.state}")
            return self.obslog_model.state
        else:
            if any([file.state!=ModelState.UNCHANGED for file in self.pypeit_files.values()]):
                msgs.info("GUI state is CHANGED")
                return ModelState.CHANGED
            else:
                msgs.info("GUI state is UNCHANGED")
                return ModelState.UNCHANGED

    def reset(self):
        """Reset the model to an empty state."""

        msgs.info(f"Resetting to NEW state.")
        self.closeAllFiles()
        self.obslog_model.reset()
        self.stateChanged.emit()

    @property
    def clipboard(self):
        """PypeItMetadataModel: The current metadata that has been copied or pasted within the Setup GUI."""
        return self._clipboard
    
    @clipboard.setter
    def clipboard(self,value):
        # Don't overwrite the clipboard instance, so that views can keep their signals attached
        if value is None or value.metadata is None:
            self._clipboard.setMetadata(None)
        else:
            self._clipboard.setMetadata(value.metadata)

    def run_setup(self):
        """Run setup on the raw data in the obslog."""

        raw_data_files =  self.obslog_model.scan_raw_data_directories()
        
        if len(raw_data_files) == 0:
            raise ValueError("Could not find any files in any of the raw data paths.")

        self._pypeit_setup = PypeItSetup.from_rawfiles(raw_data_files, self.obslog_model.spec_name)

        # These were taken from the default parameters in pypeit_obslog
        self._pypeit_setup.run(setup_only=True,
                               groupings=True,
                               clean_config=False)

        self.obslog_model.setMetadata(self._pypeit_setup.fitstbl)
        self.createFilesForConfigs()
        self.stateChanged.emit()

    def open_pypeit_file(self, pypeit_file):
        """Open an existing pypeit file and load it into this setup.

        Args:
            pypeit_file (str): The pypeit file to open.
        """
        pf = PypeItFile.from_file(pypeit_file, preserve_comments=True, vet=False)

        if pf.config is not None:
            try:
                spec_name = pf.config['rdx']['spectrograph']
            except:
                raise ValueError(f"PypeIt input file {pypeit_file} is missing a spectrograph.")
        else:
            raise ValueError(f"PypeIt input file {pypeit_file} is missing a configuration section with a spectrograph specified.")

        if spec_name not in available_spectrographs:
            raise ValueError(f"PypeIt input file {pypeit_file} contains an unknown spectrograph name: '{spec_name}'")


        if pf.data is not None:
            # If the data block is fully setup, it won't have all of the metadata fields originally created by
            # PypeItSetup, so re-create it

            filenames = pf.filenames
            directories = [str(Path(f).parent) for f in filenames ]

            # Add in extra columns to the metadata that PypeItMetadata needs
            pf.data.add_column(directories, name='directory')
            pf.data.add_column([pf.setup_name] * len(pf.data), name='setup')
            spec, par, config_specific_file = pf.get_pypeitpar()

            # Set configuration key columns not already in the metadata
            for key in spec.configuration_keys():
                if key not in pf.data.colnames:
                    # Some keys may allow None as a valid value, so we don't treat a missing
                    # key as an error
                    value = pf.setup.get(key,None)
                    pf.data.add_column([value]*len(pf.data),name=key)

            # Now set types based on core metadata
            core_metadata_model = pypeit.core.meta.get_meta_data_model()
            for colname in pf.data.colnames:
                # Treat columns not in the core metadata as strings, (these can be things like 'filename', 'frametype', etc)
                if colname not in core_metadata_model:
                    col_dtype = str
                else:
                    col_dtype = core_metadata_model[colname]['dtype']
                if col_dtype is not str:
                    # input files will in non-existent values as blank strings, convert these to None
                    string_values = [True if isinstance(v,str) else False for v in pf.data[colname]]
                    pf.data[colname][string_values] = None
                # Convert the object dtypes to the correct type
                pf.data[colname] = astropy.table.Column(data=pf.data[colname],dtype=col_dtype)

            # Build the PypeItMetaData object with the converted data from the input file
            metadata = PypeItMetaData(spec, par, data=pf.data, strict=False)
            # Make sure the frame type bits get set
            user_frametypes = {filename: frametype for filename, frametype in pf.data['filename','frametype']}
            metadata.get_frame_types(user=user_frametypes)

            self._pypeit_setup = PypeItSetup(filenames,
                                             setups=[pf.setup_name],
                                             cfg_lines=pf.cfg_lines,
                                             pypeit_file=pypeit_file)
            self._pypeit_setup.fitstbl = metadata

            self.obslog_model.setMetadata(metadata)        
            self.createFilesForConfigs(state=ModelState.UNCHANGED)
            self.stateChanged.emit()
        else:
            # Create an empty file and add what information (paths, name, parameters) we do have
            self.obslog_model.reset()
            self.obslog_model.set_spectrograph(spec_name)
            if len(pf.file_paths) > 0:
                for path in len(pf.file_paths):
                    self.obslog_model.add_raw_data_directory(path)
                
            params = PypeItParamsProxy(self.obslog_model.spectrograph.default_pypeit_par(),pf.cfg_lines)
            if pf.setup is not None and len(pf.setup) != 0:
                setup_name = pf.setup_name
            else:
                setup_name = "A"
            empty_metadata = self.obslog_model.metadata_model.createCopyForRows([])
            pf_model = PypeItFileModel(empty_metadata,
                                       ModelState.NEW,
                                       setup_name,
                                       params_model=params)

            pf_model.stateChanged.connect(self.stateChanged)

            self.pypeit_files[setup_name] = pf_model            

            self.filesAdded.emit([pf_model])
            self.stateChanged.emit()


    def removeFile(self, name):
        del self.pypeit_files[name]
        self.filesDeleted.emit([name])

    def createEmptyPypeItFile(self, new_name):
        # Create a new empty configuration.
        msgs.info(f"Creating new pypeit file model for {new_name}")

        # Create an empty copy of the obslog metadata for the new file
        empty_metadata = self.obslog_model.metadata_model.createCopyForRows([])
        default_params = PypeItParamsProxy(self.obslog_model.spectrograph.default_pypeit_par())
        
        pf_model = PypeItFileModel(empty_metadata,
                                   ModelState.NEW,
                                   new_name,
                                   params_model=default_params)

        pf_model.stateChanged.connect(self.stateChanged)

        self.pypeit_files[new_name] = pf_model            

        self.filesAdded.emit([pf_model])
        self.stateChanged.emit()
        return pf_model

    def closeAllFiles(self):
        # Delete all files
        deleted_files = list(self.pypeit_files.keys())
        self.pypeit_files = dict()
        if len(deleted_files) > 0:
            self.filesDeleted.emit(deleted_files)

    def createFilesForConfigs(self, configs=None, state=ModelState.CHANGED):
        """
        Private method to create PypeIt files for new configurations.

        Args:
            unique_config (dict): A the new configurations to set.
            pypeit_setup (:obj:`PypeItSetup`): The PypeItSetup object with parameter info for the new files.
            state (ModelState, Optional): The state of the new model. Defaults to CHANGED, meaning it has not been saved.

        """
        if configs is None:
            configs = self.obslog_model.metadata_model.metadata.unique_configurations()
            msgs.info(f"Creating files for all unique configurations: {configs}")
        else:
            msgs.info(f"Creating files for configs {configs}")

        # Create a new PypeItFileModel for each unique configuration
        config_names = list(configs.keys())
        for config_name in config_names:
            pf_model = PypeItFileModel(metadata_model=self.obslog_model.metadata_model.createCopyForConfig(config_name),
                                       state=state,
                                       name_stem=config_name,
                                       config_dict=configs[config_name],
                                       params_model=PypeItParamsProxy(self._pypeit_setup.par, self._pypeit_setup.user_cfg),
                                       )

            # Any change to the file models can change the over all state, so connect the signals
            pf_model.stateChanged.connect(self.stateChanged)
            self.pypeit_files[config_name] = pf_model

        msgs.info(f"Current files: {self.pypeit_files}")
        if len(config_names) > 0:
            self.filesAdded.emit(list(self.pypeit_files.values()))

