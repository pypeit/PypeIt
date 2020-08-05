# Interactive Table Editor
import tkinter as tk
from tkinter.ttk import Treeview

from astropy.table import Table
import matplotlib.pyplot as plt

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt

from astropy.coordinates import Angle, SkyCoord, EarthLocation, AltAz, ICRS
from astropy.time import Time
import astropy
from astropy import units as u

# Some Thoughts on how this could work
# have modified table and path variable
# on exit, for each modified row, edit the corresponding FITS file


# TODO: allow for deleting row. low priority

loc = EarthLocation.of_site('Palomar')

def update_airmass(row: astropy.table.Row):
    skycoord = SkyCoord(row['ra'], row['dec'], unit=u.deg)
    time = Time(row['mjd'], format='mjd')
    altaz = AltAz(obstime=time, location=loc)
    row['airmass'] = skycoord.transform_to(altaz).secz

def get_zenith_ra_dec(time) -> SkyCoord:
    time = Time(time, format='mjd')
    altaz = AltAz(alt=Angle(90, unit=u.deg), az=Angle(0, unit=u.deg), obstime=time, location=loc)
    return altaz.transform_to(ICRS)

class TableModel(QtCore.QAbstractTableModel):
    def __init__(self, data: Table, cols = None):
        super(TableModel, self).__init__()
        self._data = data #Table(data, masked = True, copy=False)
        self._cols = cols if cols is not None else data.colnames
        self._mask = self._data.mask
        self._rowCount = len(self._data)
        self._colCount = len(self._cols)
    
    def data(self, index, role):
        if role == Qt.DisplayRole or role == Qt.EditRole:
            col = self._cols[index.column()]
            value = self._data[col][index.row()]
            if col == 'ra':
                ra = Angle(value, unit=u.deg)
                value = ra.to_string(unit=u.hour)
            elif col == 'dec':
                dec = Angle(value, unit=u.deg)
                value = dec.to_string()
            elif col == 'dispangle':
                ang = Angle(value, unit=u.deg)
                value = ang.to_string()
            elif col == 'mjd' or col == 'airmass':
                value = f'{value:.4f}'
            return str(value)
        if role == Qt.BackgroundRole:
            col = self._cols[index.column()]
            masked = self._mask[col][index.row()]
            if masked:
                return QtGui.QColor(Qt.red)
    
    def setData(self, index: QtCore.QModelIndex, value: QtCore.QVariant, role = Qt.EditRole):
        try:
            col = self._cols[index.column()]
            if col == 'ra':
                ra = Angle(value, unit=u.hour)
                self._data[col][index.row()] = ra.degree
                if self._data['dec'][index.row()]:
                    update_airmass(self._data[index.row()])
                self.dataChanged.emit(index, index)
                return True
            if col == 'dec':
                dec = Angle(value, unit=u.deg)
                self._data[col][index.row()] = dec.degree
                if self._data['dec'][index.row()]:
                    update_airmass(self._data[index.row()])
                self.dataChanged.emit(index, index)
                return True
            self._data[col][index.row()] = value
            self.dataChanged.emit(index, index)
            return True
        except:
            return False
    
    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        return QtCore.QAbstractTableModel.flags(self, index) | Qt.ItemIsEditable
    
    def rowCount(self, index):
        return self._rowCount
    
    def columnCount(self, index):
        return self._colCount
    
    def headerData(self, section, orientation, role):
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return self._cols[section]
            if orientation == Qt.Vertical:
                return ''
    
    def frametype(self, index: QtCore.QModelIndex):
        return self._data['frametype'][index.row()]
    
    def _set_data_to_zenith(self, index: QtCore.QModelIndex):
        row = index.row()
        time = self._data['mjd'][row]
        ra_dec = get_zenith_ra_dec(time)
        self._data['ra'][row] = ra_dec.ra.deg
        self._data['dec'][row] = ra_dec.ra.deg
        self._data['airmass'][row] = 1.0



class TableView(QtWidgets.QTableView):
    def __init__(self, model: QtCore.QAbstractTableModel, parent=None):
        super().__init__(parent)
        self.setModel(model)

        self.delegate = Delegate(self, model._cols)
        self.setItemDelegate(self.delegate)
        
        self.verticalHeader().setVisible(False)
        self.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeToContents)
        self.setEditTriggers(QtWidgets.QAbstractItemView.DoubleClicked)
        self.setContextMenuPolicy(Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self.show_context_menu)
        self.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)

    def show_context_menu(self, point: QtCore.QPoint):
        index = self.indexAt(point)
        if index.isValid():
            self._menu = QtWidgets.QMenu(self)
            frametype = self.model().frametype(index)
            if ('bias' in frametype) or ('arc' in frametype) or ('flat' in frametype):
                action = self._menu.addAction('Set RA/Dec and Airmass to Zenith')
                action.triggered.connect(lambda: self.model()._set_data_to_zenith(index))
                self._menu.popup(self.viewport().mapToGlobal(point))


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, data: Table, cols = None):
        super().__init__()

        self.table = TableView(TableModel(data, cols))
        
        self.setCentralWidget(self.table)
        
        
        self.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)


class Delegate(QtWidgets.QStyledItemDelegate):
    def __init__(self, parent, cols: tuple):
        super().__init__(parent)
        self._cols = cols
    
    def createEditor(self, parent, option, index):
        editor = QtWidgets.QLineEdit(parent)
        editor.setFrame(False)
        return editor
    
    def setEditorData(self, editor: QtWidgets.QLineEdit, index: QtCore.QModelIndex):
        value = index.model().data(index, Qt.EditRole)
        editor.setText(value)
    
    def setModelData(self, editor: QtWidgets.QLineEdit, model: QtCore.QAbstractTableModel, index: QtCore.QModelIndex):
        value = editor.text()
        model.setData(index, value)
    
    def updateEditorGeometry(self, editor: QtWidgets.QLineEdit, option: QtWidgets.QStyleOptionViewItem, index: QtCore.QModelIndex):
        editor.setGeometry(option.rect)
        

if __name__ == '__main__':
    blue_table = Table.read('/Users/milan/Documents/GitHub/PypeIt/table_red.dat', format='ascii')
    
    cols = ('filename', 'frametype', 'ra', 'dec', 'target', 'dispname', 'binning', 'mjd', 'airmass', 'exptime', 'dispangle', 'dichroic', 'slitwid')
    app = QtWidgets.QApplication([])
    window = MainWindow(blue_table, cols)
    window.show()
    app.exec_()
    # TODO:
    #   - save on quit
    #   - default size / column sizing