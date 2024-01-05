import numpy as np
from PyQt6.QtWidgets import QApplication, QWidget, QTableWidget, QCheckBox, QComboBox, \
    QVBoxLayout, QSizePolicy, QMainWindow, QHBoxLayout
from PyQt6.QtGui import QFont
from PyQt6.QtCore import Qt, pyqtSignal, pyqtSlot

class FitTypes(QComboBox):
    """
    QCombobox selector for type of fit to use for cross calibration
    'Mean'
    'Linear'
    '2nd Order'
    '3rd Order'
    'Cubic Spline': Not currently active

    """
    def __init__(self, row, col, parent):
        super().__init__(parent)
        self.setStyleSheet('font-size: 10px;'
                           'selection - color: rgb(0, 0, 0);'
                           'selection - background - color: rgb(230, 230, 230);')
        self.addItems(['Mean', 'Linear', '2nd Order', '3rd Order', 'Cubic Spline'])
        self.setProperty('row', row)
        self.setProperty('col', col)

        # Todo: Cubic spline fit in SpectrumFit class
        self.model().item(4).setEnabled(False)

    def getValue(self):
        return self.currentIndex()


class CustomCheck(QWidget):
    def __init__(self, row, col, parent):
        super().__init__(parent)
        self.cbox = QCheckBox()
        self.cbox.setStyleSheet('background-color: transparent')
        self.cbox.setProperty('row', row)
        self.cbox.setProperty('col', col)
        layout = QHBoxLayout(self)
        layout.addWidget(self.cbox)
        layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.setContentsMargins(0,0,0,0)
        self.setLayout(layout)

    def isChecked(self):
        return int(self.cbox.isChecked())


class ModelHeaderCodes():
    def __init__(self):
        self.codes = {
            'ACF (a1)': b"\x41\x43\x46\x20\x28\x61\xE2\x82\x81\x29",
            'Drift (a2)': b"\x44\x72\x69\x66\x74\x20\x28\x61\xE2\x82\x82\x29",
            'Tau': b"\x54\x61\x75\x20\x28\xCF\x84\x29",
            }

    def getHeaderLabel(self, key:str):
        hdr = self.codes[key]
        if isinstance(hdr, bytes):
            label = hdr.decode('UTF-8')
        else:
            label = hdr
        return label


class ModelDesignTable(QTableWidget):
    tableModified = pyqtSignal()
    rowIDs = ['a1', 'a2', 'tau']
    colIDs = ['order', 'weighted', 'robust']

    def __init__(self):
        super().__init__(3,3)
        fnt = QFont()
        fnt.setPixelSize(12)
        fnt.setBold(True)
        self.rowNames = ['ACF (a1)', 'Drift (a2)', 'Tau']
        self.rowIDs = ['a1', 'a2', 'tau']
        self.setHorizontalHeaderLabels(['Function', 'Weight', 'Robust'])
        self.setVerticalHeaderLabels(self.rowNames)
        vHdr = ModelHeaderCodes()
        for row, rowName in enumerate(self.rowNames):
            label = vHdr.getHeaderLabel(rowName)
            for col in range(3):
                if col == 0:
                    selector = FitTypes(row, col, self)
                    selector.setCurrentIndex(1)
                    selector.currentIndexChanged.connect(self.emitTableModified)
                    self.setCellWidget(row, 0, selector)
                    self.setColumnWidth(col, 110)
                else:
                    ccbox = CustomCheck(row, col, self)
                    #Todo: Have to debug robust fitting first
                    ccbox.cbox.setEnabled(col<2)
                    ccbox.cbox.released.connect(self.emitTableModified)
                    self.setCellWidget(row, col, ccbox)
                    self.setColumnWidth(col, 65)
                    if col == 1:
                        ccbox.cbox.setChecked(True)
                self.horizontalHeaderItem(col).setFont(fnt)
            self.verticalHeaderItem(row).setText(label)
            self.verticalHeaderItem(row).setFont(fnt)

    def emitTableModified(self):
        source = self.sender()
        row = source.property('row')
        col = source.property('col')
        design = self.getFitPars()
        self.tableModified.emit()

    def getFitPars(self):
        design = {'a1': {'order': None, 'weighted': None, 'robust': None},
                  'a2': {'order': None, 'weighted': None, 'robust': None},
                  'tau': {'order': None, 'weighted': None, 'robust': None}}
        for row, (rowkey, rowval) in enumerate(design.items()):
            for col, colkey in enumerate(rowval.keys()):
                widget = self.cellWidget(row, col)
                if colkey == 'order':
                    design[rowkey][colkey] = widget.currentIndex()
                else:
                    design[rowkey][colkey] = widget.isChecked()
        return design


class AppDemo(QWidget):
    def __init__(self):
        super().__init__()
        mainLayout = QVBoxLayout()
        self.table = ModelDesignTable()
        self.table.getFitPars()
        mainLayout.addWidget(self.table)
        self.setMinimumSize(305, 120)
        self.setMaximumSize(305, 120)
        self.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Maximum)
        mainLayout.setContentsMargins(0,0,0,0)
        self.setLayout(mainLayout)



if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    demo = AppDemo()
    demo.show()
    app.exit(app.exec())