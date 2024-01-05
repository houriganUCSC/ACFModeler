import sys
import matplotlib.pyplot as plt
import numpy as np
import re
from PyQt6.QtWidgets import QApplication, QWidget, QTableWidget,\
    QVBoxLayout, QTableWidgetItem, QSizePolicy
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont
from sklearn.preprocessing import PolynomialFeatures
import statsmodels.api as sm
from src.records.Session import Session


plusminus = u"\u00B1"

class TextItem(QTableWidgetItem):
    def __init__(self, row, col, value):
        super().__init__()
        if type(value) == str:
            self.setText(str(value))
        else:
            self.setText(value)
        self.setProperty('row', row)
        self.setProperty('col', col)
        self.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
        fnt = QFont()
        fnt.setPixelSize(12)
        self.setFont(fnt)
        self.setFlags(self.flags() & ~Qt.ItemFlag.ItemIsEditable)
        self.row = -1
        self.col = -1

    def setProperty(self, name, value: int):
        if name == 'row':
            self.row = value
        elif name == 'col':
            self.col = value

    def property(self, name):
        if name == 'row':
            value = self.row
        elif name == 'col':
            value = self.col
        return value

class IsotopeFitHeaderCodes(dict):
    '''
    Class object that contains formatting information for the QTableWidget
    "hexCode":  is
    '''
    def __init__(self):
        self.codes = {
            'Mass': {'hexCode': 'Mass', "numFormat": '0.2f', "width": 60},
            'n Used': {'hexCode': 'n Used', "numFormat": 'd', "width": 75},
            'n': {'hexCode': 'n', "numFormat": 'd', "width": 75},
            'an Only':{'hexCode': 'an', "numFormat": 'd', "width": 60},
            'a1': {'hexCode': b"\x61\xE2\x82\x81", "numFormat": '.4G', "width": 75},
            'se_a1': {'hexCode': b"\xCF\x83\x28\x61\xE2\x82\x81\x29", "numFormat": '.4G', "width": 75},
            'a2': {'hexCode': b"\x61\xE2\x82\x82", "numFormat": '.4G', "width": 75},
            'se_a2': {'hexCode': b"\xCF\x83\x28\x61\xE2\x82\x82\x29", "numFormat": '.4G', "width": 75},
            'tau': {'hexCode': b"\xCF\x84", "numFormat": '.4G', "width": 75},
            'se_tau': {'hexCode': b"\xCF\x83\x28\xCF\x84\x29", "numFormat": '.4G', "width": 75},
            'rSqr': {'hexCode': b"\x72\xC2\xB2", "numFormat": '0.3f', "width": 50},
            "redChi2": {'hexCode': b"\xCF\x87\xE1\xB5\xA3\xC2\xB2", "numFormat": '0.2f', "width": 50},
            "ACF": {'hexCode': "ACF", "numFormat": '0.1f', "width": 75},
            "Drift": {'hexCode': "Drift", "numFormat": '0.1%', "width": 75},
            "dtCorrPct": {'hexCode': 'tau %', "numFormat": '0.1%', "width": 75},
            "max P": {'hexCode': 'max P', "numFormat": '0.2G', "width": 75},
            }

    def getHeaderLabel(self, key:str):
        hdr = self.codes[key]['hexCode']
        if isinstance(hdr, bytes):
            label = hdr.decode('UTF-8')
        else:
            label = hdr
        return label

    def getDefaultWidth(self, key:str):
        return self.codes[key]["width"]

    def formatValue(self, key:str, val: float):
        valFormat = self.codes[key]['numFormat']
        if np.isnan(val) or val is None:
            valText = ''
        else:
            valText = f'{val:{valFormat}}'
            if 'G' in valFormat:
                valText = valText.replace("E-0", "E-").replace("E+0", "E+")

        return valText

class IsotopeFitTableWidget(QTableWidget):
    """
    Populate table
    :param ccdf:
    :return:

    """
    def __init__(self, session: Session):
        self.colHeaders = ['Mass', 'n', 'n Used', 'an Only','max P','a1', 'se_a1', 'a2', 'se_a2', 'tau', 'se_tau',
                           'rSqr', 'redChi2', 'ACF', 'Drift','dtCorrPct']
        self.session = session
        rows = len(self.session.isotopes)
        cols = len(self.colHeaders)
        super().__init__(rows, cols)
        self.hdr = IsotopeFitHeaderCodes()
        fnt = QFont()
        fnt.setPixelSize(12)
        fnt.setBold(True)
        self.pars = ['a1', 'a2', 'tau']
        # if self.session.isotopes == []:
        #     self.isotopes = list(self.session.isotopes)
        #     self.setVerticalHeaderLabels(self.session.isotopes)
        self.setVerticalHeaderLabels(self.session.isotopes)
        self.setHorizontalHeaderLabels(self.colHeaders)
        for i,colName in enumerate(self.colHeaders):
            label = self.hdr.getHeaderLabel(colName)
            self.setHorizontalHeaderItem(i, QTableWidgetItem(label))
            self.horizontalHeaderItem(i).setFont(fnt)
            self.setColumnWidth(i, self.hdr.getDefaultWidth(colName))

        self.verticalHeader().setDefaultSectionSize(25)

        self.setStyleSheet('selection - color: rgb(0, 0, 0);'
                           'selection - background - color: rgb(230, 230, 230);'
                           'background - color: rgb(255, 255, 255)'
                           'color: rgb(0, 0, 0);')
        self.populateFilterTable(self.session)

    def initFitDictionary(self):
        pass

    def highlightChangedRow(self):
        self.selectRow(self.sender().property('row'))
        self.extractFitVals()

    def setComboBoxCell(self):
        source = self.sender()
        row = source.property('row')
        col = source.property('col')
        self.setCurrentCell(row, col)

    def populateFilterTable(self, session: Session):
        print(__name__, session.isotopes)
        self.setVerticalHeaderLabels(session.isotopes)
        for massName, massRecord in session.masses.items():
            print(__name__, massName, )
            row = np.where(session.isotopes == massName)[0][0]
            for col, colName in enumerate(self.colHeaders):
                if colName == 'Mass':
                    item = QTableWidgetItem(self.hdr.formatValue(colName, massRecord.aveMass))
                    item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                    self.setItem(row, col, item)
                elif colName in ['n', 'n Used', 'an Only', 'max P']:
                    if colName == 'n': valText = f'{massRecord.nObs:d}'
                    elif colName == 'n Used':valText = f'{massRecord.nIn:d}'
                    elif colName == 'an Only': valText = f'{massRecord.anOnly:d}'
                    elif colName == 'max P': valText = f'{massRecord.maxP: 0.3G}'.replace("E+0","E+")
                    item = QTableWidgetItem(valText)
                    item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                    self.setItem(row, col, item)
                else: self.setItem(row, col, QTableWidgetItem(''))


    def unpopulateFilterTable(self):
        for row in range(self.rowCount()):
            for col in range(self.columnCount()):
                self.setItem(row, col, QTableWidgetItem(''))

    def populateFitTable(self, session: Session):
        for massName, massRecord in session.masses.items():
            for colName in ['a1', 'se_a1', 'a2', 'se_a2', 'tau', 'se_tau', 'rSqr', 'redChi2', 'ACF', 'Drift', 'dtCorrPct']:
                col = self.colHeaders.index(colName)
                row = np.where(session.isotopes == massName)[0][0]
                val = massRecord.fits["self"][colName]
                valText = self.hdr.formatValue(colName, val)
                item = QTableWidgetItem(valText)
                item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                self.setItem(row, col, item)

    def unpopulateFitTable(self):
        for colName in ['a1', 'se_a1', 'a2', 'se_a2', 'tau', 'se_tau', 'rSqr', 'redChi2', 'ACF', 'Drift', 'dtCorrPct']:
            col = self.colHeaders.index(colName)
            for row in range(self.rowCount()):
                self.setItem(row, col, QTableWidgetItem(''))

class AppDemo(QWidget):
    def __init__(self):
        super().__init__()
        isotopes = ['Hg202', 'Pb204', 'Pb206', 'Pb207', 'Pb208', 'Th232', 'U235', 'U238']
        mainLayout = QVBoxLayout()
        session = Session()
        session.isotopes = isotopes
        self.table = IsotopeFitTableWidget(session)
        self.table.initFitDictionary()
        self.setMinimumSize(1375, 26*(len(isotopes)+1))
        self.setSizePolicy(QSizePolicy.Policy.Minimum, QSizePolicy.Policy.MinimumExpanding)
        mainLayout.addWidget(self.table)
        mainLayout.setContentsMargins(5,5,5,5)
        self.setLayout(mainLayout)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    demo = AppDemo()
    demo.show()
    app.exit(app.exec())



