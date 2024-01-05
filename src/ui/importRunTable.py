from PyQt6.QtWidgets import (QApplication, QDialog, QTableWidget, QTableWidgetItem, QDoubleSpinBox, QWidget,
                             QGridLayout, QComboBox, QLabel, QPushButton)
from PyQt6.QtCore import Qt, QMutex
from src.records.Session import Session
import numpy as np


class IsotopeIDQuery(QDialog):

    def __init__(self, session: Session, isotopes:list = None, truncMasses = None, n = None, parent = None):
        super(IsotopeIDQuery, self).__init__(parent)
        self.session = session
        self.setWindowTitle("Analysis Configuration")
        self.setFixedWidth(400)
        grid = QGridLayout()
        if isotopes is None:
            isotopes = list(self.session.masses.keys())
            truncMasses = [mass.aveMass for mass in self.session.masses.values()]
            n = int(len(self.session.scanTime) / len(self.session.samples))
        self.table = IsotopeTable(isotopes, truncMasses)
        grid.addWidget(self.table, 0, 0, 1, 3)

        self.dtSpinner = QDoubleSpinBox()
        self.dtSpinner.setRange(0,200)
        self.dtSpinner.setValue(20.0)
        grid.addWidget(self.dtSpinner, 1, 1, 1, 1)

        dtLabel = QLabel()
        dtLabel.setText("Detector Deadtime:")
        grid.addWidget(dtLabel, 1, 0, 1, 1)
        dtUnits = QLabel()
        dtUnits.setText("nanoseconds")
        grid.addWidget(dtUnits, 1, 2, 1, 1)

        runsLabel = QLabel()
        runsLabel.setText("Runs:")
        grid.addWidget(runsLabel, 2, 0, 1, 1)
        self.runsValue = QLabel()
        self.runsValue.setText(f'{n:d}')
        grid.addWidget(self.runsValue, 2, 1, 1, 1)
        passesLabel = QLabel()
        passesLabel.setText("Passes:")
        grid.addWidget(passesLabel, 3, 0, 1, 1)
        self.passesValue = QLabel()
        if n == 0:
            passesText = "0"
        else:
            passesText = "1"
        self.passesValue.setText(passesText)
        grid.addWidget(self.passesValue, 3, 1, 1, 1)
        cyclesLabel = QLabel()
        cyclesLabel.setText("Cycles: ")
        grid.addWidget(cyclesLabel, 4, 0, 1, 1)
        self.cyclesValue = QLabel()
        self.cyclesValue.setText(f'{n:d}')
        grid.addWidget(self.cyclesValue, 4, 1, 1, 1)

        Ok = QPushButton()
        Ok.setText("PROCEED")
        Ok.clicked.connect(self.getMetaData)
        grid.addWidget(Ok, 5, 1, 1,1)

        self.setLayout(grid)
        self.setWindowModality(Qt.WindowModality.ApplicationModal)

    def getMetaData(self):
        isotopes = self.table.getIsotopes()
        self.session.method = {'isotopes': isotopes,
                               'runs': int(self.runsValue.text()),
                               'passes': int(self.passesValue.text()),
                               'cycle': int(self.cyclesValue.text()),
                               'masses': len(isotopes),
                               'deadTime': self.dtSpinner.value()}
        self.close()

class IsotopeTable(QTableWidget):
    def __init__(self, isotopes:list, truncMasses):
        super().__init__()
        if isinstance(isotopes[0], int):
            isotopes = [str(isotope) for isotope in isotopes]
        self.setVerticalHeaderLabels(isotopes)
        self.setRowCount(len(isotopes))
        self.setColumnCount(3)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        for row, mass in enumerate(truncMasses):
            self.setRowHeight(row, 15)
            label = QLabel()
            label.setText(f'{mass:0.2f}')
            label.setAlignment(Qt.AlignmentFlag.AlignCenter)
            self.setCellWidget(row,0, label)
            drop = ElementDrop()
            drop.massIndex(mass)
            drop.currentIndexChanged.connect(self.updateID)
            ID = QLabel()
            ID.setText(f'{mass:0.0f}{drop.currentText()}')
            ID.setAlignment(Qt.AlignmentFlag.AlignCenter)
            self.setCellWidget(row, 1, drop)
            self.setCellWidget(row, 2, ID)

        self.setColumnWidth(0, 60)
        self.setColumnWidth(1, 60)
        self.setHorizontalHeaderLabels(['Mass', 'Element', 'ID'])

    def updateID(self):
        for row in range(self.rowCount()):
            mass = float(self.cellWidget(row, 0).text())
            el = self.cellWidget(row, 1).currentText()
            self.cellWidget(row,2).setText(f'{mass:0.0f}{el}')

    def getIsotopes(self):
        isotopes = []
        for row in range(self.rowCount()):
            isotopes.append(self.cellWidget(row, 2).text())
        return isotopes

class ElementDrop(QComboBox):
    Elements = {
        "H": 1,
        "He": 4,
        "Li": 7,
        "Be": 9,
        "B": 11,
        "C": 12,
        "N": 14,
        "O": 16,
        "F": 19,
        "Ne": 20,
        "Na": 23,
        "Mg": 24,
        "Al": 27,
        "Si": 28,
        "P": 31,
        "S": 32,
        "Cl": 35,
        "K": 39,
        "Ar": 40,
        "Ca": 40,
        "Sc": 45,
        "Ti": 48,
        "V": 51,
        "Cr": 52,
        "Mn": 55,
        "Fe": 56,
        "Ni": 59,
        "Co": 59,
        "Cu": 64,
        "Zn": 65,
        "Ga": 70,
        "Ge": 73,
        "As": 75,
        "Se": 79,
        "Br": 80,
        "Kr": 84,
        "Rb": 85,
        "Sr": 88,
        "Y": 89,
        "Zr": 91,
        "Nb": 93,
        "Mo": 96,
        "Tc": 98,
        "Ru": 101,
        "Rh": 103,
        "Pd": 106,
        "Ag": 108,
        "Cd": 112,
        "In": 115,
        "Sn": 119,
        "Sb": 122,
        "I": 127,
        "Te": 128,
        "Xe": 131,
        "Cs": 133,
        "Ba": 137,
        "La": 139,
        "Ce": 140,
        "Pr": 141,
        "Nd": 144,
        "Pm": 145,
        "Sm": 150,
        "Eu": 152,
        "Gd": 157,
        "Tb": 159,
        "Dy": 163,
        "Ho": 165,
        "Er": 167,
        "Tm": 169,
        "Yb": 173,
        "Lu": 175,
        "Hf": 178,
        "Ta": 181,
        "W": 184,
        "Re": 186,
        "Os": 190,
        "Ir": 192,
        "Pt": 195,
        "Au": 197,
        "Hg": 201,
        "Tl": 204,
        "Pb": 207,
        "Bi": 209,
        "Po": 209,
        "At": 210,
        "Rn": 222,
        "Fr": 223,
        "Ra": 226,
        "Ac": 227,
        "Pa": 231,
        "Th": 232,
        "Np": 237,
        "U": 238,
        "Pu": 242
    }
    def __init__(self):
        super().__init__()
        self.addItems(list(self.Elements.keys()))

    def massIndex(self, mass):
        massSeries = list(self.Elements.values())
        massSeries = np.array(massSeries, dtype=int)
        massKeys = list(self.Elements.keys())
        if mass not in massSeries:
            dev = abs(massSeries - mass)
            key = massKeys[np.argmin(dev)]
            idx = np.argmin(dev)
        else:
            idx = list(self.Elements.values()).index(mass)
            key = massKeys[list(self.Elements.values()).index(mass)]
        self.setCurrentIndex(idx)


if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    sesh = Session()
    isotopes = ['1', '2', '3']
    masses = [203.97, 206.98, 237.8]
    ui = IsotopeIDQuery(sesh, isotopes, masses, 320)
    ui.show()
    sys.exit(app.exec())