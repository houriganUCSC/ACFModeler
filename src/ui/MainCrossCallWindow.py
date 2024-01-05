import sys
import numpy
from PyQt6.QtWidgets import QApplication, QMainWindow, QTabWidget, QGridLayout, QComboBox,\
    QPushButton, QLabel, QSizePolicy, QTableWidget, QWidget

from PyQt6.QtCore import pyqtSlot

import src.ui.importWidget as importer
import src.ui.isotopeFitWidget as fitter
import src.ui.spectrumFitWidget as spectrum
import src.records.Session as session

from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas


class ACFModelMain(QMainWindow):

    def __init__(self):
        super().__init__()
        self.setFixedWidth(1220)
        self.setFixedHeight(820)
        self.tabs = QTabWidget()
        self.tabs.setFixedWidth(1210)
        self.tabs.setFixedHeight(810)
        self.tabs.setContentsMargins(5, 5, 5, 5)

        self.session = session.Session()
        self.session.pickleFile = None

        self.importer = importer.Ui_ImportWidget(self.session)
        self.importer.setupUi()
        # self.importer.samplesImported.connect(self.dataImported)
        self.importer.dataPickled.connect(self.dataImported)
        self.tabs.addTab(self.importer, "Import Raw Data")

        self.tabs.addTab(QWidget(), "Filter & Fit Raw Data")
        self.tabs.setTabEnabled(1, False)

        self.tabs.addTab(QWidget(), "Model Mass Spectrum")
        self.tabs.setTabEnabled(2, False)

        self.setCentralWidget(self.tabs)

        self.tabs.currentChanged.connect(self.tabChanged)

    def tabChanged(self, current: int):
        print(f'{__name__} tab changed: {self.session}')
        if current == 1 and self.session.status['new']:
            self.rawFitWidget.dataImported()
            self.session.status['new'] = False
            self.session.status['imported'] = True

    @pyqtSlot(str)
    def dataImported(self, pickleFilePath: str):
        print(__name__, "Data Imported")
        self.session = self.importer.session
        self.tabs.removeTab(1)
        self.rawFitWidget = fitter.FilterFitWidget(self.session)
        self.rawFitWidget.setupUi()
        self.rawFitWidget.regressionPickled.connect(self.dataFit)
        self.tabs.insertTab(1, self.rawFitWidget,"Filter & Fit Raw Data")
        self.session.pickleFile = pickleFilePath
        print("Data Imported", pickleFilePath)

    @pyqtSlot(str)
    def dataFit(self, pickleFilePath:str):
        self.massSpectrumWidget = spectrum.SpectrumWidget(self.session)
        self.tabs.removeTab(2)
        self.tabs.insertTab(2, self.massSpectrumWidget, "Model Mass Spectrum")
        self.session.pickleFile = pickleFilePath
        self.massSpectrumWidget.spectrumFitSaved.connect(self.spectrumFitPickled)
        print("Data Regressed", pickleFilePath)


    @pyqtSlot(str)
    def spectrumFitPickled(self, pickleFilePath:str):
        self.session.pickleFile = pickleFilePath
        print("Spectrum Regression Saved", pickleFilePath)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    demo = ACFModelMain()
    demo.show()
    app.exit(app.exec())