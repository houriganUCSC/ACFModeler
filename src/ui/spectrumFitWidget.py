import sys
import pickle
import os

from PyQt6.QtWidgets import QApplication, QWidget, QGridLayout, QComboBox, QPushButton,\
    QFileDialog, QLabel, QDoubleSpinBox, QCheckBox, QHBoxLayout
from PyQt6.QtCore import QRect, Qt, pyqtSlot
from PyQt6.QtCore import pyqtSlot, pyqtSignal

from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas

# Project imports
import src.ui.spectrumFitTable as acfTable
import src.ui.spectrumModelDesignTable as modelDesign
import src.records.Session as Session
from src.fileIO.ThermoE2XR import ThermoFIN2


class SpectrumWidget(QWidget):
    spectrumFitSaved = pyqtSignal(str)
    def __init__(self, session: Session):
        super().__init__()
        self.session = session
        self.setFixedWidth(1190)
        self.setFixedHeight(800)
        self.setContentsMargins(10,5,10,5)


        axlayout = QGridLayout()
        axlayout.setVerticalSpacing(3)


        # Fixed width
        fw = 90
        combo_fw = 190
        spb_fw = 78

        # Top Left
        row = 0
        label = QLabel()
        label.setText("1. Configure Spectrum Fit Model")
        label.setContentsMargins(0, 10, 0, 0)
        label.setFixedWidth(300)
        label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        label.setStyleSheet("font-weight: bold")
        axlayout.addWidget(label, row, 0, 1, 3)

        row += 1
        self.modelTable = modelDesign.ModelDesignTable()
        self.modelTable.setFixedHeight(115)
        self.modelTable.setFixedWidth(300)
        self.modelTable.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        axlayout.addWidget(self.modelTable, row, 0, 1, 3)

        row += 1
        label = QLabel()
        label.setText("2. Explore Spectrum Fit PLots")
        label.setContentsMargins(0, 10, 0, 0)
        label.setFixedWidth(300)
        label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        label.setStyleSheet("font-weight: bold")
        axlayout.addWidget(label, row, 0, 1, 3)

        row += 1
        label = QLabel()
        label.setText("Plot:")
        label.setContentsMargins(0, 0, 0, 0)
        label.setFixedWidth(fw)
        axlayout.addWidget(label, row, 0, 1, 1)
        self.plotSelectorCombo = QComboBox()
        self.plotSelectorCombo.addItems(['All', 'a1', 'a2', 'tau', 'ACF', 'Drift'])
        self.plotSelectorCombo.setProperty('ID', 'plotSelector')
        self.plotSelectorCombo.setFixedWidth(combo_fw)
        axlayout.addWidget(self.plotSelectorCombo, row, 1, 1, 2)

        row +=1
        for i, t in enumerate(['', 'Min', 'Max']):
            label = QLabel()
            label.setText(t)
            label.setContentsMargins(0, 0, 0, 0)
            label.setFixedWidth(fw)
            label.setAlignment(Qt.AlignmentFlag.AlignCenter)
            axlayout.addWidget(label, row, i)

        row+=1
        label = QLabel()
        label.setText('X-axis limits:')
        label.setContentsMargins(0, 0, 0, 0)
        label.setFixedWidth(fw)
        self.xLimMin = QDoubleSpinBox()
        self.xLimMin.setRange(1, 250)
        self.xLimMin.setContentsMargins(0, 0, 0, 0)
        self.xLimMin.setFixedWidth(spb_fw)
        self.xLimMin.setKeyboardTracking(False)
        self.xLimMin.setProperty('ID', 'xLimMin')
        self.xLimMin.setValue(1)
        self.xLimMin.valueChanged.connect(self.redrawPlots)
        self.xLimMax = QDoubleSpinBox()
        self.xLimMax.setRange(1, 250)
        self.xLimMax.setContentsMargins(0, 0, 0, 0)
        self.xLimMax.setFixedWidth(spb_fw)
        self.xLimMax.setKeyboardTracking(False)
        self.xLimMax.setProperty('ID', 'xLimMax')
        self.xLimMax.setValue(240)
        self.xLimMax.valueChanged.connect(self.redrawPlots)
        axlayout.addWidget(label, row, 0)
        axlayout.addWidget(self.xLimMin, row, 1)
        axlayout.addWidget(self.xLimMax, row, 2)

        row += 1
        label = QLabel()
        label.setText('Y-axis limits:')
        label.setContentsMargins(0, 0, 0, 0)
        label.setFixedWidth(fw)
        self.yLimMin = QDoubleSpinBox()
        self.yLimMin.setRange(-1000, 1000)
        self.yLimMin.setContentsMargins(0, 0, 0, 0)
        self.yLimMin.setFixedWidth(spb_fw)
        self.yLimMin.setEnabled(self.plotSelectorCombo.currentIndex() > 0)
        self.yLimMin.setKeyboardTracking(False)
        self.yLimMin.setProperty('ID','yLimMin')
        self.yLimMin.valueChanged.connect(self.redrawPlots)
        self.yLimMax = QDoubleSpinBox()
        self.yLimMax.setRange(-1000, 1000)
        self.yLimMax.setContentsMargins(0, 0, 0, 0)
        self.yLimMax.setFixedWidth(spb_fw)
        self.yLimMax.setEnabled(self.plotSelectorCombo.currentIndex() > 0)
        self.yLimMax.setKeyboardTracking(False)
        self.yLimMax.setProperty('ID', 'yLimMax')
        self.yLimMax.valueChanged.connect(self.redrawPlots)
        axlayout.addWidget(label, row, 0)
        axlayout.addWidget(self.yLimMin, row, 1)
        axlayout.addWidget(self.yLimMax, row, 2)

        row += 1
        label = QLabel()
        label.setText('Annotation:')
        label.setContentsMargins(0, 0, 0, 0)
        label.setFixedWidth(fw)
        axlayout.addWidget(label, row, 0, 1, 1)

        # self.annoCheckBox.released.connect(self.redrawPlots)
        self.annoCombo = QComboBox()
        self.annoCombo.addItems(["None", "Upper left", "Upper center", "Upper right",
                                 "Middle left", "Middle right",
                                 "Lower left", "Lower center", "Lower right"])
        self.annoCombo.setContentsMargins(0,0,0,0)
        self.annoCombo.setFixedWidth(combo_fw)
        self.annoCombo.setCurrentIndex(0)
        axlayout.addWidget(self.annoCombo, row, 1, 1, 2)
        self.annoCombo.currentIndexChanged.connect(self.redrawPlots)

        row+=1
        label = QLabel()
        label.setText("Save Plot(s)")
        label.setContentsMargins(0, 0, 0, 0)
        label.setFixedWidth(fw)
        self.savePlotBtn = QPushButton("Save")
        self.savePlotBtn.setFixedWidth(2*fw)
        self.savePlotBtn.setContentsMargins(0, 0, 0, 0)
        axlayout.addWidget(label, row, 0, 1, 1)
        axlayout.addWidget(self.savePlotBtn, row, 1, 1, 2)

        row += 1
        label = QLabel()
        label.setText("3. Commit Spectrum Fits")
        label.setContentsMargins(0, 10, 0, 0)
        label.setFixedWidth(300)
        label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        label.setStyleSheet("font-weight: bold")
        axlayout.addWidget(label, row, 0, 1, 3)

        row += 1
        label = QLabel()
        label.setText("Save Files(s)")
        label.setContentsMargins(0, 0, 0, 0)
        label.setFixedWidth(fw)
        self.saveBinBtn = QPushButton("Save Binary")
        self.saveBinBtn.setFixedWidth(2*fw)
        self.saveBinBtn.setContentsMargins(0, 0, 0, 0)
        axlayout.addWidget(label, row, 0, 1, 1)
        axlayout.addWidget(self.saveBinBtn, row, 1, 1, 2)

        row+=1
        self.saveSumBtn = QPushButton("Save Summary")
        self.saveSumBtn.setFixedWidth(2*fw)
        self.saveSumBtn.setContentsMargins(0, 0, 0, 0)
        # Todo: Create report generator
        self.saveSumBtn.setEnabled(False)
        axlayout.addWidget(self.saveSumBtn, row, 1, 1, 2)

        row += 1
        label = QLabel()
        label.setText("4. Export Chromatogram Files")
        label.setContentsMargins(0, 10, 0, 0)
        label.setFixedWidth(300)
        label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        label.setStyleSheet("font-weight: bold")
        axlayout.addWidget(label, row, 0, 1, 3)

        row += 1
        hbox = QHBoxLayout()
        hbox.setSpacing(0)
        hbox.setContentsMargins(0, 0, 0, 0)
        label = QLabel()
        label.setText("Format:")
        label.setContentsMargins(0, 0, 0, 0)
        label.setStyleSheet("padding: 0 0 0 0px;; border: 0 0 0 0px")
        label.setFixedWidth(50)
        align = Qt.AlignmentFlag.AlignRight |Qt.AlignmentFlag.AlignVCenter
        label.setAlignment(align)
        hbox.addWidget(label)
        self.expTypeCombo = QComboBox()
        self.expTypeCombo.addItems(['*FIN2'])
        self.expTypeCombo.setFixedWidth(120)
        self.expTypeCombo.setContentsMargins(0, 0, 0, 0)
        self.expTypeCombo.setEnabled(False)
        hbox.addWidget(self.expTypeCombo)

        self.expStartBtn = QPushButton("Export")
        self.expStartBtn.setFixedWidth(80)
        self.expStartBtn.setContentsMargins(0, 0, 0, 0)
        self.expStartBtn.setEnabled(False)
        hbox.addWidget(self.expStartBtn)
        hbox.setContentsMargins(0, 0, 0, 0)

        axlayout.addLayout(hbox, row, 0, 1, 3)

        row+=1
        label = QLabel("Equations:")
        label.setContentsMargins(0, 0, 0, 0)
        label.setFixedWidth(fw)
        label.setAlignment(Qt.AlignmentFlag.AlignLeft)
        axlayout.addWidget(label, row, 0, 1, 1)

        row += 1
        self.fitEquations = QLabel()
        self.fitEquations.setFixedHeight(50)
        axlayout.addWidget(self.fitEquations, row, 0, 1, 3)

        # Top Right
        self.fig = Figure(figsize=(5,4))
        self.canvas = FigureCanvas(self.fig)
        self.fig.set_canvas(self.canvas)

        # Bottom
        self.spectrumTable = acfTable.CrossCalTableWidget(self)

        # self.updatePlotLims()
        self.spectrumTable.setContentsMargins(0, 0, 0, 0)

        '''Create and populate Widget layout grid'''

        self.grid = QGridLayout()
        self.grid.setContentsMargins(0, 0, 0, 0)

        # Add Left-hand side widgets

        # Add right-hand side widget
        self.grid.addLayout(axlayout, 0, 0, 1, 1)
        self.grid.addWidget(self.canvas, 0, 1, 1, 1)

        # Add bottom widget
        self.grid.addWidget(self.spectrumTable, 1, 0, 1, 2)
        self.setLayout(self.grid)

        self.connectEvents()

    def connectEvents(self):
        self.modelTable.tableModified.connect(self.spectrumTable.createPlots)
        self.plotSelectorCombo.currentIndexChanged.connect(self.redrawPlots)
        self.spectrumTable.plotsDrawn.connect(self.updatePlotLims)
        self.spectrumTable.dataPostProcessed.connect(self.enableExport)
        # self.saveBinBtn.released.connect(self.pickleData)
        self.saveBinBtn.released.connect(self.pickleData)
        self.savePlotBtn.released.connect(self.savePlot)
        self.expStartBtn.released.connect(self.exportChromText)


    @pyqtSlot()
    def recalculateSpectrumFit(self):
        self.spectrumTable.createPlots()

    def updatePlotLims(self):
        xlim = self.spectrumTable.axs[0].get_xlim()
        ylim = self.spectrumTable.axs[0].get_ylim()
        self.yLimMin.setValue(ylim[0])
        self.yLimMax.setValue(ylim[1])
        self.xLimMin.setValue(xlim[0])
        self.xLimMax.setValue(xlim[1])
        self.yLimMin.setEnabled(self.plotSelectorCombo.currentIndex() > 0)
        self.yLimMax.setEnabled(self.plotSelectorCombo.currentIndex() > 0)

    def redrawPlots(self):
        sendy = self.sender()
        sendyID = sendy.property('ID')
        kwargs = {}
        kwargs['xlim'] = (self.xLimMin.value(), self.xLimMax.value())
        kwargs['anno'] = self.annoCombo.currentIndex()
        if not sendyID == "plotSelector":
            if self.plotSelectorCombo.currentIndex() > 0:
                kwargs['ylim'] = (self.yLimMin.value(), self.yLimMax.value())
        self.spectrumTable.createPlots(**kwargs)

    def pickleData(self):
        self.spectrumTable.commitFits()
        self.startDir = None
        dlg = QFileDialog()
        dlg.setFileMode(QFileDialog.FileMode.AnyFile)
        pickleFile, filter = dlg.getSaveFileName(
            None,
            "Select Data Files",
            self.session.picklePath,
            "Python Pickle File (*.p)",
            options=QFileDialog.Option.DontUseNativeDialog
        )
        if pickleFile:
            with open(pickleFile, 'wb') as fid:
                pickle.dump(self.session, fid)
            fid.close()
            self.spectrumFitSaved.emit(pickleFile)
    def enableExport(self):
        self.expTypeCombo.setEnabled(True)
        self.expStartBtn.setEnabled(True)

    def savePlot(self):
        plotType = self.plotSelectorCombo.currentText()
        plotName = f'spectrum_fit_{plotType}'
        startingPath = os.path.join(self.session.startDir, plotName)
        fileName, _ = QFileDialog.getSaveFileName(self,
                                                  "Save Figure",
                                                  startingPath,
                                                  "PDF Files(*.pdf);; PNG Files(*.png);; JPG Files(*.jpg)",
                                                  options=QFileDialog.Option.DontUseNativeDialog)
        if fileName:
            self.fig.savefig(fileName, dpi=300)

    def exportChromText(self):
        form = self.expTypeCombo.currentText()
        if form == "*FIN2":
            exporter = ThermoFIN2(self.session)
            finDir = self.session.startDir.joinpath('FIN2')
            finDir.mkdir(exist_ok=True)
            dlg = QFileDialog()
            dlqOptions = QFileDialog.Option.DontUseNativeDialog | QFileDialog.Option.DontUseNativeDialog
            fin2Dir = dlg.getExistingDirectory(None, "", str(finDir), options=dlqOptions)
            finName = os.path.basename(self.session.startDir)+'.FIN'
            exporter.write(fin2Dir, finName, True)



if __name__ == "__main__":
    app = QApplication(sys.argv)
    session = Session.Session()

    # pPath ="/Users/jeremyhourigan/PycharmProjects/ACFModeler/test/SESSION.p"
    pPath = "/Users/jeremyhourigan/PycharmProjects/ACFModeler/test/SESSION_SEQ1_SPECTRUM.p"
    # pPath = "/Users/jeremyhourigan/PycharmProjects/ACFModeler/test/pennState30um40pct_FIT.p"
    # pPath = "/Users/jeremyhourigan/My Drive/LA-ICP-MS/External_User_Data/PennState/pennState30um40pct.p"
    with open(pPath, 'rb') as fid:
        session = pickle.load(fid)
    fid.close()
    demo = SpectrumWidget(session)
    demo.show()
    app.exit(app.exec())
