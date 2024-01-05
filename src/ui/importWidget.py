import os
import time
import numpy as np
import pathlib
from src.records.Session import Session, Mass
from src.fileIO.ThermoE2XR import ThermoDAT
import pickle


from PyQt6.QtWidgets import QApplication, QWidget, QLabel, QComboBox, QPushButton, \
    QSizePolicy, QTreeWidget, QFileDialog, QTreeWidgetItem, QFrame, QProgressBar, QMessageBox, \
    QGridLayout, QInputDialog
from PyQt6.QtCore import QRunnable, QObject, QThreadPool, pyqtSignal as Signal, pyqtSlot as Slot, QMutex, \
    QRect,  QCoreApplication, QMetaObject, Qt

import matplotlib

from datetime import datetime
# from src.fileIO.ThermoElement import *

# Matplotlib imports
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas

matplotlib.use("QTAgg")

""" 
IMPORT WORKER CLASSES 
Provide functionality to multi-thread the import process.  
QRunnable class allows for multi-thread and emission of signals to the main application instance
Connections to these signals allows real-time updating of progress and trees. 

see: https://www.pythontutorial.net/pyqt/qthreadpool/
"""


class ImporterSignals(QObject):
    progressMax = Signal(int)
    progressInc = Signal()
    progressMsg = Signal(str)
    addToImport = Signal(dict)
    removeFromQueue = Signal()
    importStatus = Signal(str)
    queueStatus = Signal(str)
    completed = Signal()


class ImportWorker(QRunnable):
    def __init__(self, parent, qTree: QTreeWidget,
                 queueDict: dict,
                 session: Session,
                 startTime: dict,
                 mutex: QMutex):
        super().__init__()
        self.parent = parent
        self.mutex = mutex
        self.qTree = qTree
        # self.impTree = impTree
        self.nQ = 0
        self.nImp = 0
        self.queueDict = queueDict
        self.session = session
        # self.dataDict = dataDictionary
        self.startTime = startTime
        self.treeData = {'name': '',
                         'fileName': '',
                         'fullFileName': '',
                         'dirName': '',
                         'takeParent': -1,
                         'takeChild': -1,
                         'removeDir': False}
        self.signals = ImporterSignals()


    def treeCounts(self):
        temp = 0
        for i in range(self.qTree.invisibleRootItem().childCount()):
            temp += self.qTree.invisibleRootItem().child(i).childCount()
        self.nQ = temp

    @Slot()
    def run(self):
        m = 0
        self.treeCounts()
        self.signals.progressMax.emit(self.nQ)
        # Iterate Sequence Directories
        nSeq = self.qTree.invisibleRootItem().childCount()
        self.mutex.tryLock(1000)
        self.session.startTime = time.localtime()
        for i in range(nSeq):
            self.treeData['removeDir'] = False
            dirItem = self.qTree.invisibleRootItem().child(m)
            dirName = dirItem.text(0)
            # Check for checkbox status
            if dirItem.checkState(0) == Qt.CheckState.Checked or Qt.CheckState.PartiallyChecked:
                # Import children of checked or partially
                # Iterate through analysis in Sequence
                k = 0
                for j in range(dirItem.childCount()):
                    fileItem = dirItem.child(k)
                    fileName = fileItem.text(0)
                    name = fileName.strip('.dat')
                    self.signals.progressMsg.emit(f'Importing {fileName}')
                    if fileItem.checkState(0) == Qt.CheckState.Checked:
                        fpath = self.queueDict[dirName][fileName]["path"]
                        self.treeData['fullFileName'] = fpath
                        if '.dat' in fpath:
                            ThermoDAT(fpath, self.session).parseDAT()
                        else:
                            print("No parser for this file type")
                            break
                        self.session.status['imported'] = True
                        self.treeData["name"] = self.session.samples[name].group
                        self.treeData["fileName"] = self.session.samples[name].name
                        self.treeData['sampleNum'] = self.session.samples[name].smpNum
                        self.treeData["dirName"] = dirName
                        self.treeData["fullFileName"] = fileName
                        self.treeData['takeParent'] = m
                        self.treeData['takeChild'] = k
                        self.treeData['mutex'] = self.mutex
                        self.signals.addToImport.emit(self.treeData)
                        self.signals.progressInc.emit()
                        # time.sleep(0.002)
                        while not self.mutex.tryLock(10):
                            print("locked")
                    else:
                        k += 1
                    if dirItem.childCount() == 0:
                        self.treeData['removeDir'] = True
                        self.signals.addToImport.emit(self.treeData)
                        while not self.mutex.tryLock(10):
                            print("locked")
            else:
                pass
                # m += 1
        self.session.status['imported'] = True
        self.session.timeOffsets()
        self.session.startTime = np.min(self.session.scanTime)
        for massName, massRecord in self.session.masses.items():
            massRecord.count()
        self.session.updateDeadTime()
        self.session.status['imported'] = True
        self.session.status['new'] = True
        self.signals.completed.emit()


class ExportWorker(QRunnable):
    def __init__(self, session: Session, mutex: QMutex):
        super().__init__()
        self.session = session
        self.mutex = QMutex
        self.signals = ImporterSignals()

    @Slot()
    def run(self):
        m = 0
        smpNames = list(self.session.samples.keys())
        self.signals.progressMax.emit(len(smpNames))
        # Iterate Sequence Directories
        self.mutex.tryLock(1000)
        for smpName in smpNames:
            primaryKey = self.session.samples[smpName].ID

        self.session.startTime = time.localtime()
        for i in range(nSeq):
            self.treeData['removeDir'] = False
            dirItem = self.qTree.invisibleRootItem().child(m)
            dirName = dirItem.text(0)
            # Check for checkbox status
            if dirItem.checkState(0) == Qt.Checked or Qt.PartiallyChecked:
                # Import children of checked or partially
                # Iterate through analysis in Sequence
                k = 0
                for j in range(dirItem.childCount()):
                    fileItem = dirItem.child(k)
                    fileName = fileItem.text(0)
                    name = fileName.strip('.dat')
                    self.signals.progressMsg.emit(f'Importing {fileName}')
                    if fileItem.checkState(0) == Qt.CheckState.Checked:
                        fpath = self.queueDict[dirName][fileName]["path"]
                        self.treeData['fullFileName'] = fpath
                        # self.smp = self.session.createSample(name)
                        # datReader.parseDAT(fpath, self.session, self.mutex)
                        self.session.status['imported'] = True
                        # print(type(s.fileTimes["DAT"])
                        self.treeData["name"] = self.session.samples[name].group
                        self.treeData["fileName"] = self.session.samples[name].name
                        self.treeData['sampleNum'] = self.session.samples[name].smpNum
                        self.treeData["dirName"] = dirName
                        self.treeData["fullFileName"] = fileName
                        self.treeData['takeParent'] = m
                        self.treeData['takeChild'] = k
                        self.treeData['mutex'] = self.mutex
                        print(self.treeData)
                        self.signals.addToImport.emit(self.treeData)
                        self.signals.progressInc.emit()
                        while not self.mutex.tryLock(10):
                            print("locked")
                    else:
                        k += 1
                    if dirItem.childCount() == 0:
                        self.treeData['removeDir'] = True
                        self.signals.addToImport.emit(self.treeData)
                        while not self.mutex.tryLock(10):
                            print("locked")
            else:
                pass
                # m += 1
        self.session.status['imported'] = True
        self.session.timeOffsets()
        self.session.startTime = np.min(self.session.scanTime)
        self.session.calculateTimeSeries()
        for massName, massRecord in self.session.masses.items():
            massRecord.totalObservations()
        self.signals.completed.emit()


class Ui_ImportWidget(QWidget):
    samplesImported = Signal()
    dataPickled = Signal(str)
    def __init__(self, sessionData: Session):
        super().__init__()
        self.session = sessionData
        print("import init", self.session)
        self.queueDict = {}
        self.startTime = {}
        self.startDir = r'G:\My Drive\LA-ICP-MS\External_User_Data\PennState'
        # self.startDir = r'G:\My Drive\LA-ICP-MS\External_User_Data\Ingersoll\RVI22-01\Data\RVI_22_01_SEQ1'
        # self.startDir = '/Volumes/GoogleDrive/My Drive/LA-ICP-MS/External_User_Data/Ingersoll/RVI22-01/Data/RVI_22_01_SEQ1'
        self.queueGroups = {}
        self.importGroups = {}
        self.uiMutex = QMutex()
        self.session.status['new'] = False
        self.session.status['imported'] = False
        # print(f'import: {self.session}')

    def setupUi(self):
        sizePolicy = QSizePolicy(QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Minimum)
        self.setObjectName("Widget")
        self.setFixedWidth(1200)
        self.setFixedHeight(800)
        self.setSizePolicy(sizePolicy)

        self.grid = QGridLayout()
        self.grid.setGeometry(QRect(0, 0, 1190, 790))
        self.grid.setContentsMargins(5, 5, 5, 200)
        for i, width in enumerate([200, 200, 690]):
            self.grid.setColumnMinimumWidth(i, width)
            self.grid.setColumnStretch(i, 0)

        self.setLayout(self.grid)

        """ QUEUE QComboBox: File Extension Filter Combobox"""
        self.queueFilter = QComboBox()
        self.queueFilter.setObjectName("queueFilter")
        self.queueFilter.addItem("")
        self.queueFilter.addItem("")
        self.grid.addWidget(self.queueFilter, 0, 0, 1, 1)

        """ QUEUE Status: QLabel to store status string"""
        self.queueStatus = QLabel()
        self.queueStatus.setObjectName("queueStatus")
        # Todo header to queue QTreeWidget?

        """QUEUE TREE:  Tree Widget to show directories and files"""
        self.queueTreeWidget = QTreeWidget()
        self.queueTreeWidget.setObjectName("queueTreeWidget")
        self.queueTreeWidget.headerItem().setText(0, "1")
        self.queueTreeWidget.setGeometry(QRect(0, 0, 190, 350))
        self.queueTreeWidget.setMaximumSize(190, 350)
        self.queueTreeWidget.setMinimumSize(190, 350)
        self.queueTreeWidget.setSizePolicy(sizePolicy)
        self.grid.addWidget(self.queueTreeWidget, 1, 0, 1, 1)

        """QUEUE ADD FILES:  QPushButton to add new files"""
        self.queueAddFiles = QPushButton()
        self.queueAddFiles.setObjectName("queueAddFiles")
        self.grid.addWidget(self.queueAddFiles, 2, 0, 1, 1)

        """QUEUE ADD DIR:  QPushButton to add new directories"""
        self.queueAddDir = QPushButton()
        self.queueAddDir.setObjectName("queueAddDir")
        self.grid.addWidget(self.queueAddDir, 3, 0, 1, 1)

        """QUEUE REMOVE:  QPushButton to remove selected files"""
        self.queueRemove = QPushButton()
        self.queueRemove.setObjectName("queueRemove")
        self.grid.addWidget(self.queueRemove, 4, 0, 1, 1)

        """ QUEUE SPACER:  Spacer to make things look better"""
        # spacerItem = QSpacerItem(140, 42, QSizePolicy.Minimum, QSizePolicy.Fixed)
        # self.grid.addItem(spacerItem, 5, 0, 1, 1)

        """ IMPORTED:  QComboBox for Chrom Text File Format"""
        self.rawExportFileTypeCombo = QComboBox()
        self.rawExportFileTypeCombo.setObjectName("rawExportFileTypeCombo")
        self.rawExportFileTypeCombo.addItem("")
        self.rawExportFileTypeCombo.addItem("")
        self.grid.addWidget(self.rawExportFileTypeCombo, 0, 1, 1, 1)

        """ IMPORTED Status: QLabel to store status string"""
        self.importedStatus = QLabel()
        self.importedStatus.setObjectName("importedStatus")
        # TODO HeaderItem to ImportTreeWidget?

        """IMPORTED TREE:  Tree Widget to show directories and files"""
        self.importedTreeWidget = QTreeWidget(self)
        self.importedTreeWidget.setObjectName("importedTreeWidget")
        self.importedTreeWidget.setGeometry(QRect(0, 0, 190, 350))
        self.importedTreeWidget.setMaximumSize(190, 350)
        self.importedTreeWidget.setMinimumSize(190, 350)
        self.importedTreeWidget.setSizePolicy(sizePolicy)
        self.importedTreeWidget.headerItem().setText(0, "1")
        self.importedTreeWidget.setSortingEnabled(False)
        # self.importedTreeWidget.sortByColumn(1, Qt.SortOrder.AscendingOrder)
        self.grid.addWidget(self.importedTreeWidget, 1, 1, 1, 1)
        self.importedTreeWidget.sortByColumn(0, Qt.SortOrder.AscendingOrder)
        self.importedTreeWidget.setSortingEnabled(True)

        """IMPORTED IMPORT QUEUED:  QPushButton to import selected files"""
        self.importQueuedBtn = QPushButton()
        self.importQueuedBtn.setObjectName("importQueuedBtn")
        self.grid.addWidget(self.importQueuedBtn, 2, 1, 1, 1)

        """IMPORTED IMPORT BINARY:  QPushButton to import pickled files"""
        self.importBinary = QPushButton()
        self.importBinary.setObjectName("importBinary")
        self.grid.addWidget(self.importBinary, 3, 1, 1, 1)

        """IMPORTED SAVE:  QPushButton to save (Pickle) selected files"""
        self.saveImportedBtn = QPushButton(self)
        self.saveImportedBtn.setObjectName("saveImportedBtn")
        self.grid.addWidget(self.saveImportedBtn, 4, 1, 1, 1)

        """IMPORTED EXPORT:  QPushButton to save as Chrom Text according to FileTypeComboBox"""
        self.rawExportChromBtn = QPushButton()
        self.rawExportChromBtn.setObjectName("rawExportChromBtn")
        self.grid.addWidget(self.rawExportChromBtn, 5, 1, 1, 1)

        """TIME-SERIES LABEL:  QLabel for Time-series"""
        self.timeSeriesLabel = QLabel()
        self.timeSeriesLabel.setObjectName("timeSeriesLabel")
        self.timeSeriesLabel.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.timeSeriesLabel.setStyleSheet('font-weight: bold')
        self.grid.addWidget(self.timeSeriesLabel, 0, 2, 1, 1)

        """TIME-SERIES PLOT:  Configure MatPlotLib Figure Canvas"""
        # Create Matlplotlib figure
        self.timesSeriesFigure = Figure(figsize=(5, 3))
        # FigureCanvas is the Matplotlib Widget for Pyplots
        self.canvas = FigureCanvas(self.timesSeriesFigure)
        self.canvas.setObjectName("timeSeriesCanvas")
        # Add FigureCanvas to grid
        self.grid.addWidget(self.canvas, 1, 2, 5, 1)
        # Create Axes for rendering
        self.dynamic_ax = self.canvas.figure.subplots()

        """PROGRESS Text:  QLabel to display Import process status messages"""
        self.progressText = QLabel(self)
        self.progressText.setObjectName("progressText")
        self.grid.addWidget(self.progressText, 6, 0, 1, 2)

        """PROGRESS Bar:  QProgress to display Import process progress"""
        self.progressBar = QProgressBar(self)
        self.progressBar.setObjectName("progressBar")
        self.progressBar.setProperty("value", 24)
        self.grid.addWidget(self.progressBar, 6, 2, 1, 1)

        """ LOCALIZE:  Possible local language capabilities"""
        self.retranslateUi()

        """WIDGET EVENTS:  ESTABLISHED CONNECTIONS FOR UI EVENTS"""
        self.queueAddFiles.clicked.connect(self.addFiles)
        self.queueAddDir.clicked.connect(self.addDir)
        self.queueRemove.clicked.connect(self.removeFiles)
        self.importQueuedBtn.clicked.connect(self.importQueued)
        self.rawExportChromBtn.clicked.connect(self.generateChromText)
        self.saveImportedBtn.clicked.connect(self.pickleImported)
        self.importBinary.clicked.connect(self.unpickleImported)
        self.importedTreeWidget.itemDoubleClicked.connect(self.plotSelectedTimeSeries)
        self.graphSelectionModel = self.importedTreeWidget.selectionModel()
        self.graphSelectionModel.selectionChanged.connect(self.plotSelectedTimeSeries)
        QMetaObject.connectSlotsByName(self)

        # Intialize Progress
        self.progressBar.setRange(0, 1)
        self.progressBar.setValue(0)
        self.progressText.setText("Add directories or files to queue")
        self.importedTreeWidget.setHeaderLabels(["Sample"])
        self.queueTreeWidget.setHeaderLabels(["Directory"])
        self.importedStatus.setText("")
        self.queueStatus.setText("")
        self.inited = True

    def setTreeColumnWidth(self):
        width = max(self.queueFilter.width(), self.rawExportFileTypeCombo.width())
        objects = [self.queueFilter, self.queueTreeWidget, self.queueAddDir, self.queueAddFiles, self.queueRemove,
                   self.rawExportFileTypeCombo, self.importedTreeWidget, self.importQueuedBtn, self.importBinary,
                   self.saveImportedBtn, self.rawExportChromBtn]
        for obj in objects:
            obj.setMaximumWidth(width)

    def retranslateUi(self):
        _translate = QCoreApplication.translate
        self.setWindowTitle(_translate("Widget", "Widget"))
        self.queueFilter.setItemText(0, _translate("Widget", "Thermo Element (*.dat)"))
        self.queueFilter.setItemText(1, _translate("Widget", "Perkin-Elmer"))
        self.queueStatus.setText(_translate("Widget", "TextLabel"))
        self.queueAddFiles.setText(_translate("Widget", "Add Files"))
        self.queueAddDir.setText(_translate("Widget", "Add Directory"))
        self.queueRemove.setText(_translate("Widget", "Remove"))
        self.rawExportFileTypeCombo.setItemText(0, _translate("Widget", "Thermo Element (*.FIN2)"))
        self.rawExportFileTypeCombo.setItemText(1, _translate("Widget", "Perkin-Elmer"))
        self.importedStatus.setText(_translate("Widget", "TextLabel"))
        self.importQueuedBtn.setText(_translate("Widget", "Import Queued Files"))
        self.importBinary.setText(_translate("Widget", "Import Binary"))
        self.saveImportedBtn.setText(_translate("Widget", "Save Session"))
        self.rawExportChromBtn.setText(_translate("Widget", "Export Chromatogram"))
        self.timeSeriesLabel.setText(_translate("Widget", "Time Series Graph"))
        self.progressText.setText(_translate("Widget", "Text"))
        self.setTreeColumnWidth()

    def importQueued(self):
        self.initializeProgress()
        pool = QThreadPool.globalInstance()
        importWorker = ImportWorker(self,
            self.queueTreeWidget,
            self.queueDict,
            self.session,
            self.startTime,
            self.uiMutex
        )

        importWorker.signals.progressMsg.connect(self.progressText.setText)
        importWorker.signals.progressInc.connect(self.incrementProgress)
        importWorker.signals.progressMax.connect(self.progressRange)
        importWorker.signals.addToImport.connect(self.transferToImportTree)
        importWorker.signals.queueStatus.connect(self.queueStatus.setText)
        importWorker.signals.importStatus.connect(self.importedStatus.setText)
        importWorker.signals.completed.connect(self.setProgressComplete)
        pool.start(importWorker)

    def generateChromText(self):
        expType = self.rawExportFileTypeCombo.currentText()
        if expType == 'Thermo Element (*.FIN2)':
            writeFIN2(self.dataDir)
        else:
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Critical)
            msg.setText("Error")
            msg.setInformativeText('No Exporter for this file type')
            msg.setWindowTitle("Error")
            msg.exec_()

    def removeFiles(self):
        pass

    def addFiles(self):
        self.getFilterText()
        dlg = QFileDialog()
        dlg.setFileMode(QFileDialog.FileMode.AnyFile)
        dlg.setDirectory(self.startDir)
        dataFiles, filter = dlg.getOpenFileNames(None, "Select Data Files", self.startDir, self.fileFilter)
        new = 0
        # print(dataFiles)
        for file in dataFiles:
            head, fileName = os.path.split(str(file))
            head, dirName = os.path.split(head)
            if dirName not in self.queueDict:
                self.queueDict[dirName] = {}
                self.startTime[dirName] = time.localtime(datetime.now().timestamp())
            if fileName.endswith(self.extension):
                self.queueDict[dirName][fileName] = {
                    "imported": False,
                    "path": (str(file))
                }
            new += 1
        if new:
            self.tree_from_dict(self.queueDict, self.queueTreeWidget)

        self.queueStatus.setText(f'In Queue: {self.getQueueCount(self.queueTreeWidget)} files')
        self.importQueuedBtn.setEnabled(self.getQueueCount(self.queueTreeWidget))
        self.queueRemove.setEnabled(self.getQueueCount(self.queueTreeWidget))

    def addDir(self):
        self.getFilterText()
        dlg = QFileDialog()
        dirPath = dlg.getExistingDirectory(None, 'Select Data Directory', self.startDir,
                                           QFileDialog.Option.ShowDirsOnly)
        head, dirName = os.path.split(dirPath)
        if dirName not in self.queueDict:
            self.queueDict[dirName] = {}
            self.startTime[dirName] = time.localtime(datetime.now().timestamp())
        new = 0
        for file in os.listdir(dirPath):
            if file.endswith(self.extension):
                self.queueDict[dirName][file] = {
                    "imported": False,
                    "path": os.path.join(dirPath, file)
                }
                new += 1
        if new:
            self.tree_from_dict(self.queueDict, self.queueTreeWidget)
        else:
            del self.queueDict[dirName]
            dlg = QMessageBox()
            dlg.setText(f"No {self.extension} files found")
            dlg.exec()

        paths = []
        for dirName, dirEntries in self.queueDict.items():
            firstFile = list(dirEntries.keys())[0]
            paths.append(dirEntries[firstFile]['path'])
        if len(paths) == 1:
            self.commonPath = pathlib.Path(paths[0]).parent
        else:
            self.commonPath = pathlib.Path(os.path.commonpath(paths))
        self.queueStatus.setText(f'In Queue: {self.getQueueCount(self.queueTreeWidget)} files')
        self.importQueuedBtn.setEnabled(self.getQueueCount(self.queueTreeWidget))
        self.queueRemove.setEnabled(self.getQueueCount(self.queueTreeWidget))

    def getQueueCount(self, tree):
        inTree = 0
        for i in range(tree.invisibleRootItem().childCount()):
            inTree += tree.invisibleRootItem().child(i).childCount()
        return inTree

    """ PROGRESS UPDATE FUNCTIONS"""

    def initializeProgress(self):
        self.progressText.setText("Initializing import....")
        self.progressBar.setRange(0, 1)
        self.progressBar.setValue(0)
        self.progressBar.setTextVisible(False)

    def progressRange(self, max):
        self.progressBar.setRange(0, max)
        self.progressBar.setFormat("%p%")
        self.progressBar.setTextVisible(True)

    def incrementProgress(self):
        self.progressBar.setValue(self.progressBar.value() + 1)

    def setProgressComplete(self):
        self.progressText.setText("Import of queued files complete")
        self.progressBar.setRange(0, 1)
        self.progressBar.setValue(1)
        self.importedTreeWidget.sortByColumn(1, Qt.SortOrder.AscendingOrder)
        self.importedTreeWidget.setSortingEnabled(True)
        if self.session.method is None:
            # TODO:  fix missing setup information
            from src.ui.importRunTable import IsotopeIDQuery
            methodGetter = IsotopeIDQuery(self.session)
            methodGetter.exec()
            oldMassKeys = list(self.session.masses.keys())
            newMassKeys = self.session.method["isotopes"]
            if oldMassKeys != newMassKeys:
                for old, new in zip(oldMassKeys, newMassKeys):
                    self.session.masses[new] = self.session.masses.pop(old)
                self.session.isotopes = self.session.method["isotopes"]
            self.session.machineDeadTime = self.session.method["deadTime"]

        else:
            self.session.updateDeadTime()
            self.updateSessionIsotopes()
        self.samplesImported.emit()

    def updateSessionIsotopes(self):
        for smpName, smpRecord in self.session.samples.items():
            self.session.isotopes = np.union1d(self.session.isotopes, smpRecord.isotopes)

    """ TREE WIDGET FUNCTIONS """

    def tree_from_dict(self, data=None, parent=None):
        root = parent.invisibleRootItem()
        for folder, files in data.items():
            if not folder in self.queueGroups:
                item = QTreeWidgetItem(root)
                item.setText(0, folder)
                item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserTristate | Qt.ItemFlag.ItemIsUserCheckable)
                item.setCheckState(0, Qt.CheckState.Checked)
                self.queueGroups[folder] = root.indexOfChild(item)
                for file in files:
                    child = QTreeWidgetItem(item)
                    child.setFlags(child.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                    child.setText(0, file)
                    child.setCheckState(0, Qt.CheckState.Checked)

    def transferToImportTree(self, smp: dict):
        sampleName = smp['name']
        fileName = smp['fileName']
        fullFileName = smp['fullFileName']
        dirName = smp['dirName']
        smpNum = smp['sampleNum']
        mutex = smp['mutex']
        queued = self.queueTreeWidget.invisibleRootItem()
        if not smp['removeDir']:
            imported = self.importedTreeWidget.invisibleRootItem()

            # If sample not in tree
            if sampleName not in list(self.importGroups.keys()):
                #create tree item
                QTreeWidgetItem(imported).setText(0, sampleName)
                # update sample item indices
                for idx in range(imported.childCount()):
                    key = imported.child(idx).text(0)
                    self.importGroups[key] = idx
            sample = imported.child(self.importGroups[sampleName])
            fileItem = QTreeWidgetItem(sample)
            fileItem.setText(0, fileName)
            fileItem.setText(1, str(smpNum))
            fileItem.setText(2, fileName)
            #self.importedTreeWidget.setColumnHidden(1, False)
            #self.importedTreeWidget.setColumnHidden(2, False)
            queued.child(smp['takeParent']).takeChild(smp['takeChild'])
        else:
            queued.takeChild(smp['takeParent'])
        msg = f'Queued: {self.getQueueCount(self.queueTreeWidget)} files'
        self.queueStatus.setText(msg)
        msg = f'Imported: {self.getQueueCount(self.importedTreeWidget)} files'
        self.importedStatus.setText(msg)
        mutex.unlock()

    """ PLOTTING FUNCTIONS"""

    def plotSelectedTimeSeries(self, selected, deselected):
        item = self.importedTreeWidget.selectedItems()[0]
        file = item.text(2)
        print(file)
        self.session.plotTimeSeries(file, self.dynamic_ax)

    def pickleImported(self):
        self.session.commonPath = self.commonPath
        startingName = os.path.basename(self.commonPath)
        self.session.startDir = self.commonPath.joinpath('ACF_Model')
        self.session.startDir.mkdir(exist_ok=True)
        startingPath = self.session.startDir
        startingPath = startingPath.joinpath(startingName).with_suffix('.p')
        print(startingPath)
        dlg = QFileDialog()
        dlg.setFileMode(QFileDialog.FileMode.AnyFile)
        dlg.setDirectory(str(self.session.startDir))
        pickleFile, filter = dlg.getSaveFileName(
            None,
            "Select Data Files",
            str(startingPath),
            "Python Pickle File (*.p)",
            options=QFileDialog.Option.DontUseNativeDialog)
        if pickleFile:
            with open(pickleFile, 'wb') as fid:
                pickle.dump(self.session, fid)
            fid.close()
            self.session.picklePath = str(pickleFile)
            self.dataPickled.emit(pickleFile)

    def unpickleImported(self):
        print(f'{__name__}: unpickle {self.session}')
        dlg = QFileDialog()
        dlg.setFileMode(QFileDialog.FileMode.AnyFile)
        pickleFile, filter = dlg.getOpenFileName(
            None,
            "Select Data Files",
            self.startDir,
            "Python Pickle File (*.p)")
        if pickleFile:
            with open(pickleFile, 'rb') as fid:
                self.session = pickle.load(fid)
            fid.close()
            binGroups = {}
            self.importedTreeWidget.clear()
            print("binary read", self.session, self.session.isotopes)
            for smpName, smpRecord in self.session.samples.items():
                group = smpRecord.group
                name = smpRecord.name
                num = smpRecord.smpNum
                fileName = name
                dirName = smpRecord.seqDir
                if not group in binGroups.keys():
                    binGroups[group] = {}
                binGroups[group][num] = name
            imported = self.importedTreeWidget.invisibleRootItem()
            for group, files in binGroups.items():
                sample = QTreeWidgetItem(imported)
                sample.setText(0, group)
                for num in sorted(files.keys()):
                    fileItem = QTreeWidgetItem(sample)
                    fileItem.setText(0, binGroups[group][num])
                    fileItem.setText(1, str(num))
                    fileItem.setText(2, binGroups[group][num])
                # fileItem.setText(2, fullFileName)
                self.importedTreeWidget.setColumnHidden(1, True)
                self.importedTreeWidget.setColumnHidden(2, False)
            self.rawExportChromBtn.setEnabled(self.getQueueCount(self.importedTreeWidget))
            self.saveImportedBtn.setEnabled(self.getQueueCount(self.importedTreeWidget))
            self.importedTreeWidget.sortByColumn(1, Qt.SortOrder.AscendingOrder)
            self.importedTreeWidget.setSortingEnabled(True)
            self.session.status["new"] = True
            self.session.status["imported"] = True
            self.dataPickled.emit(pickleFile)
            # Todo: figure out sorting
        else:
            print("No File to Unpickle")

    def getFilterText(self):
        t = self.queueFilter.currentText()
        if t == "Thermo Element (*.dat)":
            self.fileFilter = "Element Data Files (*.dat)"
            self.extension = '.dat'
        elif t == "P-E NexION":
            self.fileFilter = "Perkin-Elmer NexION (*.*)"
            self.extension = '.csv'
        else:
            self.fileFilter = "Any File (*.*)"
            self.extension = '.txt'


if __name__ == "__main__":
    import sys
    from src.records.Session import Session

    app = QApplication(sys.argv)
    session = Session()
    ui = Ui_ImportWidget(session)
    ui.setupUi()
    ui.show()
    sys.exit(app.exec())
