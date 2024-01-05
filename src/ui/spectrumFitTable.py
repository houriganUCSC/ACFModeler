import sys

import PyQt6.QtCore
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib import ticker
import numpy as np
import re
from PyQt6.QtWidgets import QApplication, QWidget, QTableWidget, QComboBox, QTabWidget, \
    QCheckBox, QVBoxLayout, QTableWidgetItem, QSizePolicy, QHBoxLayout
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QFont
from sklearn.preprocessing import PolynomialFeatures
import statsmodels.api as sm
# from src.acfRegression import SessionFits
from src.records import Session as Session
from src.records.Session import SpectrumFits
from src.ui.headerCodes import CrossCalHeaderCodes
from src.ui.spectrumModelDesignTable import ModelDesignTable

plusminus = u"\u00B1"
class CrossCalFitType(QComboBox):
    """
    QCombobox selector for type of fit to use for cross calibration
    'self' use data from isotope itself to for cross calibration
    'internal' use data from all masses to create cross-cal vs. mass
    'external' use externally calibrated value

    """
    def __init__(self, row, col):
        super().__init__()
        self.setStyleSheet('font-size: 12px')
        self.addItems(['Self', 'Internal', 'External'])
        self.setProperty('row', row)
        self.setProperty('col', col)

    def getComboValue(self):
        return self.currentText()

    def disableSelfCal(self):
        self.model().item(0).setEnabled(False)


class TauFitType(QComboBox):
    """
    QCombobox selector for type of fit to use for cross calibration
    'self' use data from isotope itself to for cross calibration
    'internal' use data from all masses to create cross-cal vs. mass
    'external' use machine stored deadtime value
    """

    def __init__(self, row, col):
        super().__init__()
        self.setStyleSheet('font-size: 12px')
        self.addItems(['Self', 'Internal', 'External'])
        self.setProperty('row', row)
        self.setProperty('col', col)

    def getComboValue(self):
        return self.currentText()

    def disableSelfCal(self):
        self.model().item(0).setEnabled(False)


class CustomCheck(QWidget):
    toggled = pyqtSignal()
    def __init__(self, row, col, name):
        super().__init__()
        self.setProperty('row', row)
        self.setProperty('col', col)
        self.setProperty('name', name)
        self.cbox = QCheckBox()
        self.cbox.setStyleSheet('background-color: transparent')
        self.cbox.released.connect(self.emitToggled)
        layout = QHBoxLayout(self)
        layout.addWidget(self.cbox)
        layout.setAlignment(Qt.AlignmentFlag.AlignCenter)
        layout.setContentsMargins(0,0,0,0)
        self.setLayout(layout)

    def isChecked(self):
        return int(self.cbox.isChecked())

    def setChecked(self, state: bool = True):
        self.cbox.setChecked(state)

    def emitToggled(self):
        self.toggled.emit()


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
        # Disable editing
        self.setFlags(self.flags() & ~Qt.ItemFlag.ItemIsEditable)
        # self.row = -1
        # self.col = -1

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


class CrossCalHeaderCodes(dict):
    # Tau Type b"\xCE\xB1\x20\x53\x6F\x75\x72\x63\x65"
    def __init__(self):
        self.codes = {
            'Mass': {'hexCode': 'Mass', "numFormat":'0.2f', "width": 75},
            'n': {'hexCode': 'n', "numFormat":'0.1f', "width": 50},
            'an': {'hexCode': 'an', "numFormat": '0.1f', "width": 50},
            'pMax': {'hexCode': 'pMax', "numFormat": '0.2f', "width": 50},
            'ACF Type': {'hexCode':b"\xCE\xB1\x20\x53\x6F\x75\x72\x63\x65", "numFormat":'acfCombo', "width": 100},
            'use ACF': {'hexCode': b"\x75\x73\x65\x20\xCE\xB1", "numFormat": 'checkBox', "width": 50},
            'a1': {'hexCode':b"\x61\xE2\x82\x81", "numFormat": 'val+error', "width": 100},
            'a2': {'hexCode':b"\x61\xE2\x82\x82", "numFormat": 'val+error', "width": 100},
            'Tau Type': {'hexCode': b"\xCF\x84\x20\x53\x6F\x75\x72\x63\x65", "numFormat": 'tauCombo', "width": 100},
            'use Tau': {'hexCode': b"\x75\x73\x65\x20\xCF\x84", "numFormat": 'checkBox', "width": 50},
            'tau': {'hexCode':b"\xCF\x84", "numFormat": 'val+error', "width": 100},
            'int_a1': {'hexCode':b"\x69\x6E\x74\x20\x61\xE2\x82\x81", "numFormat": 'val+error', "width": 110},
            'int_a2': {'hexCode':b"\x69\x6E\x74\x20\x61\xE2\x82\x82", "numFormat": 'val+error', "width": 110},
            'int_tau': {'hexCode':b"\x69\x6E\x74\x20\xCF\x84", "numFormat": 'val+error', "width": 110},
            'se_a1': {'hexCode':b"\xCF\x83\x28\x61\xE2\x81\x81\x29", "numFormat":'.4G', "width": 50},
            'se_a2': {'hexCode':b"\xCF\x83\x28\x61\xE2\x81\x82\x29", "numFormat":'.4G', "width": 50},
            'se_tau': {'hexCode':b"\xCF\x83\x28\xCF\x84\x29", "numFormat":'.4G', "width": 50},
            'rSqr': {'hexCode':b"\x72\xC2\xB2", "numFormat":'0.3f', "width": 50},
            "redChi2": {'hexCode':b"\xCF\x87\xE1\xB5\xA3\xC2\xB2", "numFormat":'0.2f', "width": 50}
            }

    def getHeaderLabel(self, key: str):
        hdr = self.codes[key]["hexCode"]
        if isinstance(hdr, bytes):
            label = hdr.decode('UTF-8')
        else:
            label = hdr
        return label

    def getDefaultWidth(self, key:str):
        return self.codes[key]["width"]

    def formatValue(self, key:str, val: float):
        valText = ''
        if not np.isnan(val):
            valFormat = self.codes[key]["numFormat"]
            valText = f'{val:{valFormat}}'
            if 'G' or 'E' in valFormat:
                valText = valText.replace("E-0", "E-").replace("E+0", "E+")
        return valText


class CrossCalTableWidget(QTableWidget):
    """
    Populate table
    :param ccdf:
    :return:

    """
    tableModified = pyqtSignal()
    useModified = pyqtSignal()
    plotsDrawn = pyqtSignal()
    dataPostProcessed = pyqtSignal()
    def __init__(self, parent):
        self.colHeaders = ['Mass', 'n', 'an', 'pMax',
                      'ACF Type', 'use ACF', 'a1', 'a2',
                      'Tau Type', 'use Tau', 'tau',
                      'rSqr', 'redChi2',
                      'int_a1', 'int_a2',
                      'int_tau']
        self.hdr = CrossCalHeaderCodes()
        self.fnt = QFont()
        self.fnt.setPixelSize(12)
        self.fnt.setBold(True)
        self.pars = ['a1', 'a2', 'tau']
        self.session = parent.session
        self.fig = parent.fig
        self.canvas = parent.canvas
        self.isotopes = self.session.isotopes
        self.modelSetup = parent.modelTable
        self.whichPlots = parent.plotSelectorCombo
        super().__init__(len(self.isotopes), len(self.colHeaders))
        self.setVerticalHeaderLabels(self.isotopes)
        for i,colName in enumerate(self.colHeaders):
            label = self.hdr.getHeaderLabel(colName)
            self.setHorizontalHeaderItem(i, QTableWidgetItem(label))
            self.horizontalHeaderItem(i).setFont(self.fnt)
            self.setColumnWidth(i, self.hdr.getDefaultWidth(colName))
        for j in range(len(self.isotopes)):
            self.verticalHeaderItem(j).setFont(self.fnt)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.verticalHeader().setDefaultSectionSize(25)

        self.setStyleSheet('selection - color: rgb(0, 0, 0);'
                           'selection - background - color: rgb(230, 230, 230);'
                           'background - color: rgb(255, 255, 255)'
                           'color: rgb(0, 0, 0);')
        if self.session.status['fit']:
            self.session.spectrumFit = SpectrumFits()
            self.session.spectrumFit.getMassFits(self.session)
            self.initFitTable()
            self.recommendFitType()
            self.createPlots()

    def initFitTable(self):
        sf = self.session.spectrumFit
        for i, isotope in enumerate(self.isotopes):
            for j, colName in enumerate(self.colHeaders):
                colType = self.hdr.codes[colName]['numFormat']
                alphaNaN = np.isnan(sf['a1']['Y'][i]) or np.isnan(sf['a1']['Y'][i])
                tauNaN = np.isnan(sf['tau']['Y'][i])

                if colType in ['acfCombo', 'tauCombo']:
                    if colType == 'acfCombo':
                        cellWidget = CrossCalFitType(i, j)
                        if alphaNaN:
                            self.session.spectrumFit['alphaSource'][i] = 1
                            cellWidget.setCurrentIndex(1)
                            cellWidget.disableSelfCal()
                        else:
                            self.session.spectrumFit['alphaSource'][i] = 0
                            cellWidget.setCurrentIndex(0)
                            cellWidget.currentIndexChanged.connect(self.tableModified)
                    elif colType == 'tauCombo':
                        cellWidget = TauFitType(i, j)
                        if tauNaN:
                            self.session.spectrumFit['tauSource'][i] = 1
                            cellWidget.setCurrentIndex(1)
                            cellWidget.disableSelfCal()
                        else:
                            self.session.spectrumFit['tauSource'][i] = 0
                            cellWidget.setCurrentIndex(0)
                            cellWidget.currentIndexChanged.connect(self.tableModified)
                    self.setCellWidget(i, j, cellWidget)

                elif colType == 'checkBox':
                    cellWidget = CustomCheck(i, j, colName)
                    if colName == 'use ACF':
                        if alphaNaN:
                            self.session.spectrumFit['a1']['mask'][i] = False
                            self.session.spectrumFit['a2']['mask'][i] = False
                            cellWidget.setChecked(False)
                            cellWidget.cbox.setEnabled(False)
                        else:
                            self.session.spectrumFit['a1']['mask'][i] = True
                            self.session.spectrumFit['a2']['mask'][i] = True
                            cellWidget.setChecked(True)
                            cellWidget.toggled.connect(self.useModified)
                    elif colName == 'use Tau':
                        if tauNaN:
                            self.session.spectrumFit['a1']['mask'][i] = False
                            cellWidget.setChecked(False)
                            cellWidget.cbox.setEnabled(False)
                        else:
                            self.session.spectrumFit['tau']['mask'][i] = True
                            cellWidget.setChecked(True)
                            cellWidget.toggled.connect(self.useModified)
                    self.setCellWidget(i, j, cellWidget)

                elif colType == 'val+error':
                    if colName in ['int_a1', 'int_a2','int_tau']:
                        valText = ''
                    else:
                        val = self.session.spectrumFit[colName]['Y'][i]
                        se_val = self.session.spectrumFit[colName]['seY'][i]
                        if np.isnan(val):
                            valText = ''
                        else:
                            rse_val = 100 * se_val / val
                            if rse_val > 100:
                                rseText = f'{rse_val:0.1E}%'
                            elif rse_val > 10:
                                rseText = f'{rse_val:.0F}%'
                            elif rse_val > 1:
                                rseText = f'{rse_val:.1F}%'
                            elif rse_val > 0.1:
                                rseText = f'{rse_val:.2F}%'
                            elif rse_val> 0.01:
                                rseText = f'{rse_val:.2F}%'
                            else:
                                rseText = f'{rse_val:0.1E}%'
                            valText = f'{val:.4G} ({rseText})'.replace("E-0", "E-").replace("E+0", "E+")
                    cellWidget = TextItem(i, j, valText)
                    self.setItem(i, j, cellWidget)
                else:
                    if colName in ['rSqr', 'redChi2']:
                        val = self.session.spectrumFit[colName][i]
                    elif colName == 'Mass':
                        val = self.session.spectrumFit['mass'][i]
                    elif colName == 'n':
                        val = self.session.spectrumFit['nQual'][i]/1000
                    elif colName == 'an':
                        val = self.session.spectrumFit['anOnly'][i]/1000
                    elif colName == 'pMax':
                        val = self.session.spectrumFit['pMax'][i] / 1E+6
                    valText = self.hdr.formatValue(colName, val)
                    cellWidget = TextItem(i, j, valText)
                    self.setItem(i, j, cellWidget)
        self.resizeColumnsToContents()

    def updateTable(self):
        for i, par in enumerate(self.pars):
            col = self.colHeaders.index('int_'+par)
            for row, intText in enumerate(self.session.spectrumFit[par]['intText']):
                cellWidget = TextItem(row, col, intText)
                self.setItem(row, col, cellWidget)
        self.resizeColumnsToContents()

    def recommendFitType(self):
        sf = self.session.spectrumFit
        """ Evaluate alpha for Fit Type Suggestion """
        a1errs = np.array(sf['a1']['seY']) / np.array(sf['a1']['Y'])
        a2errs = np.array(sf['a2']['seY']) / np.array(sf['a2']['Y'])
        rSqrs = sf['rSqr']
        for i, (a1e, a2e, r) in enumerate(zip(a1errs, a2errs, rSqrs)):
            cbox = self.cellWidget(i, self.colHeaders.index('use ACF'))
            source = self.cellWidget(i, self.colHeaders.index('ACF Type'))
            ok = a1e < 0.005 and a2e < 0.01 and r > 0.95
            cbox.setChecked(ok)
            sf['a1']['mask'][i] = ok
            sf['a2']['mask'][i] = ok
            sf['alphaSource'][i] = int(~ok)
            source.setCurrentIndex(int(~ok))


        """ Evaluate Tau for Fit Type Suggestion """
        tauerrs = np.array(sf['a1']['seY']) / np.array(sf['a1']['Y'])
        pmaxs = sf['pMax']
        for i, (te, r, p) in enumerate(zip(tauerrs, rSqrs, pmaxs)):
            cbox = self.cellWidget(i, self.colHeaders.index('use Tau'))
            source = self.cellWidget(i, self.colHeaders.index('Tau Type'))
            ok =  te < 0.01 and r > 0.97 and p>1E+6
            cbox.setChecked(ok)
            sf['tau']['mask'][i] = ok
            sf['tauSource'][i] = int(~ok)
            source.setCurrentIndex(int(~ok))
            cbox.setChecked(ok)

    def createPlots(self, **kwargs):
        self.session.spectrumFit.fitByMass(self.modelSetup)
        self.updateTable()

        whichPlots = self.whichPlots.currentText()
        plots = ['a1', 'a2', 'tau', 'ACF', 'Drift']
        labels = [f'a$_{1}$ (x10$^{{-3}}$ sec)',
                  f'a$_{2}$ (x10$^{{-6}}$)',
                  f'$\\tau$ (x10$^{{-9}}$ sec)',
                  f'$\\alpha$$_{0}$ (sec$^{{-1}}$)',
                  f'Drift (%)']
        ymults =  [1E+3, 1E+6, 1E+9, 1, 1]
        self.fig.clear()
        axs = None
        if whichPlots == 'All':
            plots = plots
            ax0 = self.canvas.figure.add_subplot(3, 1, 1)
            ax1 = self.canvas.figure.add_subplot(3, 1, 2)
            ax2 = self.canvas.figure.add_subplot(3, 1, 3)
            self.axs = [ax0, ax1, ax2]
            plots = plots[:3]
            labels = labels[:3]
            ymults = ymults[:3]
        else:
            idx = plots.index(whichPlots)
            plots = [plots[idx]]
            labels = [labels[idx]]
            ymults = [ymults[idx]]
            self.axs = [self.canvas.figure.add_subplot(1, 1, 1)]

        anno = True
        annoPos = (0.1, 0.9)
        for key, value in kwargs.items():
            if key == 'xlim':
                for ax in self.axs:
                    ax.set_xlim(value)
            if key == 'ylim' and len(self.axs) == 1:
                self.axs[0].set_ylim(value)
            if key == 'anno':
                annoPositions = [None,(0.1, 0.9), (0.4, 0.9), (0.7, 0.9),
                                 (0.1, 0.5), (0.7, 0.5),
                                 (0.1, 0.1), (0.4, 0.1), (0.7, 0.1)]
                annoPos = annoPositions[value]

        sf = self.session.spectrumFit
        for i, (plot, label, ymult) in enumerate(zip(plots, labels, ymults)):
            sf = self.session.spectrumFit
            x = np.array(sf['mass'])
            if plot == 'ACF':
                mask = np.array(sf['a1']['mask'])
                y = 1/np.array(np.where(mask, sf['a1']['Y'], np.nan))
                # Todo: Add Formal Error Prop
                # yerr = 1/np.array(np.where(mask, sf['a1']['seY'], np.nan))
                yerr = np.zeros_like(y)
            elif plot == 'Drift':
                mask = np.array(sf['a1']['mask'])
                dt = np.max(self.session.scanTime)-np.min(self.session.scanTime)
                y0 = np.array(np.where(mask, sf['a1']['Y'], np.nan))
                dydt = np.array(np.where(mask, sf['a2']['Y'], np.nan))
                acf0 = 1/y0
                acff = 1/(y0+dydt*dt)
                y = 100*(1-acf0/acff)
                # Todo: Add Prop Error
                yerr = np.zeros_like(y)
            else:
                mask = np.array(sf[plot]['mask'])
                y = np.array(np.where(mask, sf[plot]['Y'], np.nan))
                yerr = np.array(np.where(mask, sf[plot]['seY'], np.nan))
            self.axs[i].errorbar(x, y*ymult, yerr*ymult, mec= 'black', mfc= 'blue', fmt='o', capthick=1, capsize=5)
            self.axs[i].set_ylabel(label)
            self.axs[i].set_xlabel('Mass (amu)')
            self.axs[i].get_xaxis().set_visible(i == len(self.axs)-1)

            if plot == 'a1':
                ax2ticks = [ymult/tick for tick in self.axs[i].get_yticks()]
                ax2tickLabels = self.axs[i].get_yticklabels()
                ax2 = self.axs[i].twinx()
                lims = np.array(self.axs[i].get_ylim())/ymult
                ax2.set_ylim((1/(lims[0]), 1/(lims[1])))
                ax2.set_yticks(ax2ticks)
                ax2.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}"))
                ax2.set_ylabel(f"$\\alpha_{0}$")


            # Todo: Add drift % scale to a2 plot. May be too complicated with variable ACF_0
            # if plot == 'a2':
            #     endTime = np.max(self.session.scanTime) - np.min(self.session.scanTime)
            #     ax2ticks = [1/endTime*tick for tick in axs[i].get_yticks()]
            #     ax2tickLabels = axs[i].get_yticklabels()
            #     ax2 = axs[i].twinx()
            #     lims = axs[i].get_ylim()
            #     ax2.set_ylim((1/endTime*lims[0], 1/endTime*lims[1]))
            #     ax2.set_yticks(ax2ticks)
            #     ax2.set_ylabel(f"$\\Delta$$\\alpha$")

            if plot == 'tau':
                self.axs[i].axhline(self.session.machineDeadTime * ymult)
                ax2ticks = np.array(self.axs[i].get_yticks())/ymult
                obs = 5E+6/(1+5E+6*self.session.machineDeadTime)
                corrTicks= 100*(obs/(1-obs*ax2ticks)/5E+6 - 1)

                ax2 = self.axs[i].twinx()
                axlims = np.array(self.axs[i].get_ylim())/ymult
                corrlims = 100*(obs/(1-obs*axlims)/5E+6 -1)
                ax2.set_ylim(corrlims)
                ax2.set_yticks(corrTicks)
                ax2.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.1f}%"))
                ax2.set_ylabel(f"$\Delta\\tau_{{machine}}$ @ 5 Mcps")
                # Todo: Add machine tau annotation

            # Outliers
            if plot == 'ACF':
                yout = 1/np.array(np.where(~mask, sf['a1']['Y'], np.nan))
                # Todo: Add Formal Error Prop
                # yerrout = 1/np.array(np.where(~mask, sf['a1']['seY'], np.nan))
                yerrout = np.zeros_like(yout)
            elif plot == 'Drift':
                dt = np.max(self.session.scanTime)-np.min(self.session.scanTime)
                y0 = np.array(np.where(~mask, sf['a1']['Y'], np.nan))
                dydt = np.array(np.where(~mask, sf['a2']['Y'], np.nan))
                acf0 = 1/y0
                acff = 1/(y0+dydt*dt)
                yout = 100*(1-acf0/acff)
                # Todo: Add Prop Error
                yerrout = np.zeros_like(y)
            else:
                yout = np.array(np.where(~mask, sf[plot]['Y'], np.nan))
                yerrout = np.array(np.where(~mask, sf[plot]['seY'], np.nan))
            self.axs[i].errorbar(x, ymult*yout, ymult*yerrout, mec= 'black', mfc= 'yellow', fmt='o', capthick=1, capsize=5)

            if not plot in ['ACF', 'Drift']:
                # Best fit line
                ypred = np.array(self.session.spectrumFit[plot]['best']['Y'])
                self.axs[i].plot(x, ymult*ypred, linestyle='solid', linewidth=1, color='red')

                # Lower confidence interval
                lci = np.array(self.session.spectrumFit[plot]['ci']['l'])
                self.axs[i].plot(x, ymult*lci, linestyle='solid', linewidth=0.5, color='darkgray')

                # Upper confidence interval
                uci = np.array(self.session.spectrumFit[plot]['ci']['u'])
                self.axs[i].plot(x, ymult*uci, linestyle='solid', linewidth=0.5, color='darkgray')

                # Fill between upper and lower CIs
                self.axs[i].fill_between(x, ymult*uci, ymult*lci, color='lightgray', alpha=0.5)

                # Annotate graph
                if annoPos is not None:
                    self.axs[i].annotate(self.buildEquationString(plot,False, False),
                                         xy=annoPos, xycoords='axes fraction', verticalalignment='top',fontsize=8)

            self.axs[i].get_xaxis().set_visible(i == len(plots) - 1)

        self.plotsDrawn.emit()
        self.canvas.draw()

    def buildEquationString(self, par, inclUnc = True, addGoF = False):
        # Build annotation string
        vals = self.session.spectrumFit[par]['coefs']
        se_vals = self.session.spectrumFit[par]['se_coefs']
        rSqr = self.session.spectrumFit[par]['rSqr']
        redChi2 = self.session.spectrumFit[par]['rChi2']
        eq = par + ' = '
        if par == 'a1':
            eq = f'a$_{1}$ = '
        elif par == 'a2':
            eq = f'a$_{2}$ = '
        elif par == 'tau':
            eq = f'$\\tau$ = '


        for j,[val, se_val] in enumerate(zip(vals, se_vals)):
            coeffs = f'{val: 0.2E}'
            if inclUnc:
                rse = 100*(se_val/val)
                coeffs = coeffs + f'($\pm${rse: 0.1f}%)'
            if j == 0:
                eq = eq + coeffs
            elif j == 1:
                eq = eq + ' + ('+ coeffs +')M'
            else:
                eq = eq + ' + (' + coeffs + f')M$^{j}$'
            # if inclUnc or j == len(vals) - 1:
            #     eq = eq + '\n'
        if addGoF:
            eq = eq + f'r$^{2}$= {rSqr: 0.3f}' + f'MSWD= {redChi2: 0.2f}'
        return eq

    @PyQt6.QtCore.pyqtSlot()
    def useModified(self):
        sendy = self.sender()
        row = sendy.property("row")
        pars = []
        if sendy.property("name") == 'use ACF':
            pars = ['a1', 'a2']
        elif sendy.property("name") == 'use Tau':
            pars = ['tau']
        for par in pars:
            self.session.spectrumFit[par]['mask'][row] = sendy.cbox.isChecked()
        self.session.spectrumFit.fitByMass(self.modelSetup)
        self.updateTable()
        self.createPlots()

    def commitFits(self):
        sf = self.session.spectrumFit
        for isotope, massRecord in self.session.masses.items():
            if not hasattr(massRecord, 'postProcessPars'):
                massRecord.__setattr__('postProcessPars', {})
            row = list(self.isotopes).index(isotope)
            # Source strings
            a = self.cellWidget(row, 4).currentText()
            t = self.cellWidget(row, 8).currentText()
            massRecord.postProcessPars["alphaSource"] = a
            massRecord.postProcessPars["tauSource"] = t
            for par, source in zip(['a1', 'a2', 'tau'],[a, a, t]):
                for var in ['','se']:
                    val = None
                    if source == 'Self':
                        val = sf[par][var+'Y'][row]
                    elif source == 'Internal':
                        val = sf[par]['best'][var+'Y'][row]
                    massRecord.postProcessPars[var+par] = val
            massRecord.postProcessTimeSeries()
            self.dataPostProcessed.emit()
            """ Debugging Plots """
            # y = massRecord.modeledTimeSeries
            # x = massRecord.timeSeries
            # plt.scatter(y,x/y, marker='.', s=0.01)
            # plt.title(isotope)
            # plt.show()

class AppDemo(QWidget):
    def __init__(self):
        super().__init__()
        isotopes = ['Hg202', 'Pb204', 'Pb206', 'Pb207', 'Pb208', 'Th232', 'U235', 'U238']
        mainLayout = QVBoxLayout()
        self.table = CrossCalTableWidget(isotopes)
        self.table.initFitDictionary()
        w = self.table.size().width()
        h = self.table.size().height()
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



