import typing
import numpy as np
import statsmodels.api as sm
from mpl_toolkits.mplot3d import axes3d
import writeThermoElement
from matplotlib import pyplot as plt
from PyQt5.QtWidgets import QApplication, QWidget, QTableWidget, QComboBox, QVBoxLayout, QTableWidgetItem
from PyQt5.QtCore import Qt
from sklearn.preprocessing import PolynomialFeatures
from statsmodels.sandbox.regression.predstd import wls_prediction_std

# MODEL
class IsotopeFit():

    def __init__(self):
        self.massCenter = 0
        self.dwell = 0.01
        self.fitParams = {
                "ACF Type": "",
                "Tau Type": "",
                "n": 0,
                "of": 0,
                "a1": np.nan,
                "a2": np.nan,
                "se_a1": np.nan,
                "se_a2": np.nan,
                "tau": np.nan,
                "se_tau": np.nan,
                "rSqr": np.nan,
                "redChi2": np.nan,
                "ea1": np.nan,
                "ea2": np.nan,
                "se_ea1": np.nan,
                "se_ea2": np.nan,
                "etau": np.nan,
                "se_etau": np.nan
            }
        self.allData = {
            "time":np.array([]),
            "P":np.array([]),
            "A":np.array([]),
            "inliers":np.array([]),
            "qualifiers":np.array([])
        }
        self.massTimeSeries = np.array([])

    def loadData(self,cycle_time, P, A, dwell):
        """
        :param cycle_time: (float) time in seconds relative to beginning of session
        :param P: 2D float array of deadtime corrected pulse count data of dimensions "cycles" x "channels"
        :param A: 2D float array of analog values data (bit value of ADC for Element)
        :param dwell: (float) scalar of isotope-specific dwell time.  Used to calculate count statistics error
        :return: populated allData dictionary with 2D arrays
        """
        self.allData["time"]= cycle_time
        self.allData["P"]= P
        self.allData["A"]= A
        self.allData["W"] = P**3*dwell
        self.allData["qualifiers"] = np.ones_like(cycle_time)
        self.allData["inliers"] = np.ones_like(cycle_time)
        self.dwell = dwell

    def qualify(self, Pmax, Pmin = 1, Amin = 1):
        """
        qualify creates a 1D logical array where True satisfies logical test
        if a datum "qualifies" the  Pmax > Pi > Pmin AND Ai > Amin
        :param Pmax: integer of minimum pulse counts
        :param Pmin: integer of maximum pulse counts
        :param Amin: integer of minimum bit value
        :return: qualifiers - 1D logical array where True is in the appropriate A and P range
        """
        self.allData["qualifiers"] = np.where(Pmin<self.allData["P"]<Pmax and self.allData["A"]>Amin)

    def filterOutlier(self, tukey = 3):
        """
        Filters outliers based on the median & interquartile range (iqr).
        This non-parametric filter performs better on spikes than mean / SD filters
        Spikes - short period high amplitude transients are unlikely to "qualify"
        With a gross outlier of 3xiqr ~1% of the data are rejected.
        :param tukey: float that is a multiplier
        :return: inliers: 1D logical array of not gross outliers
        """
        if tukey is None:
            return
        acf = self.allData["P"]/self.allData["A"]
        acf[np.isinf(acf)] = np.nan
        med = np.nanmedian(acf)
        Q3 = np.nanpercentile(acf,75,interpolation='midpoint')
        Q1 = np.nanpercentile(acf,25,interpolation='midpoint')
        iqr = Q3-Q1
        self.allData["inliers"] = np.abs(acf-med)/iqr <= tukey

    def regress(self, deadtime = 0, regressionType = "Weighted", mType = None):
        """
        Regresses 1/P vs. 1/A & t/A
        :param deadtime: float used to remove deadtime correction from corrected data
        :param regressionType: string type of linear regression model.
            "Weighted",
            "Ordinary",
            "Robust"
            "Weighted Robust"
        :param mType: string M-estimator type (default: HuberT) for robust linear regression
            https://www.statsmodels.org/stable/rlm.html
            "HuberT",
            "Hampel",
            "LeastSquares",
            "AndrewWave",
            "RamsayE",
            "RobustNorm",
            "TrimmedMean",
            "TukeyBiweight"
        :return: parameters of regression; uploaded to IsotopeFit class data
        """
        self.regressType = regressionType
        if regressionType == "Robust":
            accepts = self.allData["qualifiers"]
        else: accepts = self.allData["inliers"] & self.allData["qualifiers"]
        self.fitParams["n"] = accepts.count(True)
        self.fitParams["of"] = len(accepts)
        P = self.allData["P"][accepts]
        A = self.allData["A"][accepts]
        t = self.allData["time"][accepts]
        W = self.allData["W"][accepts]
        P = P/(1+P*deadtime)  #Remove Deadtime Correction
        y = 1/P
        X = np.ones(len(y),3)
        X[:,1] = 1/A
        X[:,2] = t/A
        if regressionType == "Weighted":
            model = sm.WLS(y, X, weights=W)
        elif regressionType == "Ordinary":
            model = sm.OLS(y, X)
        elif regressionType == "Robust" or "Weighted Robust":
            if mType == "HuberT": M =sm.robust.norms.HuberT()
            elif mType == "Hampel": M =sm.robust.norms.Hampel()
            elif mType == "LeastSquares": M =sm.robust.norms.LeastSquares()
            elif mType == "AndrewWave": M =sm.robust.norms.AndrewWave()
            elif mType == "RamsayE": M =sm.robust.norms.RamsayE()
            elif mType == "RobustNorm": M =sm.robust.norms.RobustNorm()
            elif mType == "TrimmedMean": M =sm.robust.norms.TrimmedMean()
            elif mType == "TukeyBiweight": M =sm.robust.norms.TukeyBiweight()
            else: M = sm.robust.norms.TukeyBiweight()
            if regressionType == "Weighted Robust":
                w = len(W)*np.sqrt(W)/np.sum(np.sqrt(W))
                y = w*y
                X = w*X
            model = sm.RLM(y, X, M)
        results = model.fit()
        self.fitParams["tau"] = results.params[0]
        self.fitParams["a1"] = results.params[1]
        self.fitParams["a2"] = results.params[2]
        self.fitParams["se_tau"] = results.bse[0]
        self.fitParams["se_a1"] = results.bse[1]
        self.fitParams["se_a2"] = results.bse[2]
        self.fitParams["rSqr"] = results.rsquared

    def plotIsotope(self, deadtime):
        """
        :param deadtime: (float) for removal of deadtime correction
        :return:
        """
        if self.regressType == "Robust":
            accepts = self.allData["qualifiers"]
        else: accepts = self.allData["inliers"] & self.allData["qualifiers"]
        P = self.allData["P"]
        A = self.allData["A"]
        t = self.allData["time"]
        P = P/(1+P*deadtime)
        z = 1/P
        x = 1/A
        y = t/A
        tau = self.fitParams["tau"]
        a1 = self.fitParams["a1"]
        a2 = self.fitParams["a2"]
        (xMesh, yMesh) = np.meshgrid(np.linspace(min(x),max(x),num = 10, endpoint=True),
                                     np.linspace(min(y),max(y),num = 10, endpoint=True))
        zMesh = tau + a1*xMesh + a2*yMesh
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.scatter3d(x[accepts],y[accepts],z[accepts],color ='green', marker ='.')
        ax.scatter3d(x[not accepts], y[not accepts], z[not accepts], color ='red', marker ='.')
        wf = ax.plot_surface(xMesh, yMesh, zMesh, color='blue',linewidth = 0.5)
        wf.set_facecolor((0.5,0.5,0.5,0.25)) #50% greyscale; 25% opacity
        return x,y,z

    def processTimeSeries(self, p, a, t, pMax):
        """
        Produces a 1D time series of intensity data post-processed with user selected cross-calibration
        parameters.  Tau Type and ACF Type relate to the QComboBox values of A QTableWidget

        :param p: deadtime corrected pulse count value
        :param a: analog value (bit count of ADC for Element 2/XR)
        :param t: time of observation since beginning of session
        :param pMax: Maximum pulse count value for
        :return: loads massTimeSeries into class data
        """
        if self.fitParams["Tau Type"] ==  "Internal":
            tau = self.fitParams["tau"]
        elif self.fitParams["Tau Type"] == "Model":
            tau = self.fitParams["etau"]
        elif self.fitParams["Tau Type"] == "Global":
            tau = 0 #place holder for Keep Original Deadtime
        if not np.isnan(tau):
            p = p/(1-p*tau)

        # ACF and Drift Correct
        if self.fitParams["ACF_Type"] == "Internal":
            a1 = self.fitParams["a1"]
            a2 = self.fitParams["a2"]
        elif self.fitParams["ACF_Type"] == "Model":
            a1 = self.fitParams["ea1"]
            a2 = self.fitParams["ea2"]

        acf = 1/(a1+a2*t)
        # Calculate analog equivalent pulses across all channels
        aCounts = (a.T * acf).T
        useAnalog = np.logical_or(p>pMax, np.isnan(p))
        p = np.where(useAnalog, aCounts,p)

        # Average for cycles across all channels
        self.massTimeSeries = np.nanmean(p, axis=1)

class SessionFits():
    """
    SessionFits class contains attributes and methods for modeling an entire sessions worth of data.
    This class heavily leverages the IsotopeFit class which contains attributes and methods for
    processing each isotope in session.
    """

    def __init__(self):
        self.regressType = None
        self.mType = None
        self.deadtime = 0
        self.pMax = 5E+6
        self.pMin = 1
        self.aMin = 1
        self.tukeyThreshold = 3
        self.fitOrder ={"a1":3,
                        "a2":0,
                        "tau":0}
        self.weight ={"a1":True,
                      "a2":True,
                      "tau":True}
        self.sessionFits = {}  #Dictionary of IsotopeFits from IsotopeFit class

    def extract(self, dataDir):
        """
        Extracts all data from imported data directory (dataDir) which contains a cycles x channels
        P and A values.  2D array is reshaped into a 1D arrays for

        :param dataDir: Dictionary of raw data in a hierarchy of Directory-->File
        :return: Data loaded into class data with self.loadData()
        """

        for dir,files in dataDir.items():
            for file in list(files.keys()):
                s = dataDir[dir][file]["data"]
                abs_t = dataDir[dir][file].scans["scanTime"]
                mass_del_t = 0
                for i, isotope in enumerate(s.isotopes):
                    self.sessionFits[isotope] = IsotopeFit()  #New instance of IsotopeFit class
                    dwell = s.masses[isotope].dwell
                    p2d = s.masses[isotope].raw["pulse"]
                    rows, cols = p2d.shape
                    p = np.reshape(p2d,(rows*cols,1))
                    a2d = s.masses[isotope].raw["analog"]
                    a = np.reshape(a2d,(rows*cols,1))
                    t2d = np.zeros_like(p2d)
                    for j in rows:
                        for k in cols:
                            abs_t[j]+mass_del_t+(k*dwell)
                    mass_del_t += s.masses[isotope].totaldwell
                    t = np.reshape(t2d,(rows*cols,1))
                    self.sessionFits[isotope].loadData(t,p,a,dwell)

    def filter(self, pMax, pMin, aMin, tukeyThreshold = None):
        """
        For all isotopes:
        1. populates "qualifiers" (1D logical array) based on whether data are between pMin & pMax and greater than Amin
        2. populates "inliers" (1D logical array) based on interquartile range test

        :param pMax: Maximum pulse count rate
        :param pMin: Minimum pulse count rate
        :param aMin: Minimum analog value
        :param tukeyThreshold: gross outlier multiplier (tukeyThreshold * interquartile range)
        :return: updates class data with qualifiers and inliers
        """
        for vals in self.sessionFits.values():
            self.pMax = pMax
            self.pMin = pMin
            self.aMin = aMin
            vals.qualify(pMax, pMin, aMin)
            if tukeyThreshold is not None:
                self.tukeyThreshold = tukeyThreshold
                vals.filterOutlier(tukeyThreshold)

    def regress(self, deadtime, regressionType, mType):
        """
        For all isotopes in session:
        1. performs linear regression using methods in the IsotopeFit class

        :param deadtime: float deadtime value for removal of deadtime correction imposed on machine-reported data
        :param regressionType: string regression type (ordinary, weighted or robust)
        :param mType: string M-estimator type for robust regression (see IsotopeFit class for details)
        :return: Regressed values are loaded into IsotopeFit class data
        """
        self.deadtime = deadtime
        self.regressType = regressionType
        if mType is not None:
            self.mType = mType
        for vals in self.sessionFits.values():
            vals.regress(deadtime, regressionType, mType)

    def plotMassFit(self, isotope):
        s = self.sessionFits[isotope]
        x,y,z = s.plotIsotope(self.deadtime)
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.scatter3d(x,y,z)

    def fitsByMass(self):
        """
        1. Extracts cross calibration fits for each mass into self.massFits
        2. Filters fits based on user input whether to use "internal" fits (from IsotopeFit class) or "external"
        fits (a1, q2, tau vs. mass)
        3. Plots used and unused data for fit v. mass
        4. Plots model fit and confidence intervals around fits
        5. Uploads fit and uncertainty (confidence range/2) as "external fits" in IsotopeFit class

        :return:
        """
        session = self.sessionFits
        fig,(ax1,ax2,axtau) = plt.subplots(1,3, sharex=True)
        massFits = {
            "isotope": np.array([]),
            "mass":np.array([]),
            "acfInternal":([]),
            "tauInternal":([]),
            "a1": np.array([]),
            "a1err": np.array([]),
            "a2": np.array([]),
            "a2err": np.array([]),
            "tau": np.array([]),
            "tauerr": np.array([]),
        }
        for isotope in list(session.keys()):
            s = session[isotope]
            massFits["acfInternal"] = np.append(s.fitParams["ACF TYPE"] == 0)
            massFits["tauInternal"] = np.append(s.fitParams["Tau TYPE"] == 0)
            massFits["mass"] = np.append(s.massCenter)
            massFits["isotope"] = np.append(isotope)
            for param in ["a1", "a1err", "a2", "a2err","tau", "tauerr"]:
                massFits[param] = np.append(s.fitParams[param])

        # Sort by Mass
        idx = massFits["mass"].ravel().argsort()
        massFits["mass"] = massFits["mass"].ravel()[idx]
        for param in ["a1", "a1err", "a2", "a2err","tau", "tauerr", "acfInternal", "tauInternal", "isotope"]:
            massFits[param] = massFits[param][idx]

        # Filter array
        iin = massFits["acfInternal"]  # True = Use internal fit; False use external
        tin = massFits["tauInternal"]
        mass = massFits["mass"]

        # Data not used to model acf parameters vs. mass
        ax1.errorbar(mass[iin],massFits["a1"][iin],massFits["a1err"][iin])  # add color, size
        ax2.errorbar(mass[iin],massFits["a2"][iin],massFits["a2err"][iin])  # add color, size
        axtau.errorbar(mass[tin],massFits["tau"][tin],massFits["tauerr"][tin])  # add color, size

        # Data not used to model acf parameters vs. mass
        ax1.errorbar(mass[not iin],massFits["a1"][not iin],massFits["a1err"][not iin])
        ax2.errorbar(mass[not iin],massFits["a2"][not iin],massFits["a2err"][not iin])
        axtau.errorbar(mass[not tin],massFits["tau"][not tin],massFits["tauerr"][not tin])

        # Model session-wise Fits by mass

        for axis,param in zip([ax1,ax2,axtau],["a1","a2","tau"]):
            if param == "tau":
                mask = "tauInternal"
            else: mask = "acfInternal"
            x = PolynomialFeatures(self.fitOrder[param]).fit_transform(mass)
            # Replace data not used for Cross-Cal v. Mass Fit with NaN.
            # NaNs are 'missing' data for the regression.
            y = np.where(massFits[mask], massFits[param], np.nan)
            if self.weight[param]:
                w = 1/mass[param+"err"]**2
            else: w = np.ones_like(x)
            w[not mask] = np.nan
            model = sm.WLS(y, x, w, missing = 'drop')  #'drop'  = Ignore NaN
            p, uci, lci = wls_prediction_std(model) # Prediction, upper & lower confidence intervals (95% by default)

        # Iterate through plot axes and parameters to plot modeled fit and confidence intervals
            axis.plot(mass, p, color='black', linestyle='solid', marker= None)
            axis.plot(mass, uci, color='blue', linestyle='dashed', marker= None)
            axis.plot(mass, lci, color='blue', linestyle='dashed', marker= None)

        # Update Isotope Fits with external model results
        for i,isotope in enumerate(massFits["isotope"]):
            for par in ["a1","a2","tau"]:
                session[isotope].fitParam["e"+par] = p[i]
                session[isotope].fitParam["se_e"+par] = (uci[i] - lci[i])/2  #Add standard error

#USER INTERFACE
class ComboACF(QComboBox):
    def __init__(self, parent, seshFit: typing.Type[SessionFits], row = -1):
        super().__init__(parent)
        self.setStyleSheet('font-size: 15px')
        self.addItems(['Internal', 'Model', 'External'])
        self.currentIndexChanged.connect(seshFit.fitsByMass)
        self.row = row

    def getACFType(self):
        print(f'Row:{self.row}: {self.currentText()}')
        # Trigger Recalculation of Model

class ComboTau(QComboBox):
    def __init__(self, parent, seshFit: typing.Type[SessionFits], row= -1):
        super().__init__(parent)
        self.setStyleSheet('font-size: 15px')
        self.addItems(['Internal', 'Model', 'Global'])
        self.currentIndexChanged.connect(seshFit.fitsByMass)
        self.row = row

    def getTauType(self):
        print(f'Row:{self.row}: {self.currentText()}')
        # Trigger Recalculation of Model

class FitTableWidget(QTableWidget):
    def __init__(self, isotopes:list, columns:list):
        self.rowNames = isotopes
        self.colNames = columns
        super().__init__(len(isotopes), len(columns))
        self.setHorizontalHeaderLabels(self.colNames)
        self.verticalHeader().setDefaultSectionSize(25)
        self.horizontalHeader().setDefaultSectionSize(75)
        self.setColumnWidth(1, 125)
        self.setColumnWidth(2, 125)

        for row, isotope in enumerate(self.rowNames):
            acfSelector = ComboACF(self, row)
            tauSelector = ComboTau(self, row)
            isotopeItem = QTableWidgetItem()
            isotopeItem.setText(isotope)
            isotopeItem.setFlags(~Qt.ItemIsEditable) #Values are read-only
            self.setCellWidget(row, 1, acfSelector)
            self.setCellWidget(row, 2, tauSelector)
            self.setItem(row, 0, isotopeItem)

    def populate(self, session: SessionFits()):
        fits = session.sessionFits
        for isotope in list(fits.keys()):
            fit = fits[isotope]
            row = self.rowNames.index(isotope)
            for key in list(fit.keys()):
                col = self.colNames.index(key)
                if key == "ACF Type":
                    acfSelector = ComboACF(self, row)
                    acfSelector.setCurrentIndex(fit[key])
                    self.setCellWidget(row, col, acfSelector)
                elif key == "Tau Type":
                    tauSelector = ComboTau(self, row)
                    tauSelector.setCurrentIndex(fit[key])
                    self.setCellWidget(row, col, tauSelector)
                else:
                    col = self.colNames.index(key)
                    parValue = fit[key]
                    val = QTableWidgetItem()
                    if key in ["a1", "ea1", "se_a1", "se_ea1"]:
                        val.setText(f'{parValue:0.2g}')
                        self.setItem(row, col, val)
                    elif key in ["a2", "ea2", "se_a2", "se_ea2"]:
                        val.setText(f'{parValue:0.2g}')
                        self.setItem(row, col, val)
                    elif key in ["tau", "etau", "se_tau", "se_etau"]:
                        val.setText(f'{parValue:0.2g}')
                        self.setItem(row, col, val)




