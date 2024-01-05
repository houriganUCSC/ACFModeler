import numpy as np

import statsmodels.api as sm
from sklearn.preprocessing import PolynomialFeatures

# Project imports
from src.ui.spectrumModelDesignTable import ModelDesignTable

class Session:
    def __init__(self):
        """

        :rtype: object
        """
        self.commonPath = ''
        self.startDir = ''
        self.picklePath = ''
        self.unique = 0
        self.massSpec = ''
        self.massSettle = 0.0005
        self.chSettle = 0.001
        self.isotopes = []
        self.method = None
        self.startTime = None
        self.status = {'imported': False, 'filtered': False, 'fit': False, 'new': False}
        self.pCross = 4E+6
        self.ignoreFaraday = True
        self.isotopeFit = {'algorithm': None, 'norm': None, 'pMax': 5E+6, 'pMin': 0, 'aMin': 1000, 'outlier': 0}
        self.machineDeadTime = 0
        self.inclUnc = False
        self.scanTime = np.array([])
        self.EDAC = np.array([])
        self.FCF = np.array([])
        self.sampleKeys = np.array([])
        self.samples = {}
        self.masses = {}
        self.spectrumFit = SpectrumFits()

    def getChromData(self, postProcessed = False):
        for smpName, smpRecord in self.samples.items():
            mask = np.where(self.sampleKeys == smpRecord.ID)
            time = self.scanTime[mask]
            time = time - np.amin(time)
            chromData = np.expand_dims(time, axis=1)
            chromHdr = ['Time']
            for massName, massObj in self.masses.items():
                if not postProcessed:
                    chromData = np.append(chromData, massObj.timeSeries[mask], axis=1)
                else:
                    chromData = np.append(chromData, massObj.modeledTimeSeries[mask], axis=1)
                chromHdr.append(massName)
            # Todo: Create Export formatted file.

    def plotTimeSeries(self, sampleName, dynamic_ax, postProcessed = False):

        if sampleName in self.samples.keys():
            primaryKey = self.samples[sampleName].ID
            isotopes = []
            dynamic_ax.cla()
            mask = np.where(self.sampleKeys == primaryKey)
            x = self.scanTime[mask]
            x = x - np.min(x)
            for massName, massObj in self.masses.items():
                if postProcessed:
                    y = massObj.modelTimeSeries[mask]
                else:
                    y = massObj.timeSeries[mask]
                isotopes.append(massName)
                dynamic_ax.semilogy(x, y)
            dynamic_ax.set_ylabel("CPS")
            dynamic_ax.set_xlabel("Time(sec)")
            dynamic_ax.set_title(sampleName)
            dynamic_ax.legend(isotopes)
            dynamic_ax.figure.canvas.draw()

    def calculateTimeSeries(self):
        if self.status["imported"] == True:
            for massName, massRecord in self.masses.items():
                massRecord.calculateTimeSeries(self)

    def filterRawData(self):
        if self.status["imported"] == True:
            for massName, massRecord in self.masses.items():
                massRecord.filter()
            print(massName, len(massRecord.timeSeries))
            relTime = self.scanTime - self.startTime
            acfKey = list(self.masses.keys())[0]
            acf = self.masses[acfKey].ACF
            X = np.ones((len(acf), 2))
            X[:, 1] = relTime
            model = sm.OLS(acf, X)
            results = model.fit()
            acf0 = results.params[0]
            acfF = acf0 + results.params[1] * max(relTime)
            drift = acfF / acf0 - 1
            self.machineACF0 = acf0
            self.machineDrift = drift
            # print(f'{__name__} Machine ACF0: {acf0:0.2f}, Drift: {drift:0.1%}')

    def regressRawData(self):
        if self.status["filtered"] == True:
            for massName, massRecord in self.masses.items():
                massRecord.regress(self, massName)
            self.status['fit'] = True

    def releaseRawDataRegression(self):
        for massRecord in self.masses.values():
            massRecord.fits = MassFits()
        self.status['fit'] = False

    def regressSpectrum(self):
        keys = ["a1", "a2", "tau"]
        uses = ["useACF", "useACF", "useTau"]
        """ EXTRACT DATA TO FIT """
        self.spectrumFit = SpectrumFits()
        spFit = self.spectrumFit
        for massName, massRecord in self.masses.items():
            massFits = massRecord.fits
            # Get Values that are checked "use" and append to X, Y, seY arrays
            # unchecked values are replaced by np.nan
            # Inverse operation for outY and outSeY arrays
            spFit["Mass"].append(massRecord.aveMass)
            for key, use in zip(keys, uses):
                if massRecord[use]:
                    spFit[key]['in']["Y"].append(massFits["self"][key])
                    spFit[key]['in']['seY'].append(massFits["self"]["se_" + key])
                    spFit[key]["out"]["Y"].append(np.nan)
                    spFit[key]["out"]["seY"].append(np.nan)
                else:
                    spFit[key]['out']["Y"].append(massFits["self"][key])
                    spFit[key]['out']['seY'].append(massFits["self"]["se_" + key])
                    spFit[key]["in"]["Y"].append(np.nan)
                    spFit[key]["in"]["seY"].append(np.nan)

            """ PERFORM FITTING ROUTINE """
            # Cycle through a1, a2, tau
            for key in keys:
                order = spFit[key]["order"]
                if order <= 3:
                    polyFeat = PolynomialFeatures(degree=order)
                    xP = polyFeat.fit_tansform(spFit["Mass"])
                    if self[key]["weighted"]:
                        w = 1/(spFit[key]["in"]["seY"]**2)
                    else: w = np.ones_like(spFit["Mass"])
                y = spFit[key]["in"]["Y"]
                if self[key]["robust"]:
                    model = sm.RLM(y, xP, M=sm.robust.norms.TukeyBiweight, missing='drop').fit()
                else:
                    model = sm.WLS(y, xP, weights=w, missing='drop').fit()
            else:
                pass #Todo: Implement Cubic Spline
            pred = model.get_prediction(xP)

            self[key]["best"] = pred.summary_frame()["mean"]
            self[key]["seBest"] = pred.summary_frame()["mean_se"]
            self[key]["lci"] = pred.summary_frame()["obs_ci_lower"] #Alternate "mean_ci_lower"
            self[key]["uci"] = pred.summary_frame()["obs_ci_upper"] #Alternate "mean_ci_upper"

            # Reduced Chi^2
            wrss = np.sum(w*(y-self[key]["best"])**2/self[key]["seY"])  # Sum of squares of weighted residuals
            v = np.sum(w) / ((np.sum(w)) ** 2 - np.sum(w ** 2))  #weighted normalizing factor

            self[key]["coefs"] = model.params
            self[key]["seCoefs"]= model.bse
            self[key]["rSqr"] = model.rsquared
            self[key]["rChi2"] = wrss/v


        """ UPDATE MODEL PARAMETERS FOR POST PROCESSING"""
        types = ["ACFType", "ACFType", "TauType"]
        for massName, massRecord in session.filterRecords.items():
            for i, [key, type] in enumerate(zip(keys, types)):
                """
                "model" is an list of value for each parameter. For example if key is "a1" and typeKey is "internal"
                then "model" will contain the fitted value for "a1" based on other isotopes in the mass spectrum. 
                if typeKey is "self" then "model" will contain the auto (self) cross calibrated value
                """
                typeKey = massRecord[type]
                Y = self[key]["best"][i]
                seY = self[key]["seBest"][i]
                if typeKey == 'self':
                    if not np.isnan(seY = self[key]["Y"][i]):  # If 'self' calibrated and not NaN
                        Y = self[key]["Y"][i]
                        seY = self[key]["seY"][i]
                self[key]["model"].append(Y)
                self[key]["seModel"].append(seY)
                # Todo: implement "external" type (i.e. use machine values)

    def postProcessTimeSeries(self):
        if self.status["modeled"] == True:
            for massName, massRecord in self.masses.items():
                massRecord.postProcessTimeSeries(self)

    def commitSpectrumFits(self):
        pass

    def timeOffsets(self, massSettle = 0, chSettle = 0):
        if self.status["imported"] == True:
            delta = 0.00
            for massName, massRecord in self.masses.items():
                recordShape = massRecord.pulse.shape
                massRecord.cycles = recordShape[0]
                massRecord.nObs = np.prod(recordShape)
                massRecord.timeOffset = delta + self.massSettle
                delta = delta + massRecord.channels*(massRecord.chDwell + self.chSettle) + self.massSettle

    def updateDeadTime(self):
        key = list(self.samples.keys())[0]
        try:
            dt = self.samples[key].metaData["deadTime"]
        except:
                dt = np.nan
        self.machineDeadTime = dt
        return dt


class Sample:
    def __init__(self):
        self.ID = None
        self.group = ''
        self.name = ''
        self.fileTimes = {}
        self.filePaths = {}
        self.metaData = {}
        self.isotopes = None
        self.session = None


class Mass:
    def __init__(self, session: Session):
        self.session = session
        self.timeOffset = 0
        self.chDwell = 0
        self.totalDwell = 0
        self.channels = 0
        self.cycles = 0
        self.nObs = 0
        self.nIn = 0
        self.nQual = 0
        self.anOnly = 0
        self.maxP = 0
        self.dtMaxP = 0
        self.dtCorrPct = 0
        self.fits = None
        self.magMass = []
        self.actMass = []
        self.aveMass = 0.0
        self.smpID = np.array([])
        self.ACF = np.array([])
        self.pulse = np.array([[]])
        self.analog = np.array([[]])
        self.faraday = np.array([[]])
        self.reported = np.array([[]])
        self.timeSeries = np.array([])
        self.modeledTimeSeries = np.array([])
        self.filteredTime = np.array([])
        self.filteredPulse = np.array([])
        self.filteredAnalog = np.array([])
        self.anOnlyTime = np.array([])
        self.fits = MassFits()
        self.postProcessPars = {"tauSource": None, "tau": None, "setau": None,
                                "alphaSource": None, "a1": None, "sea1": None, "a2": None, "sea2": None}


    def filter(self):
        pMax = self.session.isotopeFit["pMax"]
        pMin = self.session.isotopeFit["pMin"]
        aMin = self.session.isotopeFit["aMin"]
        if pMin == 0:
            pMin = aMin * np.nanmean(self.ACF)
        outlier = self.session.isotopeFit["outlier"]
        if outlier == 0:
            outlier == 10
        chs = self.channels
        time2d = np.tile(self.session.scanTime,(chs, 1)).T
        offsets = np.arange(0, chs) * (self.chDwell + self.session.chSettle)
        offsets = offsets + self.timeOffset
        time2d = time2d + offsets
        t = time2d.flatten()
        print(self.aveMass, t)
        p = self.pulse.flatten()
        a = self.analog.flatten()
        print(self.aveMass, len(a), p/a)
        anOnlyMask = np.logical_and(np.isnan(p), np.isreal(a))
        self.anOnlyTime = t[anOnlyMask]
        self.anOnly = np.sum(anOnlyMask)
        mask = np.where(np.logical_and(pMax > p, p > pMin, a > aMin))
        self.filteredTime = t[mask]
        self.filteredPulse = p[mask]
        print(self.aveMass, len(a), len(self.filteredPulse), self.filteredPulse)
        self.filteredAnalog = a[mask]
        self.nQual = len(self.filteredTime)     # Measurements in fitting range
        self.filteredAnalog[self.filteredAnalog == 0] = np.nan
        acf = self.filteredPulse / self.filteredAnalog
        #acf[np.isinf(acf)] = np.nan
        if outlier > 0:
            med = np.nanmedian(acf)
            Q3 = np.nanpercentile(acf, 75, interpolation='midpoint')  # Upper quartile
            Q1 = np.nanpercentile(acf, 25, interpolation='midpoint')  # Lower quartile
            iqr = Q3-Q1 # Interquartile range
            mask = []
            if iqr > 0:
                mask = np.where(np.abs((acf-med)/iqr) <= outlier)    # True if under gross outlier threshold
            self.filteredPulse = self.filteredPulse[mask]
            self.filteredAnalog= self.filteredAnalog[mask]
            self.filteredTime = self.filteredTime[mask]
        self.nIn = len(self.filteredTime)       # Measurements that pass Tukey filtered acf values
        if self.nQual > 0:
            try:
                self.maxP = max(self.filteredPulse)
            except:
                # print(f'{__name__} filter problem')
                pass

    def regress(self, session: Session, mass: str):
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
        if self.fits is None:
            self.fits = MassFits()
        regType = self.session.isotopeFit["algorithm"]
        normType = self.session.isotopeFit['norm']
        deadtime = self.session.machineDeadTime
        if self.nIn > 50:
            P = self.filteredPulse
            P = P / (1 + P * deadtime)  # Remove Deadtime Correction
            A = self.filteredAnalog
            t = self.filteredTime - self.session.startTime  # Relative to start of experiment
            W = self.totalDwell * P ** 3
            y = 1 / P
            X = np.ones((len(y), 3))
            X[:, 1] = 1 / A
            X[:, 2] = t / A
            if regType == "Robust":
                M = sm.robust.norms.TukeyBiweight()
                norm = normType.replace(" ","")
                if norm == "Whitened":
                    w = len(W) * np.sqrt(W) / np.sum(np.sqrt(W))
                    y = w * y
                    X = w * X
                elif norm == "HuberT": M = sm.robust.norms.HuberT()
                elif norm == "Hampel": M = sm.robust.norms.Hampel()
                elif norm == "LeastSquares": M = sm.robust.norms.LeastSquares()
                elif norm == "AndrewWave": M = sm.robust.norms.AndrewWave()
                elif norm == "RamsayE": M = sm.robust.norms.RamsayE()
                elif norm == "RobustNorm": M = sm.robust.norms.RobustNorm()
                elif norm == "TrimmedMean": M = sm.robust.norms.TrimmedMean()
                elif norm == "TukeyBiweight": M = sm.robust.norms.TukeyBiweight()
                else: M = sm.robust.norms.TukeyBiweight()
                model3D = sm.RLM(y, X, M)

            elif "Weighted" in regType:
                model3D = sm.WLS(y, X, weights=W)
                x2D = np.ones((len(y), 2))
                x2D = np.delete(X, 2, 1)
                model2D = sm.WLS(y, x2D, weights=W)
                results2D = model2D.fit()
                pred2D = results2D.predict(x2D)
                ci = results2D.conf_int()
                self.fits["2D"] = {}
                self.fits["2D"]["tau"] = results2D.params[0]
                self.fits["2D"]["a1"] = results2D.params[1]
                self.fits["2D"]["dtau"] = results2D.bse[0]
                self.fits["2D"]["da1"] = results2D.bse[1]

            elif "Ordinary" in regType:
                model3D = sm.OLS(y, X)


            # Calculate Reduced Chi-Square (aka MSWD)
            results = model3D.fit()
            pred = results.predict(X)
            res = y - pred
            # https://en.wikipedia.org/wiki/Reduced_chi-squared_statistic
            redChi2 = np.sum(W * res ** 2)/(len(res)-3)
            a1 = results.params[1]
            a2 = results.params[2]
            initialACF = 1/a1
            delTime = np.max(self.session.scanTime) - np.min(self.session.scanTime)
            finalACF = 1/(a1 + a2 * delTime)
            drift = (finalACF/initialACF)-1
            self.dtMaxP = self.maxP/(1-self.maxP*results.params[0])
            self.dtCorrPct = self.dtMaxP/self.maxP-1
            rse = []
            for i, [param, err] in enumerate(zip(results.params, results.bse)):
                rse.append(err/param)
            if "Robust" in regType:
                r2 = 'N/A`'
            else:
                r2 = {f'{results.rsquared:0.3f}'}
                # print(f'{mass}: ACF - {initialACF:0.2f}({100*rse[1]:0.2G}%), Drift = {drift:0.2%}({100*rse[2]:0.2G}%);'
                #   f' Tau = {results.params[0]:0.4G}({rse[0]:0.1%}); rSqr = {r2};'
                #   f' Dead Corr % at Pmax ({self.maxP}) = {self.dtCorrPct:0.2%}')

            self.fits["self"]['tau'] = results.params[0]
            self.fits["self"]["a1"] = results.params[1]
            self.fits["self"]["a2"] = results.params[2]
            self.fits["self"]["se_tau"] = results.bse[0]
            self.fits["self"]["se_a1"] = results.bse[1]
            self.fits["self"]["se_a2"] = results.bse[2]
            if "Robust" in regType:
                self.fits["self"]["rSqr"] = 'N/A'
            else:
                self.fits["self"]["rSqr"] = results.rsquared
            self.fits["self"]["redChi2"] = redChi2
            self.fits["self"]["ACF"] = initialACF
            self.fits['self']["Drift"] = drift
            self.fits["self"]["dtCorrPct"] = self.dtCorrPct
        else:
            self.fits["self"]['tau'] = np.nan
            self.fits["self"]["a1"] = np.nan
            self.fits["self"]["a2"] = np.nan
            self.fits["self"]["se_tau"] = np.nan
            self.fits["self"]["se_a1"] = np.nan
            self.fits["self"]["se_a2"] = np.nan
            self.fits["self"]["rSqr"] = np.nan
            self.fits["self"]["redChi2"] = np.nan
            self.fits["self"]["ACF"] = np.nan
            self.fits['self']["Drift"] = np.nan
            self.fits['self']["dtCorrPct"] = np.nan

    def calculateTimeSeries(self):
        # Multiply raw analog ADC values by stored ACF to get equivalent counts
        pCross = self.session.pCross # Cross-over to analog counts
        analogCounts = (self.analog.T * self.ACF).T
        # Filter pulse count array NaNs or P-greater-than-threshold values are replaced with corresponding ACF-scaled analog value
        reported = np.where(pCross>self.pulse, self.pulse, analogCounts)
        if not self.session.ignoreFaraday:
            reported = np.where(np.isnan(self.reported), (self.faraday.T * self.FCF).T, self.reported)
        self.timeSeries = np.nanmean(reported, axis=1)

    def postProcessTimeSeries(self):
        pCross = self.session.pCross
        inclUnc = self.session.inclUnc
        t = self.session.scanTime - np.amin(self.session.scanTime)
        acfType = self.fits['acfType']
        fits = self.postProcessPars
        if fits['alphaSource'] in ['Self', 'Internal']:
            a1 = fits['a1']
            a2 = fits['a2']
            se_a1 = fits['sea1']
            se_a2 = fits['sea2']
            acf = 1 / (a1 + a2 * t)
            modelAnalog = (self.analog.T * acf).T
            if inclUnc:
                B = acf ** 2
                dA = np.sqrt(se_a1 ** 2 + (t * se_a2) ** 2)
                dA = (self.analog.T * B * dA).T
            else:
                dP = None
        else:
            modelAnalog = (self.analog.T * self.ACF).T
            dP = None
        if fits['tauSource'] in ['Self', 'Internal']:
            tau = fits['tau']
            se_tau = fits['setau']
            pUncorr = self.pulse / (1 + self.pulse * self.session.machineDeadTime)
            modelPulse = pUncorr / (1 - pUncorr * tau)
            if inclUnc:
                D = 1 / (1 - self.pulse * tau)
                dP = np.sqrt(self.pulse / self.dwell + (self.pulse ** 2 * se_tau) ** 2)
                dP = D * dP
            else:
                dP = None
        else:
            modelPulse = self.pulse
            dP = None
        tmp = np.where(pCross > modelPulse, modelPulse, modelAnalog)
        self.modeledTimeSeries = np.nanmean(tmp, axis=1)
        if not dP is None and not dA is None:
            self.modeledTimeSeriesError = np.nanmean(np.where(modelPulse > pCross, dA, dP), axis = 1)

    def count(self):
        self.nObs = self.pulse.size


class MassFits(dict):

    def __init__(self):
        super().__init__()
        self['useTau'] = True
        self['useACF'] = True
        self['acfType'] = 'self'
        self['tauType'] = 'internal'
        self['self'] = {'a1': None, 'se_a1': None, 'a2': None, 'se_a2': None, 'tau': None, 'se_tau': None,
                         'rSqr': None, 'redChi2': None, 'ACF': None,'Drift': None}
        self['internal'] = {'a1': None, 'a2': None, 'tau': None,
                            'se_a1': None, 'se_a2': None, 'se_tau': None,
                            'lci': None, 'uci': None}
        self['external'] = {'a1': None, 'a2': None, 'tau': None}


class SpectrumFits(dict):
    def __init__(self):
        super().__init__()
        self['mass'] = []
        self['nObs'] = []
        self['nQual'] = []
        self['anOnly'] = []
        self['rSqr'] = []
        self['redChi2'] = []
        self['alphaSource'] = []
        self['tauSource'] = []
        self['pMax'] = []
        self['a1'] = {'order': 3,
                      'weighted': True,
                      'robust': False,
                      'mask': [],
                      'Y':[],
                      'seY':[],
                      'best': {'Y':[], 'seY':[]},
                      'ci': {'u':[], 'l':[]},
                      'coefs':[],
                      'se_coefs': [],
                      'intText':''}
        self['a2'] = {'order': 3,
                      'weighted': True,
                      'robust': False,
                      'mask': [],
                      'Y':[],
                      'seY':[],
                      'best': {'Y':[], 'seY':[]},
                      'ci': {'u':[], 'l':[]},
                      'coefs':[],
                      'se_coefs': [],
                      'intText':''}
        self['tau'] = {'order': 3,
                       'weighted': True,
                       'robust': False,
                       'mask': [],
                       'Y':[],
                       'seY':[],
                       'out': {'Y':[], 'seY':[]},
                       'best': {'Y':[], 'seY':[]},
                       'ci': {'u':[], 'l':[]},
                       'coefs':[],
                       'se_coefs': [],
                       'intText':''}

    def getMassFits(self, session: Session):
        for isotopes, mass in session.masses.items():
            self['mass'].append(mass.aveMass)
            self['nObs'].append(mass.nObs)
            self['nQual'].append(mass.nQual)
            self['anOnly'].append(mass.anOnly)
            self['rSqr'].append(mass.fits['self']['rSqr'])
            self['redChi2'].append(mass.fits['self']['redChi2'])
            self['pMax'].append(mass.dtMaxP)
            self['tauSource'].append(None)
            self['alphaSource'].append(None)
            for par in ['a1', 'a2', 'tau']:
                self[par]['mask'].append(True)
                self[par]['Y'].append(mass.fits['self'][par])
                self[par]['seY'].append(mass.fits['self']['se_'+par])


    def fitByMass(self, modelSetUp: ModelDesignTable):
        modelDesign = modelSetUp.getFitPars()
        for i, par in enumerate(['a1', 'a2', 'tau']):
            a = 0.05 # 95% CI
            order = modelDesign[par]['order']
            """ EXTRACT INLIER AND OUTLIER DATA """
            # Inliers

            x = np.array(self['mass'])
            mask = self[par]['mask']
            y = np.array(np.where(mask, self[par]['Y'], np.nan))
            yerr = np.array(np.where(mask, self[par]['seY'], np.nan))

            """ CURVE FITTING """
            # Transform to design matrix of specified order
            polyFeat = PolynomialFeatures(degree=order)
            xP = polyFeat.fit_transform(x.reshape(-1, 1))

            # Get model fitting parameters
            weighted = modelDesign[par]['weighted']
            w = {True: 1/yerr**2, False: np.ones_like(x)}[weighted]
            #Todo: Add Spline option
            #Todo: Test Robust
            robust = modelDesign[par]['robust']
            if robust:
                model = sm.RLM(y, xP, M=sm.robust.norms.HuberT(), missing='drop').fit()
            else:
                model = sm.WLS(y, xP, weights = w, missing='drop').fit()
            pred = model.get_prediction(xP)

            ypred = pred.summary_frame(a)["mean"]
            yprederr = pred.summary_frame(a)["mean_se"]
            rse = np.abs(100*(yprederr/ypred))
            numForms = ['0.5f', '0.3E', '0.3E']
            intText = [f'{yp:{numForms[i]}} (±{err:0.1f}%)' for yp, err in zip(ypred,rse)]
            intText = [it.replace("E-0", "E-").replace("E+0", "E+") for it in intText]
            self[par]['intText'] = intText
            self[par]['best']['Y'] = ypred
            self[par]['best']['seY'] = yprederr
            self[par]['ci']['l'] = pred.summary_frame()["mean_ci_lower"] #Alternate "mean_ci_lower"
            self[par]['ci']['u'] = pred.summary_frame()["mean_ci_upper"] #Alternate "mean_ci_upper"
            # Todo: use mean ci or obs ci?

            self[par]['coefs'] = model.params
            self[par]['se_coefs'] = model.bse
            self[par]['rSqr'] = model.rsquared

            # Reduced Chi^2
            wrss = np.sum(w*(y-ypred)**2/yerr)  # Sum of squares of weighted residuals
            v = np.sum(w) / ((np.sum(w)) ** 2 - np.sum(w ** 2))  #weighted normalizing factor
            self[par]['rChi2'] = wrss/v






