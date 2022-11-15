import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import pickle
import pandas as pd
import statsmodels.api as sm

class regressed():

    def __init__(self):
        self.deadtime = 6E-9
        self.aMin = 2000
        self.pMax = 5E+6
        self.pCross = 2.5E+6
        self.pMin = 50000
        self.tukeyOut = 3.5
        self.td = 0.010
        self.reports = {}
        self.summaryStats {}
        self.r = None
        self.isotopes = []
        self.startT
        self.f = {}

    def timeZero(self):
        self.startT = min([min(self.r[isotope]["time"]) for isotope in self.isotopes])

    def filterRange(self):
        self.timeZero()
        for isotope in self.isotopes
            p = self.r[isotope]["pulse"]
            a = self.r[isotope]["analog"]
            t = self.r[isotope]["time"]
            inRange = np.logical_and(p < self.pMax, p > self.pMin)
            inRange = np.logical_and(inRange, a > aMin)
            self.f[isotope].p = p[inRange]
            self.f[isotope].a = a[inRange]
            self.f[isotope].t = a[inRange]- self.startT
            self.f[isotope].acf = self.f[isotope].p/self.f[isotope].a

    def filterOutlier(self, tukey = None):
        if tukey == None:
            tukey = self.tukeyOut
        else: self.tukeyOut = tukey
        for isotope in self.isotopes:
            med = np.nanmedian(self.f[isotope].acf)
            Q3 = np.nanpercentile(self.f[isotope].acf,75,interpolation='midpoint')
            Q1 = np.nanpercentile(self.f[isotope].acf,25,interpolation='midpoint')
            iqr = Q3-Q1
            self.f[isotope].inlier = np.abs(self.f[isotope].acf-med)/iqr <= tukey

    def regress(self, drift = True, reject = True):
        session = pickle.load(open( "save.p", "rb" ))
        self.r = session.record
        self.isotopes = self.r.keys()
        self.filterRange()
        self.filterOutlier(3)
        for isotope in self.isotopes:
            iso = isotope.strip('\"')
            print(f'Analyzing {iso}...')
            inlier = self.f[isotope].inlier
            p = self.f[isotope].p
            a = self.f[isotope].a
            t = self.f[isotope].t
            inRange = self.f[isotope].inlier
            uncP = p/(1+p*self.deadtime)
            w = uncP**3*self.td
            if drift: colsX = 3
            else: colsX = 2
            regy = 1/uncP
            regX = np.ones([len(regy),3])
            regX[:,1] = 1/a
            regX[:,2] = (t/a)
            rmodel = sm.RLM(regy,regX,M=)
            if not reject:
                inRange = np.full_like(inRange,True)
            regy = 1/uncP[inRange]
            regX = np.ones([len(regy),colsX])
            regX[:,1] = 1/a[inRange]
            if drift:
                regX[:,2] = (t/a)[inRange]
            wmodel = sm.WLS(regy, regX,weights=w)


            results = model.fit()
            #print(results.summary())
            maxP = max(unCorrP)
            maxCorrP = pCross/(1-pCross*results.params[0])
            pCorrPct = (maxCorrP/pCross -1)*100
            acf_initial = 1/results.params[1]
            acf_final = 1/(results.params[1]+delT*results.params[2])
            del_acf = (acf_final/acf_initial -1)*100
            del_acf_rate = del_acf/delT
            print(f'{iso} Deadtime: {results.params[0]:0.02e}')
            print(f'{iso} Max Correction: {pCorrPct:0.1f}%')
            print(f'{iso} ACF initial: {acf_initial:0.02f}')
            print(f'{iso} ACF final: {acf_final:0.02f}')
            print(f'{iso} delta ACF: {del_acf:0.03f}%')
            print(f'{iso} delta ACF rate: {del_acf_rate:0.03e}%/sec')

            regX = np.ones([len(regy),2])
            regX[:,1] = 1/fA[inRange]
            model = sm.WLS(regy, regX,weights=w)
            results = model.fit()
            #print(results.summary())
            acf_initial = 1/results.params[1]
            maxCorrP = pCross/(1-pCross*results.params[0])
            pCorrPct = (maxCorrP/pCross -1)*100
            print(f'{iso} Deadtime: {results.params[0]:0.02e}')
            print(f'{iso} Max Correction: {pCorrPct:0.1f}%')
            print(f'{iso} ACF initial: {acf_initial:0.02f}')
            #plt.show()
