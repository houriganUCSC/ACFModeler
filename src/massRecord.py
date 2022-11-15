import numpy as np

class massRecord:
    """
    Mass record should be sufficiently generic to work on P & A values from any mass spectrometer.
    """
    def __init__(self):
        self.magMass = []
        self.actMass = []
        self.aveMass = 0.0
        self.truncMass = 0.0
        self.channels = 0
        self.dwell = 0
        self.totaldwell = 0
        self.raw = {
            "pulse": [[]],
            "analog": [[]],
            "faraday": [[]]
        }
        self.filtered = {
            "pulse": np.array([]),
            "analog":np.array([]),
            "time":np.array([])
        }
        self.reported = None
        self.massTimeSeries = None
        self.modeledTimeSeries = None

    def processTimeSeries(self, pMax, ACF, FCF=[], ignoreFaraday=True):
        # Multiply raw analog ADC values by stored ACF to get equivalent counts
        analogCounts = (self.raw["analog"].T * ACF).T
        # Filter pulse count array NaNs or P-greater-than-threshold values are replaced with corresponding ACF-scaled analog value
        useAnalog = np.logical_or(self.raw["pulse"] > pMax, np.isnan(self.raw["pulse"]))
        self.reported = np.where(useAnalog, analogCounts, self.raw["pulse"])
        if not ignoreFaraday:
            self.raw["faraday"] = (self.raw["faraday"].T * FCF).T
            self.reported = np.where(np.isnan(self.reported), self.raw["faraday"], self.reported)
        # Average across channels in a cycle ignoring NaNs to produce time-series
        self.massTimeSeries = np.nanmean(self.reported, axis=1)

    def filterRaw(self, pMax, aMin, absT):
        """
        Not currently used. Filtering is done at the session
        :param pMax:
        :param aMin:
        :param absT:
        :return:
        """
        a = self.raw["analog"]
        p = self.raw["pulse"]
        # r, c = np.where(np.logical_and(p<pMax, a>aMin))
        r, c = np.where(np.logical_not(np.isnan(p),np.isnan(a))) # 3x slower
        for i,j in zip(r,c):
            self.filtered["time"] = np.append(self.filtered["time"],absT[i]+self.dwell*j+1)
            self.filtered["analog"] = np.append(self.filtered["analog"],a[i][j])
            self.filtered["pulse"] = np.append(self.filtered["pulse"],p[i][j])

