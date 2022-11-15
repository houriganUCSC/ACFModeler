import pandas as pd
import numpy as np

class CrossCalFitModel():
    def __init__(self, masses):
        cols = ["ACF Type", "Tau Type", "tau", "se_tau", "a1","se_a1", "a2", "se_a2", "n", "of", "rSqr",
                    "etau", "se_etau", "ea1","se_ea1", "ea2", "se_ea2"]
        data = np.zeros(shape = (len(masses),len(cols)))
        self.ccdf = pd.DataFrame(data, columns=cols, index=masses)
        print(self.ccdf)
        #self.ccdf.sseindex(masses)

    def setCell(self, par, mass, val=np.nan):
        self.ccdf[par][mass] = val

    def populateTree(self):
        for row in range(self.ccdf.count(axis=0)):
            for col in range(self.ccdf.count(axis=1)):
                self.ccdf.iloc[row,col]



